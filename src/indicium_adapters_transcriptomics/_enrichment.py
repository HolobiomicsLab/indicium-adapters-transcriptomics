"""Shared enrichment back-end for transcriptomics adapters.

In v0.1.1 this module loads the bundled gene-name index TSV at import time
and exposes a ``GENE_INDEX`` singleton with a case-insensitive ``lookup()``
method. The index is small (~20 entries) and hand-curated from the two
benchmark reference papers (Warner et al. 2025, Cheung et al. 2024); a
maintainer script for full HGNC + Nemve1 regeneration arrives in v0.1.2.

The lookup accepts any of three input shapes and returns the matching row
as a dict (or ``None``):

- A bare gene symbol (case-insensitive): ``"WNT4"`` or ``"wnt4"``.
- A Nemve1 / JGI gene-model id (case-sensitive): ``"gw.51.96.1"``.
- A UniProt accession (case-sensitive): ``"P49338"``.
"""
from __future__ import annotations

import csv
import re
from importlib.resources import files
from typing import Final

# Columns of `data/gene_index.tsv` — declared here so the loader can be
# strict about the bundled schema. Adding columns means bumping this list
# and the lookup logic in tandem.
_GENE_INDEX_COLUMNS: Final[tuple[str, ...]] = (
    "symbol",
    "hgnc",
    "mgi",
    "flybase",
    "wormbase",
    "zfin",
    "ensembl",
    "ncbigene",
    "nemve1",
    "nvertx",
    "uniprot",
    "species",
    "description",
)


# Recognised species tags. Order encodes the fallback priority used by
# ``GeneIndex.lookup`` when neither an explicit hint nor the casing
# heuristic can disambiguate a collision. Human first because it
# dominates transcriptomics-paper claims; mouse second; then the
# v0.6.2 model organisms in rough corpus-prevalence order
# (zebrafish > Drosophila > C. elegans for transcriptomics data sets).
_KNOWN_SPECIES: Final[tuple[str, ...]] = (
    "human", "mouse", "zebrafish", "drosophila", "celegans", "nematostella",
)


def _load_gene_index_rows() -> list[dict[str, str]]:
    tsv = files("indicium_adapters_transcriptomics.data").joinpath("gene_index.tsv")
    with tsv.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if tuple(reader.fieldnames or ()) != _GENE_INDEX_COLUMNS:
            raise RuntimeError(
                "gene_index.tsv header does not match the declared schema. "
                f"Expected {_GENE_INDEX_COLUMNS!r}; got "
                f"{tuple(reader.fieldnames or ())!r}."
            )
        return [{k: (v or "").strip() for k, v in row.items()} for row in reader]


# C. elegans symbol pattern: three lowercase letters + hyphen + 1+
# alphanumeric chars. Examples: ``unc-22``, ``daf-16``, ``let-7``,
# ``gld-1``. Very distinctive — no other organism uses this layout.
_C_ELEGANS_PATTERN: Final[re.Pattern[str]] = re.compile(
    r"^[a-z]{3,4}-[a-zA-Z0-9]+$"
)


def _infer_species_from_casing(token: str) -> str | None:
    """Casing heuristic used when no explicit species hint is given.

    Conventions in the field (per-organism):

    * **Human (HGNC)** — ALL UPPERCASE: ``MAPK1``, ``WNT4``, ``CD3E``.
    * **Mouse (MGI)** — Title-cased: ``Mapk1``, ``Wnt4``, ``Cd3e``.
    * **C. elegans (WormBase)** — three (occasionally four) lowercase
      letters + hyphen + alphanumeric: ``unc-22``, ``daf-16``,
      ``let-7``. Very distinctive — no other organism uses this
      layout.
    * **Zebrafish (ZFIN)** and **Drosophila (FlyBase)** — both use
      lowercase / mixed-case conventions that are NOT
      machine-distinguishable from a casual lowercase / Title-case
      typo of a human or mouse symbol. They are intentionally
      **not** matched by the heuristic — claims targeting these
      organisms must pass ``species="zebrafish"`` or
      ``species="drosophila"`` explicitly. Without a hint a
      lowercase ``mmp9`` falls through to the load-order fallback
      and resolves to human via the cross-listed rule, preserving
      v0.1.x case-insensitive forgiveness.
    * **Nematostella** (Nemve1) ids carry dots (``gw.51.96.1``) and
      never come through this path — they hit the by-nemve1 lookup
      branch instead.

    Returns ``None`` when the casing matches no high-confidence
    convention.
    """
    if not token or not token[0].isalpha():
        return None
    # C. elegans rule first — most distinctive (mandatory hyphen).
    if _C_ELEGANS_PATTERN.match(token):
        return "celegans"
    if token.isupper():
        return "human"
    if token[0].isupper() and not token.isupper():
        return "mouse"
    # All-lowercase is ambiguous (zebrafish canonical / drosophila
    # canonical / lowercase typo of a human or mouse symbol) — defer
    # to the species hint or load-order fallback. See class docstring.
    return None


class GeneIndex:
    """Case-insensitive gene-symbol lookup with v0.6.0 species disambiguation.

    For each symbol we store the *list* of rows that share the
    lower-cased key — typically a single human row, but ``MAPK1`` /
    ``Mapk1`` resolve to the same key and so produce a 2-row list
    (one human, one mouse). ``lookup`` then picks among them using:

    1. An explicit ``species`` parameter from the caller (the
      "hint" path — see ``BulkTranscriptomicsAdapter.enrich_claim``
       passing ``claim.get("species")``).
    2. A casing heuristic (UPPER → human, Title-case → mouse) —
       safe because the conventions are universal in the field.
    3. Load order, which puts human first
       (see ``_KNOWN_SPECIES``) so collisions resolve to human by
       default when neither the hint nor the casing can disambiguate.

    Accession-style lookups (Nemve1, UniProt, MGI id) are still
    single-row dicts — accessions are unambiguous by construction.
    """

    def __init__(self, rows: list[dict[str, str]]):
        self._rows = rows
        self._by_symbol: dict[str, list[dict[str, str]]] = {}
        for row in rows:
            symbol = row.get("symbol")
            if symbol:
                self._by_symbol.setdefault(symbol.lower(), []).append(row)
        self._by_nemve1: dict[str, dict[str, str]] = {
            row["nemve1"]: row for row in rows if row.get("nemve1")
        }
        self._by_uniprot: dict[str, dict[str, str]] = {
            row["uniprot"]: row for row in rows if row.get("uniprot")
        }
        # MGI accession lookup — added in v0.6.0. CURIE form ``MGI:NNN``.
        self._by_mgi: dict[str, dict[str, str]] = {
            row["mgi"]: row for row in rows if row.get("mgi")
        }
        # v0.6.2 organism-specific accession indices. Stored under the
        # bare id (FBgn0003079) and also under the prefixed CURIE form
        # (FlyBase:FBgn0003079) so callers can pass either spelling.
        self._by_flybase: dict[str, dict[str, str]] = {}
        for row in rows:
            if v := row.get("flybase"):
                self._by_flybase[v] = row
                self._by_flybase[f"FlyBase:{v}"] = row
        self._by_wormbase: dict[str, dict[str, str]] = {}
        for row in rows:
            if v := row.get("wormbase"):
                self._by_wormbase[v] = row
                self._by_wormbase[f"WormBase:{v}"] = row
        self._by_zfin: dict[str, dict[str, str]] = {}
        for row in rows:
            if v := row.get("zfin"):
                self._by_zfin[v] = row
                self._by_zfin[f"ZFIN:{v}"] = row

    def lookup(
        self, name: str, *, species: str | None = None
    ) -> dict[str, str] | None:
        """Return the matching row as a dict, or ``None``.

        ``species`` is an optional hint: ``"human"`` / ``"mouse"`` /
        ``"nematostella"``. When given, it disambiguates symbol-level
        collisions (``MAPK1`` ↔ ``Mapk1``). When not given, a casing
        heuristic infers it; if still ambiguous, falls back to the
        load order (human first).

        Lookup order: gene symbol (case-insensitive, with species
        disambiguation) → Nemve1 / JGI id → UniProt accession →
        MGI accession. The accession paths are single-row hits and
        bypass the disambiguation logic.
        """
        if not isinstance(name, str):
            return None
        token = name.strip()
        if not token:
            return None

        # Symbol path — may return multi-row lists; disambiguate.
        candidates = self._by_symbol.get(token.lower(), [])
        if candidates:
            return self._pick_candidate(candidates, token, species)

        # Accession paths — exact, single-row, case-sensitive.
        if row := self._by_nemve1.get(token):
            return row
        if row := self._by_uniprot.get(token):
            return row
        if row := self._by_mgi.get(token):
            return row
        # v0.6.2 model-organism accession paths.
        if row := self._by_flybase.get(token):
            return row
        if row := self._by_wormbase.get(token):
            return row
        if row := self._by_zfin.get(token):
            return row
        return None

    @staticmethod
    def _pick_for_species(
        candidates: list[dict[str, str]], wanted: str
    ) -> dict[str, str] | None:
        """Return the first candidate that serves the wanted species.

        A row "serves" species ``X`` when either its ``species`` tag is
        ``X``, or it carries the canonical CURIE column for ``X``
        (HGNC for human, MGI for mouse, Nemve1 for nematostella).
        The latter catches Nematostella-anchored rows that ship a
        human ortholog under ``hgnc`` — a lookup of ``MMP9`` should
        return such a row when the human species is wanted, even
        though the row's primary tag is ``nematostella``."""
        canonical_col = {
            "human":        "hgnc",
            "mouse":        "mgi",
            "drosophila":   "flybase",
            "celegans":     "wormbase",
            "zebrafish":    "zfin",
            "nematostella": "nemve1",
        }.get(wanted)
        # First pass: prefer an exact species-tag match (rows whose
        # primary species is the wanted one).
        for row in candidates:
            if row.get("species") == wanted:
                return row
        # Second pass: any row that carries the wanted species' canonical
        # CURIE column — covers cross-listed rows (Nematostella anchor
        # + human ortholog under hgnc).
        if canonical_col:
            for row in candidates:
                if row.get(canonical_col):
                    return row
        return None

    @staticmethod
    def _pick_candidate(
        candidates: list[dict[str, str]],
        token: str,
        species: str | None,
    ) -> dict[str, str]:
        """Apply the v0.6.0 disambiguation rules to a symbol-collision list."""
        if len(candidates) == 1:
            return candidates[0]
        # Rule 1: explicit hint wins (with cross-listed fallback).
        if species:
            row = GeneIndex._pick_for_species(candidates, species)
            if row is not None:
                return row
            # Hint didn't match anything — fall through to heuristic.
        # Rule 2: casing heuristic (human / mouse / None).
        heuristic = _infer_species_from_casing(token)
        if heuristic:
            row = GeneIndex._pick_for_species(candidates, heuristic)
            if row is not None:
                return row
        # Rule 3: load order — first known-species row wins.
        for sp in _KNOWN_SPECIES:
            row = GeneIndex._pick_for_species(candidates, sp)
            if row is not None:
                return row
        # Last resort: first candidate (defensive — shouldn't happen
        # with current build, all rows carry one of _KNOWN_SPECIES).
        return candidates[0]

    def __len__(self) -> int:
        return len(self._rows)


GENE_INDEX: Final[GeneIndex] = GeneIndex(_load_gene_index_rows())


# --- Anatomy index (v0.2.0) -------------------------------------------------

_ANATOMY_INDEX_COLUMNS: Final[tuple[str, ...]] = (
    "term", "uberon", "label", "description", "source",
)


def _load_anatomy_index_rows() -> list[dict[str, str]]:
    tsv = files("indicium_adapters_transcriptomics.data").joinpath(
        "anatomy_index.tsv"
    )
    with tsv.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if tuple(reader.fieldnames or ()) != _ANATOMY_INDEX_COLUMNS:
            raise RuntimeError(
                "anatomy_index.tsv header does not match the declared schema. "
                f"Expected {_ANATOMY_INDEX_COLUMNS!r}; got "
                f"{tuple(reader.fieldnames or ())!r}."
            )
        return [{k: (v or "").strip() for k, v in row.items()} for row in reader]


class AnatomyIndex:
    """Case-insensitive anatomical-term → UBERON CURIE lookup.

    The bundled set is small on purpose (v0.2.0): hand-curated from a single
    reference paper, every entry verified against the EBI OLS4 UBERON index
    by ``data/build_anatomy_index.py`` at build time. Cnidarian-specific
    concepts (physa, oral / aboral pole, mesoglea, gastrodermis, …) are not
    represented in UBERON as of 2026 and so are deliberately absent.
    """

    def __init__(self, rows: list[dict[str, str]]):
        self._rows = rows
        self._by_term: dict[str, dict[str, str]] = {
            row["term"].lower(): row for row in rows if row["term"]
        }
        self._by_label: dict[str, dict[str, str]] = {
            row["label"].lower(): row for row in rows if row["label"]
        }

    def lookup(self, name: str) -> dict[str, str] | None:
        """Return the matching row dict or ``None``.

        Tries (in order): the lower-cased ``term`` column, then the UBERON
        ``label`` column. Both are case-insensitive.
        """
        if not isinstance(name, str):
            return None
        token = name.strip().lower()
        if not token:
            return None
        return self._by_term.get(token) or self._by_label.get(token)

    def __len__(self) -> int:
        return len(self._rows)


ANATOMY_INDEX: Final[AnatomyIndex] = AnatomyIndex(_load_anatomy_index_rows())


# --- Concept index (v0.3.0) -------------------------------------------------

_CONCEPT_INDEX_COLUMNS: Final[tuple[str, ...]] = (
    "term", "curie", "label", "description", "ontology", "source",
)


def _load_one_concept_tsv(filename: str) -> list[dict[str, str]]:
    """Strict-schema TSV loader shared by both concept-index sources."""
    tsv = files("indicium_adapters_transcriptomics.data").joinpath(filename)
    with tsv.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if tuple(reader.fieldnames or ()) != _CONCEPT_INDEX_COLUMNS:
            raise RuntimeError(
                f"{filename} header does not match the declared schema. "
                f"Expected {_CONCEPT_INDEX_COLUMNS!r}; got "
                f"{tuple(reader.fieldnames or ())!r}."
            )
        return [{k: (v or "").strip() for k, v in row.items()} for row in reader]


def _load_concept_index_rows() -> list[dict[str, str]]:
    """Load all bundled concept-index TSVs in priority order.

    Eleven sources merged into the same ``ConceptIndex`` singleton:

    1. ``concept_index.tsv`` — TranscriptoKG modelling vocab (v0.3.0).
    2. ``corpus_concept_index.tsv`` — OLS4 corpus harvest (v0.3.1).
    3. ``biotools_concept_index.tsv`` — bio.tools registry (v0.3.1).
    4. ``cl_concept_index.tsv`` — Cell Ontology (v0.5.0).
    5. ``chebi_concept_index.tsv`` — ChEBI chemicals / drugs (v0.8.0).
    6. ``hp_concept_index.tsv`` — Human Phenotype Ontology (v0.8.0).
    7. ``mp_concept_index.tsv`` — Mammalian Phenotype Ontology (v0.8.0).
    8. ``doid_concept_index.tsv`` — Human Disease Ontology (v0.8.0).
    9. ``pr_concept_index.tsv`` — Protein Ontology (v0.8.0).
    10. ``oba_concept_index.tsv`` — Ontology of Biological Attributes
        (v0.8.0).
    11. ``envo_concept_index.tsv`` — Environment Ontology (v0.8.0).

    Load order encodes a priority: a hypothetical collision (same term
    in multiple files) resolves to the **first** occurrence — paper-
    modelling-grounded (1) > OLS4 method-mention (2) > bio.tools-
    registry (3) > the specialised single-ontology indexes (4-11).
    The build scripts deduplicate against earlier files in practice,
    so collisions are unexpected.
    """
    sources = (
        "concept_index.tsv",
        "corpus_concept_index.tsv",
        "biotools_concept_index.tsv",
        "cl_concept_index.tsv",
        # v0.8.0 ontology-specific indexes (top-7 from the OLS4 audit).
        "chebi_concept_index.tsv",
        "hp_concept_index.tsv",
        "mp_concept_index.tsv",
        "doid_concept_index.tsv",
        "pr_concept_index.tsv",
        "oba_concept_index.tsv",
        "envo_concept_index.tsv",
    )
    rows: list[dict[str, str]] = []
    for filename in sources:
        rows.extend(_load_one_concept_tsv(filename))
    return rows


class ConceptIndex:
    """Case-insensitive concept-term → ontology CURIE lookup.

    Covers analysis concepts (NCIT), regulation processes (INO), assay /
    protocol / data set classes (OBI), statistical concepts (STATO), GO
    roots, plus information-artifact and relation types (IAO, RO, SIO)
    drawn from the TranscriptoKG modelling vocabulary.

    Each row carries the source ontology in the ``ontology`` column so
    consumers can filter the resolved CURIE by axis (e.g. only INO/NCIT
    when the qualifier is direction-of-regulation; only OBI when the
    qualifier is study-protocol metadata).
    """

    def __init__(self, rows: list[dict[str, str]]):
        self._rows = rows
        self._by_term: dict[str, dict[str, str]] = {
            row["term"].lower(): row for row in rows if row["term"]
        }
        self._by_label: dict[str, dict[str, str]] = {
            row["label"].lower(): row for row in rows if row["label"]
        }

    def lookup(self, name: str) -> dict[str, str] | None:
        """Return the matching row dict or ``None``.

        Tries (in order): the lower-cased ``term`` column, then the canonical
        ``label`` column. Both are case-insensitive.
        """
        if not isinstance(name, str):
            return None
        token = name.strip().lower()
        if not token:
            return None
        return self._by_term.get(token) or self._by_label.get(token)

    def __len__(self) -> int:
        return len(self._rows)


CONCEPT_INDEX: Final[ConceptIndex] = ConceptIndex(_load_concept_index_rows())


# --- Resolver helpers (v0.4.0) ----------------------------------------------
#
# Each resolver takes a free-text slot value (a gene name, an anatomy term,
# a concept literal) and returns ``(principal_curie, xrefs_dict) | None``.
# The adapter ``enrich_claim`` methods dispatch through these — one shared
# implementation per index, so bulk and single-cell stay in sync on the
# row → claim-shape rendering logic.

# Order in which gene-row columns are tried for the principal CURIE.
# HGNC and MGI values are stored already in CURIE form; Ensembl / NCBI
# Gene are bare ids and get prefixed here. The order matters for
# multi-id rows: a Nematostella row paired with a human ortholog has
# both ``hgnc`` and ``nemve1`` set — HGNC wins, which is the
# documented v0.1.2 contract. Mouse rows have ``mgi`` set and ``hgnc``
# empty so MGI naturally wins for them.
_GENE_PRINCIPAL_CURIE_PREFIX: Final[tuple[tuple[str, str | None], ...]] = (
    ("hgnc", None),
    ("mgi", None),
    ("flybase", "FlyBase"),
    ("wormbase", "WormBase"),
    ("zfin", "ZFIN"),
    ("ensembl", "ENSEMBL"),
    ("ncbigene", "NCBIGene"),
)

# Xref columns surfaced under ``xrefs[slot]``. Bare ids (no prefixing).
_GENE_XREF_COLUMNS: Final[tuple[str, ...]] = (
    "hgnc", "mgi", "flybase", "wormbase", "zfin",
    "ensembl", "ncbigene", "nemve1", "nvertx", "uniprot",
)


def resolve_gene(
    name: str, *, species: str | None = None
) -> tuple[str, dict[str, str]] | None:
    """Resolve ``name`` to a gene CURIE + xrefs via ``GENE_INDEX``.

    ``species`` is forwarded to ``GENE_INDEX.lookup`` and disambiguates
    symbol-level collisions (``MAPK1`` ↔ ``Mapk1``). Principal CURIE
    follows HGNC > MGI > Ensembl > NCBIGene priority; xrefs surface
    the bare HGNC / MGI / Ensembl / NCBIGene / Nemve1 / NvERTx /
    UniProt ids so downstream consumers can render any of them."""
    row = GENE_INDEX.lookup(name, species=species)
    if row is None:
        return None
    curie: str | None = None
    for col, prefix in _GENE_PRINCIPAL_CURIE_PREFIX:
        val = row.get(col)
        if val:
            curie = val if prefix is None else f"{prefix}:{val}"
            break
    if curie is None:
        return None
    xrefs = {col: row[col] for col in _GENE_XREF_COLUMNS if row.get(col)}
    return curie, xrefs


def resolve_anatomy(
    name: str, *, species: str | None = None
) -> tuple[str, dict[str, str]] | None:
    """Resolve ``name`` to a UBERON anatomy CURIE + xrefs via ``ANATOMY_INDEX``.

    ``species`` is accepted for resolver-signature uniformity (the
    bulk/single-cell dispatchers pass it to every resolver) but
    anatomy entries are not species-specific in the current bundle —
    UBERON is a multi-species ontology — so the hint is ignored."""
    del species  # uniform signature; not used here
    row = ANATOMY_INDEX.lookup(name)
    if row is None:
        return None
    xref = {"uberon": row["uberon"], "label": row["label"]}
    if row.get("description"):
        xref["description"] = row["description"]
    return row["uberon"], xref


def resolve_concept(
    name: str, *, species: str | None = None
) -> tuple[str, dict[str, str]] | None:
    """Resolve ``name`` to a concept CURIE + xrefs via ``CONCEPT_INDEX``.

    ``species`` is accepted for resolver-signature uniformity; concept
    entries (methodology / protocols / cell types) are not
    species-specific so the hint is ignored."""
    del species  # uniform signature; not used here
    row = CONCEPT_INDEX.lookup(name)
    if row is None:
        return None
    xref = {
        "curie": row["curie"],
        "label": row["label"],
        "ontology": row["ontology"],
    }
    if row.get("description"):
        xref["description"] = row["description"]
    return row["curie"], xref
