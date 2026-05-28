"""Shared enrichment back-end for transcriptomics adapters.

This module owns the three resolver indexes (genes, anatomy, concepts) and
exposes the three functions the adapters call from their ``enrich_claim``
methods:

* ``resolve_gene(name, species=None)``     uses ``GENE_INDEX``
* ``resolve_anatomy(name, species=None)``  uses ``ANATOMY_INDEX``
* ``resolve_concept(name, species=None)``  uses ``CONCEPT_INDEX``

All three return either ``(principal_curie, xrefs_dict)`` or ``None``. They
share the same signature so the adapters can dispatch through a dict mapping
resolver name to function.

The indexes are loaded once at import time from bundled TSV resources under
``data/``. Maintainer scripts under ``data/build_*.py`` regenerate those TSVs
from public sources (HGNC, MGI, OLS4, bio.tools, etc.).
"""
from __future__ import annotations

import csv
import re
from importlib.resources import files
from typing import Final

# Columns of ``data/gene_index.tsv``. Declared here so the loader can reject
# a TSV whose header drifts. Add a column means updating this tuple AND the
# lookup logic below.
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


# Species fallback order used by ``GeneIndex.lookup`` when neither an explicit
# hint nor the casing heuristic can disambiguate a symbol collision. Human
# first because it dominates transcriptomics-paper claims; mouse second; then
# the v0.6.2 model organisms in rough corpus-prevalence order.
_KNOWN_SPECIES: Final[tuple[str, ...]] = (
    "human", "mouse", "zebrafish", "drosophila", "celegans", "nematostella",
)


def _load_gene_index_rows() -> list[dict[str, str]]:
    """Read the bundled gene_index.tsv and return its rows as dicts.

    Raises ``RuntimeError`` if the header does not match the declared schema.
    """
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


# C. elegans symbol pattern: three (or four) lowercase letters, hyphen, then
# one or more alphanumeric characters. Examples: ``unc-22``, ``daf-16``,
# ``let-7``, ``gld-1``. No other organism uses this layout, which makes the
# pattern a safe discriminator.
_C_ELEGANS_PATTERN: Final[re.Pattern[str]] = re.compile(
    r"^[a-z]{3,4}-[a-zA-Z0-9]+$"
)


def _infer_species_from_casing(token: str) -> str | None:
    """Casing heuristic used when no explicit species hint is given.

    Field conventions:

    * Human (HGNC): ALL UPPERCASE, e.g. ``MAPK1``, ``WNT4``, ``CD3E``.
    * Mouse (MGI): Title-cased, e.g. ``Mapk1``, ``Wnt4``, ``Cd3e``.
    * C. elegans (WormBase): three or four lowercase letters, hyphen, then
      alphanumeric, e.g. ``unc-22``, ``daf-16``, ``let-7``.
    * Zebrafish (ZFIN) and Drosophila (FlyBase) both use lowercase or
      mixed-case conventions that look like a lowercase typo of a human or
      mouse symbol. The heuristic does not try to guess them. Callers
      targeting those organisms must pass ``species=`` explicitly. Without
      a hint, a lowercase ``mmp9`` falls through to the load-order fallback
      and resolves to human via the cross-listed rule.
    * Nematostella (Nemve1) ids carry dots (``gw.51.96.1``) and never come
      through this path; they hit the by-nemve1 lookup branch instead.

    Returns ``None`` when the casing matches no high-confidence convention.
    """
    if not token or not token[0].isalpha():
        return None
    # C. elegans first: most distinctive thanks to the mandatory hyphen.
    if _C_ELEGANS_PATTERN.match(token):
        return "celegans"
    if token.isupper():
        return "human"
    if token[0].isupper() and not token.isupper():
        return "mouse"
    # All-lowercase is ambiguous. Defer to the species hint or load-order
    # fallback. See the class docstring.
    return None


class GeneIndex:
    """Case-insensitive gene-symbol lookup with species disambiguation.

    Internal storage keeps one *list* of rows per lower-cased symbol. Most
    symbols map to a single row, but a few collide across organisms (``MAPK1``
    and ``Mapk1`` share the lower-cased key ``mapk1`` and produce a two-row
    list, one human and one mouse). ``lookup`` then picks one row using three
    rules in order:

    1. Explicit ``species`` argument (the hint path).
    2. Casing heuristic (UPPER then human, Title-case then mouse).
    3. Load order: human first, see ``_KNOWN_SPECIES``.

    Accession-style lookups (Nemve1, UniProt, MGI id, FlyBase, WormBase,
    ZFIN) are stored separately and always single-row, since accessions are
    unambiguous by construction.
    """

    def __init__(self, rows: list[dict[str, str]]):
        self._rows = rows
        # Symbol index: one list per lower-cased symbol so multi-organism
        # collisions can coexist and be disambiguated at lookup time.
        self._by_symbol: dict[str, list[dict[str, str]]] = {}
        for row in rows:
            symbol = row.get("symbol")
            if symbol:
                self._by_symbol.setdefault(symbol.lower(), []).append(row)
        # Accession indexes: single-row, exact-match. These are unambiguous
        # because each accession namespace is organism-scoped.
        self._by_nemve1: dict[str, dict[str, str]] = {
            row["nemve1"]: row for row in rows if row.get("nemve1")
        }
        self._by_uniprot: dict[str, dict[str, str]] = {
            row["uniprot"]: row for row in rows if row.get("uniprot")
        }
        # MGI accessions are stored as ``MGI:NNN`` already.
        self._by_mgi: dict[str, dict[str, str]] = {
            row["mgi"]: row for row in rows if row.get("mgi")
        }
        # Model-organism accession indexes (added in v0.6.2). Each
        # accession is registered under both the bare id (``FBgn0003079``)
        # and the prefixed CURIE form (``FlyBase:FBgn0003079``) so callers
        # can pass either spelling.
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
        """Return the matching row, or ``None``.

        ``species`` is an optional hint such as ``"human"`` or ``"mouse"``.
        When given, it disambiguates symbol collisions. When not given, the
        casing heuristic infers a species; if still ambiguous, load order
        decides (human wins).

        Lookup order: gene symbol (case-insensitive, with species
        disambiguation), then Nemve1, UniProt, MGI, FlyBase, WormBase, ZFIN
        accessions. The accession paths are single-row and bypass the
        disambiguation logic.
        """
        if not isinstance(name, str):
            return None
        token = name.strip()
        if not token:
            return None

        # Symbol path: may return several rows, so call the picker.
        candidates = self._by_symbol.get(token.lower(), [])
        if candidates:
            return self._pick_candidate(candidates, token, species)

        # Accession paths: single row, exact match, case-sensitive.
        if row := self._by_nemve1.get(token):
            return row
        if row := self._by_uniprot.get(token):
            return row
        if row := self._by_mgi.get(token):
            return row
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

        A row "serves" species ``X`` when either its ``species`` tag equals
        ``X``, or it carries the canonical CURIE column for ``X`` (HGNC for
        human, MGI for mouse, Nemve1 for Nematostella, FlyBase for
        Drosophila, WormBase for C. elegans, ZFIN for zebrafish).

        The second clause matters for cross-listed rows: a Nematostella
        row that also ships a human ortholog under ``hgnc`` should answer
        a human lookup of ``MMP9``, even though its primary ``species``
        tag is ``nematostella``.
        """
        canonical_col = {
            "human":        "hgnc",
            "mouse":        "mgi",
            "drosophila":   "flybase",
            "celegans":     "wormbase",
            "zebrafish":    "zfin",
            "nematostella": "nemve1",
        }.get(wanted)
        # First pass: rows whose primary species tag is the wanted one.
        for row in candidates:
            if row.get("species") == wanted:
                return row
        # Second pass: cross-listed rows that carry the wanted species'
        # canonical CURIE column.
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
        """Apply the three disambiguation rules to a symbol-collision list."""
        if len(candidates) == 1:
            return candidates[0]
        # Rule 1: explicit hint wins, with cross-listed fallback.
        if species:
            row = GeneIndex._pick_for_species(candidates, species)
            if row is not None:
                return row
            # Hint didn't match anything in the candidate set, fall through
            # to the heuristic.
        # Rule 2: casing heuristic.
        heuristic = _infer_species_from_casing(token)
        if heuristic:
            row = GeneIndex._pick_for_species(candidates, heuristic)
            if row is not None:
                return row
        # Rule 3: load order. The first known-species row in ``_KNOWN_SPECIES``
        # wins, which means human is the default tie-breaker.
        for sp in _KNOWN_SPECIES:
            row = GeneIndex._pick_for_species(candidates, sp)
            if row is not None:
                return row
        # Defensive fallback. Should not be reached under the current build
        # because every row carries one of ``_KNOWN_SPECIES``.
        return candidates[0]

    def __len__(self) -> int:
        return len(self._rows)


# Module-level singleton. Loaded once at import time.
GENE_INDEX: Final[GeneIndex] = GeneIndex(_load_gene_index_rows())


# --- Anatomy index ---------------------------------------------------------

_ANATOMY_INDEX_COLUMNS: Final[tuple[str, ...]] = (
    "term", "uberon", "label", "description", "source",
)


def _load_anatomy_index_rows() -> list[dict[str, str]]:
    """Read the bundled anatomy_index.tsv and return its rows as dicts."""
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
    """Case-insensitive anatomical-term to UBERON CURIE lookup.

    The bundled set is intentionally small: hand-curated from a reference
    paper and verified against the EBI OLS4 UBERON index by
    ``data/build_anatomy_index.py`` at build time. Cnidarian-specific
    structures (physa, oral / aboral pole, mesoglea, gastrodermis) are not
    represented in UBERON as of 2026, so they are not in the table.
    """

    def __init__(self, rows: list[dict[str, str]]):
        self._rows = rows
        # Two lower-cased indexes: the user-facing ``term`` column and the
        # canonical UBERON ``label`` column. Lookup tries both.
        self._by_term: dict[str, dict[str, str]] = {
            row["term"].lower(): row for row in rows if row["term"]
        }
        self._by_label: dict[str, dict[str, str]] = {
            row["label"].lower(): row for row in rows if row["label"]
        }

    def lookup(self, name: str) -> dict[str, str] | None:
        """Return the matching row, or ``None``.

        Tries the ``term`` column first, then falls back to the canonical
        UBERON ``label``. Both are case-insensitive.
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


# --- Concept index ---------------------------------------------------------

_CONCEPT_INDEX_COLUMNS: Final[tuple[str, ...]] = (
    "term", "curie", "label", "description", "ontology", "source",
)


def _load_one_concept_tsv(filename: str) -> list[dict[str, str]]:
    """Strict-schema TSV loader shared by every concept-index source.

    Tolerates a *missing* file (returns ``[]``): the index may not have been
    built yet, or may have been left out on purpose. Stays strict on a
    *malformed* header (raises): silent corruption would be worse than a
    loud failure. This lets the package import cleanly during maintainer
    build steps that regenerate a subset of the TSVs sequentially.
    """
    tsv = files("indicium_adapters_transcriptomics.data").joinpath(filename)
    try:
        fh = tsv.open("r", encoding="utf-8", newline="")
    except FileNotFoundError:
        return []
    with fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if tuple(reader.fieldnames or ()) != _CONCEPT_INDEX_COLUMNS:
            raise RuntimeError(
                f"{filename} header does not match the declared schema. "
                f"Expected {_CONCEPT_INDEX_COLUMNS!r}; got "
                f"{tuple(reader.fieldnames or ())!r}."
            )
        return [{k: (v or "").strip() for k, v in row.items()} for row in reader]


def _load_concept_index_rows() -> list[dict[str, str]]:
    """Load every bundled concept-index TSV in priority order.

    Nineteen sources merged into the same ``ConceptIndex`` singleton. The
    load order encodes a priority: when the same term appears in several
    files, the first occurrence wins. Paper-modelling vocabulary ranks
    above OLS4 method mentions, which rank above bio.tools, which rank
    above the specialised single-ontology indexes. In practice the build
    scripts deduplicate against earlier files, so collisions are rare.

    Sources, in load order:

    1.  ``concept_index.tsv``           TranscriptoKG modelling vocab.
    2.  ``corpus_concept_index.tsv``    OLS4 corpus harvest.
    3.  ``biotools_concept_index.tsv``  bio.tools registry.
    4.  ``cl_concept_index.tsv``        Cell Ontology.
    5.  ``chebi_concept_index.tsv``     ChEBI chemicals / drugs.
    6.  ``hp_concept_index.tsv``        Human Phenotype Ontology.
    7.  ``mp_concept_index.tsv``        Mammalian Phenotype Ontology.
    8.  ``doid_concept_index.tsv``      Human Disease Ontology.
    9.  ``pr_concept_index.tsv``        Protein Ontology.
    10. ``oba_concept_index.tsv``       Ontology of Biological Attributes.
    11. ``envo_concept_index.tsv``      Environment Ontology.
    12. ``ms_concept_index.tsv``        Mass-spectrometry methods (Tier-1).
    13. ``uo_concept_index.tsv``        Units of measurement (Tier-1).
    14. ``mod_concept_index.tsv``       Protein modifications (Tier-1).
    15. ``fbbt_concept_index.tsv``      Drosophila anatomy (Tier-2).
    16. ``wbbt_concept_index.tsv``      C. elegans anatomy (Tier-2).
    17. ``zfa_concept_index.tsv``       Zebrafish anatomy (Tier-2).
    18. ``emapa_concept_index.tsv``     Mouse anatomy (Tier-2).
    19. ``xpo_concept_index.tsv``       Xenopus phenotype (Tier-2).

    Four candidate ontologies (EC, RXNO, HCAO, ZP) were evaluated and
    deferred to a later release. See ``data/build_extra_ontology_indexes.py``
    for the per-ontology rationale.
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
        # v0.9.0 Tier-1 universal expansion.
        "ms_concept_index.tsv",
        "uo_concept_index.tsv",
        "mod_concept_index.tsv",
        # v0.9.0 Tier-2 organism-paired anatomy and phenotype.
        "fbbt_concept_index.tsv",
        "wbbt_concept_index.tsv",
        "zfa_concept_index.tsv",
        "emapa_concept_index.tsv",
        "xpo_concept_index.tsv",
    )
    rows: list[dict[str, str]] = []
    for filename in sources:
        rows.extend(_load_one_concept_tsv(filename))
    return rows


class ConceptIndex:
    """Case-insensitive concept-term to ontology CURIE lookup.

    Covers analysis concepts (NCIT), regulation processes (INO), assay /
    protocol / data set classes (OBI), statistical concepts (STATO), GO
    roots, plus information-artifact and relation types (IAO, RO, SIO)
    from the TranscriptoKG modelling vocabulary. Later releases added cell
    types (CL), chemicals (ChEBI), phenotypes (HP, MP, OBA), diseases
    (DOID), proteins (PR), environments (ENVO), units (UO), mass-spec
    methods (MS), protein modifications (MOD), and organism-paired
    anatomy / phenotype ontologies (FBbt, WBbt, ZFA, EMAPA, XPO).

    Every row carries the source ontology in the ``ontology`` column so
    consumers can filter by axis (for example only INO/NCIT when the
    qualifier is direction-of-regulation, or only OBI when the qualifier
    is study-protocol metadata).
    """

    def __init__(self, rows: list[dict[str, str]]):
        self._rows = rows
        # Same two-index structure as ``AnatomyIndex``: user-facing
        # ``term`` column plus canonical ``label`` column, both lower-cased.
        self._by_term: dict[str, dict[str, str]] = {
            row["term"].lower(): row for row in rows if row["term"]
        }
        self._by_label: dict[str, dict[str, str]] = {
            row["label"].lower(): row for row in rows if row["label"]
        }

    def lookup(self, name: str) -> dict[str, str] | None:
        """Return the matching row, or ``None``.

        Tries the ``term`` column first, then falls back to the canonical
        ``label``. Both are case-insensitive.
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


# --- Resolver helpers ------------------------------------------------------
#
# Each resolver takes a free-text slot value (a gene name, an anatomy term,
# a concept literal) and returns ``(principal_curie, xrefs_dict)`` on a hit
# or ``None`` on a miss. The adapter ``enrich_claim`` methods dispatch
# through these, so one shared implementation per index keeps bulk and
# single-cell in sync on the row-to-claim rendering logic.

# Order in which gene-row columns are tried for the principal CURIE.
# HGNC and MGI values are stored already in CURIE form; Ensembl and NCBI
# Gene are bare ids and get prefixed here. Order matters for multi-id
# rows: a Nematostella row paired with a human ortholog has both ``hgnc``
# and ``nemve1`` set, and HGNC wins. Mouse rows have ``mgi`` set with
# ``hgnc`` empty, so MGI naturally wins for them.
_GENE_PRINCIPAL_CURIE_PREFIX: Final[tuple[tuple[str, str | None], ...]] = (
    ("hgnc", None),
    ("mgi", None),
    ("flybase", "FlyBase"),
    ("wormbase", "WormBase"),
    ("zfin", "ZFIN"),
    ("ensembl", "ENSEMBL"),
    ("ncbigene", "NCBIGene"),
)

# Xref columns surfaced under ``xrefs[slot]``. Bare ids, no prefixing.
_GENE_XREF_COLUMNS: Final[tuple[str, ...]] = (
    "hgnc", "mgi", "flybase", "wormbase", "zfin",
    "ensembl", "ncbigene", "nemve1", "nvertx", "uniprot",
)


def resolve_gene(
    name: str, *, species: str | None = None
) -> tuple[str, dict[str, str]] | None:
    """Resolve ``name`` to a gene CURIE and xrefs via ``GENE_INDEX``.

    ``species`` is forwarded to ``GENE_INDEX.lookup`` and disambiguates
    symbol-level collisions such as ``MAPK1`` vs ``Mapk1``. The principal
    CURIE follows the HGNC then MGI then FlyBase then WormBase then ZFIN
    then Ensembl then NCBIGene priority defined in
    ``_GENE_PRINCIPAL_CURIE_PREFIX``. The ``xrefs`` dict carries every
    cross-reference column the row provides (bare ids).
    """
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
    """Resolve ``name`` to a UBERON anatomy CURIE and xrefs.

    ``species`` is accepted for signature uniformity (the bulk and
    single-cell dispatchers pass it to every resolver) but UBERON is a
    multi-species ontology, so the hint is ignored here.
    """
    del species
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
    """Resolve ``name`` to a concept CURIE and xrefs via ``CONCEPT_INDEX``.

    ``species`` is accepted for signature uniformity; concept entries
    (methodology, protocols, cell types, chemicals, units, and so on) are
    not species-specific, so the hint is ignored.
    """
    del species
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
