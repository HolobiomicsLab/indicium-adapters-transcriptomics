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

This module also exposes ``build_sssom_mapping``, the helper both adapters
use to lift each resolver hit to a first-class SSSOMMapping fragment
(indicium >= 1.11.0). The fragment is a plain dict matching SSSOMMapping's
constructor fields so downstream Perspicacite can fold it directly into
``IndiciumDocument.sssom_mappings`` without round-tripping through pydantic.
"""
from __future__ import annotations

import csv
import re
from importlib.resources import files
from typing import Final

# Per-resolver confidence floor for emitted SSSOMMapping nodes. Tighter
# resolvers (regex-anchored accession; exact gene-symbol lookup) get higher
# scores; heuristic fallbacks (embedded-pathway substring scan; concept
# trim-tail) get lower scores so downstream consumers can filter on
# ``confidence >= X`` to drop the soft matches if they prefer.
#
# These are FLOORS, not point estimates. The values are intentionally
# conservative — a future refinement can refine per-row (e.g. drop confidence
# when the gene-symbol lookup needed a species-hint disambiguation), but for
# v0.14 the kind-level granularity matches what the resolvers already report
# via ``ontology_term_justifications``.
_SSSOM_CONFIDENCE_BY_KIND: Final[dict[str, float]] = {
    "accession":         0.99,  # regex-anchored, near-perfect when matched
    "gene":              0.95,
    "anatomy":           0.95,
    "concept":           0.90,
    "embedded_pathway":  0.85,  # substring fallback, less anchored
    "trimmed_concept":   0.80,  # heuristic trim-tail fallback
}


# SSSOM predicate picked per resolver kind. ``skos:exactMatch`` for the
# strict resolutions where the raw text is interchangeable with the CURIE
# (HGNC symbol, accession regex hit, canonical lexical lookup);
# ``skos:closeMatch`` for the LLM-tolerant fallbacks where we modified
# the input (trimmed a trailing noun) or extracted the id from a longer
# wrapping string. The skos terminology gain matters at the SPARQL
# boundary: a consumer querying for "interchangeable identities" can
# include only ``exactMatch`` rows; a consumer wanting broader recall
# accepts both.
_SSSOM_PREDICATE_BY_KIND: Final[dict[str, str]] = {
    "accession":         "skos:exactMatch",
    "gene":              "skos:exactMatch",
    "anatomy":           "skos:exactMatch",
    "concept":           "skos:exactMatch",
    "embedded_pathway":  "skos:closeMatch",
    "trimmed_concept":   "skos:closeMatch",
}


def _slug_for_sssom_subject(text: str) -> str:
    """Sanitize a free-text mention into a URL-safe SSSOM subject suffix.

    Lowercase, replace runs of non-alphanumeric with a single underscore,
    strip leading / trailing underscores. Stable across runs so the same
    raw mention produces the same SSSOMMapping id, enabling dedup at the
    document level.

    >>> _slug_for_sssom_subject("WNT4")
    'wnt4'
    >>> _slug_for_sssom_subject("Illumina NovaSeq 6000")
    'illumina_novaseq_6000'
    >>> _slug_for_sssom_subject("R-HSA-9012999")
    'r_hsa_9012999'
    """
    slug = re.sub(r"[^a-z0-9]+", "_", text.casefold()).strip("_")
    return slug or "unknown"


def build_sssom_mapping(
    *,
    raw_text: str,
    object_curie: str,
    object_label: str | None,
    justification: str,
    resolver_kind: str,
) -> dict:
    """Build an SSSOMMapping fragment for one resolver hit.

    Returns a dict shaped like the indicium ``SSSOMMapping`` pydantic
    constructor's keyword arguments, with no fields populated beyond what
    a single lexical / composite match can know. The caller (an adapter
    ``enrich_claim``) accumulates these fragments in
    ``claim["sssom_mappings"]`` and downstream Perspicacite stages fold
    them into ``IndiciumDocument.sssom_mappings``.

    * ``raw_text``      - the free-text mention from the LLM extractor.
    * ``object_curie``  - the resolved CURIE.
    * ``object_label``  - canonical label of the resolved term, if known
      (None for accession resolutions, present for concept / anatomy /
      most gene rows).
    * ``justification`` - ``semapv:LexicalMatching`` or
      ``semapv:CompositeMatching`` (passed through unchanged from the
      resolver's justification report).
    * ``resolver_kind`` - one of ``accession`` / ``gene`` / ``anatomy`` /
      ``concept`` / ``embedded_pathway`` / ``trimmed_concept``; drives
      the default confidence floor.

    The fragment uses the synthetic ``asb.lex:`` namespace for
    ``subject_id`` so that free-text mentions stay deduplicable across a
    document without needing a pre-existing identifier registry. The
    ``id`` field follows the ``{subject_id}--{object_id}`` convention
    shipped in the indicium v1.11.0 test suite.

    The ``predicate_id`` is ``skos:exactMatch`` for strict resolutions
    where the raw text is interchangeable with the CURIE (HGNC symbol
    lookup, accession-regex hit, canonical concept lookup) and
    ``skos:closeMatch`` for fallbacks where the input was modified
    (trim-tail) or the id was extracted from a longer wrapping string
    (embedded-pathway scan).
    """
    slug = _slug_for_sssom_subject(raw_text)
    subject_id = f"asb.lex:{slug}"
    fragment: dict = {
        "id": f"{subject_id}--{object_curie}",
        "subject_id": subject_id,
        "subject_label": raw_text,
        "predicate_id": _SSSOM_PREDICATE_BY_KIND.get(
            resolver_kind, "skos:closeMatch"
        ),
        "object_id": object_curie,
        "mapping_justification": justification,
        "confidence": _SSSOM_CONFIDENCE_BY_KIND.get(resolver_kind, 0.5),
    }
    if object_label:
        fragment["object_label"] = object_label
    return fragment

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


# --- Shared concept-index TSV reader --------------------------------------
#
# The concept index proper is loaded further down in this module, but both
# the anatomy loader (organism-paired ontologies) and the concept loader
# read TSVs in the same format. The reader is hoisted up here so the
# anatomy loader can call it without a forward reference.

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


# --- Anatomy index ---------------------------------------------------------

# The on-disk anatomy_index.tsv (UBERON, species-agnostic) keeps its
# original schema with the ``uberon`` column. Organism-paired ontologies
# (FBbt, WBbt, ZFA, EMAPA) live in concept-index-shaped TSVs already
# bundled for the concept index. The loader below normalises both shapes
# into a single unified row dict carrying ``term``, ``curie``, ``label``,
# ``description``, ``ontology``, ``species``, ``source``.
_UBERON_ANATOMY_COLUMNS: Final[tuple[str, ...]] = (
    "term", "uberon", "label", "description", "source",
)

# Organism-paired anatomy ontologies that already ship as concept-index
# TSVs and can be re-used as anatomy lookups when a species hint is
# provided. Each tuple is ``(filename, species_tag, ontology_label)``.
# XPO (Xenopus phenotype) is left in the concept index only because it
# describes phenotypes, not anatomy proper.
_ORGANISM_ANATOMY_SOURCES: Final[tuple[tuple[str, str, str], ...]] = (
    ("fbbt_concept_index.tsv",  "drosophila", "FBbt"),
    ("wbbt_concept_index.tsv",  "celegans",   "WBbt"),
    ("zfa_concept_index.tsv",   "zebrafish",  "ZFA"),
    ("emapa_concept_index.tsv", "mouse",      "EMAPA"),
)


def _load_uberon_anatomy_rows() -> list[dict[str, str]]:
    """Read the bundled anatomy_index.tsv (UBERON, species-agnostic) and
    normalise each row into the unified anatomy schema."""
    tsv = files("indicium_adapters_transcriptomics.data").joinpath(
        "anatomy_index.tsv"
    )
    with tsv.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if tuple(reader.fieldnames or ()) != _UBERON_ANATOMY_COLUMNS:
            raise RuntimeError(
                "anatomy_index.tsv header does not match the declared schema. "
                f"Expected {_UBERON_ANATOMY_COLUMNS!r}; got "
                f"{tuple(reader.fieldnames or ())!r}."
            )
        out: list[dict[str, str]] = []
        for raw in reader:
            row = {k: (v or "").strip() for k, v in raw.items()}
            out.append({
                "term":        row["term"],
                "curie":       row["uberon"],
                "label":       row["label"],
                "description": row["description"],
                "ontology":    "UBERON",
                "species":     "",  # species-agnostic
                "source":      row["source"],
            })
        return out


def _load_organism_anatomy_rows(
    filename: str, species: str, ontology: str,
) -> list[dict[str, str]]:
    """Read an organism-paired concept-index TSV (FBbt/WBbt/ZFA/EMAPA)
    and normalise each row into the unified anatomy schema.

    Tolerant of a missing file (returns ``[]``): an adapter build that
    has not yet regenerated the organism-paired ontology index should
    still import cleanly.
    """
    rows = _load_one_concept_tsv(filename)
    out: list[dict[str, str]] = []
    for row in rows:
        # Honour the ontology column when present, fall back to the
        # caller-supplied default. Both should agree in practice.
        ont = row.get("ontology") or ontology
        out.append({
            "term":        row["term"],
            "curie":       row["curie"],
            "label":       row["label"],
            "description": row.get("description", ""),
            "ontology":    ont,
            "species":     species,
            "source":      row.get("source", ""),
        })
    return out


def _load_anatomy_index_rows() -> list[dict[str, str]]:
    """Load every bundled anatomy source into one unified row list.

    Sources, in load order:

    1. ``anatomy_index.tsv`` (UBERON, species-agnostic).
    2. ``fbbt_concept_index.tsv`` (Drosophila anatomy).
    3. ``wbbt_concept_index.tsv`` (C. elegans anatomy).
    4. ``zfa_concept_index.tsv`` (zebrafish anatomy).
    5. ``emapa_concept_index.tsv`` (mouse developmental anatomy).

    A given lookup key (e.g. ``"embryo"``) can appear across multiple
    ontologies. The picker in ``AnatomyIndex._pick_candidate`` then
    chooses between them using the species hint and a UBERON-first
    fallback so collisions resolve deterministically.

    The same organism-paired TSVs are also loaded by the concept index;
    that duplication is intentional. The anatomy resolver gets first
    crack at organism-specific spatial terms when a species hint is
    given, and the concept resolver still catches the same terms when
    no hint is supplied.
    """
    rows: list[dict[str, str]] = _load_uberon_anatomy_rows()
    for filename, species, ontology in _ORGANISM_ANATOMY_SOURCES:
        rows.extend(_load_organism_anatomy_rows(filename, species, ontology))
    return rows


class AnatomyIndex:
    """Case-insensitive multi-species anatomical-term to CURIE lookup.

    Each bundled anatomy source contributes rows with a common shape:
    ``{term, curie, label, description, ontology, species, source}``.
    Species values are ``""`` (UBERON, species-agnostic) or the
    organism tag matching one of the bulk-adapter's recognised species
    (``drosophila``, ``celegans``, ``zebrafish``, ``mouse``).

    Lookup is two-stage:

    1. Find every row whose ``term`` or canonical ``label`` matches the
       input string (case-insensitive).
    2. Pick one row using the species hint plus a deterministic fallback:
       species-tag exact match, then species-agnostic UBERON, then load
       order.

    Cnidarian-specific structures (physa, oral / aboral pole, mesoglea,
    gastrodermis) are not represented in UBERON as of 2026 and are
    documented as triaged out of the curated bundle.
    """

    def __init__(self, rows: list[dict[str, str]]):
        self._rows = rows
        # Multi-row indexes so the same lookup key can resolve to several
        # organism-specific ontologies. The picker then chooses based on
        # species hint or load order.
        self._by_term: dict[str, list[dict[str, str]]] = {}
        self._by_label: dict[str, list[dict[str, str]]] = {}
        for row in rows:
            if row.get("term"):
                self._by_term.setdefault(row["term"].lower(), []).append(row)
            if row.get("label"):
                self._by_label.setdefault(row["label"].lower(), []).append(row)

    def lookup(
        self, name: str, *, species: str | None = None,
    ) -> dict[str, str] | None:
        """Return the matching row, or ``None``.

        Tries the ``term`` column first, then falls back to the canonical
        ``label``. Both are case-insensitive.

        ``species`` is an optional hint such as ``"mouse"`` /
        ``"drosophila"`` / ``"celegans"`` / ``"zebrafish"``. When given,
        species-specific entries (EMAPA / FBbt / WBbt / ZFA) are
        preferred over species-agnostic UBERON entries for the same
        term. Without a hint, UBERON wins.
        """
        if not isinstance(name, str):
            return None
        token = name.strip().lower()
        if not token:
            return None
        candidates = self._by_term.get(token) or self._by_label.get(token) or []
        if not candidates:
            return None
        return self._pick_candidate(candidates, species)

    @staticmethod
    def _pick_candidate(
        candidates: list[dict[str, str]], species: str | None,
    ) -> dict[str, str]:
        """Apply the two-rule disambiguation to a list of anatomy candidates."""
        if len(candidates) == 1:
            return candidates[0]
        # Rule 1: explicit species hint wins when an entry tags that species.
        if species:
            for row in candidates:
                if row.get("species") == species:
                    return row
        # Rule 2: species-agnostic UBERON entries win when no hint matches.
        for row in candidates:
            if not row.get("species"):
                return row
        # Rule 3: defensive load-order fallback. Should not be reached when
        # at least one UBERON entry exists for the term.
        return candidates[0]

    def __len__(self) -> int:
        return len(self._rows)


ANATOMY_INDEX: Final[AnatomyIndex] = AnatomyIndex(_load_anatomy_index_rows())


# --- Concept index ---------------------------------------------------------
#
# ``_CONCEPT_INDEX_COLUMNS`` and ``_load_one_concept_tsv`` are defined earlier
# in the file (just before the anatomy section) because the anatomy loader
# reuses the same TSV reader to consume the organism-paired ontology indexes
# (FBbt / WBbt / ZFA / EMAPA). Keeping the helper near the top avoids a
# forward reference.


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
    20. ``go_concept_index.tsv``        GO biological process / molecular
        function concrete terms (v0.12).
    21. ``zp_concept_index.tsv``        Zebrafish Phenotype (v0.12,
        re-tried with canonical-form candidates).

    Two candidate ontologies (RXNO, HCAO) remain deferred because they
    re-publish CHEBI and CL terms respectively, with low effective
    transcripto yield. See ``data/build_extra_ontology_indexes.py``
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
        # v0.12 expansion: GO biological-process / molecular-function
        # concrete terms (the static concept_index.tsv only shipped the
        # three GO roots), and ZP zebrafish phenotype re-tried with
        # canonical-form candidates.
        "go_concept_index.tsv",
        "zp_concept_index.tsv",
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
    """Resolve ``name`` to an anatomy CURIE and xrefs via ``ANATOMY_INDEX``.

    ``species`` is forwarded to ``ANATOMY_INDEX.lookup`` so organism-paired
    ontologies (FBbt for Drosophila, WBbt for C. elegans, ZFA for
    zebrafish, EMAPA for mouse) can be picked when an explicit hint is
    given. Without a hint, species-agnostic UBERON entries win.

    The returned ``xrefs`` dict carries ``curie``, ``label``,
    ``ontology`` (e.g. ``"UBERON"``, ``"FBbt"``, ``"EMAPA"``), and the
    optional ``species`` tag and ``description``. The pre-v0.12 shape
    that used a dedicated ``uberon`` key has been replaced because the
    resolver is no longer UBERON-only.
    """
    row = ANATOMY_INDEX.lookup(name, species=species)
    if row is None:
        return None
    xref = {
        "curie":    row["curie"],
        "label":    row["label"],
        "ontology": row["ontology"],
    }
    if row.get("species"):
        xref["species"] = row["species"]
    if row.get("description"):
        xref["description"] = row["description"]
    return row["curie"], xref


# Trim-rules for ``resolve_concept``'s LLM-tolerant fallback. Each entry
# is ``(suffix_to_match_lowercase, replacement)``. When the original
# lookup misses, we try each rule in order: if the lowercased input ends
# with ``suffix_to_match``, replace that suffix with ``replacement`` and
# re-lookup. The first rule whose re-lookup hits wins; otherwise we fall
# through to the next rule.
#
# Two complementary shapes for plural cell-type mentions:
# ``"cnidocyte cells"`` -> ``"cnidocyte"`` (suffix removed entirely),
# ``"T cells"`` -> ``"T cell"`` (plural collapsed to canonical singular).
# Longer suffixes precede shorter overlaps so ``"signaling pathway"``
# is not mis-trimmed by the bare ``" pathway"`` rule.
#
# Conservative on purpose: each rule requires a leading whitespace so a
# trim never bisects a word (``"mesencephalon"`` is not chopped by the
# ``" cell"`` rule because there is no preceding space). One trim per
# attempt; no recursion.
_CONCEPT_TRIM_RULES: Final[tuple[tuple[str, str], ...]] = (
    (" signaling pathway", ""),
    (" signalling pathway", ""),
    (" pathway", ""),
    (" process", ""),
    (" tissue", ""),
    (" cells", ""),       # ``cnidocyte cells`` -> ``cnidocyte``
    (" cells", " cell"),  # ``T cells`` -> ``T cell``
    (" cell", ""),
    (" genes", ""),
    (" gene", ""),
)


def resolve_concept(
    name: str, *, species: str | None = None
) -> tuple[str, dict[str, str]] | None:
    """Resolve ``name`` to a concept CURIE and xrefs via ``CONCEPT_INDEX``.

    ``species`` is accepted for signature uniformity; concept entries
    (methodology, protocols, cell types, chemicals, units, and so on) are
    not species-specific, so the hint is ignored.

    If the primary lookup misses, retries with one common trailing noun
    trimmed (e.g. ``"cnidocyte cells"`` -> ``"cnidocyte"``, ``"T cells"``
    -> ``"T cell"``, ``"Wnt signaling pathway"`` -> ``"Wnt"``). The
    fallback fires only when the original term is not in the index, so
    the success path is unchanged for terms that already match.

    The returned ``xrefs`` dict carries an internal ``"_resolver_kind"``
    sentinel — ``"concept"`` for a primary lookup hit, ``"trimmed_concept"``
    when the trim-tail fallback was used. The dispatcher in
    ``enrich_claim`` pops this sentinel before storing ``xrefs`` on the
    claim, so it never reaches consumers; it carries the kind into
    ``build_sssom_mapping`` so the SSSOM fragment gets the right
    confidence floor and ``skos:closeMatch`` predicate for fallback hits.
    """
    del species
    used_fallback = False
    row = CONCEPT_INDEX.lookup(name)
    if row is None and isinstance(name, str):
        # Tolerant fallback: real LLM extractors (claude-sonnet-4-5
        # verified 2026-06-04) sometimes append common trailing nouns
        # to a canonical term. Try trimming them once and re-lookup.
        stripped = name.rstrip()
        lowered = stripped.lower()
        for match_suffix, replacement in _CONCEPT_TRIM_RULES:
            if lowered.endswith(match_suffix):
                trimmed = stripped[: -len(match_suffix)] + replacement
                trimmed = trimmed.rstrip()
                if trimmed:
                    row = CONCEPT_INDEX.lookup(trimmed)
                    if row is not None:
                        used_fallback = True
                        break
        if row is None:
            return None
    if row is None:
        return None
    xref = {
        "_resolver_kind": "trimmed_concept" if used_fallback else "concept",
        "curie": row["curie"],
        "label": row["label"],
        "ontology": row["ontology"],
    }
    if row.get("description"):
        xref["description"] = row["description"]
    return row["curie"], xref
