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
from importlib.resources import files
from typing import Final

# Columns of `data/gene_index.tsv` — declared here so the loader can be
# strict about the bundled schema. Adding columns means bumping this list
# and the lookup logic in tandem.
_GENE_INDEX_COLUMNS: Final[tuple[str, ...]] = (
    "symbol",
    "hgnc",
    "ensembl",
    "ncbigene",
    "nemve1",
    "nvertx",
    "uniprot",
    "description",
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


class GeneIndex:
    """Case-insensitive gene-symbol lookup with a few aliasing strategies.

    Pure data; loaded once at import time. The class only exposes a single
    ``lookup()`` method so adapters can stay decoupled from the on-disk
    format and the alias logic.
    """

    def __init__(self, rows: list[dict[str, str]]):
        self._rows = rows
        self._by_symbol: dict[str, dict[str, str]] = {
            row["symbol"].lower(): row for row in rows if row["symbol"]
        }
        self._by_nemve1: dict[str, dict[str, str]] = {
            row["nemve1"]: row for row in rows if row["nemve1"]
        }
        self._by_uniprot: dict[str, dict[str, str]] = {
            row["uniprot"]: row for row in rows if row["uniprot"]
        }

    def lookup(self, name: str) -> dict[str, str] | None:
        """Return the matching row as a dict, or ``None``.

        Tries (in order): gene symbol (case-insensitive) → Nemve1 / JGI id
        (case-sensitive) → UniProt accession (case-sensitive).
        """
        if not isinstance(name, str):
            return None
        token = name.strip()
        if not token:
            return None
        if row := self._by_symbol.get(token.lower()):
            return row
        if row := self._by_nemve1.get(token):
            return row
        if row := self._by_uniprot.get(token):
            return row
        return None

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
    """Load the static + OLS4-corpus + bio.tools + CL concept indexes.

    Four source TSVs are merged into the same ``ConceptIndex`` singleton:

    1. ``concept_index.tsv`` — TranscriptoKG modelling vocabulary
       (NCIT / INO / OBI / STATO / GO concepts, v0.3.0).
    2. ``corpus_concept_index.tsv`` — OLS4 label/synonym + prefix +
       substring matches across the 6-paper mini-corpus
       (EFO / OBI / NCIT / STATO platforms / protocols / stats, v0.3.1).
    3. ``biotools_concept_index.tsv`` — bio.tools registry hits for the
       bioinformatics software referenced in the corpus (v0.3.1).
    4. ``cl_concept_index.tsv`` — Cell Ontology hits for the cell-type
       vocabulary mentioned across the cell-atlas papers in the corpus
       (Sebé-Pedrós 2018 Nematostella, Tabula Muris/Sapiens, v0.5.0).

    Load order encodes a priority: a hypothetical collision (same term
    in multiple files) resolves to the **first** occurrence — paper-
    modelling-grounded (1) > OLS4 method-mention (2) > bio.tools-
    registry (3) > CL cell-type (4). In practice the build scripts
    deduplicate against earlier files so collisions are unexpected.
    """
    rows = _load_one_concept_tsv("concept_index.tsv")
    rows.extend(_load_one_concept_tsv("corpus_concept_index.tsv"))
    rows.extend(_load_one_concept_tsv("biotools_concept_index.tsv"))
    rows.extend(_load_one_concept_tsv("cl_concept_index.tsv"))
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
# HGNC values are stored already in CURIE form; the others are bare ids
# (`ENSG…`, `54361`, …) and get prefixed here.
_GENE_PRINCIPAL_CURIE_PREFIX: Final[tuple[tuple[str, str | None], ...]] = (
    ("hgnc", None),
    ("ensembl", "ENSEMBL"),
    ("ncbigene", "NCBIGene"),
)

# Xref columns surfaced under ``xrefs[slot]``. Bare ids (no prefixing).
_GENE_XREF_COLUMNS: Final[tuple[str, ...]] = (
    "hgnc", "ensembl", "ncbigene", "nemve1", "nvertx", "uniprot",
)


def resolve_gene(name: str) -> tuple[str, dict[str, str]] | None:
    """Resolve ``name`` to a gene CURIE + xrefs via ``GENE_INDEX``.

    Principal CURIE follows HGNC > Ensembl > NCBIGene priority; xrefs
    surface the bare HGNC / Ensembl / NCBIGene / Nemve1 / NvERTx /
    UniProt ids so downstream consumers can render any of them."""
    row = GENE_INDEX.lookup(name)
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


def resolve_anatomy(name: str) -> tuple[str, dict[str, str]] | None:
    """Resolve ``name`` to a UBERON anatomy CURIE + xrefs via ``ANATOMY_INDEX``."""
    row = ANATOMY_INDEX.lookup(name)
    if row is None:
        return None
    xref = {"uberon": row["uberon"], "label": row["label"]}
    if row.get("description"):
        xref["description"] = row["description"]
    return row["uberon"], xref


def resolve_concept(name: str) -> tuple[str, dict[str, str]] | None:
    """Resolve ``name`` to a concept CURIE + xrefs via ``CONCEPT_INDEX``."""
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
