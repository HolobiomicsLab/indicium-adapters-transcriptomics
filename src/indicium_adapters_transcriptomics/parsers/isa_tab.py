"""ISA-Tab parser.

`ISA-Tab <https://isa-specs.readthedocs.io/>`__ is a tab-separated
serialisation of an Investigation / Study / Assay tree. Each submission
lives in a directory containing:

* exactly one ``i_*.txt`` (the investigation file): top-level metadata
  such as title, authors, contacts, publications.
* one or more ``s_*.txt`` (study files, sample sheets): ``Source Name``,
  ``Sample Name``, ``Characteristics[organism]``,
  ``Characteristics[organism part]``, etc.
* one or more ``a_*.txt`` (assay files): per-assay measurement metadata
  such as ``Protocol REF``, ``Parameter Value[instrument model]``,
  ``Data Files``.

The investigation file has its own quirky shape. Each line is
``Section Header\\tcolumn0\\tcolumn1\\t...``, organised into blocks
named ``ONTOLOGY SOURCE REFERENCE``, ``INVESTIGATION``, ``STUDY``,
``PUBLICATION``, etc. Each block has a fixed set of field rows.

This parser stays narrow:

* From ``i_*.txt``: extracts ``Study Identifier`` and ``Study PubMed ID``.
* From any ``s_*.txt``: extracts distinct ``Characteristics[organism]``
  values, used to derive the species hint.
* From any ``a_*.txt``: extracts distinct
  ``Parameter Value[instrument model]`` and ``Parameter Value[platform]``
  values.

That is enough to emit ``consolidates_dataset``, ``sequenced_on`` and
``uses_database`` claims for a typical ISA submission.
"""
from __future__ import annotations

import csv
from pathlib import Path


def parse_isa_tab(directory: str | Path) -> list[dict]:
    """Parse an ISA-Tab submission directory and emit study-level claims.

    ``directory`` must contain at least one ``i_*.txt`` investigation
    file. Study and assay files are picked up by their ``s_*.txt`` and
    ``a_*.txt`` prefixes. Returns an empty list when the directory has
    no investigation file.
    """
    root = Path(directory)
    if not root.is_dir():
        return []

    investigation_files = sorted(root.glob("i_*.txt"))
    if not investigation_files:
        return []

    inv_fields = _parse_investigation(investigation_files[0])
    # Prefer ``Study Identifier`` but fall back to ``Investigation
    # Identifier`` for submissions that only fill in the top level.
    study_id = (
        inv_fields.get("Study Identifier", [""])[0]
        or inv_fields.get("Investigation Identifier", [""])[0]
    )
    if not study_id:
        return []

    pubmeds = [v for v in inv_fields.get("Study PubMed ID", []) if v]

    # Sample organisms from any ``s_*.txt`` file.
    organisms: list[str] = []
    for s_path in sorted(root.glob("s_*.txt")):
        organisms.extend(_collect_column(s_path, "Characteristics[organism]"))
    species = _detect_species(organisms)

    # Instrument models from any ``a_*.txt`` file. We try the canonical
    # ``Parameter Value[instrument model]`` column first, then fall back
    # to ``Parameter Value[platform]`` for submissions that file the
    # instrument under the platform column.
    instruments: list[str] = []
    for a_path in sorted(root.glob("a_*.txt")):
        for col in (
            "Parameter Value[instrument model]",
            "Parameter Value[platform]",
        ):
            instruments.extend(_collect_column(a_path, col))

    claims: list[dict] = []

    def _emit(qualifier: str, obj: str) -> None:
        """Append one claim, skipping empty objects and attaching the
        species hint when available."""
        if not obj:
            return
        claim: dict = {
            "subject": study_id,
            "qualifier": qualifier,
            "object": obj,
        }
        if species:
            claim["species"] = species
        claims.append(claim)

    for pm in pubmeds:
        _emit("uses_database", f"PubMed:{pm}")
    for inst in _dedupe(instruments):
        _emit("sequenced_on", inst)

    return claims


# --- investigation-file parser --------------------------------------------


def _parse_investigation(path: Path) -> dict[str, list[str]]:
    """Parse an ISA investigation TSV into ``{field: list of values}``.

    Investigation files use TSV where column 0 is the field name and the
    remaining columns hold one value per study (often a single value per
    row). Empty cells are dropped and duplicate field rows are merged.
    """
    out: dict[str, list[str]] = {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row:
                continue
            key = (row[0] or "").strip()
            if not key or key.startswith("#"):
                continue
            values = [c.strip() for c in row[1:] if c and c.strip()]
            if not values:
                continue
            out.setdefault(key, []).extend(values)
    return out


# --- study / assay column collector ---------------------------------------


def _collect_column(path: Path, column_name: str) -> list[str]:
    """Return all non-empty values from the given column of a TSV.

    Matches the column header case-sensitively, since ISA-Tab columns
    are conventionally capitalised exactly. Returns an empty list when
    the file is missing, empty, or has no such column.
    """
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return []
        try:
            idx = header.index(column_name)
        except ValueError:
            return []
        values: list[str] = []
        for row in reader:
            if idx < len(row):
                v = row[idx].strip()
                if v:
                    values.append(v)
    return values


def _dedupe(seq: list[str]) -> list[str]:
    """Return ``seq`` with duplicates removed, preserving first-seen order."""
    seen: set[str] = set()
    out: list[str] = []
    for s in seq:
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


# Substring-to-tag mapping for organism detection. Same recognised set as
# the other parsers in this package.
_SPECIES_FROM_NAME: tuple[tuple[str, str], ...] = (
    ("homo sapiens",              "human"),
    ("mus musculus",              "mouse"),
    ("drosophila melanogaster",   "drosophila"),
    ("caenorhabditis elegans",    "celegans"),
    ("danio rerio",               "zebrafish"),
    ("nematostella vectensis",    "nematostella"),
)


def _detect_species(organisms: list[str]) -> str | None:
    """Return the first recognised species tag found in ``organisms``."""
    for organism in organisms:
        lc = organism.strip().lower()
        for needle, tag in _SPECIES_FROM_NAME:
            if needle in lc:
                return tag
    return None
