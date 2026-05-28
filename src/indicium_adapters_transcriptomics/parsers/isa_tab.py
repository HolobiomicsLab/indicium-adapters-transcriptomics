"""ISA-Tab parser (v0.7.0).

`ISA-Tab <https://isa-specs.readthedocs.io/>`__ is a tab-separated
serialisation of an Investigation / Study / Assay tree. Each
submission lives in a directory containing:

* exactly one ``i_*.txt`` — the **investigation** file (top-level
  metadata: title, authors, contacts, publications, etc.).
* one or more ``s_*.txt`` — **study** files (sample sheets:
  ``Source Name``, ``Sample Name``, ``Characteristics[organism]``,
  ``Characteristics[organism part]``, …).
* one or more ``a_*.txt`` — **assay** files (per-assay measurement
  metadata: ``Protocol REF``, ``Parameter Value[instrument model]``,
  ``Data Files``, …).

The investigation file has its own quirky shape: each line is
``Section Header\\tcolumn0\\tcolumn1\\t...``, organised by ``ONTOLOGY
SOURCE REFERENCE``, ``INVESTIGATION``, ``STUDY``, ``PUBLICATION``,
... blocks. Each block has a fixed set of field rows.

For v0.7.0 we keep the parser narrow:

* From ``i_*.txt`` — extract ``Study Identifier`` and
  ``Study PubMed ID`` only.
* From any ``s_*.txt`` — extract distinct ``Characteristics[organism]``
  values for the species hint.
* From any ``a_*.txt`` — extract distinct
  ``Parameter Value[instrument model]`` values (or
  ``Parameter Value[platform]``).

That's sufficient to emit ``consolidates_dataset`` + ``sequenced_on``
+ ``uses_database`` claims for a typical ISA submission.
"""
from __future__ import annotations

import csv
from pathlib import Path


def parse_isa_tab(directory: str | Path) -> list[dict]:
    """Parse an ISA-Tab submission directory and emit study-level claim dicts.

    ``directory`` must contain at least one ``i_*.txt`` investigation
    file. Study and assay files are picked up by their ``s_*.txt`` /
    ``a_*.txt`` prefix.

    Returns an empty list if the directory has no ``i_*.txt`` file.
    """
    root = Path(directory)
    if not root.is_dir():
        return []

    investigation_files = sorted(root.glob("i_*.txt"))
    if not investigation_files:
        return []

    inv_fields = _parse_investigation(investigation_files[0])
    study_id = (
        inv_fields.get("Study Identifier", [""])[0]
        or inv_fields.get("Investigation Identifier", [""])[0]
    )
    if not study_id:
        return []

    pubmeds = [v for v in inv_fields.get("Study PubMed ID", []) if v]

    # Sample organism(s) from any s_*.txt
    organisms: list[str] = []
    for s_path in sorted(root.glob("s_*.txt")):
        organisms.extend(_collect_column(s_path, "Characteristics[organism]"))
    species = _detect_species(organisms)

    # Instrument model(s) from any a_*.txt
    instruments: list[str] = []
    for a_path in sorted(root.glob("a_*.txt")):
        for col in (
            "Parameter Value[instrument model]",
            "Parameter Value[platform]",
        ):
            instruments.extend(_collect_column(a_path, col))

    claims: list[dict] = []

    def _emit(qualifier: str, obj: str) -> None:
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
    """ISA investigation files are TSV where row 0 of each block is the
    field name and the remaining columns are values across "studies"
    (often a single value per row). We collapse to ``{field → list of
    non-empty values}``."""
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
    """Read a TSV with a header row and return all non-empty values in
    the column whose header matches ``column_name`` (case-sensitive,
    as ISA-Tab columns are conventionally capitalised exactly)."""
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
    seen: set[str] = set()
    out: list[str] = []
    for s in seq:
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


_SPECIES_FROM_NAME: tuple[tuple[str, str], ...] = (
    ("homo sapiens",              "human"),
    ("mus musculus",              "mouse"),
    ("drosophila melanogaster",   "drosophila"),
    ("caenorhabditis elegans",    "celegans"),
    ("danio rerio",               "zebrafish"),
    ("nematostella vectensis",    "nematostella"),
)


def _detect_species(organisms: list[str]) -> str | None:
    for organism in organisms:
        lc = organism.strip().lower()
        for needle, tag in _SPECIES_FROM_NAME:
            if needle in lc:
                return tag
    return None
