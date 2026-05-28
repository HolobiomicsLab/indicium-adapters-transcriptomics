"""GEO SOFT format parser (v0.7.0).

`GEO SOFT <https://www.ncbi.nlm.nih.gov/geo/info/soft.html>`__ is the
Gene Expression Omnibus's native text serialisation. A SOFT document
contains one or more records, each opened by ``^TYPE = id`` and
followed by ``!field = value`` lines, until the next ``^`` or EOF.

Top-level record types we care about:

* ``^SERIES`` — the study (one per document). Carries GSE accession,
  title, platform reference(s), pubmed id, sample references, etc.
* ``^SAMPLE`` — one per biological sample. Carries GSM accession,
  organism, characteristics (tissue / time-point / treatment),
  instrument model, library protocol.
* ``^PLATFORM`` — sequencing platform / array. Carries GPL accession,
  title, technology, distribution organism.

This parser is **minimum viable**: it extracts the fields that almost
every transcriptomics GEO series populates and emits one claim per
distinct (series → field → value) triple. Multi-valued fields (like
multiple sample ids) emit one claim per value.

What we don't try to do in v0.7.0:

* Quoted / multi-line values (rare in practice; SOFT lines are usually
  single-line).
* Parsing ``!Sample_characteristics_ch1`` substructure ("tissue: X;
  age: Y") into separate fields — we keep the raw string.
* Cross-referencing GSE → GSM relationships beyond what the headers
  declare via ``!Series_sample_id``.
"""
from __future__ import annotations

import re
from typing import Final

# Record-opener pattern: ``^TYPE = id``. SOFT permits the ``=`` to be
# missing for some legacy variants — we tolerate either.
_RECORD_OPENER: Final[re.Pattern[str]] = re.compile(r"^\^([A-Z][A-Z_]*)\s*(?:=\s*(.*))?$")

# Field line pattern: ``!FieldName = value`` (whitespace tolerant).
_FIELD_LINE: Final[re.Pattern[str]] = re.compile(r"^!([A-Za-z][\w]*)\s*=\s*(.*)$")


def _tokenise(text: str) -> list[tuple[str, str, dict[str, list[str]]]]:
    """Split a SOFT document into ``(record_type, record_id, fields)``
    tuples in document order. ``fields`` is a ``{field_name → list of
    values}`` dict — multi-valued fields like ``!Series_sample_id``
    keep all their values."""
    records: list[tuple[str, str, dict[str, list[str]]]] = []
    cur_type: str | None = None
    cur_id: str = ""
    cur_fields: dict[str, list[str]] = {}

    def _flush() -> None:
        if cur_type is not None:
            records.append((cur_type, cur_id, cur_fields))

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if m := _RECORD_OPENER.match(line):
            _flush()
            cur_type = m.group(1)
            cur_id = (m.group(2) or "").strip()
            cur_fields = {}
            continue
        if cur_type is None:
            continue  # skip pre-record garbage
        if m := _FIELD_LINE.match(line):
            field = m.group(1)
            value = (m.group(2) or "").strip()
            cur_fields.setdefault(field, []).append(value)

    _flush()
    return records


def parse_geo_soft(text: str) -> list[dict]:
    """Parse GEO SOFT ``text`` and emit study-level claim dicts.

    Per ``^SERIES`` record, emits claims like:

    * ``{subject: <GSE>, qualifier: "consolidates_dataset", object: <GSM>}``
      for each ``!Series_sample_id`` value (so the bulk adapter
      surfaces a ``GEO:<GSM>`` CURIE on the object slot).
    * ``{subject: <GSE>, qualifier: "sequenced_on", object: <instrument>}``
      for each distinct ``!Sample_instrument_model`` value found
      among the document's sample records.
    * ``{subject: <GSE>, qualifier: "uses_database", object: <pubmed_id>}``
      when ``!Series_pubmed_id`` is set.

    Claims attach a ``species`` hint when ``!Sample_organism_ch1``
    contains a recognised model-organism name — so a downstream gene
    resolver can pick the right organism without re-deriving it from
    the symbol. Currently recognised: ``human`` / ``mouse`` /
    ``drosophila`` / ``celegans`` / ``zebrafish`` / ``nematostella``.
    """
    records = _tokenise(text)
    series_records = [(rid, fields) for rt, rid, fields in records if rt == "SERIES"]
    sample_records = [fields for rt, _, fields in records if rt == "SAMPLE"]

    claims: list[dict] = []
    for series_id, series_fields in series_records:
        gse = _first(series_fields, "Series_geo_accession") or series_id
        if not gse:
            continue
        species_hint = _detect_species(sample_records)

        # Build (qualifier, object) tuples first, then materialise
        # claims in one pass — avoids a closure capturing the loop
        # variables ``gse`` / ``species_hint`` (ruff B023).
        pending: list[tuple[str, str]] = []

        for gsm in series_fields.get("Series_sample_id", []):
            if gsm:
                pending.append(("consolidates_dataset", gsm))

        pubmed = _first(series_fields, "Series_pubmed_id")
        if pubmed:
            pending.append(("uses_database", f"PubMed:{pubmed}"))

        instruments: set[str] = set()
        for sf in sample_records:
            for inst in sf.get("Sample_instrument_model", []):
                if inst:
                    instruments.add(inst)
        for inst in sorted(instruments):
            pending.append(("sequenced_on", inst))

        for qualifier, obj in pending:
            claim: dict = {
                "subject": gse,
                "qualifier": qualifier,
                "object": obj,
            }
            if species_hint:
                claim["species"] = species_hint
            claims.append(claim)

    return claims


# --- helpers ---------------------------------------------------------------


def _first(fields: dict[str, list[str]], key: str) -> str:
    vs = fields.get(key) or []
    return vs[0] if vs else ""


_SPECIES_FROM_ORGANISM: Final[tuple[tuple[str, str], ...]] = (
    ("homo sapiens",              "human"),
    ("mus musculus",              "mouse"),
    ("drosophila melanogaster",   "drosophila"),
    ("caenorhabditis elegans",    "celegans"),
    ("danio rerio",               "zebrafish"),
    ("nematostella vectensis",    "nematostella"),
)


def _detect_species(sample_records: list[dict[str, list[str]]]) -> str | None:
    """Return the first recognised species tag found in any
    ``!Sample_organism_ch1`` value across the sample records."""
    for sf in sample_records:
        for organism in sf.get("Sample_organism_ch1", []):
            lc = organism.strip().lower()
            for needle, tag in _SPECIES_FROM_ORGANISM:
                if needle in lc:
                    return tag
    return None
