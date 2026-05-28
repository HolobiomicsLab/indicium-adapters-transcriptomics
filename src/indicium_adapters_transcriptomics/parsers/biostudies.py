"""BioStudies JSON parser.

`BioStudies <https://www.ebi.ac.uk/biostudies>`__ exposes each submission
as a JSON document with a top-level ``accno`` and a ``section`` tree.
Studies typically follow this shape::

    {
      "accno": "S-BIAD123",
      "section": {
        "type": "Study",
        "attributes": [
          {"name": "Title", "value": "..."},
          {"name": "Organism", "value": "Mus musculus"}
        ],
        "subsections": [
          {
            "type": "Author",
            "attributes": [{"name": "Name", "value": "..."}]
          },
          {
            "type": "Files",
            "links": [{"url": "..."}]
          }
        ]
      }
    }

The parser walks ``section`` and its ``subsections`` recursively,
collecting all ``attributes`` into a single flat ``{name: value}`` dict,
then emits study-level claims from it.
"""
from __future__ import annotations

import json
from typing import Any


def parse_biostudies_json(text_or_obj: str | dict) -> list[dict]:
    """Parse a BioStudies JSON document and return study-level claims.

    ``text_or_obj`` can be either the raw JSON string or an already
    parsed dict (convenient when the caller already has the record in
    memory).

    Recognised attributes (case-insensitive):

    * ``Organism``: adds a ``species`` hint on every emitted claim.
    * ``Pubmed`` or ``PubMed``: emits ``uses_database PubMed:<id>``.
    * ``Instrument``, ``Instrument model``, ``Platform``: emits
      ``sequenced_on <value>`` (one claim, first matching key wins).

    The top-level ``accno`` becomes the subject of every emitted claim.
    Nested ``links`` are not surfaced as claims: they are file
    references, not assertions.
    """
    if isinstance(text_or_obj, str):
        data = json.loads(text_or_obj)
    else:
        data = text_or_obj

    accno = (data.get("accno") or "").strip()
    if not accno:
        return []

    attrs = _collect_attributes(data.get("section") or {})
    species = _detect_species_from_attr(attrs.get("organism"))
    claims: list[dict] = []

    def _emit(qualifier: str, obj: str) -> None:
        """Append one claim, skipping empty objects and attaching the
        species hint when one was detected."""
        if not obj:
            return
        claim: dict = {
            "subject": accno,
            "qualifier": qualifier,
            "object": obj,
        }
        if species:
            claim["species"] = species
        claims.append(claim)

    pubmed = attrs.get("pubmed") or attrs.get("pmid")
    if pubmed:
        _emit("uses_database", f"PubMed:{pubmed}")

    # Only one instrument claim per study. Pick the first attribute key
    # that resolves to a non-empty value.
    for instrument_key in ("instrument", "instrument model", "platform"):
        inst = attrs.get(instrument_key)
        if inst:
            _emit("sequenced_on", inst)
            break

    return claims


# --- helpers --------------------------------------------------------------


def _collect_attributes(section: dict[str, Any]) -> dict[str, str]:
    """Flatten a section's attributes (and those of its sub-sections) into
    a single ``{name.lower(): value}`` dict.

    First-write wins on duplicates. Top-level attributes are treated as
    the most authoritative, so child sections never overwrite a
    parent-level value.
    """
    out: dict[str, str] = {}

    def _walk(sec: dict[str, Any]) -> None:
        for attr in sec.get("attributes") or []:
            name = (attr.get("name") or "").strip().lower()
            value = (attr.get("value") or "").strip()
            if name and value and name not in out:
                out[name] = value
        for sub in sec.get("subsections") or []:
            _walk(sub)

    _walk(section)
    return out


# Substring-to-tag mapping for organism detection. See the equivalent
# table in ``geo_soft.py``; the two parsers share the same set of
# recognised model organisms.
_SPECIES_FROM_NAME: tuple[tuple[str, str], ...] = (
    ("homo sapiens",              "human"),
    ("mus musculus",              "mouse"),
    ("drosophila melanogaster",   "drosophila"),
    ("caenorhabditis elegans",    "celegans"),
    ("danio rerio",               "zebrafish"),
    ("nematostella vectensis",    "nematostella"),
)


def _detect_species_from_attr(organism: str | None) -> str | None:
    """Map a BioStudies ``Organism`` attribute value to a species tag.

    Returns ``None`` when the value is missing or unrecognised.
    """
    if not organism:
        return None
    lc = organism.strip().lower()
    for needle, tag in _SPECIES_FROM_NAME:
        if needle in lc:
            return tag
    return None
