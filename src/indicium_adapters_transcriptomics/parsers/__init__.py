"""Repository-metadata parsers (v0.7.0).

Each parser takes a raw metadata source (a GEO SOFT text dump, an ENA
XML response, a BioStudies JSON document, an ISA-Tab directory) and
returns a list of **claim dicts** in the shape consumed by
``BulkTranscriptomicsAdapter.enrich_claim`` /
``SingleCellTranscriptomicsAdapter.enrich_claim``.

The parsers extract **study-level** assertions only (which platform was
used, which dataset is consolidated, which organism the samples came
from, etc.). Feature-level claims (``differentially_expressed_in`` /
``marker_of`` / ...) require the actual count tables and DEG outputs —
those stay the responsibility of the upstream LLM extractor.

Each parser is intentionally minimum-viable in v0.7.0: it covers the
fields most commonly populated in real metadata, ignores the rest, and
emits a small handful of claims per record. v0.7.x follow-ups can
deepen field coverage when a real use-case demands it.

Usage::

    from indicium_adapters_transcriptomics.parsers import parse_geo_soft
    from indicium_adapters_transcriptomics import BulkTranscriptomicsAdapter

    adapter = BulkTranscriptomicsAdapter()
    for claim in parse_geo_soft(soft_text):
        enriched = adapter.enrich_claim(claim)
        print(enriched.get("ontology_terms"))
"""
from __future__ import annotations

from indicium_adapters_transcriptomics.parsers.biostudies import (
    parse_biostudies_json,
)
from indicium_adapters_transcriptomics.parsers.ena_xml import parse_ena_xml
from indicium_adapters_transcriptomics.parsers.geo_soft import parse_geo_soft
from indicium_adapters_transcriptomics.parsers.isa_tab import parse_isa_tab

__all__ = [
    "parse_geo_soft",
    "parse_ena_xml",
    "parse_biostudies_json",
    "parse_isa_tab",
]
