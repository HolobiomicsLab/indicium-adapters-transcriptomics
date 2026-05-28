"""Repository-metadata parsers (added in v0.7.0).

Each parser takes a raw metadata source (a GEO SOFT text dump, an ENA XML
response, a BioStudies JSON document, an ISA-Tab directory) and returns a
list of claim dicts ready for ``BulkTranscriptomicsAdapter.enrich_claim``
or ``SingleCellTranscriptomicsAdapter.enrich_claim``.

The parsers only extract study-level assertions: which platform was used,
which dataset is consolidated, which organism the samples came from, and
so on. Feature-level claims (``differentially_expressed_in``,
``marker_of``, etc.) need the actual count tables and DEG outputs, which
stay the responsibility of the upstream LLM extractor.

Each parser is minimum-viable: it covers the fields most commonly populated
in real metadata, ignores the rest, and emits a small handful of claims
per record. Later releases can deepen field coverage when a real use-case
calls for it.

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
