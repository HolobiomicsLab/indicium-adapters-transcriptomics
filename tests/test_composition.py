"""Cross-adapter composition invariants.

Mirrors the metabolomics package's ``test_composition.py`` contract: sibling
adapters in the same wheel MUST have disjoint qualifier sets so a composed
adapter's qualifier vocabulary is a clean union with no shadowing.
"""
from __future__ import annotations

from indicium_adapters import compose_adapters, discover_adapters

from indicium_adapters_transcriptomics import (
    BulkTranscriptomicsAdapter,
    SingleCellTranscriptomicsAdapter,
)


def test_both_adapters_discoverable_via_entry_points():
    found = discover_adapters()
    assert "bulk-transcriptomics" in found
    assert "single-cell-transcriptomics" in found


def test_domain_ids_are_distinct():
    assert (
        BulkTranscriptomicsAdapter.domain_id
        != SingleCellTranscriptomicsAdapter.domain_id
    )


def test_qualifier_sets_are_disjoint():
    """Required by the indicium-adapters composition contract. Overlap would
    break ``compose_adapters`` ordering guarantees and confuse Perspicacité's
    extraction prompt (it consumes ``extraction_context()`` to constrain
    LLM claim extraction)."""
    overlap = (
        BulkTranscriptomicsAdapter.qualifiers
        & SingleCellTranscriptomicsAdapter.qualifiers
    )
    assert overlap == frozenset(), (
        f"Qualifiers overlap between Bulk and SingleCell adapters: {overlap}"
    )


def test_compose_adapters_unions_qualifiers():
    composite = compose_adapters([
        BulkTranscriptomicsAdapter(),
        SingleCellTranscriptomicsAdapter(),
    ])
    expected_qualifiers = (
        BulkTranscriptomicsAdapter.qualifiers
        | SingleCellTranscriptomicsAdapter.qualifiers
    )
    assert composite.qualifiers == expected_qualifiers
