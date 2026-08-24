"""Cross-adapter composition invariants.

Mirrors the metabolomics package's ``test_composition.py`` contract: sibling
adapters in the same wheel MUST have disjoint qualifier sets so a composed
adapter's qualifier vocabulary is a clean union with no shadowing. The
classification vocabulary is the opposite case: it is package-level, so both
siblings must declare the same terms under the same ``primary_domain``.
"""
from __future__ import annotations

from indicium_adapters import compose_adapters, discover_adapters

from indicium_adapters_transcriptomics import (
    BulkTranscriptomicsAdapter,
    SingleCellTranscriptomicsAdapter,
)
from indicium_adapters_transcriptomics._vocabulary import (
    TRANSCRIPTOMICS_SUBDOMAINS,
    TRANSCRIPTOMICS_TECHNIQUES,
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


def test_both_adapters_declare_the_transcriptomics_primary_domain():
    """ASB groups sibling adapters into one field by ``primary_domain``. An
    adapter that declares none contributes no domain, and its papers
    classify as "other"."""
    assert BulkTranscriptomicsAdapter.primary_domain == "transcriptomics"
    assert SingleCellTranscriptomicsAdapter.primary_domain == "transcriptomics"


def test_classification_vocabulary_is_shared_by_both_adapters():
    """Unlike qualifiers, subdomains and techniques are package-level: the
    consumer unions them into one set, so a divergent sibling would silently
    widen the declared vocabulary."""
    assert BulkTranscriptomicsAdapter.subdomains == TRANSCRIPTOMICS_SUBDOMAINS
    assert SingleCellTranscriptomicsAdapter.subdomains == TRANSCRIPTOMICS_SUBDOMAINS
    assert BulkTranscriptomicsAdapter.techniques == TRANSCRIPTOMICS_TECHNIQUES
    assert SingleCellTranscriptomicsAdapter.techniques == TRANSCRIPTOMICS_TECHNIQUES
