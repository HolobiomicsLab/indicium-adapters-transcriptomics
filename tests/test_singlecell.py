"""SingleCellTranscriptomicsAdapter conformance + behaviour tests.

v0.1.0 ships the adapter as a structural stub: discoverable, Protocol-compliant,
emits SHACL, but ``enrich_claim()`` is a no-op pending the v0.2.0 reference
single-cell paper. These tests pin that v0.1.0 behaviour so the v0.2.0 work has
a clear contract to extend.
"""
from __future__ import annotations

import rdflib
from indicium_adapters import CURRENT_API_VERSION, DomainAdapter, SHACLProvider

from indicium_adapters_transcriptomics import SingleCellTranscriptomicsAdapter

# --- Protocol conformance ---------------------------------------------------


def test_implements_domain_adapter_protocol():
    a = SingleCellTranscriptomicsAdapter()
    assert isinstance(a, DomainAdapter)


def test_implements_shacl_provider_protocol():
    a = SingleCellTranscriptomicsAdapter()
    assert isinstance(a, SHACLProvider)


def test_api_version_matches_core():
    assert SingleCellTranscriptomicsAdapter.api_version == CURRENT_API_VERSION


def test_domain_id_is_canonical():
    assert SingleCellTranscriptomicsAdapter.domain_id == "single-cell-transcriptomics"


# --- Static attributes ------------------------------------------------------


def test_qualifiers_are_cell_resolution_specific():
    qs = SingleCellTranscriptomicsAdapter.qualifiers
    assert isinstance(qs, frozenset)
    # Cell-resolution qualifiers that don't apply to bulk experiments.
    assert "marker_of_cell_type" in qs
    assert "enriched_in_cluster" in qs
    assert "localized_in_tissue" in qs
    assert "localized_in_spatial_domain" in qs


def test_ontology_prefixes_include_cell_and_anatomy_axes():
    p = SingleCellTranscriptomicsAdapter.ontology_prefixes
    # The two axes that single-cell most depends on.
    assert "CL" in p
    assert "UBERON" in p


# --- extraction_context -----------------------------------------------------


def test_extraction_context_names_the_domain_and_qualifiers():
    ctx = SingleCellTranscriptomicsAdapter().extraction_context()
    assert "single-cell" in ctx.lower()
    for q in ("marker_of_cell_type", "enriched_in_cluster", "localized_in_tissue"):
        assert q in ctx


# --- enrich_claim (no-op in v0.1.0) -----------------------------------------


def test_enrich_claim_is_noop_in_v0_1_0():
    """v0.1.0 stub: no enrichment yet, CL + UBERON lookups land in v0.2.0."""
    a = SingleCellTranscriptomicsAdapter()
    claim = {
        "subject": "GZMB",
        "qualifier": "marker_of_cell_type",
        "object": "cytotoxic T cell",
    }
    out = a.enrich_claim(claim)
    assert out == claim
    assert out is not claim, "must return a defensive copy, not the input"


def test_enrich_claim_does_not_mutate_input():
    a = SingleCellTranscriptomicsAdapter()
    original = {
        "subject": "CD8A",
        "qualifier": "marker_of_cell_type",
        "object": "T cell",
    }
    snapshot = dict(original)
    a.enrich_claim(original)
    assert original == snapshot


# --- SHACL ------------------------------------------------------------------


def test_shacl_shapes_returns_a_parsed_graph():
    g = SingleCellTranscriptomicsAdapter().shacl_shapes()
    assert isinstance(g, rdflib.Graph)
    assert len(g) > 0
