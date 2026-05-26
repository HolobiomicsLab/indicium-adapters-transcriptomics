"""SingleCellTranscriptomicsAdapter conformance + behaviour tests.

v0.2.0 wires UBERON anatomy enrichment into ``enrich_claim()`` for the
spatial-localisation qualifiers. Cell-type (CL) enrichment is on the
v0.2.1 roadmap.
"""
from __future__ import annotations

import rdflib
from indicium_adapters import CURRENT_API_VERSION, DomainAdapter, SHACLProvider

from indicium_adapters_transcriptomics import SingleCellTranscriptomicsAdapter
from indicium_adapters_transcriptomics._enrichment import ANATOMY_INDEX

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


# --- enrich_claim — anatomy lookup (v0.2.0) ---------------------------------


def test_anatomy_index_has_a_meaningful_number_of_entries():
    """The v0.2.0 starter set is small but every entry is verified against
    UBERON. The floor is set at 3 so dropping below signals a regression
    in the build pipeline."""
    assert len(ANATOMY_INDEX) >= 3


def test_enrich_claim_does_not_mutate_input():
    a = SingleCellTranscriptomicsAdapter()
    original = {
        "subject": "MMP9",
        "qualifier": "localized_in_tissue",
        "object": "mesentery",
    }
    snapshot = dict(original)
    a.enrich_claim(original)
    assert original == snapshot


def test_enrich_claim_returns_new_dict():
    a = SingleCellTranscriptomicsAdapter()
    claim = {"subject": "x", "qualifier": "localized_in_tissue", "object": "mesentery"}
    out = a.enrich_claim(claim)
    assert out is not claim


def test_enrich_claim_resolves_uberon_anatomy_in_object_slot():
    """`localized_in_tissue` typically uses a gene as subject and an
    anatomical term as object — the object slot should pick up UBERON."""
    a = SingleCellTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "MMP9",
        "qualifier": "localized_in_tissue",
        "object": "mesentery",
    })
    assert out["ontology_terms"]["object"] == "UBERON:0002095"
    xrefs = out["xrefs"]["object"]
    assert xrefs["uberon"] == "UBERON:0002095"
    assert xrefs["label"].lower() == "mesentery"


def test_enrich_claim_is_case_insensitive_on_anatomy_terms():
    a = SingleCellTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "x",
        "qualifier": "localized_in_tissue",
        "object": "PHARYNX",
    })
    assert out["ontology_terms"]["object"] == "UBERON:0006562"


def test_enrich_claim_leaves_cnidarian_specific_terms_alone():
    """Cnidarian-specific anatomy (physa, mesoglea, gastrodermis, etc.) has
    no UBERON entry as of 2026 and is deliberately absent from the bundle.
    The adapter must not invent a CURIE for such terms."""
    a = SingleCellTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "MMP14",
        "qualifier": "localized_in_spatial_domain",
        "object": "physa",
    })
    assert "ontology_terms" not in out


def test_enrich_claim_resolves_both_subject_and_object_slots():
    """An adapter that only enriches one of the two slots would be a
    regression — both must be tried."""
    a = SingleCellTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "mesentery",
        "qualifier": "localized_in_tissue",
        "object": "pharynx",
    })
    assert out["ontology_terms"]["subject"] == "UBERON:0002095"
    assert out["ontology_terms"]["object"] == "UBERON:0006562"


# --- SHACL ------------------------------------------------------------------


def test_shacl_shapes_returns_a_parsed_graph():
    g = SingleCellTranscriptomicsAdapter().shacl_shapes()
    assert isinstance(g, rdflib.Graph)
    assert len(g) > 0
