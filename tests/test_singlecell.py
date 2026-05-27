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
    The adapter must not invent a CURIE for such terms.

    Pre-v0.4.0 this test could assert ``"ontology_terms" not in out`` because
    single-cell had no gene resolver and so MMP14 (the subject, a gene)
    also went unresolved. v0.4.0 wires the gene resolver in — MMP14 now
    correctly resolves to HGNC:7160 — so the test narrows to the object
    slot, which is what the cnidarian-anatomy guarantee actually covers.
    """
    a = SingleCellTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "MMP14",
        "qualifier": "localized_in_spatial_domain",
        "object": "physa",
    })
    assert "object" not in out.get("ontology_terms", {})


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


# --- gene enrichment for single-cell (v0.4.0) -------------------------------


def test_singlecell_resolves_gene_subject_for_marker_of_cell_type():
    """v0.4.0: single-cell now wires the shared gene resolver. A claim
    like ``marker_of_cell_type: WNT4 cnidocyte`` resolves ``WNT4`` to
    its HGNC CURIE — pre-v0.4.0 it was unenriched because single-cell
    only had anatomy + concept resolvers."""
    a = SingleCellTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "WNT4",
        "qualifier": "marker_of_cell_type",
        "object": "cnidocyte",
    })
    assert out["ontology_terms"]["subject"] == "HGNC:12783"
    xrefs = out["xrefs"]["subject"]
    assert xrefs.get("hgnc") == "HGNC:12783"
    # cnidocyte is a Nematostella-specific cell type; no CL CURIE yet (v0.5.0).
    assert "object" not in out.get("ontology_terms", {})


def test_singlecell_gene_resolver_does_not_override_anatomy():
    """Single-cell resolver priority is gene → anatomy → concept. For
    a slot whose value is an anatomy term (not a gene), anatomy still
    wins because gene-lookup returns ``None`` and the cascade falls
    through. ``localized_in_tissue: <gene> mesentery`` keeps the
    anatomy CURIE on the object slot."""
    a = SingleCellTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "MMP9",
        "qualifier": "localized_in_tissue",
        "object": "mesentery",
    })
    assert out["ontology_terms"]["subject"] == "HGNC:7176"
    assert out["ontology_terms"]["object"] == "UBERON:0002095"


# --- concept enrichment scaffolding (v0.3.0) --------------------------------


def test_concept_lookup_path_is_wired_into_singlecell_adapter(monkeypatch):
    """Wiring test for the concept-lookup branch in the SingleCell
    `enrich_claim`. Monkey-patches a single-row index in isolation from
    the bundled TSV: a term that does not match anatomy must still
    resolve via the concept index when an entry exists."""
    from indicium_adapters_transcriptomics import _enrichment, singlecell

    fake_row = {
        "term": "rna-seq",
        "curie": "OBI:0001177",
        "label": "RNA-seq assay",
        "description": "An assay that sequences expressed RNA.",
        "ontology": "OBI",
        "source": "monkeypatch fixture",
    }
    fake_index = _enrichment.ConceptIndex([fake_row])
    # v0.4.0: single-cell no longer imports CONCEPT_INDEX directly — it
    # dispatches through the shared resolve_concept helper, which reads
    # _enrichment.CONCEPT_INDEX. Patch target is the source of truth.
    monkeypatch.setattr(_enrichment, "CONCEPT_INDEX", fake_index)

    a = singlecell.SingleCellTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "study_42",
        "qualifier": "localized_in_spatial_domain",
        "object": "rna-seq",
    })
    assert out["ontology_terms"]["object"] == "OBI:0001177"
    assert out["xrefs"]["object"]["ontology"] == "OBI"


# --- SHACL ------------------------------------------------------------------


def test_shacl_shapes_returns_a_parsed_graph():
    g = SingleCellTranscriptomicsAdapter().shacl_shapes()
    assert isinstance(g, rdflib.Graph)
    assert len(g) > 0
