"""BulkTranscriptomicsAdapter conformance + behaviour tests."""
from __future__ import annotations

import pytest
import rdflib
from indicium_adapters import CURRENT_API_VERSION, DomainAdapter, SHACLProvider

from indicium_adapters_transcriptomics import BulkTranscriptomicsAdapter, accession_to_curie

# --- Protocol conformance ---------------------------------------------------


def test_implements_domain_adapter_protocol():
    a = BulkTranscriptomicsAdapter()
    assert isinstance(a, DomainAdapter)


def test_implements_shacl_provider_protocol():
    a = BulkTranscriptomicsAdapter()
    assert isinstance(a, SHACLProvider)


def test_api_version_matches_core():
    assert BulkTranscriptomicsAdapter.api_version == CURRENT_API_VERSION


def test_domain_id_is_canonical():
    assert BulkTranscriptomicsAdapter.domain_id == "bulk-transcriptomics"


# --- Static attributes ------------------------------------------------------


def test_qualifiers_are_a_nonempty_frozenset():
    qs = BulkTranscriptomicsAdapter.qualifiers
    assert isinstance(qs, frozenset)
    assert len(qs) > 0
    assert all(isinstance(q, str) and q == q.lower() for q in qs)


def test_qualifiers_include_core_claim_and_study_terms():
    qs = BulkTranscriptomicsAdapter.qualifiers
    # Claim-level
    assert "differentially_expressed_in" in qs
    assert "co_expressed_with" in qs
    assert "enriched_for" in qs
    # Study-level
    assert "sequenced_on" in qs
    assert "aligned_with" in qs
    assert "quantified_with" in qs
    assert "normalized_with" in qs


def test_ontology_prefixes_cover_required_axes():
    p = BulkTranscriptomicsAdapter.ontology_prefixes
    # Experimental / statistical
    assert {"EFO", "OBI", "STATO"} <= p.keys()
    # Bio / taxonomy / anatomy
    assert {"NCBITaxon", "CL", "UBERON"} <= p.keys()
    # Function / sequence / disease
    assert {"GO", "SO", "MONDO", "PATO"} <= p.keys()
    # Tooling
    assert "EDAM" in p
    # Gene identifiers
    assert {"HGNC", "MGI", "ENSEMBL", "NCBIGene"} <= p.keys()
    # Repositories
    assert {"GEO", "SRA", "ENA", "ArrayExpress", "ENCODE", "BioProject"} <= p.keys()


def test_ontology_prefix_values_are_iri_bases():
    for prefix, iri in BulkTranscriptomicsAdapter.ontology_prefixes.items():
        assert isinstance(iri, str)
        assert iri.startswith(("http://", "https://"))
        assert iri.endswith(("_", ":", "/")), (
            f"prefix {prefix!r} IRI must end with the boundary char used "
            f"by CURIE expansion (got {iri!r})"
        )


# --- extraction_context -----------------------------------------------------


def test_extraction_context_names_the_domain_and_qualifiers():
    a = BulkTranscriptomicsAdapter()
    ctx = a.extraction_context()
    assert "bulk transcriptomics" in ctx.lower()
    for q in ("differentially_expressed_in", "sequenced_on", "normalized_with"):
        assert q in ctx


# --- enrich_claim -----------------------------------------------------------


def test_enrich_claim_is_pure_does_not_mutate_input():
    a = BulkTranscriptomicsAdapter()
    original = {
        "id": "c1",
        "subject": "GSE168938",
        "qualifier": "consolidates_dataset",
        "object": "study",
    }
    snapshot = dict(original)
    a.enrich_claim(original)
    assert original == snapshot


def test_enrich_claim_returns_new_dict():
    a = BulkTranscriptomicsAdapter()
    claim = {"subject": "BRCA1", "qualifier": "expressed_in", "object": "breast"}
    out = a.enrich_claim(claim)
    assert out is not claim


def test_enrich_claim_resolves_geo_accession_to_curie():
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "GSE168938",
        "qualifier": "consolidates_dataset",
        "object": "paper",
    })
    assert out["ontology_terms"]["subject"] == "GEO:GSE168938"


def test_enrich_claim_resolves_sra_accession_in_object():
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "study",
        "qualifier": "consolidates_dataset",
        "object": "PRJEB104129",
    })
    assert out["ontology_terms"]["object"] == "BioProject:PRJEB104129"


def test_enrich_claim_leaves_unrecognized_terms_alone():
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "an unfamiliar token",
        "qualifier": "expressed_in",
        "object": "tissue",
    })
    assert "ontology_terms" not in out


@pytest.mark.parametrize(
    "accession, expected",
    [
        ("GSE168938", "GEO:GSE168938"),
        ("GSM4567890", "GEO:GSM4567890"),
        ("SRR15934891", "SRA:SRR15934891"),
        ("ERR15934891", "ENA:ERR15934891"),
        ("PRJEB104129", "BioProject:PRJEB104129"),
        ("E-MTAB-1234", "ArrayExpress:E-MTAB-1234"),
        ("ENCSR123ABC", "ENCODE:ENCSR123ABC"),
        ("phs000123.v1.p1", "dbGaP:phs000123.v1.p1"),
        ("S-BIAD123", "BioStudies:S-BIAD123"),
        ("not-an-accession", None),
        ("", None),
    ],
)
def test_accession_to_curie_table(accession, expected):
    assert accession_to_curie(accession) == expected


# --- SHACL ------------------------------------------------------------------


def test_shacl_shapes_returns_a_parsed_graph():
    a = BulkTranscriptomicsAdapter()
    g = a.shacl_shapes()
    assert isinstance(g, rdflib.Graph)
    assert len(g) > 0
