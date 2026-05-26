"""BulkTranscriptomicsAdapter conformance + behaviour tests."""
from __future__ import annotations

import pytest
import rdflib
from indicium_adapters import CURRENT_API_VERSION, DomainAdapter, SHACLProvider

from indicium_adapters_transcriptomics import BulkTranscriptomicsAdapter, accession_to_curie
from indicium_adapters_transcriptomics._enrichment import GENE_INDEX

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


# --- gene-symbol enrichment (v0.1.1) ----------------------------------------


def test_gene_index_has_a_meaningful_number_of_entries():
    """v0.1.2 bundled index covers the full HGNC complete set plus the G2V2
    Nematostella cross-reference — ~50k rows. Floor is set well below the
    actual count so the test signals a regression rather than tracks the
    exact number."""
    assert len(GENE_INDEX) >= 30_000


def test_enrich_claim_resolves_gene_symbol_to_hgnc():
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "WNT4",
        "qualifier": "regulated_by",
        "object": "cWnt pathway",
    })
    assert out["ontology_terms"]["subject"] == "HGNC:12783"


def test_enrich_claim_resolves_gene_symbol_case_insensitive():
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "mmp9",
        "qualifier": "expressed_in",
        "object": "mesenteries",
    })
    assert out["ontology_terms"]["subject"] == "HGNC:7176"


def test_enrich_claim_resolves_nemve1_jgi_id():
    """A Nemve1 / JGI gene-model id in `subject` should resolve to the same
    HGNC CURIE as the gene symbol (via the same TSV row)."""
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "gw.51.96.1",
        "qualifier": "expressed_in",
        "object": "physa",
    })
    # gw.51.96.1 is FGF8 in the bundle
    assert out["ontology_terms"]["subject"] == "HGNC:3686"


def test_enrich_claim_resolves_uniprot_accession():
    """A canonical human UniProt accession should resolve to its HGNC CURIE
    via the bundled HGNC complete set."""
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "P56705",   # canonical human WNT4
        "qualifier": "expressed_in",
        "object": "oral pole",
    })
    assert out["ontology_terms"]["subject"] == "HGNC:12783"


def test_enrich_claim_attaches_full_xrefs_for_a_known_gene():
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "WNT4",
        "qualifier": "regulated_by",
        "object": "cWnt pathway",
    })
    xrefs = out["xrefs"]["subject"]
    # Spot-check three independent authorities: HGNC (the principal), the
    # paired Ensembl human gene and the canonical UniProt accession that
    # HGNC ships first under WNT4's `uniprot_ids`.
    assert xrefs.get("hgnc") == "HGNC:12783"
    assert xrefs.get("ensembl") == "ENSG00000162552"
    assert xrefs.get("uniprot") == "P56705"


def test_nemve1_id_resolves_to_human_ortholog_via_symbol_fallback():
    """A G2V2 Nematostella row whose UniProt did not match an HGNC row should
    still pick up the HGNC ortholog via the Preferred_name → HGNC symbol
    fallback (added in v0.1.2). FGF8 is one such case: G2V2's first UniProt
    accession is non-human, but Preferred_name='fgf8' matches HGNC's FGF8."""
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "gw.51.96.1",     # FGF8 JGI Nemve1 id
        "qualifier": "expressed_in",
        "object": "physa",
    })
    assert out["ontology_terms"]["subject"] == "HGNC:3686"
    xrefs = out["xrefs"]["subject"]
    assert xrefs.get("nemve1") == "gw.51.96.1"
    assert xrefs.get("hgnc") == "HGNC:3686"


def test_pure_human_gene_resolves_even_without_nematostella_ortholog():
    """The v0.1.2 build script also emits one row per HGNC gene that has no
    Nematostella mapping in G2V2, so any plain human-gene mention still
    resolves. MAPK1 is one such case — it's in HGNC, has no clear Nemve1
    ortholog in G2V2."""
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "MAPK1",
        "qualifier": "regulated_by",
        "object": "MEK1",
    })
    assert out["ontology_terms"]["subject"] == "HGNC:6871"
    xrefs = out["xrefs"]["subject"]
    assert xrefs.get("hgnc") == "HGNC:6871"
    # Pure-HGNC rows carry no Nematostella xrefs.
    assert "nemve1" not in xrefs


def test_enrich_claim_object_slot_also_resolved():
    """Same lookup logic must apply to both subject and object slots."""
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "U0126 treatment",
        "qualifier": "differentially_expressed_in",
        "object": "WNT4",
    })
    assert out["ontology_terms"]["object"] == "HGNC:12783"


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
