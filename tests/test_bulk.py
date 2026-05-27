"""BulkTranscriptomicsAdapter conformance + behaviour tests."""
from __future__ import annotations

import pytest
import rdflib
from indicium_adapters import CURRENT_API_VERSION, DomainAdapter, SHACLProvider

from indicium_adapters_transcriptomics import BulkTranscriptomicsAdapter, accession_to_curie
from indicium_adapters_transcriptomics._enrichment import CONCEPT_INDEX, GENE_INDEX

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


# --- concept enrichment scaffolding (v0.3.0) --------------------------------


def test_concept_index_is_loadable():
    """The CONCEPT_INDEX singleton must load at import time without raising,
    even when the TSV is empty (e.g. pending OLS4 verification)."""
    assert isinstance(len(CONCEPT_INDEX), int)
    assert len(CONCEPT_INDEX) >= 0


def test_concept_lookup_returns_none_for_unknown_term():
    assert CONCEPT_INDEX.lookup("an entirely unfamiliar term") is None


def test_concept_lookup_handles_empty_and_non_str_input():
    assert CONCEPT_INDEX.lookup("") is None
    assert CONCEPT_INDEX.lookup(None) is None  # type: ignore[arg-type]


def test_enrich_claim_is_safe_noop_for_unknown_concept_terms():
    """When neither subject nor object resolves through any of the
    enrichment indices (gene / accession / anatomy / concept), the claim
    must pass through untouched — the concept fallback must not invent
    a CURIE."""
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "completely-bogus-term-zzz",
        "qualifier": "differentially_expressed_in",
        "object": "U0126",
    })
    assert "ontology_terms" not in out


# --- corpus concept index (v0.3.1) ------------------------------------------


def test_corpus_concept_terms_resolve_via_concept_index():
    """v0.3.1 ships a corpus-driven concept index harvested from the
    6-paper mini-corpus. 16 terms resolved via OLS4 across three match
    modes (exact / prefix / substring). Pin each one so a regression
    in either the harvester or the dual-TSV loader is caught
    immediately."""
    expected = {
        # term              (curie,           ontology, axis,         mode)
        # --- exact-match (label or synonym == corpus term) ---
        "illumina":         ("NCIT:C215820",  "NCIT",   "platform",   "exact"),
        "mars-seq":         ("EFO:0008796",   "EFO",    "protocol",   "exact"),
        "smart-seq2":       ("EFO:0008931",   "EFO",    "protocol",   "exact"),
        "drop-seq":         ("EFO:0008722",   "EFO",    "protocol",   "exact"),
        "slide-seq":        ("EFO:0009920",   "EFO",    "protocol",   "exact"),
        "t-test":           ("NCIT:C53231",   "NCIT",   "stat_test",  "exact"),
        # --- prefix-match (corpus term = whole-word prefix of label) ---
        "pacbio":           ("EFO:0920006",   "EFO",    "platform",   "prefix"),
        "oxford nanopore":  ("OBI:0002750",   "OBI",    "platform",   "prefix"),
        "bulk rna-seq":     ("OBI:0003090",   "OBI",    "protocol",   "prefix"),
        "visium":           ("EFO:0010961",   "EFO",    "protocol",   "prefix"),
        "wilcoxon":         ("NCIT:C53246",   "NCIT",   "stat_test",  "prefix"),
        "benjamini-hochberg": ("NCIT:C61596", "NCIT",   "stat_test",  "prefix"),
        # --- substring-match (corpus term = whole-word substring of label;
        # used for instrument families whose canonical EFO label
        # prepends a vendor name) ---
        "nextseq":          ("EFO:0008566",   "EFO",    "platform",   "substring"),
        "novaseq":          ("EFO:0008637",   "EFO",    "platform",   "substring"),
        "hiseq":            ("EFO:0004203",   "EFO",    "platform",   "substring"),
        "sequel":           ("EFO:0008630",   "EFO",    "platform",   "substring"),
    }
    for term, (curie, ontology, axis, mode) in expected.items():
        row = CONCEPT_INDEX.lookup(term)
        assert row is not None, f"corpus term {term!r} should resolve"
        assert row["curie"] == curie, (
            f"{term!r}: expected {curie}, got {row['curie']}"
        )
        assert row["ontology"] == ontology
        # source column: "corpus extraction (<axis>, <mode>-match)"
        assert axis in row["source"]
        assert f"{mode}-match" in row["source"]


def test_corpus_concept_resolution_is_end_to_end_through_enrich_claim():
    """End-to-end: a claim mentioning a corpus-resolved term must get
    the CURIE attached by `enrich_claim()` without monkey-patching."""
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "study_42",
        "qualifier": "sequenced_on",
        "object": "Illumina",
    })
    assert out["ontology_terms"]["object"] == "NCIT:C215820"
    xref = out["xrefs"]["object"]
    assert xref["curie"] == "NCIT:C215820"
    assert xref["ontology"] == "NCIT"
    assert xref["label"] == "Illumina"


def test_concept_index_carries_static_corpus_and_biotools_rows():
    """The triple-TSV loader must surface all three sources through the
    same singleton — static (v0.3.0) + OLS4 corpus (v0.3.1) + bio.tools
    (v0.3.1). Lower bound only, so adding rows in any file does not
    break the test."""
    n = len(CONCEPT_INDEX)
    # 31 static + 16 OLS4-corpus + 9 bio.tools = 56 minimum.
    assert n >= 31 + 16 + 9, f"expected ≥56 concept rows, got {n}"


# --- bio.tools concept index (v0.3.1) ---------------------------------------


def test_biotools_software_terms_resolve_via_concept_index():
    """v0.3.1 also ships a bio.tools-driven concept index covering the
    bioinformatics software referenced in the corpus. Nine terms
    resolved via the bio.tools registry (direct-ID lookup + search
    fallback with case + whitespace + hyphen normalisation —
    ``Bowtie2`` (corpus) ↔ ``Bowtie 2`` (bio.tools name) for example)."""
    expected = {
        "star":       ("biotools:star",         "aligner_quant"),
        "bowtie":     ("biotools:bowtie",       "aligner_quant"),
        "bowtie2":    ("biotools:bowtie2",      "aligner_quant"),
        "minimap2":   ("biotools:minimap2",     "aligner_quant"),
        "rsem":       ("biotools:rsem",         "aligner_quant"),
        "stringtie":  ("biotools:stringtie",    "aligner_quant"),
        "deseq2":     ("biotools:deseq2",       "de_tool"),
        "seurat":     ("biotools:seurat",       "sc_tool"),
        "picard":     ("biotools:picard_tools", "qc"),
    }
    for term, (curie, axis) in expected.items():
        row = CONCEPT_INDEX.lookup(term)
        assert row is not None, f"bio.tools term {term!r} should resolve"
        assert row["curie"] == curie
        assert row["ontology"] == "biotools"
        assert axis in row["source"]


def test_biotools_resolution_is_end_to_end_through_enrich_claim():
    """A claim mentioning a bio.tools-resolved software name gets the
    biotools CURIE attached via the standard enrich_claim path.

    ``minimap2`` is used here because it has no HGNC gene-symbol
    collision — so it resolves to ``biotools:minimap2`` regardless of
    whether the v0.4.0 qualifier-aware priority order or the default
    order is in effect. The qualifier-aware path is exercised
    separately by ``test_aligned_with_*`` below."""
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "study_42",
        "qualifier": "aligned_with",
        "object": "minimap2",
    })
    assert out["ontology_terms"]["object"] == "biotools:minimap2"
    xref = out["xrefs"]["object"]
    assert xref["curie"] == "biotools:minimap2"
    assert xref["ontology"] == "biotools"
    assert xref["label"] == "Minimap2"


# --- qualifier-aware resolution (v0.4.0) ------------------------------------


def test_aligned_with_qualifier_resolves_star_as_aligner_not_gene():
    """v0.4.0 contract: study-level qualifiers (``aligned_with`` here)
    use ``concept`` → ``accession`` → ``gene`` resolver priority. So
    a claim ``aligned_with: <study> STAR`` correctly resolves ``STAR``
    to ``biotools:star`` (the aligner) — even though ``STAR`` is also
    a valid HGNC gene symbol (HGNC:11359, steroidogenic acute
    regulatory protein). Pre-v0.4.0 this resolved to the gene because
    gene-lookup always ran first."""
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "study_42",
        "qualifier": "aligned_with",
        "object": "STAR",
    })
    assert out["ontology_terms"]["object"] == "biotools:star"
    assert out["xrefs"]["object"]["ontology"] == "biotools"


def test_expressed_in_qualifier_keeps_gene_priority_for_star():
    """v0.4.0 contract (mirror of the above): claim-level feature
    qualifiers (``expressed_in`` here) keep the default
    ``accession`` → ``gene`` → ``concept`` order, so the same name
    ``STAR`` resolves to ``HGNC:11359`` — the steroidogenic acute
    regulatory protein gene — because feature claims are dominated by
    gene mentions."""
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "STAR",
        "qualifier": "expressed_in",
        "object": "mesentery",
    })
    assert out["ontology_terms"]["subject"] == "HGNC:11359"


def test_missing_qualifier_falls_back_to_default_priority():
    """If a claim arrives without a qualifier (or with an unknown
    qualifier), the resolver chain defaults to
    ``accession`` → ``gene`` → ``concept`` — the safe choice for
    feature claims, which dominate the transcriptomics adapter's
    use cases."""
    a = BulkTranscriptomicsAdapter()
    # No qualifier at all.
    out = a.enrich_claim({"subject": "STAR"})
    assert out["ontology_terms"]["subject"] == "HGNC:11359"

    # Unknown qualifier.
    out = a.enrich_claim({
        "subject": "STAR",
        "qualifier": "some_qualifier_we_dont_recognize",
    })
    assert out["ontology_terms"]["subject"] == "HGNC:11359"


def test_study_level_priority_still_falls_through_to_gene():
    """When a study-level qualifier names a term that's *not* in the
    concept index and *not* an accession, the resolver chain falls
    through to gene-lookup — so claim shapes like
    ``aligned_with: <study> WNT4`` (nonsensical but possible from a
    sloppy extractor) still attach a CURIE rather than dropping the
    slot entirely. The principle: relaxed priority never *removes*
    resolution, only re-orders it."""
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "study_42",
        "qualifier": "aligned_with",
        "object": "WNT4",  # not in concept index, but is a gene
    })
    assert out["ontology_terms"]["object"] == "HGNC:12783"


def test_concept_lookup_attaches_curie_and_xrefs(monkeypatch):
    """Wiring contract for the concept-lookup branch in `enrich_claim`:
    monkey-patches a single-row index in to exercise the path in
    isolation from the bundled TSV and pin the output shape — principal
    CURIE under `ontology_terms[slot]`, label + ontology + description
    under `xrefs[slot]`."""
    from indicium_adapters_transcriptomics import _enrichment, bulk

    fake_row = {
        "term": "fold change",
        "curie": "NCIT:C60693",
        "label": "Fold Change",
        "description": "Ratio between two measurements.",
        "ontology": "NCIT",
        "source": "monkeypatch fixture",
    }
    fake_index = _enrichment.ConceptIndex([fake_row])
    monkeypatch.setattr(bulk, "CONCEPT_INDEX", fake_index)

    a = bulk.BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "fold change",
        "qualifier": "differentially_expressed_in",
        "object": "U0126",
    })
    assert out["ontology_terms"]["subject"] == "NCIT:C60693"
    xref = out["xrefs"]["subject"]
    assert xref["curie"] == "NCIT:C60693"
    assert xref["label"] == "Fold Change"
    assert xref["ontology"] == "NCIT"
    assert "description" in xref


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
