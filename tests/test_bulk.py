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


# --- mouse (MGI) gene resolution (v0.6.0) -----------------------------------


def test_mouse_symbol_resolves_to_mgi_via_casing_heuristic():
    """v0.6.0: a Title-cased gene symbol (``Mapk1``) triggers the
    casing heuristic — the gene resolver picks the mouse (MGI) row,
    not the human (HGNC) row, even though both exist under the same
    case-insensitive key in the index."""
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "Mapk1",
        "qualifier": "differentially_expressed_in",
        "object": "regenerating_tissue",
    })
    assert out["ontology_terms"]["subject"] == "MGI:1346858"
    assert out["xrefs"]["subject"].get("mgi") == "MGI:1346858"


def test_human_uppercase_symbol_still_resolves_to_hgnc():
    """Mirror test: an UPPER-cased symbol picks human via the casing
    heuristic. Verifies the heuristic is two-sided and that adding
    mouse to the index didn't break human resolution for the same
    symbol family."""
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "MAPK1",
        "qualifier": "differentially_expressed_in",
        "object": "regenerating_tissue",
    })
    assert out["ontology_terms"]["subject"] == "HGNC:6871"


def test_species_hint_overrides_casing_heuristic():
    """v0.6.0: an explicit ``species`` field on the claim overrides
    the casing heuristic. A weirdly-cased symbol ``MAPK1`` with
    ``species="mouse"`` resolves to the MGI row — useful when an LLM
    extractor normalises casing but knows the paper is about mouse."""
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "MAPK1",
        "qualifier": "differentially_expressed_in",
        "object": "regenerating_tissue",
        "species": "mouse",
    })
    assert out["ontology_terms"]["subject"] == "MGI:1346858"


def test_mgi_accession_resolves_directly():
    """An MGI accession (``MGI:1346858``) in a claim slot resolves
    via the gene resolver's accession-lookup path, bypassing the
    symbol-disambiguation logic entirely. Useful for callers who
    want unambiguous mouse-gene references."""
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "MGI:1346858",
        "qualifier": "differentially_expressed_in",
        "object": "regenerating_tissue",
    })
    assert out["ontology_terms"]["subject"] == "MGI:1346858"
    # MGI value is in the principal CURIE; symbol comes from the row.
    assert out["xrefs"]["subject"].get("mgi") == "MGI:1346858"


def test_mouse_row_carries_full_xref_map():
    """v0.6.1: mouse rows now ship Ensembl / NCBI Gene / UniProt
    cross-references in addition to the MGI accession (parsed from
    MRK_ENSEMBL.rpt, MGI_EntrezGene.rpt, MRK_SwissProt_TrEMBL.rpt).
    Pin Mapk1's full xref map so a regression in any of the three
    MGI report parsers is caught."""
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "Mapk1",
        "qualifier": "differentially_expressed_in",
        "object": "regenerating_tissue",
    })
    xrefs = out["xrefs"]["subject"]
    assert xrefs.get("mgi") == "MGI:1346858"
    assert xrefs.get("ensembl") == "ENSMUSG00000063358"
    assert xrefs.get("ncbigene") == "26413"
    # UniProt is the SwissProt-or-TrEMBL accession from MGI's report;
    # the canonical SwissProt one for mouse Mapk1 is P63085 but MGI
    # may list a TrEMBL accession first — accept any non-empty value.
    assert xrefs.get("uniprot")


# --- v0.6.2: Drosophila / C. elegans / Zebrafish ---------------------------


def test_drosophila_unique_symbol_resolves_to_flybase():
    """v0.6.2: a drosophila-unique gene name (``dpp`` = decapentaplegic,
    no overlap with other organisms) resolves to its FlyBase CURIE
    without needing a species hint — the single-candidate shortcut in
    GeneIndex.lookup picks the only row.

    v0.6.3: the row's xrefs now carry the UniProt accession too
    (P07713) thanks to the UniProt per-organism idmapping join.
    """
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "dpp",
        "qualifier": "expressed_in",
        "object": "wing_disc",
    })
    assert out["ontology_terms"]["subject"] == "FlyBase:FBgn0000490"
    xrefs = out["xrefs"]["subject"]
    assert xrefs.get("flybase") == "FBgn0000490"
    assert xrefs.get("uniprot") == "P07713"  # v0.6.3 UniProt join


def test_drosophila_via_species_hint_disambiguates_collision():
    """``Egfr`` is Title-cased so the heuristic picks mouse; a
    species="drosophila" hint overrides and routes to FlyBase."""
    a = BulkTranscriptomicsAdapter()
    # Without hint — Title-case → mouse.
    out = a.enrich_claim({
        "subject": "Egfr",
        "qualifier": "expressed_in",
        "object": "wing_disc",
    })
    assert out["ontology_terms"]["subject"] == "MGI:95294"

    # With hint — overrides to drosophila.
    out = a.enrich_claim({
        "subject": "Egfr",
        "qualifier": "expressed_in",
        "object": "wing_disc",
        "species": "drosophila",
    })
    assert out["ontology_terms"]["subject"] == "FlyBase:FBgn0003731"


def test_celegans_xxx_dash_n_pattern_triggers_heuristic():
    """v0.6.2 casing heuristic recognises the C. elegans symbol
    pattern (3-4 lowercase letters + hyphen + alphanumeric) — most
    distinctive convention across the supported organisms.
    ``unc-22`` → WormBase:WBGene00006759.

    v0.6.3: the xrefs map now also carries the UniProt accession
    Q23551 (twitchin) thanks to the UniProt per-organism join.
    """
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "unc-22",
        "qualifier": "expressed_in",
        "object": "muscle",
    })
    assert out["ontology_terms"]["subject"] == "WormBase:WBGene00006759"
    assert out["xrefs"]["subject"].get("uniprot") == "Q23551"  # v0.6.3


def test_zebrafish_unique_symbol_resolves_via_single_candidate():
    """``tncb`` (tenascin Cb) is zebrafish-unique in the bundled
    index. Single-candidate lookup short-circuits the heuristic and
    returns the only row."""
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "tncb",
        "qualifier": "expressed_in",
        "object": "myoseptum",
    })
    assert out["ontology_terms"]["subject"] == "ZFIN:ZDB-GENE-980526-104"


def test_zebrafish_via_species_hint_for_collision():
    """``nanog`` exists as a Nematostella-anchored row (with human
    HGNC ortholog), a mouse Title-case row, and a zebrafish row.
    Without a hint the lowercase casing is ambiguous so the
    load-order fallback picks human (cross-listed via the
    Nematostella row's hgnc ortholog). With ``species="zebrafish"``
    the resolver routes to ZFIN."""
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "nanog",
        "qualifier": "marker_of",
        "object": "embryonic_stem_cell",
        "species": "zebrafish",
    })
    assert out["ontology_terms"]["subject"] == "ZFIN:ZDB-GENE-030131-5486"


def test_organism_accessions_resolve_directly():
    """FlyBase / WormBase / ZFIN accessions in claim slots resolve
    via the direct-accession lookup paths, bypassing the symbol
    disambiguation entirely. The lookup accepts both the bare id
    (``FBgn0000490``) and the prefixed form (``FlyBase:FBgn0000490``)."""
    a = BulkTranscriptomicsAdapter()
    for token, expected in (
        ("FBgn0000490",         "FlyBase:FBgn0000490"),
        ("FlyBase:FBgn0000490", "FlyBase:FBgn0000490"),
        ("WBGene00006759",      "WormBase:WBGene00006759"),
        ("ZDB-GENE-980526-104", "ZFIN:ZDB-GENE-980526-104"),
    ):
        out = a.enrich_claim({
            "subject": token,
            "qualifier": "expressed_in",
            "object": "x",
        })
        assert out["ontology_terms"]["subject"] == expected, (
            f"{token!r}: expected {expected}, got {out['ontology_terms']['subject']}"
        )


def test_mouse_only_symbol_resolves_without_collision():
    """A mouse-specific gene with no human equivalent of the same
    name (the gene ``Xist`` exists in both species but its mouse
    capitalisation is ``Xist``; the human one is ``XIST`` so the
    collision is benign). Pick something clearer: ``Foxp3``
    (uppercase ``FOXP3`` is the human ortholog with HGNC; the
    Title-cased ``Foxp3`` is the mouse Treg marker)."""
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "Foxp3",
        "qualifier": "marker_of",
        "object": "regulatory_T_cell",
    })
    assert out["ontology_terms"]["subject"].startswith("MGI:")


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
    a CURIE.

    Note (v0.8.0): pre-v0.8.0 this test used ``U0126`` as the object,
    but U0126 (the MEK inhibitor used in Warner 2025) now resolves to
    ``CHEBI:90693`` via the new ChEBI index — which is the v0.8.0
    feature working correctly. Both slots now need to be unresolvable
    terms."""
    a = BulkTranscriptomicsAdapter()
    out = a.enrich_claim({
        "subject": "completely-bogus-term-zzz",
        "qualifier": "differentially_expressed_in",
        "object": "another-bogus-term-qqq",
    })
    assert "ontology_terms" not in out


# --- corpus concept index (v0.3.1) ------------------------------------------


def test_corpus_concept_terms_resolve_via_concept_index():
    """v0.3.1 ships a corpus-driven concept index harvested from the
    6-paper mini-corpus. 17 terms resolved via OLS4 across three match
    modes (exact / prefix / substring). Pin each one so a regression
    in either the harvester or the dual-TSV loader is caught
    immediately.

    v0.5.0 also fixed an OLS4-API field-name bug: the search endpoint
    exposes synonyms under ``exact_synonyms`` rather than ``synonym``.
    Several v0.3.1 resolutions silently fell through to prefix or
    substring matching before the fix; with the fix those terms now
    take the more canonical synonym-exact match instead. The pins
    below reflect the v0.5.0 (corrected) results:

    - ``pacbio`` → NCIT:C146819 SMRT Sequencing (the technology,
      previously EFO:0920006 PacBio Revio via prefix).
    - ``oxford nanopore`` → NCIT:C146818 Nanopore Sequencing
      (previously OBI:0002750 ONT MinION via prefix).
    - ``bulk rna-seq`` → NCIT:C219552 Bulk RNA Sequencing
      (previously OBI:0003090 bulk RNA-seq assay via prefix).
    - ``10x genomics`` → EFO:0008995 10x technology (new resolution;
      unresolved in v0.3.1 because the prefix-match required the
      label to start with "10x Genomics" but EFO uses "10x" alone).
    """
    expected = {
        # term              (curie,           ontology, axis,         mode)
        # --- exact-match (label or synonym == corpus term) ---
        "illumina":         ("NCIT:C215820",  "NCIT",   "platform",   "exact"),
        "pacbio":           ("NCIT:C146819",  "NCIT",   "platform",   "exact"),
        "oxford nanopore":  ("NCIT:C146818",  "NCIT",   "platform",   "exact"),
        "mars-seq":         ("EFO:0008796",   "EFO",    "protocol",   "exact"),
        "smart-seq2":       ("EFO:0008931",   "EFO",    "protocol",   "exact"),
        "drop-seq":         ("EFO:0008722",   "EFO",    "protocol",   "exact"),
        "slide-seq":        ("EFO:0009920",   "EFO",    "protocol",   "exact"),
        "bulk rna-seq":     ("NCIT:C219552",  "NCIT",   "protocol",   "exact"),
        "10x genomics":     ("EFO:0008995",   "EFO",    "protocol",   "exact"),
        "visium":           ("EFO:0010961",   "EFO",    "protocol",   "exact"),
        "t-test":           ("NCIT:C53231",   "NCIT",   "stat_test",  "exact"),
        # --- prefix-match (corpus term = whole-word prefix of label) ---
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


def test_concept_index_carries_all_11_concept_sources():
    """The quad-TSV loader of v0.5.0 expanded to **eleven** sources at
    v0.8.0 — the original four (static / OLS4 corpus / bio.tools / CL)
    plus seven ontology-specific indexes (CHEBI / HP / MP / DOID / PR /
    OBA / ENVO) harvested from the OLS4 audit. Lower bound only, so
    adding rows in any file does not break the test."""
    n = len(CONCEPT_INDEX)
    # v0.3.0 31 + v0.3.1 17 corpus + v0.3.1 9 biotools + v0.5.0 35 CL
    # + v0.8.0: CHEBI 36 + HP 30 + MP 19 + DOID 29 + PR 18 + OBA 7
    # + ENVO 15 = 246 rows minimum.
    assert n >= 31 + 17 + 9 + 35 + 36 + 30 + 19 + 29 + 18 + 7 + 15, (
        f"expected ≥246 concept rows, got {n}"
    )


# --- v0.8.0 ontology-specific harvested indexes -----------------------------


def test_v080_chebi_chemicals_resolve():
    """v0.8.0 ChEBI harvest pins. The U0126 hit is the Warner 2025
    MEK inhibitor — this is the resolution v0.8.0 was *for*."""
    expected = {
        "u0126":         "CHEBI:90693",
        "cisplatin":     "CHEBI:27899",
        "paclitaxel":    "CHEBI:45863",
        "estradiol":     "CHEBI:23965",
        "wortmannin":    "CHEBI:52289",
    }
    for term, curie in expected.items():
        row = CONCEPT_INDEX.lookup(term)
        assert row is not None, f"CHEBI term {term!r} should resolve"
        assert row["curie"] == curie
        assert row["ontology"] == "CHEBI"


def test_v080_hp_human_phenotype_resolve():
    """v0.8.0 Human Phenotype harvest pins. HPO is the standard
    clinical phenotype vocabulary."""
    expected = {
        "intellectual disability": "HP:0001249",
        "microcephaly":            "HP:0000252",
        "seizure":                 "HP:0001250",
    }
    for term, curie in expected.items():
        row = CONCEPT_INDEX.lookup(term)
        assert row is not None, f"HP term {term!r} should resolve"
        assert row["curie"] == curie
        assert row["ontology"] == "HP"


def test_v080_mp_mouse_phenotype_resolve():
    """v0.8.0 Mammalian Phenotype harvest pins (mouse studies)."""
    expected = {
        "abnormal body weight":   "MP:0001259",
        "decreased body weight":  "MP:0001262",
    }
    for term, curie in expected.items():
        row = CONCEPT_INDEX.lookup(term)
        assert row is not None, f"MP term {term!r} should resolve"
        assert row["curie"] == curie
        assert row["ontology"] == "MP"


def test_v080_doid_disease_resolve():
    """v0.8.0 Human Disease Ontology pins (clinical / cancer context)."""
    expected = {
        "cancer":         "DOID:162",
        "breast cancer":  "DOID:1612",
        "lung cancer":    "DOID:1324",
    }
    for term, curie in expected.items():
        row = CONCEPT_INDEX.lookup(term)
        assert row is not None, f"DOID term {term!r} should resolve"
        assert row["curie"] == curie
        assert row["ontology"] == "DOID"


def test_v080_pr_protein_resolve():
    """v0.8.0 Protein Ontology pins. Note: PR is heavily nested —
    generic candidate names land on specific sub-class entries via
    synonym matching (``transcription factor`` → ``transcription
    factor Sp``). The CURIEs are valid PR entries but not always the
    most canonical; v0.8.x can refine candidates."""
    expected = {
        "alpha-synuclein":   "PR:000015322",
        "tubulin":           "PR:000028799",
        "ubiquitin":         "PR:000035728",
    }
    for term, curie in expected.items():
        row = CONCEPT_INDEX.lookup(term)
        assert row is not None, f"PR term {term!r} should resolve"
        assert row["curie"] == curie
        assert row["ontology"] == "PR"


def test_v080_oba_attribute_resolve():
    """v0.8.0 Ontology of Biological Attributes pins. OBA imports
    terms from PATO and GO; the ``curie`` may use those prefixes
    while the ``ontology`` column reflects the harvesting source
    (OBA) as a provenance hint."""
    row = CONCEPT_INDEX.lookup("heart rate")
    assert row is not None
    assert row["curie"] == "OBA:1001087"
    assert row["ontology"] == "OBA"

    # PATO-via-OBA: harvested via the OBA index but the canonical
    # term is in PATO.
    row = CONCEPT_INDEX.lookup("fertility")
    assert row is not None
    assert row["curie"] == "PATO:0000274"


def test_v080_envo_environment_resolve():
    """v0.8.0 Environment Ontology pins — holobiont / sample-context."""
    expected = {
        "marine biome":     "ENVO:00000447",
        "freshwater biome": "ENVO:00000873",
        "coral reef":       "ENVO:00000150",
    }
    for term, curie in expected.items():
        row = CONCEPT_INDEX.lookup(term)
        assert row is not None, f"ENVO term {term!r} should resolve"
        assert row["curie"] == curie
        assert row["ontology"] == "ENVO"


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
    # v0.4.0: bulk no longer imports CONCEPT_INDEX directly — it uses the
    # shared resolve_concept helper from _enrichment, which reads the
    # module-level singleton. So the patch target is _enrichment.
    monkeypatch.setattr(_enrichment, "CONCEPT_INDEX", fake_index)

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
