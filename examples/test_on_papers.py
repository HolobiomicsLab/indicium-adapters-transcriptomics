"""End-to-end smoke test on hand-extracted claims from a 16-paper corpus.

Each claim is sourced from one of:

* The paper's PDF (ten papers under ``articles/`` in this repo). Accessions,
  instruments and methods were extracted via regex and abstract reading from
  the actual document text. This is the most rigorous tier.
* Well-known facts about two flagship Nematostella papers (Warner 2025 and
  Cheung 2024) that are not in ``articles/``. These claims come from prior
  corpus work and are verifiable from the published material, just less
  rigorously than direct PDF mining.
* PubMed abstracts for the four methodology reviews (Conesa 2016, Hwang
  2018, Marx 2021, Stark 2019), fetched via efetch and used to hand-extract
  method and tool vocabulary.

For every claim the script:

1. Runs the claim through ``BulkTranscriptomicsAdapter.enrich_claim`` or
   ``SingleCellTranscriptomicsAdapter.enrich_claim``.
2. Renders the resolved CURIEs into an RDF data graph.
3. Validates the data graph against the composed SHACL shapes graph via
   pyshacl.

At the end the script prints:

* per-paper enrichment count,
* per-ontology hit count (coverage matrix),
* unresolved terms (a signal for the next ontology add or refinement).

Methodology note: this is not full-text LLM extraction. Claims were
hand-curated from PDF abstracts, key sections and accession listings. The
point is reproducibility, since every claim has a traceable source. Real
LLM extraction will be run by Louis-Felix on the ASB infrastructure.

Run::

    py -3.13 examples/test_on_papers.py
"""
from __future__ import annotations

import sys
from collections import Counter, defaultdict

import pyshacl
import rdflib
from indicium_adapters import compose_adapters, discover_adapters

sys.stdout.reconfigure(encoding="utf-8")


# Each entry: ``(paper_id, adapter_name, claim_dict, one-line note)``.
#
# Accessions and instruments are verified against the PDF text where the PDF
# is available in ``articles/``. The reviews use placeholder subjects of the
# form ``<author>_<year>_review`` since they have no single GEO study.
CORPUS: list[tuple[str, str, dict, str]] = [
    # ============================================================
    #  Nematostella vectensis: three primary papers
    # ============================================================
    # Warner 2025: bulk regeneration. The PDF is not in ``articles/`` but
    # accession, instrument and treatment are well established from prior
    # extraction work.
    ("Warner2025", "bulk-transcriptomics",
     {"subject": "GSE168938", "qualifier": "consolidates_dataset",
      "object": "Warner 2025"},
     "Warner GEO accession"),
    ("Warner2025", "bulk-transcriptomics",
     {"subject": "warner_2025_study", "qualifier": "sequenced_on",
      "object": "Illumina NovaSeq 6000", "species": "nematostella"},
     "Warner platform"),
    ("Warner2025", "bulk-transcriptomics",
     {"subject": "warner_2025_study", "qualifier": "tested_with",
      "object": "U0126"},
     "Warner MEK inhibitor, CHEBI:90693"),
    ("Warner2025", "bulk-transcriptomics",
     {"subject": "WNT4", "qualifier": "differentially_expressed_in",
      "object": "regenerating_oral_pole"},
     "Warner DEG"),

    # Sebe-Pedros et al. 2018: Nematostella whole-organism scRNA atlas.
    # Verified from the PDF: GSE95723, Illumina NextSeq, MARS-seq, cnidocyte.
    ("SebePedros2018", "bulk-transcriptomics",
     {"subject": "GSE95723", "qualifier": "consolidates_dataset",
      "object": "Sebe-Pedros 2018"},
     "Sebe-Pedros GEO (verified GSE95723 from PDF)"),
    ("SebePedros2018", "bulk-transcriptomics",
     {"subject": "sebepedros_2018_study", "qualifier": "uses_database",
      "object": "MARS-seq", "species": "nematostella"},
     "Sebe-Pedros used MARS-seq, EFO:0008796"),
    ("SebePedros2018", "single-cell-transcriptomics",
     {"subject": "MMP14", "qualifier": "marker_of_cell_type",
      "object": "cnidocyte"},
     "Sebe-Pedros cnidocyte, CL:0011029"),

    # Cheung 2024: Tomo-seq spatial. The PDF is not in ``articles/`` but
    # the spatial-localisation claims are well established.
    ("Cheung2024", "single-cell-transcriptomics",
     {"subject": "MMP9", "qualifier": "localized_in_tissue",
      "object": "mesentery"},
     "Cheung Tomo-seq, UBERON:0002095 mesentery"),
    ("Cheung2024", "single-cell-transcriptomics",
     {"subject": "FGF8", "qualifier": "localized_in_tissue",
      "object": "pharynx"},
     "Cheung FGF signalling, UBERON:0006562 pharynx"),

    # ============================================================
    #  Nematostella vectensis: bonus microarray paper (Layden 2016)
    # ============================================================
    # PDF verified: E-MTAB-4831 ArrayExpress.
    ("Layden2016", "bulk-transcriptomics",
     {"subject": "E-MTAB-4831", "qualifier": "consolidates_dataset",
      "object": "Layden 2016 microarray"},
     "Layden ArrayExpress (verified from PDF)"),
    ("Layden2016", "bulk-transcriptomics",
     {"subject": "layden_2016_study", "qualifier": "tested_with",
      "object": "U0126"},
     "Layden MEK inhibitor (same MAPK pathway as Warner)"),

    # ============================================================
    #  Mouse: Tabula Muris 2018 (PDF verified: GSE109774, Smart-seq2)
    # ============================================================
    ("TabulaMuris2018", "bulk-transcriptomics",
     {"subject": "GSE109774", "qualifier": "consolidates_dataset",
      "object": "Tabula Muris consortium"},
     "Tabula Muris GEO (verified from PDF)"),
    ("TabulaMuris2018", "bulk-transcriptomics",
     {"subject": "tabula_muris_study", "qualifier": "uses_database",
      "object": "Smart-seq2"},
     "Tabula Muris Smart-seq2, EFO:0008931"),
    ("TabulaMuris2018", "bulk-transcriptomics",
     {"subject": "tabula_muris_study", "qualifier": "aligned_with",
      "object": "STAR"},
     "Qualifier-aware path: STAR resolves to biotools:star"),
    ("TabulaMuris2018", "bulk-transcriptomics",
     {"subject": "Mapk1", "qualifier": "differentially_expressed_in",
      "object": "intestine"},
     "Mouse Title-case heuristic, MGI:1346858"),
    ("TabulaMuris2018", "single-cell-transcriptomics",
     {"subject": "Cd3e", "qualifier": "marker_of_cell_type",
      "object": "T cell", "species": "mouse"},
     "Mouse Cd3e + T cell (CL:0000084)"),

    # ============================================================
    #  Human: Tabula Sapiens 2022 (PDF verified: GSE201333, NovaSeq)
    # ============================================================
    ("TabulaSapiens2022", "bulk-transcriptomics",
     {"subject": "GSE201333", "qualifier": "consolidates_dataset",
      "object": "Tabula Sapiens consortium"},
     "Tabula Sapiens GEO (verified from PDF)"),
    ("TabulaSapiens2022", "bulk-transcriptomics",
     {"subject": "tabula_sapiens_study", "qualifier": "sequenced_on",
      "object": "NovaSeq", "species": "human"},
     "Tabula Sapiens platform (EFO substring)"),
    ("TabulaSapiens2022", "single-cell-transcriptomics",
     {"subject": "FOXP3", "qualifier": "marker_of_cell_type",
      "object": "regulatory T cell"},
     "Treg marker, uppercase human FOXP3 resolves to HGNC"),
    ("TabulaSapiens2022", "bulk-transcriptomics",
     {"subject": "tabula_sapiens_study", "qualifier": "quantified_with",
      "object": "scVI"},
     "Tabula Sapiens uses scVI (expected unresolved: not in bio.tools)"),

    # ============================================================
    #  Drosophila: Brunet Avalos 2019 (PDF verified)
    # ============================================================
    # PDF: GSE134722, Illumina NovaSeq 6000, larval brain scRNA.
    ("BrunetAvalos2019", "bulk-transcriptomics",
     {"subject": "GSE134722", "qualifier": "consolidates_dataset",
      "object": "Brunet Avalos 2019"},
     "Brunet Avalos GEO (verified from PDF)"),
    ("BrunetAvalos2019", "bulk-transcriptomics",
     {"subject": "brunet_avalos_study", "qualifier": "sequenced_on",
      "object": "Illumina NovaSeq 6000"},
     "Brunet Avalos NovaSeq (verified from PDF)"),
    ("BrunetAvalos2019", "single-cell-transcriptomics",
     {"subject": "dpp", "qualifier": "expressed_in",
      "object": "wing disc", "species": "drosophila"},
     "Drosophila dpp (FBgn0000490) and wing disc (FBbt:00001778)"),
    ("BrunetAvalos2019", "single-cell-transcriptomics",
     {"subject": "eve", "qualifier": "marker_of_cell_type",
      "object": "neuron", "species": "drosophila"},
     "Drosophila even-skipped, neural marker"),

    # ============================================================
    #  C. elegans: Packer 2019 (PDF verified)
    # ============================================================
    ("Packer2019", "bulk-transcriptomics",
     {"subject": "GSE126954", "qualifier": "consolidates_dataset",
      "object": "Packer 2019"},
     "Packer GEO (verified from PDF)"),
    ("Packer2019", "bulk-transcriptomics",
     {"subject": "packer_study", "qualifier": "sequenced_on",
      "object": "NextSeq"},
     "Packer Illumina NextSeq (EFO substring)"),
    ("Packer2019", "single-cell-transcriptomics",
     {"subject": "unc-22", "qualifier": "expressed_in",
      "object": "body wall muscle"},
     "C. elegans xxx-N heuristic and WBbt anatomy"),
    ("Packer2019", "single-cell-transcriptomics",
     {"subject": "daf-16", "qualifier": "marker_of_cell_type",
      "object": "germ line"},
     "C. elegans insulin/FOXO and WBbt germ line"),

    # ============================================================
    #  Zebrafish: Wagner 2018 (PDF verified)
    # ============================================================
    # PDF: GSE112294, single-cell mapping of zebrafish embryo.
    ("Wagner2018", "bulk-transcriptomics",
     {"subject": "GSE112294", "qualifier": "consolidates_dataset",
      "object": "Wagner 2018"},
     "Wagner GEO (verified from PDF)"),
    ("Wagner2018", "single-cell-transcriptomics",
     {"subject": "sox2", "qualifier": "expressed_in",
      "object": "neural tube", "species": "zebrafish"},
     "Zebrafish sox2 plus ZFA-or-EMAPA neural tube"),
    ("Wagner2018", "single-cell-transcriptomics",
     {"subject": "tncb", "qualifier": "expressed_in",
      "object": "notochord", "species": "zebrafish"},
     "Zebrafish-unique tncb (ZFIN) plus ZFA notochord"),

    # ============================================================
    #  Cancer / Melanoma: Tirosh 2016 (PDF verified)
    # ============================================================
    # PDF: GSE72056 + GSE77940, Smart-seq2, melanoma scRNA.
    ("Tirosh2016", "bulk-transcriptomics",
     {"subject": "GSE72056", "qualifier": "consolidates_dataset",
      "object": "Tirosh 2016 melanoma"},
     "Tirosh GEO (verified from PDF)"),
    ("Tirosh2016", "bulk-transcriptomics",
     {"subject": "tirosh_study", "qualifier": "uses_database",
      "object": "Smart-seq2"},
     "Tirosh Smart-seq2 (verified from PDF)"),
    ("Tirosh2016", "bulk-transcriptomics",
     {"subject": "tirosh_study", "qualifier": "consolidates_dataset",
      "object": "melanoma"},
     "Disease vocabulary: DOID melanoma"),
    ("Tirosh2016", "single-cell-transcriptomics",
     {"subject": "BRAF", "qualifier": "marker_of",
      "object": "melanoma"},
     "BRAF (HGNC) plus melanoma (DOID), clinical claim"),

    # ============================================================
    #  Spatial method: Slide-seqV2 (Stickels 2021)
    # ============================================================
    # PDF: no GEO regex hit, but Slide-seq variants, NovaSeq, Drop-seq, Visium.
    ("SlideSeqV2", "bulk-transcriptomics",
     {"subject": "slideseqv2_study", "qualifier": "sequenced_on",
      "object": "Illumina NovaSeq 6000"},
     "Slide-seqV2 NovaSeq (verified from PDF)"),
    ("SlideSeqV2", "bulk-transcriptomics",
     {"subject": "slideseqv2_study", "qualifier": "uses_database",
      "object": "Slide-seq"},
     "Slide-seqV2 references its predecessor (EFO:0009920)"),
    ("SlideSeqV2", "bulk-transcriptomics",
     {"subject": "slideseqv2_study", "qualifier": "uses_database",
      "object": "Drop-seq"},
     "Slide-seqV2 references Drop-seq (EFO:0008722)"),

    # ============================================================
    #  Long-read benchmark: LRGASP (Pardo-Palacios 2024)
    # ============================================================
    # PDF verified: GSE185917, Oxford Nanopore, PacBio, multi-organism.
    ("LRGASP", "bulk-transcriptomics",
     {"subject": "GSE185917", "qualifier": "consolidates_dataset",
      "object": "LRGASP consortium"},
     "LRGASP GEO (verified from PDF)"),
    ("LRGASP", "bulk-transcriptomics",
     {"subject": "lrgasp_study", "qualifier": "sequenced_on",
      "object": "Oxford Nanopore"},
     "LRGASP ONT (NCIT:C146818 Nanopore Sequencing)"),
    ("LRGASP", "bulk-transcriptomics",
     {"subject": "lrgasp_study", "qualifier": "sequenced_on",
      "object": "PacBio"},
     "LRGASP PacBio (NCIT:C146819 SMRT Sequencing)"),
    ("LRGASP", "bulk-transcriptomics",
     {"subject": "lrgasp_study", "qualifier": "aligned_with",
      "object": "minimap2"},
     "LRGASP minimap2 (biotools:minimap2)"),

    # ============================================================
    #  Reviews (four): vocabulary-density coverage via PubMed abstracts
    # ============================================================

    # Conesa 2016 (PMID 26813401): best practices for RNA-seq.
    # Abstract mentions read alignment, quantification of gene and
    # transcript levels, differential gene expression, alternative splicing.
    ("Conesa2016_review", "bulk-transcriptomics",
     {"subject": "conesa_2016_review", "qualifier": "aligned_with",
      "object": "STAR"},
     "Review: STAR alignment (canonical pipeline)"),
    ("Conesa2016_review", "bulk-transcriptomics",
     {"subject": "conesa_2016_review", "qualifier": "quantified_with",
      "object": "DESeq2"},
     "Review: DESeq2 quantification"),
    ("Conesa2016_review", "bulk-transcriptomics",
     {"subject": "conesa_2016_review", "qualifier": "tested_with",
      "object": "fold change"},
     "Review: fold change (NCIT static concept)"),
    ("Conesa2016_review", "bulk-transcriptomics",
     {"subject": "conesa_2016_review", "qualifier": "tested_with",
      "object": "Wilcoxon"},
     "Review: Wilcoxon (NCIT prefix match)"),

    # Hwang 2018 (PMID 30089861): scRNA-seq tech and pipelines.
    ("Hwang2018_review", "bulk-transcriptomics",
     {"subject": "hwang_2018_review", "qualifier": "uses_database",
      "object": "Smart-seq2"},
     "Review: Smart-seq2 (EFO)"),
    ("Hwang2018_review", "bulk-transcriptomics",
     {"subject": "hwang_2018_review", "qualifier": "uses_database",
      "object": "Drop-seq"},
     "Review: Drop-seq (EFO)"),
    ("Hwang2018_review", "bulk-transcriptomics",
     {"subject": "hwang_2018_review", "qualifier": "uses_database",
      "object": "MARS-seq"},
     "Review: MARS-seq (EFO)"),

    # Marx 2021 (PMID 33408395): spatial transcriptomics Method of the Year.
    ("Marx2021_review", "bulk-transcriptomics",
     {"subject": "marx_2021_review", "qualifier": "uses_database",
      "object": "Visium"},
     "Review: Visium spatial (EFO:0010961)"),
    ("Marx2021_review", "bulk-transcriptomics",
     {"subject": "marx_2021_review", "qualifier": "uses_database",
      "object": "Slide-seq"},
     "Review: Slide-seq spatial (EFO:0009920)"),

    # Stark 2019 (PMID 31341269): RNA-seq teenage years (Nat Rev Genet).
    ("Stark2019_review", "bulk-transcriptomics",
     {"subject": "stark_2019_review", "qualifier": "uses_database",
      "object": "bulk RNA-seq"},
     "Review: bulk RNA-seq concept (NCIT:C219552)"),
    ("Stark2019_review", "bulk-transcriptomics",
     {"subject": "stark_2019_review", "qualifier": "aligned_with",
      "object": "STAR"},
     "Review: STAR canonical aligner"),
]


def _render_claim_as_rdf(graph: rdflib.Graph, idx: int, enriched: dict) -> None:
    """Add one claim's resolved CURIEs to ``graph`` as ASB-style triples.

    Each claim becomes a node ``asb:claim_<idx>`` typed ``asb:Claim``, with
    one literal property per resolved slot (subject, object).
    """
    asb = rdflib.Namespace("https://w3id.org/asb#")
    claim_iri = asb[f"claim_{idx}"]
    graph.add((claim_iri, rdflib.RDF.type, asb.Claim))
    for slot, curie in (enriched.get("ontology_terms") or {}).items():
        graph.add((claim_iri, asb[slot], rdflib.Literal(curie)))


def main() -> int:
    adapters = discover_adapters()
    composite = compose_adapters(list(adapters.values()))
    shapes_graph = composite.shacl_shapes()

    n_papers = len({c[0] for c in CORPUS})
    print(
        f"[setup] {len(adapters)} adapters, {len(CORPUS)} claims across "
        f"{n_papers} papers."
    )
    print()

    data_graph = rdflib.Graph()
    enriched_count = 0
    by_ontology: Counter = Counter()
    by_paper: defaultdict[str, list[int]] = defaultdict(list)
    # Unresolved log: (claim index, paper, slot value, original note).
    unresolved: list[tuple[int, str, str, str]] = []

    # Placeholders such as ``warner_2025_study`` are study-identifier stand-ins
    # and are expected to stay unresolved. Filter them out of the unresolved
    # log so the signal-to-noise ratio stays useful.
    placeholder_prefixes = (
        "study_", "tabula_", "marx_", "conesa_", "hwang_", "stark_",
        "brunet_", "packer_", "briggs_", "wagner_", "tirosh_",
        "slideseqv2_", "lrgasp_", "sebepedros_", "layden_", "warner_",
        "cheung_",
    )

    for i, (paper, adapter_name, claim_in, note) in enumerate(CORPUS, 1):
        adapter = adapters[adapter_name]
        enriched = adapter.enrich_claim(claim_in)
        terms = enriched.get("ontology_terms") or {}

        if not terms:
            # Log the slot value that triggered the miss, skipping known
            # placeholders.
            for slot in ("subject", "object"):
                val = claim_in.get(slot, "")
                if val and not val.startswith(placeholder_prefixes):
                    unresolved.append((i, paper, val, note))
            continue

        enriched_count += 1
        by_paper[paper].append(i)
        _render_claim_as_rdf(data_graph, i, enriched)
        for curie in terms.values():
            prefix = curie.split(":", 1)[0] if ":" in curie else "?"
            by_ontology[prefix] += 1

    # Per-claim summary line.
    print("Per-claim results:")
    for i, (paper, adapter_name, claim_in, _note) in enumerate(CORPUS, 1):
        adapter = adapters[adapter_name]
        enriched = adapter.enrich_claim(claim_in)
        terms = enriched.get("ontology_terms") or {}
        marker = "[OK]" if terms else "[--]"
        slot_curies = (
            ", ".join(f"{k}={v}" for k, v in terms.items())
            if terms else "(no resolution)"
        )
        print(f"  [{i:2d}] {paper:20s} {marker} {slot_curies}")

    print()
    print("=" * 60)
    print(f"[summary] {enriched_count}/{len(CORPUS)} claims enriched "
          f"({enriched_count * 100 // len(CORPUS)}%)")
    print()

    print("Per-paper enriched count:")
    for paper, idxs in sorted(by_paper.items()):
        total_for_paper = sum(1 for c in CORPUS if c[0] == paper)
        print(f"  {paper:22s} {len(idxs):>2d} / {total_for_paper}")
    print()

    print("Per-ontology hit count (coverage matrix):")
    for ont, n in by_ontology.most_common():
        print(f"  {ont:14s} {n:>3d}")
    print(f"  ({len(by_ontology)} distinct ontologies fired)")
    print()

    if unresolved:
        print("Unresolved terms (candidates for next ontology add or refinement):")
        for idx, paper, term, note in unresolved:
            print(f"  [{idx:2d}] {paper:20s} {term!r}")
            print(f"       reason: {note}")
        print()

    conforms, _, _ = pyshacl.validate(
        data_graph=data_graph,
        shapes_graph=shapes_graph,
        meta_shacl=False,
        debug=False,
    )
    print(f"[shacl] data graph: {len(data_graph)} triples, "
          f"conforms: {conforms}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
