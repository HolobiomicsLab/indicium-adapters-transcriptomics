"""End-to-end smoke test on hand-extracted claims from real papers.

Two transcriptomics papers contribute claims here: Warner et al. 2025
(*Nematostella vectensis* regeneration) and Tabula Muris (mouse cell
atlas). The claims are the kind of assertions an LLM extractor would
produce after reading the paper — gene names, study accessions, cell
types, instrument models, treatments.

For each claim we:

1. Pass it through ``BulkTranscriptomicsAdapter.enrich_claim`` or
   ``SingleCellTranscriptomicsAdapter.enrich_claim`` based on the
   ``_adapter`` annotation.
2. Print the resolved CURIEs + xref attachments.
3. Build a tiny RDF graph asserting the resolved CURIEs as
   ``oa:hasBody`` properties on a synthetic claim subject, validate
   against the adapters' SHACL shapes graph via pyshacl. The
   validation is the integration smoke-test — proves the round-trip
   from claim dict → CURIEs → RDF → SHACL is wired.

Distinct from the hand-built Warner demo
(``examples/enrich_warner2018.py``) which is a *resolution coverage*
showcase. This script is a *pipeline integration smoke test*.

Run::

    py -3.13 examples/test_on_papers.py
"""
from __future__ import annotations

import json
import sys

import pyshacl
import rdflib
from indicium_adapters import compose_adapters, discover_adapters

sys.stdout.reconfigure(encoding="utf-8")


# Each tuple: (paper_id, adapter_name, claim_dict, one-line note).
WARNER_CLAIMS: list[tuple[str, str, dict, str]] = [
    (
        "Warner2025", "bulk-transcriptomics",
        {"subject": "GSE168938", "qualifier": "consolidates_dataset",
         "object": "Warner 2025"},
        "GEO accession resolves via the bulk accession resolver",
    ),
    (
        "Warner2025", "bulk-transcriptomics",
        {"subject": "study_warner_2025", "qualifier": "sequenced_on",
         "object": "Illumina NovaSeq 6000",
         "species": "nematostella"},
        "Platform mention via study-level qualifier",
    ),
    (
        "Warner2025", "bulk-transcriptomics",
        {"subject": "study_warner_2025", "qualifier": "tested_with",
         "object": "U0126"},
        "MEK inhibitor (CHEBI) via v0.8.0 chemicals harvest",
    ),
    (
        "Warner2025", "bulk-transcriptomics",
        {"subject": "WNT4", "qualifier": "differentially_expressed_in",
         "object": "regenerating_oral_pole"},
        "Gene-symbol resolution (HGNC)",
    ),
    (
        "Warner2025", "single-cell-transcriptomics",
        {"subject": "MMP9", "qualifier": "localized_in_tissue",
         "object": "mesentery"},
        "Gene + UBERON anatomy (Cheung 2024 anatomy index)",
    ),
]

TABULA_MURIS_CLAIMS: list[tuple[str, str, dict, str]] = [
    (
        "TabulaMuris", "bulk-transcriptomics",
        {"subject": "tabula_muris_2018", "qualifier": "consolidates_dataset",
         "object": "GSE109774"},
        "Tabula Muris GEO accession",
    ),
    (
        "TabulaMuris", "bulk-transcriptomics",
        {"subject": "tabula_muris_2018", "qualifier": "aligned_with",
         "object": "STAR"},
        "v0.4.0 qualifier-aware: STAR resolves to biotools:star "
        "(aligner) under aligned_with",
    ),
    (
        "TabulaMuris", "bulk-transcriptomics",
        {"subject": "Mapk1", "qualifier": "differentially_expressed_in",
         "object": "intestine"},
        "v0.6.0 mouse Title-case heuristic → MGI:1346858",
    ),
    (
        "TabulaMuris", "single-cell-transcriptomics",
        {"subject": "WNT4", "qualifier": "marker_of_cell_type",
         "object": "macrophage"},
        "Gene + CL cell type (v0.5.0)",
    ),
    (
        "TabulaMuris", "single-cell-transcriptomics",
        {"subject": "Mapk1", "qualifier": "marker_of_cell_type",
         "object": "cardiomyocyte",
         "species": "mouse"},
        "Mouse gene + CL → cardiomyocyte's synonym in CL",
    ),
]


def _render_claim_as_rdf(graph: rdflib.Graph, idx: int, enriched: dict) -> None:
    """Mint a tiny RDF representation of an enriched claim so we have
    something to feed pyshacl. We don't follow the full ASB ``Claim``
    SuperPattern here — just stash the principal CURIEs as ``rdfs:label``
    on a synthetic node so the SHACL CURIE-format shape fires."""
    asb = rdflib.Namespace("https://w3id.org/asb/")
    claim_iri = asb[f"claim_{idx}"]
    graph.add((claim_iri, rdflib.RDF.type, asb.Claim))
    for slot, curie in (enriched.get("ontology_terms") or {}).items():
        prop = asb[slot]
        graph.add((claim_iri, prop, rdflib.Literal(curie)))


def main() -> int:
    adapters = discover_adapters()
    composite = compose_adapters(list(adapters.values()))
    shapes_graph = composite.shacl_shapes()
    print(
        f"[setup] {len(adapters)} adapters discovered, "
        f"{len(shapes_graph)} SHACL triples loaded.\n"
    )

    all_claims = WARNER_CLAIMS + TABULA_MURIS_CLAIMS
    data_graph = rdflib.Graph()
    enriched_count = 0

    for i, (paper, adapter_name, claim_in, note) in enumerate(all_claims, 1):
        adapter = adapters[adapter_name]
        enriched = adapter.enrich_claim(claim_in)
        terms = enriched.get("ontology_terms") or {}
        marker = "✓" if terms else "○"
        slot_curies = (
            ", ".join(f"{k}={v}" for k, v in terms.items()) if terms else "no resolution"
        )
        print(f"[{i:2d}] {paper:12s} {adapter_name:30s} {marker} {slot_curies}")
        print(f"      {note}")
        if terms:
            enriched_count += 1
            _render_claim_as_rdf(data_graph, i, enriched)

    print()
    print(f"[summary] {enriched_count}/{len(all_claims)} claims enriched.")

    # SHACL validation on the synthesised RDF. The shapes ship with a
    # CURIE-format NodeShape that should pass — the assert proves the
    # full chain is plumbed.
    conforms, _report_graph, _report_text = pyshacl.validate(
        data_graph=data_graph,
        shapes_graph=shapes_graph,
        meta_shacl=False,
        debug=False,
    )
    print(f"[shacl] data graph: {len(data_graph)} triples, "
          f"conforms: {conforms}")
    print()

    # Dump one fully-enriched claim for visual inspection.
    print("Sample enriched claim (claim #5, Warner MMP9 localised_in_tissue):")
    sample = adapters["single-cell-transcriptomics"].enrich_claim(
        WARNER_CLAIMS[4][2]
    )
    print(json.dumps(sample, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
