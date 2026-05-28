"""End-to-end demonstration of the transcriptomics adapter on a hand-built
set of Warner2018-flavoured claims.

The 12 claims below mimic what an LLM extractor might emit after reading
Warner et al. 2025 (or the surrounding Nematostella regeneration
literature). They cover every resolution axis the adapter currently
supports, so a single run gives a visual smoke test of:

- gene-symbol resolution (HGNC + Nematostella cross-references)
- repository-accession resolution (GEO / SRA / ENA / ...)
- corpus-harvested concept CURIEs (OLS4 platforms / protocols / stats)
- bio.tools software CURIEs
- v0.4.0 qualifier-aware priority (STAR-as-aligner vs STAR-as-gene)
- single-cell UBERON anatomy resolution
- v0.4.0 single-cell gene resolution (new in v0.4.0)
- terms deliberately left unenriched (cnidarian anatomy, CL cell types)

No network calls. No LLM. Self-contained — useful as a confidence check
before integrating the adapter into the ASB pipeline.

Run::

    py -3.13 examples/enrich_warner2018.py
"""
from __future__ import annotations

import json
import sys
from typing import Final

from indicium_adapters_transcriptomics import (
    BulkTranscriptomicsAdapter,
    SingleCellTranscriptomicsAdapter,
)

sys.stdout.reconfigure(encoding="utf-8")


# (adapter_name, claim_dict, one-line rationale shown to the reader)
ClaimSpec = tuple[str, dict, str]


CLAIMS: Final[list[ClaimSpec]] = [
    # --- Bulk: gene resolution ----------------------------------------------
    (
        "bulk",
        {"subject": "WNT4", "qualifier": "differentially_expressed_in",
         "object": "regenerating_oral_pole"},
        "Gene symbol (subject) → HGNC + xref map",
    ),
    (
        "bulk",
        {"subject": "MAPK1", "qualifier": "regulated_by",
         "object": "MAP2K1"},
        "Two gene symbols → both resolve via HGNC",
    ),
    # --- Bulk: accession resolution -----------------------------------------
    (
        "bulk",
        {"subject": "study_warner_2025", "qualifier": "consolidates_dataset",
         "object": "GSE168938"},
        "GEO accession → GEO:GSE168938",
    ),
    # --- Bulk: corpus-harvested concept (platform) --------------------------
    (
        "bulk",
        {"subject": "study_warner_2025", "qualifier": "sequenced_on",
         "object": "Illumina"},
        "Platform name → NCIT:C215820 (v0.3.1 corpus harvest)",
    ),
    (
        "bulk",
        {"subject": "study_warner_2025", "qualifier": "sequenced_on",
         "object": "NovaSeq"},
        "Instrument series → EFO:0008637 Illumina NovaSeq 6000 (substring match)",
    ),
    # --- Bulk: qualifier-aware (the STAR collision) -------------------------
    (
        "bulk",
        {"subject": "study_warner_2025", "qualifier": "aligned_with",
         "object": "STAR"},
        "v0.4.0 STAR collision: study-level qualifier → biotools:star (aligner)",
    ),
    (
        "bulk",
        {"subject": "STAR", "qualifier": "expressed_in",
         "object": "regenerating_blastema"},
        "v0.4.0 STAR collision: feature qualifier → HGNC:11359 (gene)",
    ),
    # --- Bulk: bio.tools + concept index ------------------------------------
    (
        "bulk",
        {"subject": "study_warner_2025", "qualifier": "quantified_with",
         "object": "DESeq2"},
        "DE tool → biotools:deseq2",
    ),
    (
        "bulk",
        {"subject": "fold change", "qualifier": "tested_with",
         "object": "Wilcoxon"},
        "Static concept (subject) + corpus stat (object, prefix match)",
    ),
    # --- Single-cell: anatomy + new gene resolver (v0.4.0) -----------------
    (
        "singlecell",
        {"subject": "MMP9", "qualifier": "localized_in_tissue",
         "object": "mesentery"},
        "v0.4.0 single-cell: gene subject (new) + UBERON anatomy object",
    ),
    (
        "singlecell",
        {"subject": "WNT4", "qualifier": "marker_of_cell_type",
         "object": "cnidocyte"},
        "v0.4.0 single-cell: gene subject (new); cnidocyte (cell type) unresolved — CL roadmap",
    ),
    # --- Single-cell: deliberately unenriched ------------------------------
    (
        "singlecell",
        {"subject": "MMP14", "qualifier": "localized_in_spatial_domain",
         "object": "physa"},
        "Subject (gene) resolves; physa stays bare — cnidarian-specific, no UBERON",
    ),
]


def _format_dict(d: dict) -> str:
    return json.dumps(d, indent=2, ensure_ascii=False, sort_keys=True)


def _diff_keys(before: dict, after: dict) -> list[str]:
    """Top-level keys that are new in ``after``."""
    return sorted(set(after.keys()) - set(before.keys()))


def main() -> int:
    bulk = BulkTranscriptomicsAdapter()
    sc = SingleCellTranscriptomicsAdapter()

    adapters = {"bulk": bulk, "singlecell": sc}

    n_enriched = 0
    n_no_change = 0

    print("=" * 70)
    print(f"Transcriptomics adapter end-to-end demo — {len(CLAIMS)} claims")
    print("=" * 70)

    for i, (adapter_name, claim_in, rationale) in enumerate(CLAIMS, start=1):
        adapter = adapters[adapter_name]
        claim_out = adapter.enrich_claim(claim_in)

        new_keys = _diff_keys(claim_in, claim_out)
        added_terms = claim_out.get("ontology_terms") or {}
        was_enriched = bool(added_terms)

        if was_enriched:
            n_enriched += 1
        else:
            n_no_change += 1

        marker = "✓ enriched" if was_enriched else "○ no enrichment"
        print()
        print(f"[{i:2d}] {adapter_name:10s} {marker:18s} — {rationale}")
        print(f"    IN  : {_format_dict(claim_in)}")
        if was_enriched:
            curies = ", ".join(f"{slot}={c}" for slot, c in added_terms.items())
            print(f"    OUT : ontology_terms[{curies}]")
            xrefs = claim_out.get("xrefs") or {}
            for slot, xref in xrefs.items():
                # Compress xrefs to the principal id keys for readability.
                principal = (
                    xref.get("hgnc")
                    or xref.get("uberon")
                    or xref.get("curie")
                )
                label = xref.get("label", "")
                if label:
                    print(f"          xrefs[{slot}] : {principal} \"{label}\"")
                else:
                    print(f"          xrefs[{slot}] : {principal}")
        else:
            print(f"    OUT : (no resolution — slot values pass through)")
        print(f"    new keys on output dict: {new_keys}")

    print()
    print("=" * 70)
    print(f"Summary: {n_enriched}/{len(CLAIMS)} claims enriched, "
          f"{n_no_change}/{len(CLAIMS)} unchanged.")
    print("=" * 70)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
