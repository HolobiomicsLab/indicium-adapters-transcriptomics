"""End-to-end demo of the transcriptomics adapter on a hand-built corpus
of Warner2018-flavoured claims.

The claim list below mimics what an LLM extractor might emit after reading
Warner et al. 2025 (or the surrounding Nematostella regeneration
literature) plus a few cross-species expansions. The set covers every
resolution axis the adapter currently supports, so a single run gives a
visual smoke test of:

* gene-symbol resolution (HGNC plus Nematostella cross-references)
* repository-accession resolution (GEO, SRA, ENA, ArrayExpress, ...)
* corpus-harvested concept CURIEs (OLS4 platforms, protocols, stats)
* bio.tools software CURIEs
* qualifier-aware priority (STAR-as-aligner vs STAR-as-gene)
* single-cell UBERON anatomy resolution
* single-cell gene resolution
* multi-organism resolution (mouse, drosophila, C. elegans, zebrafish)
* per-ontology indexes (ChEBI, HP, DOID, EMAPA, FBbt, MS, MOD)
* terms intentionally left unenriched (cnidarian anatomy without a
  UBERON equivalent)

No network calls, no LLM, self-contained. Useful as a confidence check
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

# Force UTF-8 stdout for Windows consoles.
sys.stdout.reconfigure(encoding="utf-8")


# Each demo entry is ``(adapter_name, claim_dict, one-line rationale)``.
# The rationale is printed alongside the result so a reader scanning the
# output can tell which behaviour each line exercises.
ClaimSpec = tuple[str, dict, str]


CLAIMS: Final[list[ClaimSpec]] = [
    # --- Bulk: gene resolution -------------------------------------------
    (
        "bulk",
        {"subject": "WNT4", "qualifier": "differentially_expressed_in",
         "object": "regenerating_oral_pole"},
        "Gene symbol (subject) resolves to HGNC plus xref map",
    ),
    (
        "bulk",
        {"subject": "MAPK1", "qualifier": "regulated_by",
         "object": "MAP2K1"},
        "Two gene symbols, both resolve via HGNC",
    ),
    # --- Bulk: accession resolution --------------------------------------
    (
        "bulk",
        {"subject": "study_warner_2025", "qualifier": "consolidates_dataset",
         "object": "GSE168938"},
        "GEO accession resolves to GEO:GSE168938",
    ),
    # --- Bulk: corpus-harvested concept (platform) -----------------------
    (
        "bulk",
        {"subject": "study_warner_2025", "qualifier": "sequenced_on",
         "object": "Illumina"},
        "Platform name resolves to NCIT:C215820 (OLS4 corpus harvest)",
    ),
    (
        "bulk",
        {"subject": "study_warner_2025", "qualifier": "sequenced_on",
         "object": "NovaSeq"},
        "Instrument series resolves to EFO:0008637 Illumina NovaSeq 6000 (substring match)",
    ),
    # --- Bulk: qualifier-aware (the STAR collision) ----------------------
    (
        "bulk",
        {"subject": "study_warner_2025", "qualifier": "aligned_with",
         "object": "STAR"},
        "STAR under a study-level qualifier resolves to biotools:star (aligner)",
    ),
    (
        "bulk",
        {"subject": "STAR", "qualifier": "expressed_in",
         "object": "regenerating_blastema"},
        "STAR under a feature qualifier resolves to HGNC:11359 (gene)",
    ),
    # --- Bulk: bio.tools + concept index ---------------------------------
    (
        "bulk",
        {"subject": "study_warner_2025", "qualifier": "quantified_with",
         "object": "DESeq2"},
        "DE tool resolves to biotools:deseq2",
    ),
    (
        "bulk",
        {"subject": "fold change", "qualifier": "tested_with",
         "object": "Wilcoxon"},
        "Static concept (subject) plus corpus stat (object, prefix match)",
    ),
    # --- Mouse (MGI) via Tabula Muris-style claim ------------------------
    (
        "singlecell",
        {"subject": "Foxp3", "qualifier": "marker_of_cell_type",
         "object": "regulatory T cell"},
        "Mouse: Title-cased Foxp3 resolves to MGI via the casing heuristic; Treg resolves to CL",
    ),
    (
        "bulk",
        {"subject": "Mapk1", "qualifier": "differentially_expressed_in",
         "object": "intestine"},
        "Mouse via heuristic: Mapk1 picks MGI (and NOT MAPK1 picking HGNC)",
    ),
    (
        "bulk",
        {"subject": "MAPK1", "qualifier": "differentially_expressed_in",
         "object": "intestine", "species": "mouse"},
        "Explicit species hint overrides even when the symbol is UPPERCASE",
    ),
    # --- ChEBI / HP / DOID via per-ontology harvested indexes -----------
    (
        "bulk",
        {"subject": "study_warner_2025", "qualifier": "tested_with",
         "object": "U0126"},
        "ChEBI: U0126 (Warner's MEK inhibitor) resolves to CHEBI:90693",
    ),
    (
        "bulk",
        {"subject": "patient_42", "qualifier": "expressed_in",
         "object": "intellectual disability"},
        "HP: a clinical phenotype resolves to HP:0001249 (forces species=human)",
    ),
    (
        "bulk",
        {"subject": "tumor_sample_7", "qualifier": "consolidates_dataset",
         "object": "breast cancer"},
        "DOID: a disease label resolves to DOID:1612 (cancer subtype)",
    ),
    # --- Cross-axis multi-ontology coverage -----------------------------
    (
        "bulk",
        {"subject": "tabula_muris_study", "qualifier": "expressed_in",
         "object": "limb bud", "species": "mouse"},
        "EMAPA: mouse developmental anatomy resolves to EMAPA:35944",
    ),
    (
        "bulk",
        {"subject": "drosophila_study", "qualifier": "expressed_in",
         "object": "wing disc", "species": "drosophila"},
        "FBbt: Drosophila imaginal disc resolves to FBbt:00001778",
    ),
    (
        "bulk",
        {"subject": "ms_proteomics_study", "qualifier": "quantified_with",
         "object": "orbitrap"},
        "MS: mass-spectrometry instrument resolves to MS:1000484",
    ),
    (
        "bulk",
        {"subject": "ptm_study", "qualifier": "tested_with",
         "object": "phosphorylated residue"},
        "PSI-MOD: PTM class resolves to MOD:00696 (multi-axis coverage)",
    ),
    # --- Drosophila / C. elegans / Zebrafish ----------------------------
    (
        "bulk",
        {"subject": "dpp", "qualifier": "expressed_in",
         "object": "wing_disc"},
        "Drosophila: dpp (decapentaplegic) resolves to FlyBase, single-candidate path",
    ),
    (
        "bulk",
        {"subject": "unc-22", "qualifier": "marker_of",
         "object": "body_wall_muscle"},
        "C. elegans: unc-22 resolves to WormBase via the xxx-N heuristic",
    ),
    (
        "bulk",
        {"subject": "nanog", "qualifier": "marker_of_cell_type",
         "object": "embryonic_stem_cell", "species": "zebrafish"},
        "Zebrafish via species hint: nanog resolves to ZFIN "
        "(no hint would default to human via Nematostella cross-listing)",
    ),
    # --- Single-cell: anatomy plus gene resolver -------------------------
    (
        "singlecell",
        {"subject": "MMP9", "qualifier": "localized_in_tissue",
         "object": "mesentery"},
        "Single-cell: gene subject plus UBERON anatomy object",
    ),
    (
        "singlecell",
        {"subject": "WNT4", "qualifier": "marker_of_cell_type",
         "object": "cnidocyte"},
        "Single-cell: gene subject plus CL cell type, both slots fully resolved",
    ),
    # --- Single-cell: intentionally unenriched --------------------------
    (
        "singlecell",
        {"subject": "MMP14", "qualifier": "localized_in_spatial_domain",
         "object": "physa"},
        "Subject (gene) resolves; physa stays bare (cnidarian-specific, no UBERON term)",
    ),
]


def _format_dict(d: dict) -> str:
    """Pretty-print a claim dict with sorted keys for stable diffs."""
    return json.dumps(d, indent=2, ensure_ascii=False, sort_keys=True)


def _diff_keys(before: dict, after: dict) -> list[str]:
    """Return the top-level keys that ``after`` has but ``before`` lacks."""
    return sorted(set(after.keys()) - set(before.keys()))


def main() -> int:
    bulk = BulkTranscriptomicsAdapter()
    sc = SingleCellTranscriptomicsAdapter()

    adapters = {"bulk": bulk, "singlecell": sc}

    n_enriched = 0
    n_no_change = 0

    print("=" * 70)
    print(f"Transcriptomics adapter end-to-end demo ({len(CLAIMS)} claims)")
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

        marker = "[OK] enriched" if was_enriched else "[--] no enrichment"
        print()
        print(f"[{i:2d}] {adapter_name:10s} {marker:18s}  {rationale}")
        print(f"    IN  : {_format_dict(claim_in)}")
        if was_enriched:
            curies = ", ".join(f"{slot}={c}" for slot, c in added_terms.items())
            print(f"    OUT : ontology_terms[{curies}]")
            xrefs = claim_out.get("xrefs") or {}
            for slot, xref in xrefs.items():
                # Compress the xrefs map to the principal id key plus the
                # human-readable label, so the demo output stays scannable.
                # ``hgnc`` covers gene xrefs (which still use per-namespace
                # keys), ``curie`` covers anatomy + concept xrefs (which use
                # the generic shape since the multi-species anatomy refactor).
                principal = xref.get("hgnc") or xref.get("curie")
                label = xref.get("label", "")
                if label:
                    print(f"          xrefs[{slot}] : {principal} \"{label}\"")
                else:
                    print(f"          xrefs[{slot}] : {principal}")
        else:
            print("    OUT : (no resolution; slot values pass through)")
        print(f"    new keys on output dict: {new_keys}")

    print()
    print("=" * 70)
    print(f"Summary: {n_enriched}/{len(CLAIMS)} claims enriched, "
          f"{n_no_change}/{len(CLAIMS)} unchanged.")
    print("=" * 70)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
