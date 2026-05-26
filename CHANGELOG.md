# Changelog

## v0.1.0 — initial scaffold (2026-05-26)

### Added

- **`BulkTranscriptomicsAdapter`** (api_version=1) with 13 transcriptomics
  qualifiers (claim-level: `differentially_expressed_in`, `co_expressed_with`,
  `enriched_for`, `expressed_in`, `marker_of`, `regulated_by`; study-level:
  `sequenced_on`, `aligned_with`, `quantified_with`, `normalized_with`,
  `tested_with`, `uses_database`, `consolidates_dataset`).
- **`SingleCellTranscriptomicsAdapter`** (api_version=1) — structural stub
  with 5 cell-resolution qualifiers (`marker_of_cell_type`,
  `enriched_in_cluster`, `co_expressed_in_cell_type`, `localized_in_tissue`,
  `localized_in_spatial_domain`). `enrich_claim()` is a no-op pass-through in
  v0.1.0; real CL + UBERON lookup arrives in v0.2.0.
- `TRANSCRIPTOMICS_PREFIXES` IRI map (EFO, OBI, STATO, NCBITaxon, CL, UBERON,
  GO, SO, MONDO, PATO, EDAM plus HGNC, MGI, Ensembl, NCBI Gene, RefSeq,
  UniProt and the 8 major transcriptomics repository accession prefixes),
  shared by both adapters.
- `BulkTranscriptomicsAdapter.enrich_claim()` resolves repository accessions
  to CURIEs (GEO, SRA, ENA, ArrayExpress, ENCODE, BioProject, dbGaP,
  BioStudies) via the `accession_to_curie` helper. Gene-symbol enrichment is
  wired but no-op (real bundled index lands in v0.1.1).
- `shacl_shapes()` on both adapters returns a minimal SHACL graph validating
  CURIE format on ASB claim graphs.
- Entry-point registration under `indicium.adapters` for both
  `bulk-transcriptomics` and `single-cell-transcriptomics`.
- Pytest suite covering Protocol conformance, attribute presence, enrichment
  behaviour, accession-pattern table and the cross-adapter disjoint-qualifier
  invariant (`tests/test_composition.py`).
