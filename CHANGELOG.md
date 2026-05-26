# Changelog

## v0.1.2 — comprehensive HGNC + Nematostella gene index (2026-05-26)

### Added

- `data/build_gene_index.py` maintainer script that regenerates the bundled
  TSV from two authoritative sources: the **HGNC complete set** (downloaded
  with a fallback mirror) and the **G2V2 Nematostella cross-reference table**
  (Nv2g / Nv2t / JGI / NvERTx / NVE alignment with UniProt-derived
  orthology, GO and KEGG annotations). The script accepts both sources as
  local paths; the HGNC TSV is cached under `~/.cache/`.
- Two-pass join strategy:
  1. Nematostella-anchored rows from G2V2 — HGNC ortholog attached via
     UniProt accession when it matches HGNC, otherwise via a
     case-insensitive `Preferred_name → HGNC symbol` fallback.
  2. Pure-HGNC rows for every human gene that was *not* already paired with
     a Nematostella entry, so plain human-gene mentions (e.g. an MMP in a
     cnidarian-vs-mammal comparison) still resolve to their HGNC CURIE.
- Build-time-only dependencies (`pandas`, `openpyxl`, `requests`) declared
  under `[project.optional-dependencies] build`; the runtime keeps its
  minimal closure (indicium, indicium-adapters, rdflib, pyshacl).

### Changed

- Bundled `data/gene_index.tsv` rebuilt: **~53,000 entries** (up from 22
  hand-curated). ~7,100 Nematostella-anchored rows with HGNC ortholog
  identified, ~7,000 Nematostella-only, ~39,000 pure-HGNC human.
- Behaviour-based test refactor: gene-index tests now assert lookup
  semantics (case-insensitive symbol resolution, principal-CURIE order,
  xref attachment, both lookup paths exercised) rather than specific
  hand-curated entries. Two new regression tests cover the v0.1.2-specific
  capabilities (Nemve1 → HGNC via symbol fallback; pure-HGNC coverage for
  human-only mentions).

## v0.1.1 — bundled gene-name index (2026-05-26)

### Added

- Bundled `data/gene_index.tsv` — 22 hand-curated entries drawn from the two
  benchmark reference papers (Warner et al. 2025 supplementary DEG tables;
  Cheung et al. 2024 main text). Covers core regeneration signalling (WNT4,
  FGF8, BMP1, SPRY2, MEOX1, HMX2), MEK/ERK pathway (MAPK1, MAPK3, MAP2K1),
  cnidarian matrix and ECM (MMP9/14/17, TIMP2, COL4A1, COL24A1, VCAN, MFAP4),
  β-catenin (CTNNB1), Brachyury (TBXT) and a handful of transcription
  factors and adhesion proteins. Each row carries HGNC, Ensembl, NCBI Gene,
  Nemve1 / JGI, NvERTx and UniProt cross-references plus a description.
- `_enrichment.GeneIndex` real lookup class loading the bundled TSV at
  import time. Accepts gene symbols (case-insensitive), Nemve1 / JGI ids,
  and UniProt accessions; returns the full row dict or `None`.
- `BulkTranscriptomicsAdapter.enrich_claim()` now resolves gene names in
  `subject` / `object` to a principal CURIE (HGNC > Ensembl > NCBIGene)
  and attaches the full xref map under `xrefs[slot]`. Repository accession
  enrichment (v0.1.0) is unchanged and still takes priority over the gene
  lookup.
- Regression tests covering all three lookup paths (symbol case-insensitive,
  Nemve1 / JGI, UniProt), the principal-CURIE order, and the xref
  attachment. Coverage: 95.17%.

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
