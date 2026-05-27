# Changelog

## v0.3.1 — corpus-driven concept harvest: OLS4 (cascade modes) + bio.tools (2026-05-27)

### Added — OLS4 corpus harvester with cascade match modes

- `data/build_corpus_concept_candidates.py` maintainer script that
  ingests the per-paper vocabulary JSONs produced by
  `articles/_extract_vocab.py` over the 6-paper mini-corpus and uses
  the **OLS4 search-by-label** endpoint to propose CURIEs, then
  re-verifies each by IRI lookup. Hand-curation of CURIEs is avoided;
  the only curation is the axis → ontology preference list and the
  axis → match-mode declaration.
- Per-axis **cascade of match modes** tried in order; the first hit
  wins. Strictness beats preference: an exact match in any preferred
  ontology beats a relaxed match in the first preferred ontology.
  - `exact` — label or synonym equals the corpus term verbatim
    (case-insensitive). Default for all software axes and for high-
    precision lookups.
  - `prefix` — corpus term is a whole-word **prefix** of the OLS4
    label (lets descriptive suffixes through: "Wilcoxon" → "Wilcoxon
    Signed Rank Test"; "bulk RNA-seq" → "bulk RNA-seq assay";
    "Visium" → "Visium Spatial Gene Expression"). Used by
    `protocol` and `stat_test` axes.
  - `substring` — corpus term is a whole-word **substring** of the
    OLS4 label (handles instrument families whose canonical label
    prepends a vendor name: "NovaSeq" → "Illumina NovaSeq 6000").
    Used by the `platform` axis only; OLS4 ranking picks the
    canonical model.
- Bundled `data/corpus_concept_index.tsv` — **16 verified entries**:
  - exact-match (6): `Illumina` (NCIT:C215820), `MARS-seq`
    (EFO:0008796), `Smart-seq2` (EFO:0008931), `Drop-seq`
    (EFO:0008722), `Slide-seq` (EFO:0009920), `t-test` (NCIT:C53231,
    via synonym).
  - prefix-match (6): `PacBio` (EFO:0920006 PacBio Revio),
    `Oxford Nanopore` (OBI:0002750 Oxford Nanopore MinION),
    `bulk RNA-seq` (OBI:0003090), `Visium` (EFO:0010961),
    `Wilcoxon` (NCIT:C53246), `Benjamini-Hochberg` (NCIT:C61596).
  - substring-match (4): `NextSeq` (EFO:0008566 Illumina NextSeq 550),
    `NovaSeq` (EFO:0008637 Illumina NovaSeq 6000), `HiSeq`
    (EFO:0004203 Illumina HiSeq 2000), `Sequel` (EFO:0008630
    PacBio Sequel system).
- The match mode is tagged in the TSV's `source` column
  (`corpus extraction (<axis>, <mode>-match)`) so the provenance of
  each row is auditable.

### Added — bio.tools registry harvester (with name normalisation)

- `data/build_biotools_index.py` maintainer script picking up where
  OLS4 left off: for software-axis terms still unresolved, queries
  the **bio.tools registry**, the canonical home for bioinformatics
  tools. Two-step resolution per term:
  1. Direct API lookup `GET https://bio.tools/api/<id>/?format=json`
     with a normalised ID (lowercased, spaces → hyphens).
  2. If 404, search endpoint with paginated results, filtered
     client-side for **exact name match (case-insensitive)** — and
     additionally for a **normalised-name match** that drops
     whitespace / hyphens / underscores so corpus `Bowtie2` matches
     the canonical bio.tools name `Bowtie 2`. Without this
     normalisation step the OLS4 harvester's substring trick has no
     bio.tools equivalent.
- Bundled `data/biotools_concept_index.tsv` — **9 verified entries**:
  `STAR` (biotools:star), `Bowtie` (biotools:bowtie), **`Bowtie2`**
  (biotools:bowtie2, via name normalisation), `Minimap2`
  (biotools:minimap2), `RSEM` (biotools:rsem), `StringTie`
  (biotools:stringtie), `DESeq2` (biotools:deseq2), `Seurat`
  (biotools:seurat), `Picard` (biotools:picard_tools).
- New `biotools:` prefix in `TRANSCRIPTOMICS_PREFIXES` (resolves to
  `https://bio.tools/<id>`).
- `data/biotools_unresolved.txt` — final log of the 13 corpus terms
  not findable in OLS4 nor bio.tools. Almost all are long-read tools
  (FLAIR, TALON, IsoQuant, SQANTI, StringTie2, Bambu, FLAMES, LRGASP)
  that genuinely do not exist in either registry — they live in
  their authors' GitHub / Zenodo only.

### Added — runtime + tests

- **Triple-TSV loader** in `_enrichment._load_concept_index_rows()`:
  loads `concept_index.tsv` → `corpus_concept_index.tsv` →
  `biotools_concept_index.tsv` in that order. A hypothetical
  collision resolves to the first occurrence: paper-modelling-
  grounded (static) wins over OLS4 method-mention wins over
  bio.tools-registry.
- 6 regression tests covering: pinned resolution for each of the 16
  OLS4-corpus terms (with the match-mode tag asserted in the source
  column), pinned resolution for each of the 9 bio.tools terms, two
  end-to-end `enrich_claim()` paths (one OLS4-resolved platform via
  `Illumina`, one bio.tools-resolved aligner via `minimap2` — `STAR`
  is deliberately avoided because it also matches HGNC:11359
  steroidogenic acute regulatory protein, which gene-lookup correctly
  resolves first), and the triple-TSV size lower bound (≥56 rows).

### Notes

- Total v0.3.1 yield: **25/38 corpus terms resolved** (~66%) —
  16 via OLS4 (6 exact + 6 prefix + 4 substring) + 9 via bio.tools.
- The 13 still-unresolved terms fall into one clean bucket: niche
  long-read tools that are genuinely outside both OLS4 (no ontology
  covers them) and bio.tools (not registered). Resolution would
  require either upstream registry contributions or per-tool
  GitHub/Zenodo CURIE schemes — neither is on the v0.3.x roadmap.
- A surprise found while writing the tests: some short bio.tools
  software names (`STAR`, `Picard`) **also** exist as HGNC gene
  symbols. The current resolver runs gene-lookup before
  concept-lookup, so a bare `STAR` in a claim resolves to HGNC:11359.
  This is correct default behaviour (gene mentions are more frequent
  in transcriptomics claims), but a qualifier-aware resolver
  (`aligned_with` → prefer bio.tools; `expressed_in` → prefer HGNC)
  would be a useful v0.4.0 refinement.

## v0.3.0 — TranscriptoKG concept-index enrichment (NCIT / INO / OBI / STATO / GO) (2026-05-27)

### Added

- New ontology axes wired into `TRANSCRIPTOMICS_PREFIXES`: **INO**
  (regulation processes), **NCIT** (analysis concepts), **SIO**, **RO** and
  **IAO** — drawn from the TranscriptoKG vocabulary work so adapters and
  the upstream KG stay aligned on the same CURIE namespace set.
- `data/build_concept_index.py` maintainer script that verifies a curated
  list of 31 CURIEs against EBI OLS4 and writes a TSV with the canonical
  label and definition. The candidate set is drawn from the TranscriptoKG
  mapping & vocab (Boukhajou, 2026): 2 INO regulation terms, 20 NCIT
  analysis concepts (Fold Change, log2, Timepoint, Multiple Comparison
  Adjustment, RNA Sequencing Data File, ...), 3 GO roots, 5 OBI
  investigation / role / protocol / DE-transformation / RNA-sequencing
  terms, and 1 STATO statistical-variable anchor. The build retries on
  OLS4 5xx and throttles 0.3s between requests.
- Bundled `data/concept_index.tsv` — **31 verified entries**, every one
  re-checked against OLS4 at build time (no hand-curation in the shipped
  row). Two OBI entries (`OBI:0000249`, `OBI:0000650`) had been tagged
  incorrectly in the upstream TranscriptoKG vocab — OLS4 verification
  caught the drift and the source notes were corrected to reflect the
  canonical labels ("technical replicate role" and "differential
  expression analysis data transformation"); the terms themselves stay
  relevant to the transcriptomics domain so the rows are kept.
- `_enrichment.ConceptIndex` real lookup class loading the bundled TSV
  at import time. Accepts a term (case-insensitive) and returns the full
  row dict (CURIE / label / description / ontology) or `None`. Strict
  header check at load time: a malformed TSV raises before silent
  mis-mapping can occur.
- Concept-lookup branch wired into both adapters' `enrich_claim()` — runs
  after gene / accession / anatomy lookup, so domain-specific resolution
  always takes precedence. Adds the principal CURIE under
  `ontology_terms[slot]` and `{curie, label, ontology, description}`
  under `xrefs[slot]`.
- 7 regression tests covering: index loadability, `None` for unknown /
  empty / non-str input, safe-no-op for unknown terms, end-to-end CURIE +
  xref attachment, single-cell-side wiring via monkeypatched index.

### Changed

- Test docstrings that referenced "while the TSV is empty (pending OLS4
  verification)" updated to reflect the populated v0.3.0 state.
- `BulkTranscriptomicsAdapter` and `SingleCellTranscriptomicsAdapter` now
  expose 17 / 12 effective lookup axes respectively (gene + accession +
  anatomy + concept).

### Notes

- A second concept-harvesting pass (corpus-driven: tools / platforms /
  protocols extracted from the 6-paper mini-corpus — STAR, Seurat,
  MARS-seq, Illumina, 10x, etc.) is scheduled for v0.3.1 / v0.4.0 once a
  `build_corpus_concept_candidates.py` script using OLS4
  search-by-label is written. Hand-curation of CURIEs is explicitly
  avoided.
- The earlier ROADMAP note about contributing cnidarian-specific anatomy
  upstream to UBERON is unchanged — still on the v0.4.0 wish-list.

## v0.2.0 — UBERON anatomy enrichment for SingleCellTranscriptomicsAdapter (2026-05-26)

### Added

- `data/build_anatomy_index.py` maintainer script that resolves a curated
  list of anatomy terms (drawn from the Cheung et al. 2024 Tomo-seq paper
  on *Nematostella vectensis*) against the EBI OLS4 UBERON index. Only
  entries whose returned label exactly matches the query are written —
  hyponym promotion (e.g. *"nasal tentacle"* for query *"tentacle"*)
  cannot slip through.
- Bundled `data/anatomy_index.tsv` (4 verified UBERON anchors for v0.2.0):
  `mesentery` (UBERON:0002095), `pharynx` (UBERON:0006562), `body wall`
  (UBERON:0000309), `organism subdivision` (UBERON:0000475). The bundle
  is deliberately small — UBERON has no entry for cnidarian-specific
  concepts (physa, mesoglea, gastrodermis, polyp body form, oral / aboral
  pole, body column), so the build script documents them in
  `TRIAGED_OUT_KNOWN_CNIDARIAN` rather than emitting wrong CURIEs.
- `_enrichment.AnatomyIndex` real lookup class loading the bundled TSV at
  import time. Accepts an anatomical term (case-insensitive) and also
  matches against the UBERON `label` field; returns the full row dict
  or `None`.
- `SingleCellTranscriptomicsAdapter.enrich_claim()` now resolves anatomy
  terms in `subject` / `object` to their UBERON CURIE under
  `ontology_terms[slot]` and attaches the label + description under
  `xrefs[slot]` so consumers can render without an OLS4 round-trip.
- 7 regression tests covering: case-insensitive anatomy lookup, UBERON
  CURIE attachment, xref shape, both subject / object slots resolved,
  cnidarian-specific terms left untouched (no invented CURIEs).

### Notes

- Cell-type enrichment via the Cell Ontology (CL) for the cell-type
  qualifiers (`marker_of_cell_type`, `enriched_in_cluster`,
  `co_expressed_in_cell_type`) is scheduled for v0.2.1 once an
  scRNA-seq reference paper is selected.
- The four cnidarian-specific terms triaged out in v0.2.0 are candidates
  for a v0.3.0 effort to contribute them upstream to UBERON, or to allow
  alternative ontology prefixes (e.g. CARO, GO for ECM) in the anatomy
  index.

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
