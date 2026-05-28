# Changelog

## v0.6.1 — Mouse cross-references from MGI (Ensembl / Entrez / UniProt) (2026-05-28)

### Added

- `data/build_gene_index.py` now fetches three additional MGI
  cross-reference reports and joins them on the MGI accession:
  - `MRK_ENSEMBL.rpt` (16 MB) — mouse Ensembl gene id
    (``ENSMUSG…``). Yield: 63,787 mouse rows enriched (~85%).
  - `MGI_EntrezGene.rpt` (88 MB) — mouse NCBI Gene (Entrez) id.
    Filtered to ``Status="O" AND Marker Type="Gene"``. Yield: 45,496
    (~60%).
  - `MRK_SwissProt_TrEMBL.rpt` (2 MB) — mouse UniProt accession
    (first accession listed, may be either SwissProt-reviewed or
    TrEMBL). Yield: 23,534 (~31%) — UniProt is dominated by
    SwissProt coverage which is smaller than the full proteome.
- All three files are downloaded to the same cache directory as
  HGNC and MRK_List2 (``~/.cache/indicium-adapters-transcriptomics/``).
  ``--mgi-ensembl-cache`` / ``--mgi-entrez-cache`` /
  ``--mgi-uniprot-cache`` CLI options accept custom paths;
  ``--no-mgi-xrefs`` skips the v0.6.1 enrichment entirely
  (regression to the v0.6.0 minimal mouse rows).
- ``build_gene_index._load_mgi_xrefs`` reads the three header-less,
  positional MGI reports and merges them into a single
  ``{MGI accession → {ensembl, ncbigene, uniprot}}`` dict consumed
  by step 3 of ``build()``.

### Tests

- 1 new test pinning ``Mapk1``'s full xref map (MGI + Ensembl +
  Entrez + UniProt) — catches regressions in any of the three MGI
  report parsers.

### Notes

- Mouse rows now look semantically equivalent to human rows in
  terms of xref coverage. A typical mouse claim like
  ``Mapk1 differentially_expressed_in regenerating_tissue`` produces:
  ``ontology_terms[subject] = MGI:1346858`` +
  ``xrefs[subject] = {mgi=MGI:1346858, ensembl=ENSMUSG00000063358,
  ncbigene=26413, uniprot=Q3UF82}``.
- ``v0.6.0`` mouse rows that lacked these xrefs would now silently
  re-enrich on adapter import — no API change.

## v0.6.0 — Mouse (MGI) gene index + species-aware lookup (2026-05-28)

### Added — MGI mouse gene catalogue

- `data/build_gene_index.py` extended to fetch the **MGI marker list**
  (MRK_List2.rpt, ~78 MB). Filters to ``Marker Type == "Gene"`` to
  exclude BAC/YAC ends and syntenic markers. **75,250 mouse rows**
  added to the bundled index — every protein-coding (and most
  non-coding) mouse marker known to MGI.
- Mouse rows carry `species=mouse` and an `mgi=MGI:NNN` accession;
  ensembl / ncbigene / uniprot cross-refs are left empty in v0.6.0
  (those live in separate MGI report files — integrating them is on
  the v0.6.1 roadmap).
- Bundled `gene_index.tsv` grows from ~53k to **128,353 rows**:
  14,150 Nematostella-anchored + 38,953 pure-HGNC human + 75,250 MGI
  mouse.

### Added — species-aware lookup (B + C disambiguation)

- New schema columns ``mgi`` and ``species`` (``human`` / ``mouse`` /
  ``nematostella``). The species tag drives the v0.6.0 lookup rules.
- New optional `species` field on the claim dict (``"human"`` /
  ``"mouse"`` / ``"nematostella"``). The dispatcher forwards it to
  every resolver; the gene resolver uses it to disambiguate
  symbol-level collisions (``MAPK1`` ↔ ``Mapk1`` both lower-case
  to ``mapk1``).
- **Casing heuristic** for claims with no species hint:
  ``MAPK1`` (all uppercase) → human (HGNC); ``Mapk1`` (Title case)
  → mouse (MGI). This is the universal naming convention in the
  field and works without any LLM-side change.
- **Fallback priority** when neither hint nor heuristic disambiguates:
  human > mouse > nematostella (load order, encoded in
  ``_KNOWN_SPECIES``). Most transcriptomics claims are human-centric.
- **Cross-listed rows.** A Nematostella-anchored row that carries a
  human HGNC ortholog under ``hgnc`` is treated as "serving human"
  too. So ``MMP9`` (all-uppercase → heuristic says human) resolves
  to ``HGNC:7176`` via the Nematostella row's ortholog field, not
  to the mouse ``Mmp9`` row. Rule lives in
  ``GeneIndex._pick_for_species``.
- MGI accessions (``MGI:1346858``) in claim slots resolve directly
  via the new ``_by_mgi`` index — bypassing the symbol-collision
  logic entirely.

### Added — runtime + resolver-signature uniformity

- Principal CURIE priority extended: HGNC > **MGI** > Ensembl >
  NCBIGene. Mouse rows with empty ``hgnc`` naturally surface
  ``MGI:NNN`` as principal; Nematostella rows with a human ortholog
  still surface HGNC first.
- ``xrefs[slot]`` now also surfaces the ``mgi`` column for any row
  that carries it.
- All three resolvers (accession / gene / concept / anatomy) now
  share a uniform signature ``(name, *, species=None)`` so the
  dispatcher can forward the hint without branching. Non-gene
  resolvers explicitly ignore the hint (anatomy and concept entries
  are not species-specific in the current bundles).

### Added — tests

- 5 new tests in ``test_bulk.py``:
  - mouse via Title-case heuristic (``Mapk1`` → MGI:1346858)
  - human via UPPER-case heuristic (``MAPK1`` → HGNC:6871)
  - species hint overrides heuristic (``MAPK1``+``species=mouse``
    → MGI:1346858)
  - MGI accession resolves directly (``MGI:1346858`` → MGI:1346858)
  - mouse-only convenience pin (``Foxp3`` → MGI:1891436)

### Updated — demo

- ``examples/enrich_warner2018.py`` gains 3 v0.6.0 claims: Foxp3 →
  MGI + Treg → CL (cross-mechanism), Mapk1 → MGI via heuristic,
  MAPK1+species=mouse → MGI override. 15/15 claims enrich.

### Notes

- v0.6.1 follow-up — wire MGI Ensembl / MGI EntrezGene /
  MGI SwissProt cross-ref files into ``build_gene_index.py`` so
  mouse rows ship full xref maps (currently only ``mgi`` is
  populated).
- v0.6.2 follow-up — zebrafish (ZFIN), C. elegans (WormBase),
  Drosophila (FlyBase) following the same multi-species pattern.

## v0.5.0 — Cell Ontology enrichment + OLS4 cascade refinements (2026-05-28)

### Added — Cell Ontology (CL) harvest

- `data/build_cl_concept_index.py` maintainer script that resolves a
  curated list of 37 cell-type candidate names against the **Cell
  Ontology** via OLS4 search with the same cascade match modes used
  for the corpus harvest (exact label → exact synonym → prefix).
  Substring mode is omitted to avoid false positives — cell-type
  labels in CL are already minimal.
- Candidate list drawn from the three cell-atlas papers in the
  6-paper mini-corpus: **Sebé-Pedrós 2018** (*Nematostella vectensis*
  whole-organism cell atlas), **Tabula Muris 2018** (mouse), and
  **Tabula Sapiens 2022** (human). Plus canonical mammalian immune /
  nervous / epithelial cell types that recur in any transcriptomics
  study.
- Bundled `data/cl_concept_index.tsv` — **35 verified entries**
  (~95% yield), including: `cnidocyte` (CL:0011029, the cnidarian-
  specific cell type — CL covers it even though UBERON has no
  cnidarian anatomy), `T cell` / `B cell` / `macrophage` / `neutrophil`
  / `monocyte` (mammalian immune), `astrocyte` / `oligodendrocyte`
  (CNS), `fibroblast` / `endothelial cell` / `hepatocyte` (epithelial /
  connective), `enterocyte` / `goblet cell` / `Paneth cell` (GI), etc.
- `data/cl_concept_unresolved.txt` — 2 unresolved candidates
  (`gland cell`, `nematocyte` — the latter deduped against `cnidocyte`
  since both map to CL:0011029).
- Loader now reads four TSVs: static + OLS4 corpus + bio.tools + CL.
  Same priority ladder (first match wins) with CL placed last.

### Changed — OLS4 cascade matcher fixes

- **Label-first cascade pass in exact mode.** Previously the exact
  matcher returned the first OLS4 result whose label OR synonym
  matched the corpus term. OLS4's relevance ranking sometimes puts a
  synonym-match doc ahead of the label-match doc — e.g. searching
  ``macrophage`` against CL returns CL:0000394 plasmatocyte (which
  has macrophage as a synonym) before CL:0000235 macrophage (whose
  label *is* macrophage). The v0.5.0 cascade does label-exact across
  all docs first, then synonym-exact second; the canonical label-
  match wins.
- **Fix OLS4 search-endpoint synonym field name.** Search docs expose
  synonyms under ``exact_synonyms``, not ``synonym`` (the latter is
  used by the term-by-IRI endpoint). Pre-v0.5.0 the synonym pass was
  silently dead because the code looked at the wrong field; resolutions
  that should have used synonym matching fell through to prefix /
  substring matching instead. The fix is defensive — both
  ``exact_synonyms`` and ``synonym`` are now consulted.
- Net effect on the v0.3.1 corpus index (re-built with the fixes):
  one new resolution (``10x Genomics`` → ``EFO:0008995`` 10x technology)
  plus three canonicalisations:
  - ``pacbio`` from ``EFO:0920006`` PacBio Revio → ``NCIT:C146819``
    **SMRT Sequencing** (the technology, not a specific model).
  - ``oxford nanopore`` from ``OBI:0002750`` ONT MinION →
    ``NCIT:C146818`` **Nanopore Sequencing**.
  - ``bulk rna-seq`` from ``OBI:0003090`` bulk RNA-seq assay →
    ``NCIT:C219552`` **Bulk RNA Sequencing**.
  All three new mappings are semantically more canonical — the
  corpus term is the technology, not a specific instrument / assay
  realisation. Pinned tests updated accordingly. Corpus yield rises
  from 16 → 17 / 38 (~45%).

### Added — tests

- `test_cl_cell_types_resolve_via_concept_index` pins 12 CL
  resolutions including the macrophage-not-plasmatocyte assertion
  (proves the label-first refinement is in effect end-to-end).
- v0.4.0 test `test_singlecell_resolves_gene_subject_for_marker_of_cell_type`
  updated: pre-v0.5.0 it asserted ``cnidocyte`` stayed unenriched
  (CL not wired); v0.5.0 asserts the object slot resolves to
  ``CL:0011029``. The test now exercises both slots of the
  ``marker_of_cell_type`` qualifier — gene subject + CL object.
- `test_corpus_concept_terms_resolve_via_concept_index` updated for
  the v0.5.0 canonicalised mappings (pacbio / oxford nanopore /
  bulk rna-seq) plus the new ``10x genomics`` resolution. Pins the
  match-mode tag in the source column.
- Size invariant now ≥92 rows (was ≥56): 31 static + 17 OLS4-corpus
  + 9 bio.tools + 35 CL.

### Notes

- Single-cell now has **two-slot resolution** for all five of its
  cell-resolution qualifiers — both the gene subject (v0.4.0 wiring)
  and the cell-type object (v0.5.0 wiring) get CURIEs. The
  `marker_of_cell_type` qualifier is now fully expressive.
- The 24 corpus terms unresolved at v0.3.1 → 21 at v0.5.0. The four
  changes: +1 new (`10x Genomics`), the three canonicalisations don't
  change the unresolved count (they were already resolved, just to
  a different CURIE).
- Demo `examples/enrich_warner2018.py` updated to reflect the new
  v0.5.0 behaviour on claim 11 (cnidocyte resolves to CL:0011029
  instead of staying bare).

## v0.4.0 — qualifier-aware resolution priority (2026-05-27)

### Changed

- `BulkTranscriptomicsAdapter.enrich_claim()` is now **qualifier-aware**.
  Previously the resolver chain was hardcoded
  ``accession`` → ``gene`` → ``concept`` for all qualifiers; this always
  resolved short ambiguous names like ``STAR`` to the HGNC gene
  (HGNC:11359, steroidogenic acute regulatory protein) even when the
  claim's qualifier made it obvious the term was the aligner.
- Two priority orders are now defined and picked per claim:
  - **Default** (``accession`` → ``gene`` → ``concept``) — used for the
    claim-level feature qualifiers (``differentially_expressed_in``,
    ``co_expressed_with``, ``enriched_for``, ``expressed_in``,
    ``marker_of``, ``regulated_by``) and for ``consolidates_dataset``.
    Gene mentions dominate these claims, so gene-lookup wins for
    collision-prone names.
  - **Study-level** (``concept`` → ``accession`` → ``gene``) — used for
    the six study-metadata qualifiers ``sequenced_on``,
    ``aligned_with``, ``quantified_with``, ``normalized_with``,
    ``tested_with``, ``uses_database``. Object slot of these is a tool
    / platform / database name, so the concept index (including
    bio.tools) wins; the gene resolver only fires as a fallback.
- Internally, the three previously-inlined resolvers are now
  first-class functions (``_resolve_accession``, ``_resolve_gene``,
  ``_resolve_concept``) dispatched through a ``_RESOLVERS`` map. Each
  returns ``(curie, xrefs_dict) | None`` so the resolver-chain is now
  a simple loop. Same shape, half the code.

### Added

- 4 regression tests pinning the v0.4.0 contract:
  - ``aligned_with: <study> STAR`` → ``biotools:star`` (was HGNC:11359
    pre-v0.4.0).
  - ``expressed_in: STAR <something>`` → ``HGNC:11359`` (gene-priority
    preserved for feature qualifiers).
  - Missing or unknown qualifier → default priority (safe fallback).
  - Study-level priority still falls through to gene-lookup when the
    concept index has no hit, so a ``aligned_with`` claim whose object
    is a gene name (sloppy LLM output) still gets a CURIE rather than
    a dropped slot.

### Changed — SingleCellTranscriptomicsAdapter (symmetric refactor)

- ``SingleCellTranscriptomicsAdapter`` now uses the same
  resolver-dispatcher pattern and **gains a gene resolver**. Previously
  single-cell had only anatomy + concept lookup wired, so a slot like
  ``marker_of_cell_type: WNT4 cnidocyte`` left ``WNT4`` unenriched
  even though it's a perfectly valid HGNC symbol. Now it resolves to
  ``HGNC:12783``.
- Single-cell uses one default priority order (``gene → anatomy →
  concept``) for all its qualifiers — they're all claim-level
  (cell-resolution feature assertions), so no study-level override is
  needed. Anatomy still wins for tissue terms (gene-lookup of
  ``mesentery`` returns ``None``, the cascade falls through to
  anatomy).
- 2 new tests pin the contract: gene resolution for the cell-type
  qualifiers, and the gene-doesn't-override-anatomy fallback.

### Added — shared resolver helpers

- The three resolvers (gene / anatomy / concept) moved into
  ``_enrichment.py`` as module-level callables (``resolve_gene``,
  ``resolve_anatomy``, ``resolve_concept``), keeping the row → CURIE
  rendering logic next to the corresponding Index classes. Bulk and
  single-cell both import and dispatch through them — one source of
  truth for the rendering rules.

### Added — local CI parity helper

- ``scripts/check.ps1`` runs the same checks the GitHub Actions CI
  runs (ruff + pytest with coverage). Use it before every push;
  pytest alone misses lint failures. Direct response to the v0.3.1
  CI red caused by a ruff E741 that didn't show up locally.

### Notes

- The 13 corpus terms unresolved at v0.3.1 (Cell Ranger, scVI, FLAIR,
  TALON, IsoQuant, SQANTI, StringTie2, Bambu, FLAMES, LRGASP, 10x
  Genomics, Iso-Seq, ONT direct RNA) are unaffected by v0.4.0 — they
  were absent from both OLS4 and bio.tools, not blocked by priority
  ordering.
- A monkeypatch test that previously patched ``bulk.CONCEPT_INDEX`` /
  ``singlecell.CONCEPT_INDEX`` directly was updated to patch
  ``_enrichment.CONCEPT_INDEX`` instead — since v0.4.0, both adapters
  read the index through the shared ``resolve_concept`` helper and
  the singleton's home is ``_enrichment``.

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
