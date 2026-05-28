# Roadmap

## Design philosophy — what belongs here

Before adding anything, ask the 4-quadrant question:

1. **Constrain** the agent (qualifier set, prompt-context wording)?
2. **Resolve** identifiers (gene symbol → HGNC, accession → CURIE,
   tool name → EDAM)?
3. **Validate** output (SHACL shape, CURIE-format check)?
4. Or is it something an agent could do as well freehand?

If only the last applies, **skip it**. The adapter is a constraints +
identifier-resolution + validation layer, not a re-implementation of what
upstream LLM extractors already do.

The transcriptomics adapter has one further axis to watch: the field carries
*many* legacy vocabularies (HGNC / MGI / Ensembl / NCBI Gene / RefSeq / Nemve1
/ species-specific symbols). The strategy is to ship a *curated, paper-driven*
gene-name index rather than a comprehensive multi-species dump; expand only as
new benchmark papers reveal the need.

## Versions

### v0.1.0 — scaffold + two-adapter structure (done)
- **`BulkTranscriptomicsAdapter`** — 13 qualifiers, prefix map, extraction
  context, repository accession → CURIE enrichment (8 prefixes), minimal
  SHACL CURIE-format shape.
- **`SingleCellTranscriptomicsAdapter`** — structural stub: 5 cell-resolution
  qualifiers, prefix map, extraction context, no-op `enrich_claim()`. CL +
  UBERON enrichment behaviour scheduled for v0.2.0.
- `tests/test_composition.py` enforces disjoint qualifier sets.
- Pytest conformance + behaviour suite for both adapters.

### v0.1.1 — bundled gene-name index (Bulk)
- Maintainer script `data/build_gene_index.py` fetching HGNC complete set +
  Nematostella vectensis Nemve1 gene models.
- Shipped TSV under `data/`, loaded into `_enrichment.GENE_INDEX`.
- `BulkTranscriptomicsAdapter.enrich_claim()` resolves gene symbols in
  `subject` / `object` to the principal CURIE (HGNC for human, NemveGene for
  Nematostella) with a full xref map for downstream consumers.
- Unit tests round-trip a representative sample.

### v0.2.0 — SingleCell anatomy enrichment (UBERON, done)
- `SingleCellTranscriptomicsAdapter.enrich_claim()` resolves anatomy terms
  in `subject` / `object` to UBERON CURIEs via a small hand-curated bundled
  index (4 verified entries from Cheung et al. 2024).
- `data/build_anatomy_index.py` maintainer script verifies every candidate
  CURIE against the EBI OLS4 UBERON index; only exact-label matches are
  kept (no hyponym promotion).
- Cnidarian-specific terms (physa, mesoglea, gastrodermis, oral / aboral
  pole, body column, polyp body form) are deliberately triaged out — UBERON
  has no entry for them as of 2026. Listed in
  `data/build_anatomy_index.py:TRIAGED_OUT_KNOWN_CNIDARIAN`.

### v0.2.1 — SingleCell cell-type enrichment (CL)
- Add CL (Cell Ontology) lookup for the cell-type qualifiers
  (`marker_of_cell_type`, `enriched_in_cluster`,
  `co_expressed_in_cell_type`).
- Driven by an scRNA-seq reference paper (e.g. Sebé-Pedrós et al. 2018
  whole-organism *Nematostella* cell-type atlas).
- Same build pattern as v0.2.0: maintainer script hits OLS4 CL endpoint,
  exact-label match only.

### v0.3.0 — TranscriptoKG concept-index enrichment (done)
- `data/build_concept_index.py` verifies 31 CURIEs against EBI OLS4:
  INO regulation processes, NCIT analysis concepts, GO roots, OBI
  investigation / role / protocol / DE-transformation / RNA-sequencing,
  STATO statistical variable.
- New prefixes wired into `TRANSCRIPTOMICS_PREFIXES`: INO, NCIT, SIO,
  RO, IAO.
- `_enrichment.ConceptIndex` shared by both adapters; concept lookup
  runs after gene / accession / anatomy resolution so domain-specific
  hits always win.

### v0.3.1 — corpus-driven concept harvesting: OLS4 cascade + bio.tools normalised (done)
- `data/build_corpus_concept_candidates.py` ingests the per-paper vocab
  JSONs and runs a **per-axis cascade** of match modes against OLS4
  (`exact` → `prefix` → `substring`), re-verifying each candidate
  by IRI lookup.
  - protocol / stat_test axes get `exact` + `prefix` (catches
    "Wilcoxon" → "Wilcoxon Signed Rank Test", "Visium" → "Visium
    Spatial Gene Expression", "bulk RNA-seq" → "bulk RNA-seq assay").
  - platform axis adds `substring` (catches "NovaSeq" → "Illumina
    NovaSeq 6000", "HiSeq" → "Illumina HiSeq 2000").
- `data/build_biotools_index.py` picks up where OLS4 left off with
  direct-ID lookup + name-normalised search fallback (drops
  whitespace / hyphens / underscores so corpus `Bowtie2` matches
  bio.tools canonical name `Bowtie 2`).
- Combined yield: **25/38 terms resolved** (~66%) — 16 OLS4
  (6 exact + 6 prefix + 4 substring) + 9 bio.tools.
- 13 terms still unresolved: niche long-read tools (FLAIR, TALON,
  IsoQuant, SQANTI, StringTie2, Bambu, FLAMES, LRGASP) which do not
  exist in OLS4 ontologies or bio.tools registry; plus a few protocols
  (`10x Genomics`, `Iso-Seq`, `ONT direct RNA`) and software
  (`Cell Ranger`, `scVI`) similarly absent.
- Triple-TSV loader in `ConceptIndex`: static (v0.3.0) → OLS4
  corpus (v0.3.1) → bio.tools (v0.3.1). Collision rule: first match
  wins, encoding a "grounding priority" ladder.

### v0.6.0 — Mouse (MGI) + species-aware lookup (done)
- ``build_gene_index.py`` fetches the MGI marker list (MRK_List2.rpt)
  in parallel with HGNC + G2V2. Filters to ``Marker Type == "Gene"``.
  75k mouse rows added; bundled index grows from 53k → 128k.
- ``gene_index.tsv`` schema gains ``mgi`` and ``species`` columns.
- ``GeneIndex.lookup`` disambiguates ``MAPK1`` ↔ ``Mapk1`` collisions
  via (1) explicit ``species`` hint on the claim, (2) casing
  heuristic, (3) load-order fallback. Cross-listed rows
  (Nematostella + human HGNC ortholog) serve human lookups too.
- Principal CURIE priority: HGNC > MGI > Ensembl > NCBIGene.

### v0.6.1 — Mouse cross-refs (done)
- ``build_gene_index.py`` fetches MGI's three xref reports
  (``MRK_ENSEMBL.rpt``, ``MGI_EntrezGene.rpt``,
  ``MRK_SwissProt_TrEMBL.rpt``) and joins them on the MGI
  accession. ~85% of mouse rows ship Ensembl, ~60% ship Entrez,
  ~31% ship UniProt (SwissProt coverage limit).

### v0.6.2 — Drosophila + C. elegans + Zebrafish (done)
- ``build_gene_index.py`` fetches one NCBI gene_info file per
  organism (~1-2 MB each) and emits species-tagged rows alongside
  human / mouse / Nematostella. Total bundled index now 254,922
  rows.
- New schema columns ``flybase`` / ``wormbase`` / ``zfin`` + direct
  accession lookup paths (both bare ids and prefixed CURIE forms
  accepted).
- C. elegans casing heuristic added (xxx-N pattern). Drosophila
  and zebrafish rely on the optional ``species`` claim hint or
  the single-candidate shortcut — the lowercase-only space is too
  ambiguous to auto-route reliably.

### v0.6.3 — UniProt cross-refs for the v0.6.2 species (done)
- ``build_gene_index.py`` fetches UniProt's per-organism idmapping
  files (DROME_7227 / CAEEL_6239 / DANRE_7955) and joins them on the
  organism-specific id. Yields ~55% Drosophila, ~43% C. elegans,
  ~39% zebrafish UniProt coverage.
- ``Ensembl_Metazoa`` accepted alongside ``Ensembl`` in NCBI
  gene_info dbXrefs — defensive, currently yields 0 for the
  invertebrates because NCBI doesn't surface those in gene_info.

### v0.6.4 — Invertebrate Ensembl ids (planned)
- NCBI's gene_info doesn't carry Ensembl Genomes refs for
  Drosophila / C. elegans. Source those from **Ensembl Genomes
  FTP** (https://ftp.ensemblgenomes.org/pub/metazoa/) — TSV
  cross-reference files per release. Join on FlyBase / WormBase
  id to populate ``ensembl`` for these rows.

### v0.10.0 — Integration smoke tests + Protocol-conformance suite (done)
- ``tests/test_protocol_conformance.py`` centralises the cross-cutting
  invariants (discoverability, snake_case qualifiers, kebab-case
  domain_id, IRI prefix-map validity, SHACL NodeShape presence).
- ``examples/quickstart.py`` 30-line bring-up demo.
- ``examples/test_on_papers.py`` real-paper smoke test (10 claims
  from Warner + Tabula Muris with SHACL validation).

### v0.11.0 — Real-corpus smoke test via SciLEx (planned)
- Use SciLEx (Celian Ringwald's literature crawler) to pull a
  ~10-paper diverse corpus covering the 6 supported organisms
  (human / mouse / Nematostella / Drosophila / C. elegans / zebrafish)
  across both modalities (bulk / single-cell).
- Hand-extract ~3-5 claims per paper from titles + abstracts.
- Expand ``examples/test_on_papers.py`` to run the full pipeline
  end-to-end on this diverse corpus + report a coverage matrix.

### v1.0.0 — Stable release (planned)
- ``CITATION.cff`` (siblings don't have it — leading move).
- README badges (CI / coverage / Python / license) + ontology-
  provenance table modelled on the epigenomics adapter's table.
- Tag ``v1.0.0`` + GitHub release notes synthesising every release
  from v0.4.0.

### v0.9.0 — Broader ontology coverage for benchmark consistency (done)
- 8 ontologies added via the existing REGISTRY-driven harvester:
  Tier-1 universal (MS / UO / MOD) + Tier-2 organism-paired anatomy
  / phenotype (FBbt / WBbt / ZFA / EMAPA / XPO). 201 new CURIEs.
- 4 audit-shortlisted ontologies deferred to v0.9.x with per-
  ontology rationale documented in build_extra_ontology_indexes.py:
  - EC (re-exports GO, no enzyme-class labels)
  - RXNO (CHEBI/BFO imports, no named reactions)
  - HCAO (CL/CARO re-publish, redundant)
  - ZP ("<entity> abnormal, abnormal" label pattern needs flipping)
- Loader made tolerant of missing TSVs so maintainer build steps
  that regenerate a subset don't crash on import.

### v0.9.x — refine the deferred 4
- EC: wire ExPASy enzyme.dat directly to map enzyme names ↔ EC
  numbers, bypass OLS4's label-less EC ontology.
- RXNO: needs a reaction-name database; defer until corpus shows
  named-reaction mentions worth harvesting.
- HCAO: redundant with CL — skip unless HCAO ships project-specific
  cell types not in CL.
- ZP: write a ZP-pattern-aware harvester that flips
  "abnormal X morphology" candidates to "X morphology abnormal,
  abnormal" before searching.

### v0.8.0 — Top-7 OLS4 ontology indexes (done)
- ``data/build_extra_ontology_indexes.py`` harvests CHEBI / HP / MP /
  DOID / PR / OBA / ENVO via the same cascade pattern as v0.5.0 CL.
  154 new CURIEs total; concept index grows from ~92 to ~246 rows.
- ``articles/_ols4_audit_report.md`` documents the full 280-ontology
  OLS4 catalogue with curated picks; future ontology additions
  should consult its queue.

### v0.8.x — Incremental ontology adds (driven by benchmark misses)
- Reserve for ontologies that **show up missing on real benchmark
  output**, not preemptive. From the audit's curated queue:
  remaining 21 candidates (MAXO, EC, MS, OMIT, HCAO, EMAPA, FBbt,
  WBbt, ZFA, ZP, XPO, NBO, RXNO, MOD, UO, BTO, OGSF, CMPO, etc.).
  Pick by benchmark-miss frequency, not by alphabetical / coverage
  obsession.

### v0.7.0 — Repository metadata parsers (done)
- ``parsers/`` subpackage with four functions:
  ``parse_geo_soft(text)``, ``parse_ena_xml(xml_text)``,
  ``parse_biostudies_json(text_or_obj)``,
  ``parse_isa_tab(directory)``. Uniform contract — input →
  ``list[dict]`` of claim dicts directly fed to
  ``adapter.enrich_claim``.
- All four parsers extract a ``species`` hint from the organism
  field (when recognised) and attach it to every emitted claim —
  so downstream gene resolution picks the right organism even
  when the LLM normalised symbol casing.

### v0.7.x follow-ups
- ENA: also parse ``<SAMPLE>`` descriptors (taxon id → organism
  hint) and ``LIBRARY_DESCRIPTOR > LIBRARY_STRATEGY``.
- ISA-Tab: parse ``Comment[<concept>]`` columns for finer
  protocol / treatment metadata.
- GEO: split ``!Sample_characteristics_ch1`` "tissue: X;
  age: Y" into separate claims (currently kept as raw string).

### v0.5.0 — Cell Ontology (CL) cell-type enrichment + cascade fixes (done)
- `data/build_cl_concept_index.py` harvests CL for the cell-type
  vocabulary mentioned across the cell-atlas papers in the mini-corpus
  (Sebé-Pedrós 2018 Nematostella, Tabula Muris/Sapiens). 35/37
  candidates resolved (~95%), incl. `cnidocyte` → CL:0011029.
- Single-cell `marker_of_cell_type` / `enriched_in_cluster` /
  `co_expressed_in_cell_type` qualifiers now resolve their cell-type
  objects via CL.
- OLS4 cascade: label-first pass in exact mode (fixes the
  macrophage → plasmatocyte ranking quirk) + fix
  `exact_synonyms` field name on search docs (was silently dead).
- v0.3.1 corpus re-built with the fixes: 16 → 17 hits, plus
  canonicalised mappings for `pacbio` (→ SMRT Sequencing),
  `oxford nanopore` (→ Nanopore Sequencing), `bulk rna-seq` (→ Bulk
  RNA Sequencing).

### v0.4.0 — qualifier-aware resolution + symmetric single-cell refactor (done)
- `BulkTranscriptomicsAdapter.enrich_claim()` picks its resolver chain
  based on the claim's qualifier:
  - Study-level qualifiers (`sequenced_on`, `aligned_with`,
    `quantified_with`, `normalized_with`, `tested_with`,
    `uses_database`) use `concept → accession → gene` — fixes the
    STAR/HGNC collision: `aligned_with: <study> STAR` now resolves
    to `biotools:star` instead of HGNC:11359.
  - Feature qualifiers and `consolidates_dataset` keep the default
    `accession → gene → concept`.
- `SingleCellTranscriptomicsAdapter` adopts the same dispatcher
  pattern and **gains a gene resolver**: claims like
  `marker_of_cell_type: WNT4 cnidocyte` now resolve WNT4 to HGNC:12783
  (unenriched pre-v0.4.0). Anatomy still wins for tissue names
  through the fallback cascade.
- Resolver helpers (`resolve_gene`, `resolve_anatomy`, `resolve_concept`)
  extracted to `_enrichment.py` so bulk and single-cell share the same
  row → CURIE rendering logic.
- `scripts/check.ps1` ships a local CI parity command (ruff +
  pytest+cov) to catch lint regressions before push.

### v0.4.0+ — broader coverage
- Contribute cnidarian-specific anatomy terms upstream to UBERON / CARO
  so future versions don't need to triage them out.
- Allow alternative ontology prefixes in the anatomy index (e.g. GO for
  extracellular matrix, NCBITaxon for organism-context claims).

### v0.4.0+ — broader gene-name coverage and parsers
- Extend the gene-name index to mouse (MGI), zebrafish (ZFIN), Drosophila
  (FlyBase), C. elegans (WormBase), Arabidopsis (TAIR), yeast (SGD).
- GEO SOFT / ENA XML / BioStudies / ISA-Tab parser helpers projecting series
  metadata into claim dicts (parallel to the metabolomics mzTab parser).
