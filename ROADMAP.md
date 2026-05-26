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

### v0.2.0 — SingleCell enrichment (CL + UBERON)
- `SingleCellTranscriptomicsAdapter.enrich_claim()` gets a real
  implementation: lookup cell-type subjects/objects against a bundled CL
  subset, tissue/anatomy literals against UBERON.
- Bundled CL + UBERON subsets via `data/build_vocabs.py`.
- Driven by a reference scRNA paper supplied by the intern; the bundled
  vocabulary subsets are scoped to terms exercised by that paper, expanded
  as additional papers reveal the need.

### v0.3.0+ — broader gene-name coverage and parsers
- Extend the gene-name index to mouse (MGI), zebrafish (ZFIN), Drosophila
  (FlyBase), C. elegans (WormBase), Arabidopsis (TAIR), yeast (SGD).
- GEO SOFT / ENA XML / BioStudies / ISA-Tab parser helpers projecting series
  metadata into claim dicts (parallel to the metabolomics mzTab parser).
