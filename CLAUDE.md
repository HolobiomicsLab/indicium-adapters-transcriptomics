# CLAUDE.md — Dev Guide

## Read this before adding anything

The package's design philosophy and the 4-quadrant filter for new features
live at the top of `ROADMAP.md`. Any new contribution should justify itself
in those terms: does it **constrain** the agent, **resolve** identifiers,
**validate** output, or is it something an agent could do as well freehand?
If the last, skip it.

## Project shape

Two sibling `DomainAdapter` classes in one wheel:

- `BulkTranscriptomicsAdapter` (`bulk.py`) — bulk RNA-seq / microarray /
  time-series claims plus study-level assay metadata.
- `SingleCellTranscriptomicsAdapter` (`singlecell.py`) — cell-resolution
  claims (cell-type markers, cluster enrichment, spatial localisation).
  Structural stub in v0.1.0; CL + UBERON enrichment lands in v0.2.0.

Both target:

- `api_version = 1` matching `indicium-adapters.CURRENT_API_VERSION`.
- `DomainAdapter` and `SHACLProvider` Protocol conformance.
- Entry-point registration under `indicium.adapters`.
- Hardcoded namespace URI constants in `_namespaces.py`.
- **Disjoint qualifier sets** — enforced by `tests/test_composition.py`.

## Adding a new qualifier

1. Decide which adapter owns it (or whether it justifies a new sibling).
2. Add it to that adapter's `qualifiers` frozenset.
3. Update `extraction_context()` to mention it.
4. If it needs enrichment behaviour, add a branch in `enrich_claim()`.
5. If it needs validation, add a NodeShape in the adapter's shapes TTL.
6. Add unit tests in `tests/test_<adapter>.py`.
7. Keep qualifier sets disjoint across sibling adapters.

## Refreshing bundled data files

Maintainer-only scripts under `src/indicium_adapters_transcriptomics/data/`:

- `build_gene_index.py` (v0.1.1) — regenerates `gene_index.tsv` from HGNC's
  complete set + Nemve1 gene models. Caches to a `.json` (gitignored).

- `build_vocabs.py` (v0.2.0) — regenerates CL + UBERON + EFO inline subsets
  from their OWL releases.

Don't run either at install or import time.

## Test before commit

```bash
ruff check .
pytest --cov=indicium_adapters_transcriptomics
```

CI matrix is 3.12 / 3.13.

## Repo relationships

- Core Protocol lives in `~/git/indicium-adapters` (do not vendor it).
- The indicium claim/evidence vocabulary lives in `~/git/indicium`.
- Pipeline RDF validation for transcriptomics workflow classes (SO / UBERON /
  CL / EFO node shapes) belongs upstream in the consuming graph store — we do
  not duplicate it here.
