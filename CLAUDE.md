# CLAUDE.md — Dev Guide

## Read this before adding anything

The package's design philosophy and the 4-quadrant filter for new features
live at the top of `ROADMAP.md`. Any new contribution should justify itself
in those terms: does it **constrain** the agent, **resolve** identifiers,
**validate** output, or is it something an agent could do as well freehand?
If the last, skip it.

## Project shape

Two sibling `DomainAdapter` classes in one wheel:

- `BulkTranscriptomicsAdapter` (`bulk.py`): bulk RNA-seq, microarray and
  time-series claims plus study-level assay metadata. 13 qualifiers.
- `SingleCellTranscriptomicsAdapter` (`singlecell.py`): cell-resolution
  claims (cell-type markers, cluster enrichment, spatial localisation).
  5 qualifiers, gene + anatomy + concept enrichment all live.

Both share:

- The resolver back-end in `_enrichment.py` (`GENE_INDEX` over 254,922
  rows across six organisms; `ANATOMY_INDEX` over UBERON;
  `CONCEPT_INDEX` over 19 ontology sources for ~447 entries).
- The `TRANSCRIPTOMICS_PREFIXES` CURIE-prefix map in `_namespaces.py`.
- The classification vocabulary (`primary_domain` = `transcriptomics`,
  `subdomains`, `techniques`) in `_vocabulary.py`, declared identically
  by both adapters and unioned by ASB's classification layer.
- The repository-accession resolver in `_registries.py`.
- The SHACL shapes graph in `shapes/transcriptomics.shapes.ttl`.

Both target:

- `api_version = 1` matching `indicium-adapters.CURRENT_API_VERSION`.
- `DomainAdapter` and `SHACLProvider` Protocol conformance.
- Entry-point registration under `indicium.adapters`.
- Hardcoded namespace URI constants in `_namespaces.py`.
- **Disjoint qualifier sets**, enforced by `tests/test_composition.py`.

The `parsers/` subpackage (v0.7.0) turns repository metadata documents
(GEO SOFT, ENA XML, BioStudies JSON, ISA-Tab) into claim dicts ready
for `enrich_claim`. Each parser detects a recognised organism and
attaches a `species` hint to every emitted claim so the gene resolver
picks the right organism without re-deriving it from the symbol
casing.

## Adding a new qualifier

1. Decide which adapter owns it (or whether it justifies a new sibling).
2. Add it to that adapter's `qualifiers` frozenset.
3. Update `extraction_context()` to mention it.
4. If it needs enrichment behaviour, add a branch in `enrich_claim()`.
5. If it needs validation, add a NodeShape in the adapter's shapes TTL.
6. Add unit tests in `tests/test_<adapter>.py`.
7. Keep qualifier sets disjoint across sibling adapters.

## Refreshing bundled data files

Maintainer-only scripts under `src/indicium_adapters_transcriptomics/data/`.
None of these run at install or import time; they are invoked by hand
when the upstream source changes or when the curated candidate list
grows.

- `build_gene_index.py`: regenerates `gene_index.tsv` from HGNC, MGI,
  FlyBase (NCBI gene_info), WormBase (NCBI gene_info), ZFIN (NCBI
  gene_info), the *Nematostella* G2V2 cross-reference, and UniProt
  per-organism idmapping files. Caches to a `.json` (gitignored).
  254,922 rows as of v0.6.2.
- `build_anatomy_index.py`: regenerates `anatomy_index.tsv` from a
  curated UBERON candidate list, verified against the EBI OLS4
  UBERON index with strict label match (no hyponym promotion).
- `build_concept_index.py`: regenerates `concept_index.tsv` from
  the TranscriptoKG modelling vocabulary, re-verified via OLS4.
- `build_corpus_concept_candidates.py`: OLS4 cascade harvest
  (`exact` then `prefix` then `substring`) over corpus vocab JSONs.
- `build_biotools_index.py`: bio.tools registry harvest with
  direct-ID lookup + normalised name fallback.
- `build_cl_concept_index.py`: Cell Ontology harvest.
- `build_extra_ontology_indexes.py`: REGISTRY-driven harvester for
  the v0.8.0 and v0.9.0 ontology expansion (ChEBI, HP, MP, DOID,
  PR, OBA, ENVO, MS, UO, MOD, FBbt, WBbt, ZFA, EMAPA, XPO).

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
