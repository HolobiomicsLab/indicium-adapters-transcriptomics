# indicium-adapters-transcriptomics

Transcriptomics-domain DomainAdapters for the [indicium](https://github.com/holobiomicslab/indicium) ecosystem.

## What this does

This package lifts free-text transcriptomics mentions in ASB claims to canonical CURIEs.

A claim arrives from the upstream LLM extractor (Perspicacité) as a `(subject, qualifier, object)` triple where subject and object are typically raw strings (`"WNT4"`, `"GSE168938"`, `"STAR"`, `"mesentery"`). This package resolves those strings to ontology CURIEs (`HGNC:12783`, `GEO:GSE168938`, `biotools:star`, `UBERON:0002095`) so downstream RDF storage and SPARQL querying are stable.

The package is a pure-function claim enricher. It does no LLM extraction, makes no network calls at enrichment time, and ships every reference table (gene symbols, anatomy, concepts) as bundled package resources regenerated offline by maintainer scripts.

## When to use which adapter

| Use the bulk adapter when... | Use the single-cell adapter when... |
|---|---|
| The claim is about a gene, transcript, or set at the study level | The claim only makes sense with per-cell or per-spot resolution |
| You're describing assay metadata (platform, aligner, DE tool, accession) | You're describing a cell-type marker, cluster enrichment, or spatial localisation |
| Qualifier is `differentially_expressed_in`, `marker_of`, `sequenced_on`, etc. | Qualifier is `marker_of_cell_type`, `enriched_in_cluster`, `localized_in_tissue`, etc. |

The two qualifier sets are disjoint by construction (enforced by `tests/test_composition.py`), so a composed adapter via `compose_adapters` produces a clean union with no shadowing.

### Deliberate divergence from indicium-adapters #1

The [indicium-adapters #1 spec](https://github.com/HolobiomicsLab/indicium-adapters/issues/1) proposes a single `TranscriptomicsAdapter` with `domain_id = "transcriptomics"`. This package ships two siblings instead:

* `bulk-transcriptomics` for study-level and bulk-feature claims (the dominant case);
* `single-cell-transcriptomics` for assertions that require per-cell or per-spot resolution (`marker_of_cell_type`, `enriched_in_cluster`, `localized_in_tissue`, `localized_in_spatial_domain`).

Rationale: the two qualifier families describe semantically different data sources (a `marker_of_cell_type` claim cannot be made without scRNA / spatial data; a `differentially_expressed_in` claim usually cannot be made FROM scRNA / spatial data without aggregation). Encoding the split as two adapters lets `test_composition.py` enforce qualifier disjointness as a structural invariant — a property a single 18-qualifier frozenset would lose. The metabolomics reference repo also ships three classes (`MetabolomicsAdapter`, `MetaboKGAdapter`, `EMIAdapter`) rather than one, so the spec's singleton entry-point is a convention rather than a hard constraint.

A unifying `TranscriptomicsAdapter` umbrella class that delegates to the two siblings by qualifier may be added later if downstream consumers need the canonical `domain_id == "transcriptomics"` entry-point. For now, use `compose_adapters([bulk, sc])` to get the union behaviour.

## What it resolves

**Gene names** across six organisms (254,922 rows total) from HGNC, MGI, FlyBase, WormBase, ZFIN, plus the *Nematostella vectensis* G2V2 cross-reference (JGI Nemve1 + NvERTx + NVE + Nv2g + UniProt). Symbol-level collisions across organisms (`MAPK1` vs `Mapk1` vs `mapk1`) are disambiguated by an optional `species` claim hint or a casing heuristic (UPPERCASE for human, Title-case for mouse, `xxx-N` for *C. elegans*). Direct accession lookups (HGNC, MGI, FlyBase, WormBase, ZFIN, Nemve1, UniProt) bypass the disambiguation.

**Repository accessions** for GEO, SRA, ENA, ArrayExpress, ENCODE, BioProject, dbGaP, BioStudies.

**Anatomy** via a multi-species index of 152 entries: UBERON (35 species-agnostic terms covering nervous / digestive / respiratory / circulatory / urogenital / endocrine / lymphoid / musculoskeletal systems) plus the organism-paired ontologies FBbt (34 Drosophila), WBbt (21 C. elegans), ZFA (34 zebrafish), EMAPA (28 mouse developmental). Lookup is species-aware: a `species` hint on the claim routes to the organism-specific ontology when available, with UBERON winning by default. Cnidarian-specific terms (*physa*, *mesoglea*, *gastrodermis*) are documented as triaged out because UBERON has no entry for them.

**Concepts** from 21 ontology sources merged into a single index of ~547 verified CURIEs:

| Layer | Sources |
|---|---|
| TranscriptoKG modelling vocabulary | INO, NCIT, GO (3 roots), OBI, STATO |
| Tool registry | bio.tools (STAR, DESeq2, Seurat, minimap2, StringTie, ...) |
| OLS4 corpus harvest | platforms (Illumina, NovaSeq, ...), protocols (Smart-seq2, Visium, ...), stats (Wilcoxon, Benjamini-Hochberg) |
| Cell types | CL (cnidocyte, T cell / B cell subtypes, alveolar macrophage, Kupffer cell, pancreatic A/B/D, podocyte, keratinocyte, melanocyte, neuronal subtypes, stem cells, ... 65 entries) |
| Top-7 OLS4 audit picks | ChEBI, HP, MP, DOID, PR, OBA, ENVO |
| Universal Tier-1 | MS (mass spec), UO (units), MOD (PTMs) |
| Organism-paired Tier-2 | FBbt, WBbt, ZFA, EMAPA, XPO |
| GO biological process / molecular function | ~56 concrete terms (Wnt / MAPK / BMP / Hippo signalling, apoptotic process, cell cycle, DNA repair, gene expression, oxidative phosphorylation, immune response, regeneration, ATP binding, kinase activity, ...) |
| ZP zebrafish phenotype | swim bladder / neural tube / notochord / fin / heart / circulation morphology phenotypes (canonical "<entity> <quality>, abnormal" form) |

The load order encodes a priority ladder (paper-modelling vocab > OLS4 corpus > bio.tools > single-ontology indexes), so a hypothetical term collision resolves to the most grounded source.

## Qualifier-aware resolution

The `enrich_claim()` method picks its resolver chain based on the claim's qualifier. The motivating example is the `STAR` collision: `STAR` is both a popular aligner (bio.tools) and a gene symbol (HGNC:11359). The adapter routes them differently:

| Qualifier class | Resolver priority | Effect on `STAR` |
|---|---|---|
| Study-level (`sequenced_on`, `aligned_with`, `quantified_with`, `normalized_with`, `tested_with`, `uses_database`) | concept > accession > gene | `STAR` resolves to `biotools:star` |
| Feature-level (`differentially_expressed_in`, `upregulated_in`, `downregulated_in`, `expressed_in`, `marker_of`, `marker_for`, `regulated_by`, `enriched_for`, `enriched_in_pathway`, `co_expressed_with`, `consolidates_dataset`) | accession > gene > concept | `STAR` resolves to `HGNC:11359` |

The four qualifiers `upregulated_in`, `downregulated_in`, `marker_for` and `enriched_in_pathway` were added in v0.13 to align with the [indicium-adapters #1](https://github.com/HolobiomicsLab/indicium-adapters/issues/1) contract:

- `upregulated_in` / `downregulated_in` are signed specialisations of `differentially_expressed_in` (use them when log2FC sign is known; fall back to the generic form otherwise).
- `marker_for` is the spec-canonical form alongside the legacy `marker_of` (both stay supported; new claims should prefer `marker_for`).
- `enriched_in_pathway` narrows `enriched_for` for claims whose object resolves into REACTOME / KEGG / WP pathways (use the generic `enriched_for` for GO / phenotype / disease enrichment).

## Metadata parsers

The `parsers/` subpackage turns repository-metadata documents into claim dicts ready for `enrich_claim`. Four parsers share the same `input -> list[dict]` contract:

- `parse_geo_soft(text)` for GEO SOFT format
- `parse_ena_xml(xml_text)` for ENA / SRA XML
- `parse_biostudies_json(text_or_obj)` for BioStudies JSON
- `parse_isa_tab(directory)` for ISA-Tab investigation / study / assay bundles

Each parser detects a recognised organism in the metadata and attaches a `species` hint to every emitted claim, so the gene resolver picks the right organism without re-deriving it from the symbol casing.

## Install

This package is not published on PyPI; the repository is private for now. Install it from a checkout:

```bash
git clone <repository-url>
cd indicium-adapters-transcriptomics
```

The public dependencies are installed from PyPI with the development extra:

```bash
python -m pip install -e ".[dev]"
```

## Usage

```python
from indicium_adapters import compose_adapters, discover_adapters

found = discover_adapters()
bulk = found["bulk-transcriptomics"]
sc   = found["single-cell-transcriptomics"]

# 1) Repository accession -> CURIE
enriched = bulk.enrich_claim({
    "subject": "GSE168938",
    "qualifier": "consolidates_dataset",
    "object": "Warner et al. 2025",
})
# enriched["ontology_terms"]["subject"] == "GEO:GSE168938"

# 2) Gene symbol -> HGNC CURIE + full xref map
enriched = bulk.enrich_claim({
    "subject": "WNT4",
    "qualifier": "regulated_by",
    "object": "cWnt pathway",
})
# enriched["ontology_terms"]["subject"] == "HGNC:12783"
# enriched["xrefs"]["subject"] == {
#     "hgnc": "HGNC:12783", "ensembl": "ENSG00000162552",
#     "ncbigene": "54361", "uniprot": "P56705",
# }

# 3) Nematostella Nemve1 / JGI gene-model id resolves to the human
#    ortholog via the G2V2 cross-reference and the HGNC symbol fallback
enriched = bulk.enrich_claim({
    "subject": "gw.51.96.1",       # FGF8 in the Nemve1 v1 assembly
    "qualifier": "expressed_in",
    "object": "physa",
})
# enriched["ontology_terms"]["subject"] == "HGNC:3686"  (FGF8)
# enriched["xrefs"]["subject"] == {
#     "hgnc": "HGNC:3686", "ensembl": "ENSG00000107831",
#     "ncbigene": "2253", "nemve1": "gw.51.96.1",
#     "nvertx": "NvERTx.4.55830", "uniprot": "Q805B2",
# }

# 4) Single-cell anatomy -> UBERON or organism-specific CURIE
enriched = sc.enrich_claim({
    "subject": "MMP9",
    "qualifier": "localized_in_tissue",
    "object": "mesentery",
})
# enriched["ontology_terms"]["object"] == "UBERON:0002095"
# enriched["xrefs"]["object"] == {
#     "curie": "UBERON:0002095",
#     "label": "mesentery",
#     "ontology": "UBERON",
#     "description": "Anatomical organ component composed of a double layer ...",
# }

# 4b) Species hint routes organism-specific anatomy to the right ontology
enriched = sc.enrich_claim({
    "subject": "dpp",
    "qualifier": "localized_in_tissue",
    "object": "wing disc",
    "species": "drosophila",
})
# enriched["ontology_terms"]["object"] == "FBbt:00001778"
# enriched["xrefs"]["object"]["ontology"] == "FBbt"
# enriched["xrefs"]["object"]["species"]  == "drosophila"

# 5) Analysis-concept term -> NCIT / INO / OBI / GO / STATO CURIE
enriched = bulk.enrich_claim({
    "subject": "fold change",
    "qualifier": "tested_with",
    "object": "DESeq2",
})
# enriched["ontology_terms"]["subject"] == "NCIT:C60693"  (Fold Change)
# enriched["ontology_terms"]["object"]  == "biotools:deseq2"

# 6) Same name, different qualifier -> different CURIE.
#    "STAR" is both an aligner (biotools:star) and a gene (HGNC:11359);
#    the qualifier disambiguates.
enriched_aligner = bulk.enrich_claim({
    "subject": "study_42", "qualifier": "aligned_with", "object": "STAR",
})
# enriched_aligner["ontology_terms"]["object"] == "biotools:star"

enriched_gene = bulk.enrich_claim({
    "subject": "STAR", "qualifier": "expressed_in", "object": "mesentery",
})
# enriched_gene["ontology_terms"]["subject"] == "HGNC:11359"

# Use both adapters at once via the composition helper
composite = compose_adapters([bulk, sc])
```

The 17 bulk-adapter qualifiers cover claim-level feature assertions (`differentially_expressed_in`, `upregulated_in`, `downregulated_in`, `co_expressed_with`, `enriched_for`, `enriched_in_pathway`, `expressed_in`, `marker_of`, `marker_for`, `regulated_by`) and study-level assay metadata (`sequenced_on`, `aligned_with`, `quantified_with`, `normalized_with`, `tested_with`, `uses_database`, `consolidates_dataset`).

The 5 single-cell-adapter qualifiers (`marker_of_cell_type`, `enriched_in_cluster`, `co_expressed_in_cell_type`, `localized_in_tissue`, `localized_in_spatial_domain`) cover assertions that only make sense at per-cell or per-spot resolution.

## SSSOMMapping emission and `default_workflow_coupling` (v0.14)

Every resolver hit in `enrich_claim()` now produces a first-class [SSSOMMapping](https://w3id.org/sssom/) fragment in `claim["sssom_mappings"]` (dict keyed by mapping id). Each fragment matches `indicium.generated.model.SSSOMMapping`'s pydantic constructor kwargs (indicium ≥ 1.11.0), so downstream consumers can fold the entries straight into `IndiciumDocument.sssom_mappings` without round-tripping.

The fragment uses the synthetic `asb.lex:` namespace for the subject (a slug of the raw text) so free-text LLM mentions stay deduplicable across a document. The predicate routes via SSSOM convention:

| Resolver kind | Predicate | Confidence floor | Notes |
|---|---|---|---|
| `accession` | `skos:exactMatch` | 0.99 | regex-anchored repository accession or pathway ID |
| `gene` | `skos:exactMatch` | 0.95 | HGNC / MGI / FlyBase / WormBase / ZFIN / Ensembl symbol lookup |
| `anatomy` | `skos:exactMatch` | 0.95 | UBERON / FBbt / WBbt / ZFA / EMAPA |
| `concept` | `skos:exactMatch` | 0.90 | bio.tools, NCIT, INO, OBI, STATO, GO, CL, ... |
| `embedded_pathway` | `skos:closeMatch` | 0.85 | substring fallback for IDs wrapped in longer text ("Reactome ID R-HSA-..."). |
| `trimmed_concept` | `skos:closeMatch` | 0.80 | trim-tail fallback for "cnidocyte cells" → "cnidocyte" |

`skos:exactMatch` means the raw text is interchangeable with the CURIE (canonical identifier); `skos:closeMatch` is used for the two LLM-tolerant fallback paths where we modified the input or extracted the id from a longer wrapping string. Downstream SPARQL queries can filter on `confidence >= 0.90` or `predicate = skos:exactMatch` to drop heuristic hits.

Both adapters also expose `default_workflow_coupling = "loose"` as a class attribute (indicium ≥ 1.10.0). Bulk RNA-seq (DESeq2 / edgeR / Salmon / STAR) and single-cell (CellRanger / STARsolo / Seurat / Scanpy) pipelines are file-based and individually reproducible per Suter et al. 2025, so the spec recommends `"loose"`. Perspicacite reads this attribute when constructing `Evidence` nodes for transcriptomics claims.

## SHACL shapes

Both adapters' `shacl_shapes()` return an `rdflib.Graph` with a NodeShape that validates CURIE format on ASB `Claim` graphs. Pipeline-level RDF validation for transcriptomics workflow classes (SO, UBERON, CL, EFO shapes) lives upstream in the consuming graph store, not here.

## Validation

The package ships with two runnable demonstrations under `examples/`:

- `examples/quickstart.py`: 30-line bring-up smoke test (discover, enrich one claim per adapter, validate SHACL shapes parse).
- `examples/test_on_papers.py`: end-to-end test on a 16-paper hand-curated corpus (45 claims spanning six organisms across bulk, single-cell, spatial and long-read modalities, plus four methodology reviews). 44/45 enriched (98%), 18 distinct ontologies firing, SHACL conformance.

`REPORT.html` at the project root documents the full development history (architecture, data sources, resolution methods, release timeline, results, limitations).

## Development

```bash
git clone https://github.com/holobiomicslab/indicium-adapters-transcriptomics
cd indicium-adapters-transcriptomics
pytest
```

See `CLAUDE.md` for the developer guide and the 4-quadrant filter that decides whether a feature belongs in this adapter at all. See `CHANGELOG.md` for the release-by-release history.

## License

MIT, see `LICENSE`.
