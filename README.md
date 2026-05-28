# indicium-adapters-transcriptomics

Transcriptomics-domain DomainAdapters for the [indicium](https://github.com/holobiomicslab/indicium) ecosystem.

v0.8.0 ships two sibling adapters + a `parsers/` subpackage + a
multi-ontology concept index (11 source TSVs covering ~246 verified
CURIEs across the transcriptomics-relevant ontology landscape):

- **`BulkTranscriptomicsAdapter`** — bulk RNA-seq / microarray / time-series
  expression claims plus study-level assay metadata. Resolves repository
  accessions (GEO / SRA / ENA / ArrayExpress / ENCODE / BioProject / dbGaP /
  BioStudies) and gene names from a bundled comprehensive
  **multi-species** index (**~255,000 entries** as of v0.6.2):
  - **Human** (HGNC, ~40k) cross-referenced with the
    *Nematostella vectensis* G2V2 mapping (JGI Nemve1 + NvERTx +
    NVE + Nv2g + UniProt).
  - **Mouse** (MGI complete set, ~75k) with Ensembl / Entrez /
    UniProt cross-refs from MGI's report files (v0.6.1).
  - **Drosophila** (FlyBase, ~25k), **C. elegans** (WormBase,
    ~47k), **Zebrafish** (ZFIN, ~55k) from NCBI's per-organism
    gene_info files (v0.6.2).

  Symbol-level collisions (``MAPK1`` ↔ ``Mapk1`` ↔ ``mapk1`` ↔
  ``Mapk``) are disambiguated by an optional ``species`` claim hint
  (``"human"`` / ``"mouse"`` / ``"drosophila"`` / ``"celegans"`` /
  ``"zebrafish"`` / ``"nematostella"``) or a casing heuristic
  (UPPER → human, Title-case → mouse, ``xxx-N`` → C. elegans).
  Direct accession lookups (HGNC, MGI, FlyBase, WormBase, ZFIN,
  Nemve1, UniProt) bypass the disambiguation entirely.
- **`SingleCellTranscriptomicsAdapter`** — cell-resolution claims (cell-type
  markers, cluster enrichment, spatial localisation). Ships a small bundled
  UBERON anatomy index for the spatial-localisation qualifiers
  (`localized_in_tissue`, `localized_in_spatial_domain`), anchored on the
  Cheung et al. 2024 *Nematostella vectensis* Tomo-seq paper. Cell-type
  enrichment via the Cell Ontology (CL) shipped in v0.5.0.

The **`parsers/`** subpackage (v0.7.0) turns repository-metadata
documents into claim dicts ready for `enrich_claim`. Four parsers
share the same `input → list[dict]` contract:

- `parse_geo_soft(text)` — GEO SOFT format
- `parse_ena_xml(xml_text)` — ENA / SRA XML
- `parse_biostudies_json(text_or_obj)` — BioStudies JSON
- `parse_isa_tab(directory)` — ISA-Tab investigation/study/assay bundle

Each detects a recognised organism in the metadata and attaches a
`species` hint to every emitted claim, so the downstream gene
resolver picks the right organism without re-deriving it from the
symbol casing.

Both adapters also share a **concept index** (v0.3.0 → v0.5.0) — 92+
verified CURIEs across four sources:

- **Static (v0.3.0, 31 rows)** — drawn from the TranscriptoKG modelling
  vocabulary: INO regulation processes, NCIT analysis concepts, GO
  roots, OBI investigation / role / protocol / DE-transformation /
  RNA-sequencing, STATO statistical variables. Covers terms like
  *fold change*, *log2*, *up-regulation*, *RNA sequencing assay*,
  *multiple comparison adjustment*.
- **OLS4 corpus-harvested (v0.3.1, 16 rows)** — cascade of three
  match modes against OLS4, re-verified by IRI lookup:
  - *exact* (6): `Illumina`, `MARS-seq`, `Smart-seq2`, `Drop-seq`,
    `Slide-seq`, `t-test`.
  - *prefix* (6, lets descriptive suffixes through): `Wilcoxon` →
    Wilcoxon Signed Rank Test, `bulk RNA-seq` → bulk RNA-seq assay,
    `Visium` → Visium Spatial Gene Expression, `Benjamini-Hochberg`,
    `PacBio`, `Oxford Nanopore`.
  - *substring* (4, instrument families): `NovaSeq` → Illumina
    NovaSeq 6000, `HiSeq`, `NextSeq`, `Sequel`.
- **bio.tools registry (v0.3.1, 9 rows)** — direct-ID + normalised
  name search on https://bio.tools/, canonical registry for
  bioinformatics software: `STAR`, `Bowtie`, `Bowtie2` (matched via
  whitespace normalisation against bio.tools name "Bowtie 2"),
  `Minimap2`, `RSEM`, `StringTie`, `DESeq2`, `Seurat`, `Picard`.
- **Cell Ontology (v0.5.0, 35 rows)** — harvested from the cell-atlas
  papers in the corpus: `cnidocyte` (CL:0011029 — the Nematostella
  iconic cell type), `T cell` / `B cell` / `macrophage` / `monocyte` /
  `neutrophil` / `astrocyte` / `oligodendrocyte` / `fibroblast` /
  `endothelial cell` / `hepatocyte` / `enterocyte` / `goblet cell`
  and 23 more. Powers single-cell qualifiers
  (`marker_of_cell_type`, `enriched_in_cluster`,
  `co_expressed_in_cell_type`) end-to-end.

**v0.4.0 — qualifier-aware resolution priority.** `enrich_claim()`
picks its resolver chain based on the claim's qualifier:

- Study-level qualifiers (`sequenced_on`, `aligned_with`,
  `quantified_with`, `normalized_with`, `tested_with`, `uses_database`)
  use `concept → accession → gene` priority — so a tool name like
  `STAR` in `aligned_with: <study> STAR` resolves to `biotools:star`
  (the aligner) rather than `HGNC:11359` (the steroidogenic acute
  regulatory protein gene with the same symbol).
- Feature qualifiers (`differentially_expressed_in`, `expressed_in`,
  `marker_of`, `regulated_by`, etc.) and `consolidates_dataset` keep
  the default `accession → gene → concept` order — gene mentions
  dominate feature claims, so gene-lookup wins for collision-prone
  names.

Both register under the `indicium.adapters` entry-point group and are
discoverable via `indicium_adapters.discover_adapters()`. Their qualifier sets
are disjoint by construction (enforced by `tests/test_composition.py`).

## Install

```bash
pip install indicium-adapters-transcriptomics
```

`indicium` and `indicium-adapters` are not on PyPI yet; if installing from
source, pull them from git first.

## Usage

```python
from indicium_adapters import compose_adapters, discover_adapters

found = discover_adapters()
bulk = found["bulk-transcriptomics"]
sc   = found["single-cell-transcriptomics"]

# 1) Repository accession → CURIE
enriched = bulk.enrich_claim({
    "subject": "GSE168938",
    "qualifier": "consolidates_dataset",
    "object": "Warner et al. 2025",
})
# enriched["ontology_terms"]["subject"] == "GEO:GSE168938"

# 2) Gene symbol → HGNC CURIE + full xref map
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

# 3) Nematostella Nemve1 / JGI gene-model id → resolved to the human
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

# 4) SingleCell anatomy → UBERON CURIE
enriched = sc.enrich_claim({
    "subject": "MMP9",
    "qualifier": "localized_in_tissue",
    "object": "mesentery",
})
# enriched["ontology_terms"]["object"] == "UBERON:0002095"
# enriched["xrefs"]["object"] == {
#     "uberon": "UBERON:0002095",
#     "label": "mesentery",
#     "description": "Anatomical organ component composed of a double layer ...",
# }

# 5) Analysis-concept term → NCIT / INO / OBI / GO / STATO CURIE
enriched = bulk.enrich_claim({
    "subject": "fold change",
    "qualifier": "tested_with",
    "object": "DESeq2",
})
# enriched["ontology_terms"]["subject"] == "NCIT:C60693"  (Fold Change)
# enriched["ontology_terms"]["object"]  == "biotools:deseq2"

# 6) Same name, different qualifier → different CURIE (v0.4.0)
#    "STAR" is both an aligner (biotools:star) and a gene symbol
#    (HGNC:11359). The qualifier disambiguates.
enriched_aligner = bulk.enrich_claim({
    "subject": "study_42", "qualifier": "aligned_with", "object": "STAR",
})
# enriched_aligner["ontology_terms"]["object"] == "biotools:star"

enriched_gene = bulk.enrich_claim({
    "subject": "STAR", "qualifier": "expressed_in", "object": "mesentery",
})
# enriched_gene["ontology_terms"]["subject"] == "HGNC:11359"
# enriched["xrefs"]["object"] == {
#     "curie": "biotools:deseq2",
#     "label": "DESeq2",
#     "ontology": "biotools",
#     "description": "R/Bioconductor package for differential gene ...",
# }

# Use both adapters at once via the composition helper
composite = compose_adapters([bulk, sc])
```

The 13 bulk-adapter qualifiers cover both claim-level assertions about
features (`differentially_expressed_in`, `co_expressed_with`, `enriched_for`,
`expressed_in`, `marker_of`, `regulated_by`) and study-level assay metadata
(`sequenced_on`, `aligned_with`, `quantified_with`, `normalized_with`,
`tested_with`, `uses_database`, `consolidates_dataset`).

The 5 single-cell-adapter qualifiers cover cell-resolution claims
(`marker_of_cell_type`, `enriched_in_cluster`, `co_expressed_in_cell_type`,
`localized_in_tissue`, `localized_in_spatial_domain`) — assertions that only
make sense once you have per-cell or per-spot resolution.

## SHACL shapes

Both adapters' `shacl_shapes()` return an `rdflib.Graph` with a NodeShape that
validates CURIE format on ASB `Claim` graphs. Pipeline-level RDF validation
for transcriptomics workflow classes (SO / UBERON / CL / EFO shapes) lives
upstream in the consuming graph store, not here.

## Development

```bash
git clone https://github.com/holobiomicslab/indicium-adapters-transcriptomics
cd indicium-adapters-transcriptomics
pip install -e .[dev]
pytest
```

See `CLAUDE.md` for the dev guide and the 4-quadrant filter used to decide
whether a feature belongs in this adapter at all.

## License

MIT — see `LICENSE`.
