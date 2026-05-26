# indicium-adapters-transcriptomics

Transcriptomics-domain DomainAdapters for the [indicium](https://github.com/holobiomicslab/indicium) ecosystem.

v0.2.0 ships two sibling adapters:

- **`BulkTranscriptomicsAdapter`** — bulk RNA-seq / microarray / time-series
  expression claims plus study-level assay metadata. Resolves repository
  accessions (GEO / SRA / ENA / ArrayExpress / ENCODE / BioProject / dbGaP /
  BioStudies) and gene names from a bundled comprehensive index (~53,000
  entries) covering the full HGNC human gene catalogue cross-referenced
  with the *Nematostella vectensis* G2V2 mapping (JGI Nemve1 + NvERTx + NVE
  + Nv2g + UniProt). The lookup accepts a gene symbol, a Nemve1 / JGI id,
  or a UniProt accession; HGNC orthology is attached to Nematostella rows
  whenever the UniProt or `Preferred_name → HGNC symbol` join succeeds.
- **`SingleCellTranscriptomicsAdapter`** — cell-resolution claims (cell-type
  markers, cluster enrichment, spatial localisation). Ships a small bundled
  UBERON anatomy index for the spatial-localisation qualifiers
  (`localized_in_tissue`, `localized_in_spatial_domain`), anchored on the
  Cheung et al. 2024 *Nematostella vectensis* Tomo-seq paper. Cell-type
  enrichment via the Cell Ontology (CL) is on the v0.2.1 roadmap.

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
