# indicium-adapters-transcriptomics

Transcriptomics-domain DomainAdapters for the [indicium](https://github.com/holobiomicslab/indicium) ecosystem.

v0.1.0 ships two sibling adapters:

- **`BulkTranscriptomicsAdapter`** — bulk RNA-seq / microarray / time-series
  expression claims plus study-level assay metadata. Repository-accession
  enrichment (GEO / SRA / ENA / ArrayExpress / ENCODE / BioProject / dbGaP /
  BioStudies); gene-symbol enrichment (HGNC + Nemve1) is on the v0.1.1
  roadmap.
- **`SingleCellTranscriptomicsAdapter`** — cell-resolution claims (cell-type
  markers, cluster enrichment, spatial localisation). Structure-only in
  v0.1.0 — `enrich_claim()` becomes a real CL (Cell Ontology) + UBERON lookup
  in v0.2.0 once a reference single-cell paper anchors which terms to bundle.

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

enriched = bulk.enrich_claim({
    "subject": "GSE168938",
    "qualifier": "consolidates_dataset",
    "object": "Warner et al. 2025",
})
# enriched["ontology_terms"]["subject"] == "GEO:GSE168938"

# Use both at once
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
