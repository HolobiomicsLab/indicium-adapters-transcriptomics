"""indicium-adapters-transcriptomics — transcriptomics-domain adapters for the indicium ecosystem.

v0.1.0 ships two sibling adapters, both registered under the
``indicium.adapters`` entry-point group:

- ``BulkTranscriptomicsAdapter`` — bulk RNA-seq / microarray / time-series
  expression claims plus study-level assay metadata. Repository-accession
  enrichment in v0.1.0; gene-symbol enrichment in v0.1.1.
- ``SingleCellTranscriptomicsAdapter`` — cell-type markers, cluster
  enrichment, spatial localisation. Structure-only in v0.1.0;
  ``enrich_claim()`` becomes a real CL + UBERON lookup in v0.2.0 once a
  reference single-cell paper anchors which terms to bundle.

The ``accession_to_curie`` helper from ``_registries`` is re-exported for
bridge layers that want to normalise GEO / SRA / ENA / ArrayExpress / ENCODE /
BioProject / dbGaP / BioStudies accessions before emitting claims.
"""
from __future__ import annotations

from indicium_adapters_transcriptomics._registries import accession_to_curie
from indicium_adapters_transcriptomics.bulk import BulkTranscriptomicsAdapter
from indicium_adapters_transcriptomics.singlecell import SingleCellTranscriptomicsAdapter

__version__ = "0.3.0"
__all__ = [
    "BulkTranscriptomicsAdapter",
    "SingleCellTranscriptomicsAdapter",
    "accession_to_curie",
]
