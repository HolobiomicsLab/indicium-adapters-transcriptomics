"""indicium-adapters-transcriptomics: transcriptomics-domain adapters for the indicium ecosystem.

The package exposes two sibling adapters, both registered under the
``indicium.adapters`` entry-point group:

* ``BulkTranscriptomicsAdapter`` for bulk RNA-seq, microarray and
  time-series expression claims plus study-level assay metadata.
* ``SingleCellTranscriptomicsAdapter`` for cell-type markers, cluster
  enrichment and spatial localisation.

The ``accession_to_curie`` helper from ``_registries`` is also re-exported
for bridge layers that need to normalise GEO / SRA / ENA / ArrayExpress /
ENCODE / BioProject / dbGaP / BioStudies accessions before emitting claims.

See ``REPORT.html`` at the project root for the full development write-up
and ``CHANGELOG.md`` for the release-by-release history.
"""
from __future__ import annotations

from importlib.metadata import version

from indicium_adapters_transcriptomics._registries import accession_to_curie
from indicium_adapters_transcriptomics.bulk import BulkTranscriptomicsAdapter
from indicium_adapters_transcriptomics.singlecell import SingleCellTranscriptomicsAdapter

__version__ = version("indicium-adapters-transcriptomics")
__all__ = [
    "BulkTranscriptomicsAdapter",
    "SingleCellTranscriptomicsAdapter",
    "accession_to_curie",
]
