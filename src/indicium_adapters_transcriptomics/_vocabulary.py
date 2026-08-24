"""Classification vocabulary declared by the transcriptomics adapters.

ASB's classification layer reads ``primary_domain`` / ``subdomains`` /
``techniques`` off every installed adapter and unions the declared terms into
its controlled vocabulary
(``agentic_science_builder.classification._adapter_vocab``). An adapter that
declares no ``primary_domain`` contributes no domain at all, so a
transcriptomics paper would classify as "other" and carry no technique.

The two exported tuples are what
``BulkTranscriptomicsAdapter.subdomains`` / ``.techniques`` and
``SingleCellTranscriptomicsAdapter.subdomains`` / ``.techniques`` expose.
Both siblings expose the same tuples on purpose: the vocabulary describes the
field, not the assay resolution, and the consumer unions the siblings into one
set. Unlike ``qualifiers``, these are deliberately not disjoint.
"""
from __future__ import annotations

from typing import Final

# From the 2026-08-25 cross-domain vocabulary audit: a term two sciences
# both declare cannot tell them apart, so shared terms were dropped.
TRANSCRIPTOMICS_SUBDOMAINS: Final[tuple[str, ...]] = (
    "gene-expression",
    "alternative-splicing",
    "single-cell-transcriptomics",
    "spatial-transcriptomics",
    "non-coding-rna-expression-profiling",
)

# From the 2026-08-25 cross-domain vocabulary audit: a term two sciences
# both declare cannot tell them apart, so shared terms were dropped.
TRANSCRIPTOMICS_TECHNIQUES: Final[tuple[str, ...]] = (
    "rna-seq",
    "single-cell-rna-seq",
    "long-read-rna-sequencing",
    "small-rna-seq",
    "rt-qpcr",
    "expression-microarray",
    "transcript-quantification",
    "differential-expression-analysis",
    "transcriptome-assembly",
    "gene-co-expression-network-analysis",
    "cell-type-annotation",
    "rna-velocity",
    "trajectory-inference",
)
