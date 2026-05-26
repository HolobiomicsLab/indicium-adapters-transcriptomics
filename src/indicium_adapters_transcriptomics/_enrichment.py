"""Shared enrichment back-end for transcriptomics adapters.

In v0.1.0 this module is a placeholder. The bundled gene-index lookup (HGNC +
Nematostella Nemve1 in v0.1.0, expanded in v0.2.0) will be loaded here as a
module-level ``GENE_INDEX`` and consumed by adapter ``enrich_claim()`` methods.
"""
from __future__ import annotations

# Bundled gene name → CURIE lookup is on the v0.1.1 roadmap. The placeholder
# lookup below keeps the adapter wiring testable while the build pipeline for
# the index TSV is being put in place.


class _EmptyGeneIndex:
    """No-op gene index returned by `enrich_claim()` until v0.1.1 ships the
    bundled HGNC + Nemve1 TSV."""

    def lookup(self, name: str) -> dict | None:
        return None


GENE_INDEX = _EmptyGeneIndex()
