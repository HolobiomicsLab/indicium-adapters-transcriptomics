"""SingleCellTranscriptomicsAdapter — single-cell and spatial transcriptomics adapter.

Covers cell-resolution claims: cell-type markers, cluster-level enrichment,
spatial-domain localisation. Disjoint from ``BulkTranscriptomicsAdapter``:
generic expression assertions and study-level workflow qualifiers live in the
bulk adapter; this one adds the qualifiers that only make sense once you have
per-cell or per-spot resolution.

v0.1.0 ships the adapter structure (qualifiers, prefix map, extraction
context) but ``enrich_claim()`` is a no-op pass-through pending the v0.2.0
reference paper that will anchor the CL (Cell Ontology) and UBERON enrichment
lookups.
"""
from __future__ import annotations

from importlib.resources import files

import rdflib

from indicium_adapters_transcriptomics._namespaces import TRANSCRIPTOMICS_PREFIXES

__all__ = ["SingleCellTranscriptomicsAdapter"]


class SingleCellTranscriptomicsAdapter:
    """Single-cell / spatial transcriptomics DomainAdapter (api_version=1).

    Qualifiers are deliberately scoped to assertions that require
    cell-resolution data — bulk expression qualifiers (e.g.
    ``differentially_expressed_in``, ``marker_of``) stay in
    ``BulkTranscriptomicsAdapter`` to preserve disjointness across siblings.

    v0.1.0 is a structural stub: the adapter is discoverable, conforms to the
    Protocol, ships SHACL CURIE-format shapes, but ``enrich_claim()`` returns
    its input unchanged. The CL + UBERON lookup tables land in v0.2.0 once a
    reference single-cell paper anchors the enrichment behaviour.
    """

    domain_id = "single-cell-transcriptomics"
    api_version = 1

    qualifiers: frozenset[str] = frozenset({
        # Cell-resolution qualifiers — disjoint from the bulk set on purpose.
        "marker_of_cell_type",
        "enriched_in_cluster",
        "co_expressed_in_cell_type",
        "localized_in_tissue",
        "localized_in_spatial_domain",
    })

    ontology_prefixes: dict[str, str] = TRANSCRIPTOMICS_PREFIXES

    def extraction_context(self) -> str:
        return (
            "Domain: single-cell / spatial transcriptomics. Extra qualifiers "
            "available: marker_of_cell_type, enriched_in_cluster, "
            "co_expressed_in_cell_type, localized_in_tissue, "
            "localized_in_spatial_domain. For cell-type subjects/objects use "
            "the Cell Ontology (CL) canonical label; for tissue / anatomy use "
            "UBERON. Use bulk-transcriptomics qualifiers for assertions that "
            "do not require per-cell or per-spot resolution."
        )

    def enrich_claim(self, claim: dict) -> dict:
        """No-op pass-through in v0.1.0.

        CL (cell type) + UBERON (anatomy) enrichment lookups land in v0.2.0
        once a reference single-cell paper anchors which terms to bundle.
        Returns a defensive shallow copy so callers can safely mutate the
        result without affecting the input.
        """
        return dict(claim)

    def shacl_shapes(self) -> rdflib.Graph:
        """Return ASB-claim SHACL shapes (CURIE format validation only).

        Shares the transcriptomics shapes TTL with the bulk adapter — both
        adapters emit claims into the same ASB graph and so apply the same
        CURIE-format NodeShape.
        """
        g = rdflib.Graph()
        ttl = files("indicium_adapters_transcriptomics.shapes").joinpath(
            "transcriptomics.shapes.ttl"
        )
        g.parse(ttl.open("rb"), format="turtle")
        return g
