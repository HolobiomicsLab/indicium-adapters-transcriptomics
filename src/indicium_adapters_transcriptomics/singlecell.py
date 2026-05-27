"""SingleCellTranscriptomicsAdapter — single-cell and spatial transcriptomics adapter.

Covers cell-resolution claims: cell-type markers, cluster-level enrichment,
spatial-domain localisation. Disjoint from ``BulkTranscriptomicsAdapter``:
generic expression assertions and study-level workflow qualifiers live in the
bulk adapter; this one adds the qualifiers that only make sense once you have
per-cell or per-spot resolution.

v0.2.0 introduced UBERON anatomy enrichment for the spatial-localisation
qualifiers; v0.4.0 wires in the shared gene resolver (so a single-cell
claim like ``marker_of_cell_type: WNT4 cnidocyte`` now resolves WNT4 to
HGNC:12783) and adopts the symmetric resolver-dispatcher pattern used
by the bulk adapter. Cell-type (CL ontology) enrichment for the cell-type
qualifiers themselves is still on the v0.5.0 roadmap — requires building
a CL index from the corpus.
"""
from __future__ import annotations

from collections.abc import Callable
from importlib.resources import files
from typing import Final

import rdflib

from indicium_adapters_transcriptomics._enrichment import (
    resolve_anatomy,
    resolve_concept,
    resolve_gene,
)
from indicium_adapters_transcriptomics._namespaces import TRANSCRIPTOMICS_PREFIXES

__all__ = ["SingleCellTranscriptomicsAdapter"]


# Slot names that may carry an anatomical term and therefore benefit from
# UBERON lookup. Object-side mostly — the spatial qualifiers describe a gene
# (subject) localised in a tissue (object).
_ENRICHABLE_SLOTS: Final[tuple[str, ...]] = ("subject", "object")


# --- Resolver dispatch (v0.4.0) --------------------------------------------
#
# Single-cell has no study-level qualifiers (all five qualifiers are
# claim-level / feature assertions), so a single default priority order
# covers every claim. Gene wins on collision-prone names; anatomy and
# concept fall through.

_RESOLVERS: Final[dict[str, Callable[[str], tuple[str, dict[str, str]] | None]]] = {
    "gene":    resolve_gene,
    "anatomy": resolve_anatomy,
    "concept": resolve_concept,
}

_DEFAULT_PRIORITY: Final[tuple[str, ...]] = ("gene", "anatomy", "concept")


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
        """Enrich a single-cell claim with gene / anatomy / concept CURIEs.

        v0.4.0 contract: for each slot the resolver chain runs
        ``gene`` → ``anatomy`` → ``concept`` and the first hit wins.

        Gene first means a slot like ``MMP9`` in
        ``marker_of_cell_type: MMP9 cnidocyte`` resolves to ``HGNC:7176``
        (was unenriched pre-v0.4.0 because single-cell had no gene
        resolver wired). Anatomy still wins for tissue names because
        gene-lookup of ``mesentery`` returns ``None``; concept fires
        last as a fallback for protocol / methodology mentions.

        Cell-type (CL ontology) resolution for the cell-type qualifiers
        themselves is still on the v0.5.0 roadmap.

        Pure function — never mutates the input dict.
        """
        claim = dict(claim)
        if "ontology_terms" in claim:
            claim["ontology_terms"] = dict(claim["ontology_terms"])
        if "xrefs" in claim:
            claim["xrefs"] = dict(claim["xrefs"])

        for slot in _ENRICHABLE_SLOTS:
            name = (claim.get(slot) or "").strip()
            if not name:
                continue
            for resolver_name in _DEFAULT_PRIORITY:
                resolver = _RESOLVERS[resolver_name]
                result = resolver(name)
                if result is None:
                    continue
                curie, xrefs = result
                claim.setdefault("ontology_terms", {})[slot] = curie
                if xrefs:
                    claim.setdefault("xrefs", {})[slot] = xrefs
                break  # first hit per slot wins

        return claim

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
