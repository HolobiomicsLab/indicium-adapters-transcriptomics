"""SingleCellTranscriptomicsAdapter — single-cell and spatial transcriptomics adapter.

Covers cell-resolution claims: cell-type markers, cluster-level enrichment,
spatial-domain localisation. Disjoint from ``BulkTranscriptomicsAdapter``:
generic expression assertions and study-level workflow qualifiers live in the
bulk adapter; this one adds the qualifiers that only make sense once you have
per-cell or per-spot resolution.

v0.2.0 ships UBERON anatomy enrichment for the spatial-localisation
qualifiers (``localized_in_tissue``, ``localized_in_spatial_domain``).
Cell-type (CL ontology) enrichment for the cell-type qualifiers is on the
v0.2.1 roadmap.
"""
from __future__ import annotations

from importlib.resources import files
from typing import Final

import rdflib

from indicium_adapters_transcriptomics._enrichment import (
    ANATOMY_INDEX,
    CONCEPT_INDEX,
)
from indicium_adapters_transcriptomics._namespaces import TRANSCRIPTOMICS_PREFIXES

__all__ = ["SingleCellTranscriptomicsAdapter"]


# Slot names that may carry an anatomical term and therefore benefit from
# UBERON lookup. Object-side mostly — the spatial qualifiers describe a gene
# (subject) localised in a tissue (object).
_ENRICHABLE_SLOTS: Final[tuple[str, ...]] = ("subject", "object")


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
        """Enrich a claim with UBERON CURIEs for anatomical terms.

        v0.2.0 behaviour: recognise anatomy names in ``subject`` /
        ``object`` and attach them to ``ontology_terms`` (UBERON CURIE) and
        ``xrefs`` (label + description). The bundled set is small and
        verified against UBERON; cnidarian-specific terms (physa, mesoglea,
        gastrodermis, etc.) are absent because UBERON has no entry for
        them as of 2026 — see ``data/build_anatomy_index.py`` for the
        triaged-out list.

        Cell-type lookup (CL ontology) is on the v0.2.1 roadmap.

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
            if row := ANATOMY_INDEX.lookup(name):
                claim.setdefault("ontology_terms", {})[slot] = row["uberon"]
                # Xrefs: label + description so consumers can render the
                # CURIE without an OLS4 round-trip.
                xrefs_slot = {
                    "uberon": row["uberon"],
                    "label": row["label"],
                }
                if row.get("description"):
                    xrefs_slot["description"] = row["description"]
                claim.setdefault("xrefs", {})[slot] = xrefs_slot
                continue

            # Concept lookup (v0.3.0) — same shared index the bulk adapter
            # uses. Useful when a single-cell claim refers to a study- or
            # protocol-level concept (e.g. `acquired_in_mode` ↔ NCIT).
            if row := CONCEPT_INDEX.lookup(name):
                claim.setdefault("ontology_terms", {})[slot] = row["curie"]
                xref_slot = {
                    "curie": row["curie"],
                    "label": row["label"],
                    "ontology": row["ontology"],
                }
                if row.get("description"):
                    xref_slot["description"] = row["description"]
                claim.setdefault("xrefs", {})[slot] = xref_slot

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
