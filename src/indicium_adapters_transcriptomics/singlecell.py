"""SingleCellTranscriptomicsAdapter: single-cell and spatial transcriptomics adapter.

Covers cell-resolution claims: cell-type markers, cluster-level enrichment,
spatial-domain localisation. The qualifier set is disjoint from the bulk
adapter on purpose. Generic expression assertions and study-level workflow
metadata stay on the bulk side; this adapter only adds claims that need
per-cell or per-spot resolution to make sense.

The adapter is a pure-function claim enricher (no LLM extraction, no network
calls, no file writes) and shares the same resolver back-end as the bulk
adapter, so gene-name resolution stays consistent across the two siblings.
"""
from __future__ import annotations

from collections.abc import Callable
from importlib.resources import files
from typing import Final

import rdflib

from indicium_adapters_transcriptomics._enrichment import (
    build_sssom_mapping,
    resolve_anatomy,
    resolve_concept,
    resolve_gene,
)
from indicium_adapters_transcriptomics._namespaces import TRANSCRIPTOMICS_PREFIXES

__all__ = ["SingleCellTranscriptomicsAdapter"]


# Slot names checked by ``enrich_claim``. Single-cell claims usually have a
# gene as subject and an anatomical structure or cell type as object, but
# both slots are inspected so a claim with a gene in the object slot still
# resolves.
_ENRICHABLE_SLOTS: Final[tuple[str, ...]] = ("subject", "object")


# --- Resolver dispatch -----------------------------------------------------
#
# Single-cell has no study-level qualifiers; all five qualifiers are
# claim-level feature assertions, so a single default priority order covers
# every claim. Gene wins on collision-prone names; anatomy and concept fall
# through after it.

_RESOLVERS: Final[dict[str, Callable[[str], tuple[str, dict[str, str]] | None]]] = {
    "gene":    resolve_gene,
    "anatomy": resolve_anatomy,
    "concept": resolve_concept,
}

_DEFAULT_PRIORITY: Final[tuple[str, ...]] = ("gene", "anatomy", "concept")


class SingleCellTranscriptomicsAdapter:
    """Single-cell and spatial transcriptomics DomainAdapter (api_version=1).

    Qualifiers are scoped to assertions that require cell-resolution data.
    Bulk expression qualifiers such as ``differentially_expressed_in`` and
    ``marker_of`` stay on ``BulkTranscriptomicsAdapter`` to keep the two
    sibling qualifier sets disjoint.

    Every qualifier in the ``qualifiers`` frozenset MUST match
    ``^[a-z][a-z0-9_]*$`` (snake_case). Required by indicium >= 1.7.0:
    Perspicacite routes adapter qualifiers to the open
    ``asb:domainQualifier`` slot, whose SHACL shape enforces exactly this
    pattern. Non-conforming values would be silently rejected by
    ``indicium.validate_graph`` at the RDF boundary. The constraint is
    enforced by ``tests/test_protocol_conformance.py``.
    """

    domain_id = "single-cell-transcriptomics"
    api_version = 1
    # Single-cell pipelines (CellRanger / STARsolo / Seurat / Scanpy) are
    # file-based and independently reproducible — same recommendation as
    # the bulk sibling. Set by the spec on Evidence nodes derived from
    # these claims (indicium >= 1.10.0, Suter et al. 2025).
    default_workflow_coupling: str = "loose"

    # Every qualifier MUST match ``^[a-z][a-z0-9_]*$`` (snake_case).
    # Required by indicium >= 1.7.0: at serialisation time, Perspicacite
    # routes any qualifier from this set to ``asb:domainQualifier`` (the
    # open slot for domain vocabularies), which enforces the snake_case
    # pattern at the SHACL boundary. Non-conforming values would be
    # rejected by ``indicium.validate_graph``. Enforced for the package
    # by ``tests/test_protocol_conformance.py``.
    qualifiers: frozenset[str] = frozenset({
        # Cell-resolution qualifiers, disjoint from the bulk set.
        "marker_of_cell_type",
        "enriched_in_cluster",
        "co_expressed_in_cell_type",
        "localized_in_tissue",
        "localized_in_spatial_domain",
    })

    ontology_prefixes: dict[str, str] = TRANSCRIPTOMICS_PREFIXES

    def extraction_context(self) -> str:
        """Prompt fragment fed to the upstream LLM extractor.

        Lists the single-cell qualifiers and tells the extractor which
        ontologies to prefer for cell-type and tissue mentions.
        """
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

        For each slot the resolver chain runs gene then anatomy then
        concept, and the first hit wins.

        * Gene runs first so a slot like ``MMP9`` in
          ``marker_of_cell_type: MMP9 cnidocyte`` resolves to ``HGNC:7176``.
        * Anatomy still wins for tissue names because gene lookup of
          ``mesentery`` returns ``None``.
        * Concept fires last as a fallback for protocol or methodology
          mentions and (since v0.5.0) cell-type names that match CL labels.

        Pure function: never mutates the input dict.
        """
        # Shallow-copy the claim and its nested dicts so the caller's
        # input is preserved.
        claim = dict(claim)
        if "ontology_terms" in claim:
            claim["ontology_terms"] = dict(claim["ontology_terms"])
        if "ontology_term_justifications" in claim:
            claim["ontology_term_justifications"] = dict(claim["ontology_term_justifications"])
        if "xrefs" in claim:
            claim["xrefs"] = dict(claim["xrefs"])
        if "sssom_mappings" in claim:
            claim["sssom_mappings"] = dict(claim["sssom_mappings"])

        # Optional species hint forwarded to every resolver. Non-gene
        # resolvers ignore it.
        species_hint = claim.get("species")
        for slot in _ENRICHABLE_SLOTS:
            name = (claim.get(slot) or "").strip()
            if not name:
                continue
            for resolver_name in _DEFAULT_PRIORITY:
                resolver = _RESOLVERS[resolver_name]
                result = resolver(name, species=species_hint)
                if result is None:
                    continue
                curie, xrefs = result
                claim.setdefault("ontology_terms", {})[slot] = curie
                justification = (
                    "semapv:CompositeMatching"
                    if resolver_name == "gene"
                    else "semapv:LexicalMatching"
                )
                claim.setdefault("ontology_term_justifications", {})[slot] = justification
                if xrefs:
                    claim.setdefault("xrefs", {})[slot] = xrefs
                # Emit first-class SSSOMMapping fragment (indicium >= 1.11.0),
                # mirroring the bulk adapter. Folded into
                # IndiciumDocument.sssom_mappings downstream by Perspicacite.
                fragment = build_sssom_mapping(
                    raw_text=name,
                    object_curie=curie,
                    object_label=xrefs.get("label") if isinstance(xrefs, dict) else None,
                    justification=justification,
                    resolver_kind=resolver_name,
                )
                claim.setdefault("sssom_mappings", {})[fragment["id"]] = fragment
                break  # First hit per slot wins.

        return claim

    def shacl_shapes(self) -> rdflib.Graph:
        """Return ASB-claim SHACL shapes (CURIE-format validation only).

        Shares the transcriptomics shapes TTL with the bulk adapter: both
        adapters emit claims into the same ASB graph and apply the same
        CURIE-format NodeShape.
        """
        g = rdflib.Graph()
        ttl = files("indicium_adapters_transcriptomics.shapes").joinpath(
            "transcriptomics.shapes.ttl"
        )
        g.parse(ttl.open("rb"), format="turtle")
        return g
