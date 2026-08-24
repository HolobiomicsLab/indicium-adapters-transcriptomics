"""BulkTranscriptomicsAdapter: bulk RNA-seq, microarray and time-series adapter.

This adapter handles generic bulk-RNA-seq, microarray and time-series
expression claims plus study-level metadata (sequencing platform, alignment,
normalisation, differential-expression tool). Single-cell-specific qualifiers
(cell-type markers, cluster enrichment, spatial localisation) live in
``SingleCellTranscriptomicsAdapter`` to keep the two qualifier sets disjoint.

The adapter is a pure-function claim enricher. It does no LLM extraction,
makes no network calls at enrichment time, and does not write files. All
domain-specific data (gene lookups, ontology vocab subsets) is shipped as
bundled package resources under ``data/`` and regenerated offline by
maintainer scripts under ``data/build_*.py``.
"""
from __future__ import annotations

from collections.abc import Callable
from importlib.resources import files
from typing import Final

import rdflib

from indicium_adapters_transcriptomics._enrichment import (
    build_sssom_mapping,
    resolve_concept,
    resolve_gene,
)
from indicium_adapters_transcriptomics._namespaces import TRANSCRIPTOMICS_PREFIXES
from indicium_adapters_transcriptomics._registries import (
    accession_to_curie,
    find_embedded_pathway_curie,
)
from indicium_adapters_transcriptomics._vocabulary import (
    TRANSCRIPTOMICS_SUBDOMAINS,
    TRANSCRIPTOMICS_TECHNIQUES,
)

__all__ = ["BulkTranscriptomicsAdapter", "accession_to_curie"]


# Slot names checked by ``enrich_claim``. Kept narrow on purpose: we only try
# to resolve the subject and object of each claim, never the qualifier or
# other metadata fields.
_ENRICHABLE_SLOTS: Final[tuple[str, ...]] = ("subject", "object")


# --- Resolver dispatch -----------------------------------------------------
#
# Each resolver inspects a slot value and returns ``(curie, xrefs)`` on a hit
# or ``None`` on a miss. The gene and concept resolvers are shared with the
# single-cell adapter and live in ``_enrichment``. The accession resolver is
# bulk-only because single-cell claims do not carry repository accessions in
# the same way.


def _resolve_accession(
    name: str, *, species: str | None = None
) -> tuple[str, dict[str, str]] | None:
    """Repository / pathway-DB accession resolver.

    Signature matches the gene and concept resolvers so the dispatcher can
    treat all three uniformly. ``species`` is ignored: repository
    accessions are species-agnostic.

    Resolution is two-pass:

    1. Strict anchored match via ``accession_to_curie`` (the common case:
       ``"GSE168938"`` -> ``"GEO:GSE168938"``).
    2. Embedded pathway-DB scan via ``find_embedded_pathway_curie`` for
       Reactome / KEGG / WP IDs wrapped in longer LLM-emitted strings
       (e.g. ``"Wnt pathway (Reactome ID R-HSA-9012999)"`` ->
       ``"REACTOME:R-HSA-9012999"``). Only fires if step 1 misses, so
       legitimate inputs that already match in step 1 keep their exact
       previous behaviour.

    The returned xrefs dict carries an internal ``"_resolver_kind"``
    sentinel (``"accession"`` for path 1, ``"embedded_pathway"`` for
    path 2). The dispatcher pops it before storing xrefs on the claim,
    so downstream consumers never see the sentinel; it carries the
    kind into ``build_sssom_mapping`` so fallback hits get the right
    confidence floor and ``skos:closeMatch`` predicate.
    """
    del species
    curie = accession_to_curie(name)
    if curie is not None:
        return curie, {"_resolver_kind": "accession"}
    curie = find_embedded_pathway_curie(name)
    if curie is None:
        return None
    return curie, {"_resolver_kind": "embedded_pathway"}


_Resolver = Callable[..., tuple[str, dict[str, str]] | None]

# Dispatch table. ``enrich_claim`` looks resolvers up by name in this dict,
# in the order given by the active priority tuple below.
_RESOLVERS: Final[dict[str, _Resolver]] = {
    "accession": _resolve_accession,
    "gene":      resolve_gene,
    "concept":   resolve_concept,
}


# Default resolver priority for claim-level (feature) qualifiers. The
# subject of a feature claim is usually a gene mention, so gene-lookup wins
# over the broader concept index. ``accession`` runs first because
# repository accessions are unambiguous and should never be re-routed
# through a gene-symbol collision.
_DEFAULT_PRIORITY: Final[tuple[str, ...]] = ("accession", "gene", "concept")


# Study-level priority for metadata qualifiers, whose object is typically a
# tool / platform / database mention. ``concept`` (which contains bio.tools)
# runs first so that, for example, ``STAR`` in ``aligned_with`` resolves to
# ``biotools:star`` (the aligner) rather than ``HGNC:11359`` (the
# steroidogenic acute regulatory protein gene with the same symbol).
# ``accession`` still runs early so a GEO id in any slot stays a GEO CURIE;
# ``gene`` falls last to catch corpus mentions that are not known tools.
_STUDY_LEVEL_PRIORITY: Final[tuple[str, ...]] = ("concept", "accession", "gene")


# Qualifier-to-priority map. Listed qualifiers get the study-level order;
# every other qualifier (including ``consolidates_dataset``, whose object
# is a repository accession that the accession resolver handles first
# under either order) falls back to the default.
_QUALIFIER_PRIORITY: Final[dict[str, tuple[str, ...]]] = {
    "sequenced_on":     _STUDY_LEVEL_PRIORITY,
    "aligned_with":     _STUDY_LEVEL_PRIORITY,
    "quantified_with":  _STUDY_LEVEL_PRIORITY,
    "normalized_with":  _STUDY_LEVEL_PRIORITY,
    "tested_with":      _STUDY_LEVEL_PRIORITY,
    "uses_database":    _STUDY_LEVEL_PRIORITY,
}


class BulkTranscriptomicsAdapter:
    """Bulk transcriptomics DomainAdapter (api_version=1).

    Qualifiers cover two axes:

    * Claim-level assertions about genes / transcripts / sets:
      ``differentially_expressed_in``, ``co_expressed_with``,
      ``enriched_for``, ``expressed_in``, ``marker_of``, ``regulated_by``.
    * Study-level assay-metadata assertions: ``sequenced_on``,
      ``aligned_with``, ``quantified_with``, ``normalized_with``,
      ``tested_with``, ``uses_database``, ``consolidates_dataset``.

    Every qualifier in the ``qualifiers`` frozenset MUST match
    ``^[a-z][a-z0-9_]*$`` (snake_case). Required by indicium >= 1.7.0:
    Perspicacite routes adapter qualifiers to the open
    ``asb:domainQualifier`` slot, whose SHACL shape enforces exactly this
    pattern. Non-conforming values would be silently rejected by
    ``indicium.validate_graph`` at the RDF boundary. The constraint is
    enforced by ``tests/test_protocol_conformance.py``.

    Single-cell qualifiers (cell-type markers, cluster enrichment, spatial
    localisation) belong in ``SingleCellTranscriptomicsAdapter`` so the two
    sibling adapters stay disjoint.

    The ``default_workflow_coupling`` class attribute advertises the
    Evidence.workflow_coupling value Perspicacite should set on Evidence
    nodes derived from bulk-transcriptomics claims (indicium >= 1.10.0,
    Suter et al. 2025). Bulk RNA-seq pipelines (DESeq2 / edgeR / Salmon /
    STAR) are file-based and independently reproducible, so the spec
    recommends ``"loose"``. Consumers that know better at runtime can
    override.
    """

    domain_id = "bulk-transcriptomics"
    api_version = 1
    default_workflow_coupling: str = "loose"

    # ASB's classification layer reads ``primary_domain`` in preference to
    # ``domain_id``; it is what groups this adapter with its single-cell
    # sibling into one field. Both siblings declare the same package-level
    # vocabulary, which the consumer unions (see ``_vocabulary``).
    primary_domain: str = "transcriptomics"
    # From the 2026-08-25 cross-domain vocabulary audit: a term two sciences
    # both declare cannot tell them apart, so shared terms were dropped.
    subdomains: tuple[str, ...] = TRANSCRIPTOMICS_SUBDOMAINS
    techniques: tuple[str, ...] = TRANSCRIPTOMICS_TECHNIQUES

    # Every qualifier MUST match ``^[a-z][a-z0-9_]*$`` (snake_case).
    # Required by indicium >= 1.7.0: at serialisation time, Perspicacite
    # routes any qualifier from this set to ``asb:domainQualifier`` (the
    # open slot for domain vocabularies), which enforces the snake_case
    # pattern at the SHACL boundary. Non-conforming values would be
    # rejected by ``indicium.validate_graph``. Enforced for the package
    # by ``tests/test_protocol_conformance.py``.
    qualifiers: frozenset[str] = frozenset({
        # Claim-level assertions about features (genes / transcripts / sets).
        "differentially_expressed_in",
        "co_expressed_with",
        "enriched_for",
        "expressed_in",
        "marker_of",
        "regulated_by",
        # Spec-aligned (indicium-adapters #1) directional + pathway-specific
        # variants. ``upregulated_in`` / ``downregulated_in`` are signed
        # specialisations of ``differentially_expressed_in`` (extractor
        # picks the signed form when log2FC sign is known, otherwise
        # falls back to the generic ``differentially_expressed_in``).
        # ``marker_for`` is the spec-canonical naming alongside the
        # legacy ``marker_of`` (both stay supported; new claims should
        # prefer ``marker_for``). ``enriched_in_pathway`` narrows
        # ``enriched_for`` for claims whose object resolves into
        # REACTOME / KEGG / WP, with the generic ``enriched_for`` kept
        # for GO / phenotype / disease enrichment.
        "upregulated_in",
        "downregulated_in",
        "marker_for",
        "enriched_in_pathway",
        # Study-level assay-metadata assertions. Subject is typically a
        # study or capsule identifier; object is the platform / tool /
        # database literal. Mirrors the metabolomics adapter pattern.
        "sequenced_on",
        "aligned_with",
        "quantified_with",
        "normalized_with",
        "tested_with",
        "uses_database",
        "consolidates_dataset",
    })

    ontology_prefixes: dict[str, str] = TRANSCRIPTOMICS_PREFIXES

    def extraction_context(self) -> str:
        """Prompt fragment fed to the upstream LLM extractor.

        Lists every qualifier this adapter handles and gives canonical
        examples of subject / object shapes per qualifier family so the
        extractor can route claims to the right adapter.
        """
        return (
            "Domain: bulk transcriptomics. Extra qualifiers available: "
            "differentially_expressed_in, upregulated_in, downregulated_in, "
            "co_expressed_with, enriched_for, enriched_in_pathway, "
            "expressed_in, marker_of, marker_for, regulated_by, "
            "sequenced_on, aligned_with, quantified_with, normalized_with, "
            "tested_with, uses_database, consolidates_dataset. "
            "Use canonical gene symbols (HGNC for human, MGI for mouse, "
            "official model-organism database symbols otherwise) for gene "
            "subjects and objects where possible. "
            "Direction-of-regulation: prefer the signed specialisations "
            "``upregulated_in`` / ``downregulated_in`` when the source "
            "reports a clear direction (positive / negative log2FC, "
            "explicit 'increased' / 'decreased' wording); fall back to "
            "``differentially_expressed_in`` when direction is unknown "
            "or mixed. The signed forms are spec-aligned with the "
            "indicium-adapters contract; both shapes serialise to the "
            "same RDF class on the consuming graph store. "
            "Pathway-vs-general enrichment: use ``enriched_in_pathway`` "
            "for enrichment results whose object is a REACTOME / KEGG / "
            "WP pathway (e.g. 'REACTOME:R-HSA-9012999', 'KEGG:hsa04310', "
            "'WP:WP554'); keep ``enriched_for`` for GO terms, phenotype "
            "ontology hits, or disease-set enrichment. The enricher "
            "auto-lifts bare pathway identifiers ('R-HSA-9012999', "
            "'hsa04310', 'WP554') to the prefixed CURIE form, so emitting "
            "the raw id is fine when the source paper does not pre-prefix. "
            "Marker assertions: ``marker_for`` is the spec-canonical "
            "form (e.g. 'CD3E marker_for T cell'); ``marker_of`` is the "
            "legacy synonym kept for backward compatibility with claims "
            "emitted before the spec alignment. Single-cell claims use "
            "``marker_of_cell_type`` instead. "
            "Study-level qualifiers "
            "(sequenced_on, aligned_with, quantified_with, normalized_with, "
            "tested_with, uses_database, consolidates_dataset) take a "
            "study / capsule / paper id as subject and a platform / tool / "
            "database / accession literal as object (e.g. 'Illumina HiSeq "
            "2500', 'STAR', 'HTSeq', 'TMM', 'DESeq2', 'Ensembl', "
            "'GEO:GSE168938'). Single-cell-specific claims (cell-type "
            "markers, cluster enrichment) belong in the single-cell adapter."
        )

    def enrich_claim(self, claim: dict) -> dict:
        """Enrich a claim with transcriptomics CURIEs.

        For each slot (subject, object), the resolver chain is chosen by
        the claim's ``qualifier``:

        * Study-level qualifiers (``sequenced_on``, ``aligned_with``,
          ``quantified_with``, ``normalized_with``, ``tested_with``,
          ``uses_database``) use ``concept`` then ``accession`` then
          ``gene``, so a tool name like ``STAR`` resolves to
          ``biotools:star`` instead of ``HGNC:11359``.
        * All other qualifiers (claim-level feature assertions, plus
          ``consolidates_dataset``) use ``accession`` then ``gene`` then
          ``concept``: gene-symbol mentions dominate feature claims.

        The first resolver that returns a hit wins for that slot; later
        resolvers are not consulted. ``ontology_terms[slot]`` carries the
        principal CURIE; ``xrefs[slot]`` carries the resolver-specific
        cross-references (HGNC / MGI / Ensembl / NCBIGene / Nemve1 /
        NvERTx / UniProt for genes; ``{curie, label, ontology,
        description}`` for concepts; empty for accessions).

        Pure function: never mutates the input dict.
        """
        # Shallow-copy the claim and its nested ontology/xrefs dicts so
        # the caller's input is not mutated.
        claim = dict(claim)
        if "ontology_terms" in claim:
            claim["ontology_terms"] = dict(claim["ontology_terms"])
        if "ontology_term_justifications" in claim:
            claim["ontology_term_justifications"] = dict(claim["ontology_term_justifications"])
        if "xrefs" in claim:
            claim["xrefs"] = dict(claim["xrefs"])
        if "sssom_mappings" in claim:
            claim["sssom_mappings"] = dict(claim["sssom_mappings"])

        qualifier = (claim.get("qualifier") or "").strip()
        priority = _QUALIFIER_PRIORITY.get(qualifier, _DEFAULT_PRIORITY)
        # Optional species hint forwarded to every resolver. Non-gene
        # resolvers ignore it. When the claim has no hint, the gene
        # resolver falls back to the casing heuristic (UPPER then human,
        # Title-case then mouse) before defaulting to human via load order.
        species_hint = claim.get("species")

        for slot in _ENRICHABLE_SLOTS:
            name = (claim.get(slot) or "").strip()
            if not name:
                continue
            for resolver_name in priority:
                resolver = _RESOLVERS[resolver_name]
                result = resolver(name, species=species_hint)
                if result is None:
                    continue
                curie, xrefs = result
                # The resolver may have signalled a sub-kind via the
                # ``_resolver_kind`` sentinel inside xrefs (e.g. the
                # ``accession`` slot's pathway-substring fallback reports
                # ``embedded_pathway``; the ``concept`` slot's trim-tail
                # fallback reports ``trimmed_concept``). Pop it before
                # storing xrefs so downstream consumers never see the
                # sentinel, but use it to drive the SSSOM confidence
                # floor and ``skos:closeMatch`` predicate routing.
                if isinstance(xrefs, dict):
                    actual_kind = xrefs.pop("_resolver_kind", resolver_name)
                else:
                    actual_kind = resolver_name
                claim.setdefault("ontology_terms", {})[slot] = curie
                justification = (
                    "semapv:LexicalMatching"
                    if resolver_name == "concept"
                    else "semapv:CompositeMatching"
                )
                claim.setdefault("ontology_term_justifications", {})[slot] = justification
                if xrefs:
                    claim.setdefault("xrefs", {})[slot] = xrefs
                # Emit a first-class SSSOMMapping fragment for the hit
                # (indicium >= 1.11.0). Folded into
                # IndiciumDocument.sssom_mappings downstream by
                # Perspicacite. Label comes from xrefs when the
                # resolver carried one; accessions and gene rows
                # without a label key pass None and the fragment
                # omits ``object_label``.
                fragment = build_sssom_mapping(
                    raw_text=name,
                    object_curie=curie,
                    object_label=xrefs.get("label") if isinstance(xrefs, dict) else None,
                    justification=justification,
                    resolver_kind=actual_kind,
                )
                claim.setdefault("sssom_mappings", {})[fragment["id"]] = fragment
                break  # First hit per slot wins.

        return claim

    def shacl_shapes(self) -> rdflib.Graph:
        """Return ASB-claim SHACL shapes (CURIE-format validation only).

        Pipeline-level RDF validation for transcriptomics workflows (SO,
        UBERON, CL class shapes) lives upstream in the consuming graph
        store, not here.
        """
        g = rdflib.Graph()
        ttl = files("indicium_adapters_transcriptomics.shapes").joinpath(
            "transcriptomics.shapes.ttl"
        )
        g.parse(ttl.open("rb"), format="turtle")
        return g
