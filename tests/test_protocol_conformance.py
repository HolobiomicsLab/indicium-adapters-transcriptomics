"""Centralised Protocol-conformance contract tests (v0.10.0 / v1.0.0 prep).

The bulk and single-cell adapter test modules already assert
``isinstance(a, DomainAdapter)`` / ``isinstance(a, SHACLProvider)`` /
``api_version == CURRENT_API_VERSION`` per-adapter. This file collects
the **cross-cutting** invariants that apply to **every** adapter
exposed by the package, and adds the depth checks the sibling
epigenomics adapter does at v0.1:

- Every adapter is discoverable via the ``indicium.adapters`` entry-
  point group (not just our package's direct imports).
- Every adapter's ``domain_id`` is canonical, lowercase, and matches
  the entry-point name.
- Every adapter's ``qualifiers`` are all snake_case (lowercase + `_`
  separators) — extraction prompts and SPARQL templates depend on
  this convention.
- Every adapter's ``ontology_prefixes`` map dereferences to valid
  http(s) URIs ending with `:`, `/`, `_` or `#` (the OBO Foundry
  convention) so the ASB graph emitter can mint round-trippable
  CURIEs.
- Every adapter's ``shacl_shapes`` returns a parsed graph with at
  least one ``sh:NodeShape``.
"""
from __future__ import annotations

import re
from urllib.parse import urlparse

import pytest
import rdflib
from indicium_adapters import (
    CURRENT_API_VERSION,
    DomainAdapter,
    SHACLProvider,
    discover_adapters,
)

from indicium_adapters_transcriptomics import (
    BulkTranscriptomicsAdapter,
    SingleCellTranscriptomicsAdapter,
)

# All adapter classes the package exports. New siblings (e.g. a
# microarray-specific adapter in v1.1+) should be appended here so
# they automatically get the full conformance suite.
ALL_ADAPTER_CLASSES = (
    BulkTranscriptomicsAdapter,
    SingleCellTranscriptomicsAdapter,
)


# --- Discovery / Protocol conformance --------------------------------------


def test_every_adapter_is_discoverable_via_entry_points():
    """The indicium-adapters discovery helper walks the
    ``indicium.adapters`` entry-point group. Each adapter class
    exported by the package must be registered there — otherwise a
    downstream consumer running ``discover_adapters()`` will silently
    miss it."""
    found = discover_adapters()
    package_domain_ids = {a.domain_id for a in ALL_ADAPTER_CLASSES}
    discovered_domain_ids = set(found.keys())
    missing = package_domain_ids - discovered_domain_ids
    assert not missing, f"adapters not registered in entry points: {missing}"


@pytest.mark.parametrize(
    "adapter_cls", ALL_ADAPTER_CLASSES, ids=lambda c: c.domain_id
)
def test_implements_domain_adapter_protocol(adapter_cls):
    assert isinstance(adapter_cls(), DomainAdapter)


@pytest.mark.parametrize(
    "adapter_cls", ALL_ADAPTER_CLASSES, ids=lambda c: c.domain_id
)
def test_implements_shacl_provider_protocol(adapter_cls):
    assert isinstance(adapter_cls(), SHACLProvider)


@pytest.mark.parametrize(
    "adapter_cls", ALL_ADAPTER_CLASSES, ids=lambda c: c.domain_id
)
def test_api_version_matches_indicium_core(adapter_cls):
    """A version mismatch here means the indicium core API moved on
    and the adapter wasn't refreshed — the symptom would be subtle
    runtime breakage rather than an import error."""
    assert adapter_cls.api_version == CURRENT_API_VERSION


# --- Naming / convention invariants ----------------------------------------


_SNAKE_CASE = re.compile(r"^[a-z][a-z0-9_]*$")
_DOMAIN_ID = re.compile(r"^[a-z][a-z0-9-]*$")


@pytest.mark.parametrize(
    "adapter_cls", ALL_ADAPTER_CLASSES, ids=lambda c: c.domain_id
)
def test_domain_id_is_canonical(adapter_cls):
    """``domain_id`` must be kebab-case (lowercase + hyphens). Used as
    the entry-point key and as a dict key downstream."""
    assert _DOMAIN_ID.match(adapter_cls.domain_id), (
        f"domain_id {adapter_cls.domain_id!r} is not kebab-case"
    )


@pytest.mark.parametrize(
    "adapter_cls", ALL_ADAPTER_CLASSES, ids=lambda c: c.domain_id
)
def test_qualifiers_are_all_snake_case(adapter_cls):
    """Extraction prompts, SPARQL templates, and downstream RDF
    serialisations depend on this convention. CamelCase or kebab-case
    qualifiers would break those."""
    non_snake = [
        q for q in adapter_cls.qualifiers if not _SNAKE_CASE.match(q)
    ]
    assert not non_snake, (
        f"{adapter_cls.domain_id}: non-snake_case qualifiers: {non_snake}"
    )


# --- Ontology-prefix map validity ------------------------------------------


@pytest.mark.parametrize(
    "adapter_cls", ALL_ADAPTER_CLASSES, ids=lambda c: c.domain_id
)
def test_ontology_prefixes_resolve_to_http_uri(adapter_cls):
    """Every prefix-map entry must dereference to an http(s) URI that
    looks like a namespace base (ending in ``:`` / ``/`` / ``_`` /
    ``#``). Otherwise CURIE expansion downstream produces invalid
    IRIs."""
    bad: list[tuple[str, str]] = []
    for prefix, base in adapter_cls.ontology_prefixes.items():
        parsed = urlparse(base)
        if parsed.scheme not in ("http", "https"):
            bad.append((prefix, base))
            continue
        if not base.endswith((":", "/", "_", "#")):
            bad.append((prefix, base))
    assert not bad, (
        f"{adapter_cls.domain_id}: prefix-map entries with invalid base "
        f"URIs: {bad}"
    )


# --- SHACL shapes contract -------------------------------------------------


@pytest.mark.parametrize(
    "adapter_cls", ALL_ADAPTER_CLASSES, ids=lambda c: c.domain_id
)
def test_shacl_shapes_graph_contains_at_least_one_node_shape(adapter_cls):
    """``shacl_shapes()`` must return an ``rdflib.Graph`` with at least
    one ``sh:NodeShape``. The graph is the entry point the consumer
    feeds to pyshacl, so an empty graph would be a silent no-op
    validator."""
    g = adapter_cls().shacl_shapes()
    assert isinstance(g, rdflib.Graph)
    sh = rdflib.Namespace("http://www.w3.org/ns/shacl#")
    node_shapes = list(g.subjects(rdflib.RDF.type, sh.NodeShape))
    assert node_shapes, (
        f"{adapter_cls.domain_id}: shacl_shapes() returned a graph with "
        "zero NodeShapes"
    )
