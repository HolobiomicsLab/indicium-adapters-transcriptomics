"""Quickstart — discover, enrich, validate in 30 lines.

The shortest runnable demonstration that the adapter is wired
correctly: discover both adapters via the indicium entry-point
machinery, enrich a hand-built claim through each one, then
validate the SHACL shapes graph parses. Use this as the
"Hello World" smoke test when bringing the package into a new
environment.

The Warner2018 demo (``examples/enrich_warner2018.py``) is the
fuller showcase across 25 representative claims.

Run::

    py -3.13 examples/quickstart.py
"""
from __future__ import annotations

import json
import sys

import pyshacl
import rdflib
from indicium_adapters import compose_adapters, discover_adapters

sys.stdout.reconfigure(encoding="utf-8")


def main() -> int:
    # 1) Discover.
    adapters = discover_adapters()
    print(f"[discover] {len(adapters)} adapters: {sorted(adapters)}")
    bulk = adapters["bulk-transcriptomics"]
    sc = adapters["single-cell-transcriptomics"]

    # 2) Enrich one claim per adapter.
    bulk_claim = bulk.enrich_claim({
        "subject": "WNT4",
        "qualifier": "differentially_expressed_in",
        "object": "regenerating_oral_pole",
    })
    sc_claim = sc.enrich_claim({
        "subject": "MMP9",
        "qualifier": "localized_in_tissue",
        "object": "mesentery",
    })
    print(f"[bulk] subject CURIE: {bulk_claim['ontology_terms']['subject']}")
    print(f"[sc]   subject CURIE: {sc_claim['ontology_terms']['subject']}")
    print(f"[sc]   object  CURIE: {sc_claim['ontology_terms']['object']}")

    # 3) Verify the composed SHACL shapes graph parses cleanly under
    #    pyshacl. We validate an *empty* data graph against the shapes
    #    — the call surfaces any TTL syntax / prefix issues in the
    #    bundled shapes as a raised exception. The conformance result
    #    is uninteresting (no data → no violations) but the parse step
    #    is the actual smoke test.
    composite = compose_adapters([bulk, sc])
    shapes_graph = composite.shacl_shapes()
    empty_data = rdflib.Graph()
    conforms, _, _ = pyshacl.validate(
        data_graph=empty_data,
        shapes_graph=shapes_graph,
        meta_shacl=False,
        debug=False,
    )
    print(f"[shacl] shapes graph parsed ({len(shapes_graph)} triples), "
          f"conforms on empty data: {conforms}")

    # 4) Final sanity dump.
    print()
    print("Enriched bulk claim:")
    print(json.dumps(bulk_claim, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
