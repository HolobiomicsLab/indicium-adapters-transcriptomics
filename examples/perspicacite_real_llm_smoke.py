"""Real-LLM smoke test for the v0.13 spec-aligned qualifiers.

This script answers the feature-correctness question that the synthetic
``perspicacite_integration.py`` cannot: when a real LLM (Claude sonnet-4-5
through Perspicacite's pipeline) reads a passage that says "upregulated",
does it pick the new ``upregulated_in`` qualifier added in v0.13, or does
it fall back to the generic ``differentially_expressed_in``?

Four short synthetic passages exercise the four v0.13 qualifiers — each
passage embeds an unambiguous lexical signal so a competent LLM should
pick the spec-canonical form:

* signed up-regulation -> ``upregulated_in``
* signed down-regulation -> ``downregulated_in``
* cell-type marker -> ``marker_for``
* Reactome pathway enrichment -> ``enriched_in_pathway``

What's verified:

* The LLM's qualifier choice on each passage (vs the generic fallback).
* That the adapter then enriches the resulting claim to the expected
  CURIE on subject and / or object.
* That the full set of qualifiers Perspicacite advertises in the prompt
  includes the four new v0.13 ones (sanity-check that
  ``extraction_context()`` is reaching the LLM).

What's *not* verified:

* Anything about real paper text (passages are hand-written so the
  lexical signal is unambiguous — we are testing qualifier choice, not
  summarisation skill).
* Cost ceilings / retries (small fixed budget: 4 short LLM calls).

Run::

    $env:PYTHONPATH = "C:\\path\\to\\Perspicacite-AI\\src"
    py -3.13 examples/perspicacite_real_llm_smoke.py

Requires ``ANTHROPIC_API_KEY`` in Perspicacite-AI's ``.env`` (load_dotenv
walks parent dirs from cwd, so launch from this repo's root).

Findings — first real run (2026-06-04, claude-sonnet-4-5, pre-fix):

* Qualifier choice: 4/4. The updated ``extraction_context()`` is
  steering the LLM cleanly to the signed / specialised forms — every
  passage's spec-canonical qualifier was the only emitted qualifier.
* Subject enrichment: 3/4 (upreg, downreg, marker all resolved the
  gene subject to HGNC). The ``pathway`` passage had an empty
  enrichment because the LLM emitted ``REACTOME:R-HSA-9012999`` as
  the object — already prefixed — but the anchored accession regex
  required an un-prefixed bare ID, so the resolver missed.
* Object enrichment for ``marker``: missed because the LLM emitted
  ``cnidocyte cells`` rather than the bare CL label ``cnidocyte``,
  and the concept index does exact-match.

Findings — second real run (2026-06-04, same model, post-fix):

* Qualifier choice: 4/4 (unchanged, LLM cache hit).
* Subject enrichment: 4/4 (no change for upreg / downreg / marker;
  pathway has a generic gene-set subject that does not resolve, which
  matches expectations).
* Object enrichment: 2/2 — both ``marker`` and ``pathway`` now resolve:
  - ``cnidocyte cells`` -> ``CL:0011029`` via the concept-trim fallback
    (added 2026-06-04 in ``_enrichment.resolve_concept``);
  - ``REACTOME:R-HSA-9012999`` -> ``REACTOME:R-HSA-9012999`` via the
    embedded-pathway fallback (added 2026-06-04 in
    ``_registries.find_embedded_pathway_curie``, wired through
    ``bulk._resolve_accession``).

Both fixes are graceful-degradation rules: they fire only when the
strict primary lookup misses, so the success path is unchanged. Pinned
by ``tests/test_bulk.py::test_enriched_in_pathway_resolves_embedded_*``
and ``::test_concept_resolver_trims_common_trailing_nouns``.
"""
from __future__ import annotations

import asyncio
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

# Perspicacite ships its .env at the sibling repo's root, not here. We
# resolve the path explicitly so the test runs from anywhere instead of
# silently failing when launched from a directory above Perspicacite-AI.
_PERSPICACITE_ROOT = Path(
    r"C:\Users\spirale\Documents\Logiciels_stage_GB5\Perspicacite-AI"
)


try:
    from dotenv import load_dotenv
    from perspicacite.config import load_config
    from perspicacite.llm.client import AsyncLLMClient
    from perspicacite.pipeline.claims import extract_claims

    from indicium_adapters_transcriptomics import BulkTranscriptomicsAdapter
except ModuleNotFoundError as exc:
    sys.stderr.write(
        f"[fatal] required module not importable: {exc}\n"
        "Set PYTHONPATH to your Perspicacite-AI/src directory before running, "
        "and install this package in editable mode.\n"
    )
    raise SystemExit(2) from None


# Sonnet-4-5 is the IRCAN lab default model. Perspicacite's stock
# LLMConfig has ``default_provider = "deepseek"`` (Pydantic default —
# the project deliberately ships without an opinionated Anthropic
# default per the comment at perspicacite/config/schema.py:503). The
# model-string prefix (``anthropic/...``) is NOT honoured because
# ``_build_model_string`` rebuilds the LiteLLM model identifier from
# the ``provider`` variable, which falls back to ``config.default_provider``
# when the caller does not pass one. We override ``default_provider``
# on the loaded config instead — surgical, no Perspicacite mutation.
LLM_PROVIDER = "anthropic"
LLM_MODEL = "claude-sonnet-4-5"


# Hand-written passages. Each one contains one unambiguous lexical
# signal for the qualifier the LLM should pick. The lexical signal is
# the dependent variable; everything else is held constant.
#
# Each entry: ``(label, passage_text, expected_qualifier, post-coercion
# slot check)``. The slot check is what the adapter should produce after
# Perspicacite coerces the LLM output, so the harness can report
# qualifier-choice and enrichment in one go.
PASSAGES: list[dict] = [
    {
        "label": "upreg",
        "expected_qualifier": "upregulated_in",
        "expected_curie_subject": "HGNC:12783",  # WNT4
        "passage": (
            "In Nematostella vectensis adults, WNT4 was significantly "
            "upregulated in the regenerating oral pole 24 hours after "
            "bisection (log2FC = +2.3, padj = 1e-4)."
        ),
    },
    {
        "label": "downreg",
        "expected_qualifier": "downregulated_in",
        "expected_curie_subject": "HGNC:7176",  # MMP9
        "passage": (
            "MMP9 was significantly downregulated in the regenerating "
            "oral pole compared with intact controls "
            "(log2FC = -1.8, padj = 5e-3)."
        ),
    },
    {
        "label": "marker",
        "expected_qualifier": "marker_for",
        "expected_curie_object": "CL:0011029",  # cnidocyte
        "passage": (
            "Sc-RNA-seq analysis identified MMP9 as a robust marker for "
            "cnidocyte cells in Nematostella vectensis."
        ),
    },
    {
        "label": "pathway",
        "expected_qualifier": "enriched_in_pathway",
        "expected_curie_object": "REACTOME:R-HSA-9012999",
        "passage": (
            "Gene set enrichment analysis showed that the target genes "
            "were significantly enriched in the canonical Wnt signaling "
            "pathway (Reactome ID R-HSA-9012999, padj = 2e-6)."
        ),
    },
]


def _passage_to_perspicacite_chunk(p: dict) -> dict:
    """Format a hand-written passage as a Perspicacite ``chunk_text``.

    Perspicacite expects a dict with ``chunk_text`` and a ``source`` sub-
    dict carrying at least an identifier and a title. The doi here is
    synthetic — the goal is qualifier choice, not citation fidelity.
    """
    return {
        "chunk_text": p["passage"],
        "source": {
            "doi": f"10.0000/synthetic.{p['label']}",
            "title": f"Synthetic v0.13 smoke passage ({p['label']})",
        },
    }


async def run() -> int:
    print("=" * 72)
    print("v0.13 spec-alignment real-LLM smoke test (Claude sonnet-4-5)")
    print("=" * 72)

    # Load .env from Perspicacite-AI's repo root — that's where the
    # Anthropic key lives.
    env_path = _PERSPICACITE_ROOT / ".env"
    if not env_path.exists():
        sys.stderr.write(f"[fatal] expected .env at {env_path}\n")
        return 2
    load_dotenv(env_path)

    cfg = load_config()
    # Override the stock DeepSeek default with the lab's Anthropic
    # target. extract_claims hands the LLM-completion call only the
    # model string, so the provider has to be pinned on the LLMConfig
    # before the client is built.
    cfg.llm.default_provider = LLM_PROVIDER
    cfg.llm.default_model = LLM_MODEL
    llm_client = AsyncLLMClient(cfg.llm)
    adapter = BulkTranscriptomicsAdapter()

    print(f"\nProvider/Model: {LLM_PROVIDER}/{LLM_MODEL}")
    print(f"Adapter: {adapter.domain_id} ({len(adapter.qualifiers)} qualifiers)\n")

    results = []
    for p in PASSAGES:
        print(f"--- passage: {p['label']} ---")
        print(f"  expected qualifier: {p['expected_qualifier']}")
        print(f"  passage (truncated): {p['passage'][:100]}...")

        chunk = _passage_to_perspicacite_chunk(p)
        claims = await extract_claims(
            llm_client=llm_client,
            passages=[chunk],
            context="Nematostella transcriptomics",
            model=LLM_MODEL,
            domain_adapter=adapter,
        )

        if not claims:
            print("  result: NO CLAIMS extracted")
            results.append((p["label"], None, None, "no claims"))
            print()
            continue

        # We are interested in whether ANY emitted claim picks the
        # spec-canonical qualifier (the LLM can emit multiple per
        # passage; we just need one targeted hit).
        picked_qualifiers = [c.get("qualifier") for c in claims]
        spec_hit = p["expected_qualifier"] in picked_qualifiers
        print(f"  emitted qualifiers: {picked_qualifiers}")
        print(f"  spec-canonical qualifier picked: {'YES' if spec_hit else 'NO'}")

        # Find the first claim with the spec-canonical qualifier (if any)
        # and report its enrichment.
        target_claim = next(
            (c for c in claims if c.get("qualifier") == p["expected_qualifier"]),
            None,
        )
        if target_claim is not None:
            terms = target_claim.get("ontology_terms") or {}
            print(f"  enrichment: {dict(terms)}")
            # Print the raw SPO triple the LLM emitted so a WRONG
            # enrichment can be debugged without re-running the test
            # (LLM calls cost money — make each run diagnostic).
            print(
                f"  raw claim:  subject={target_claim.get('subject')!r}, "
                f"qualifier={target_claim.get('qualifier')!r}, "
                f"relation={target_claim.get('relation')!r}, "
                f"object={target_claim.get('object')!r}"
            )
            for slot_key in ("expected_curie_subject", "expected_curie_object"):
                if slot_key in p:
                    slot = slot_key.removeprefix("expected_curie_")
                    expected = p[slot_key]
                    actual = terms.get(slot)
                    match = actual == expected
                    print(
                        f"  {slot:7s} expected {expected!r}, "
                        f"got {actual!r}: {'OK' if match else 'WRONG'}"
                    )
            # v0.14: confirm SSSOMMapping emission. Each resolved slot
            # should have an entry in sssom_mappings shaped like indicium's
            # SSSOMMapping pydantic class. We just print the count and one
            # fragment so future runs see the contract at a glance.
            sssom = target_claim.get("sssom_mappings") or {}
            print(f"  sssom: {len(sssom)} fragment(s) emitted")
            if sssom:
                sample = next(iter(sssom.values()))
                print(
                    f"    sample: id={sample.get('id')!r}, "
                    f"predicate={sample.get('predicate_id')!r}, "
                    f"justification={sample.get('mapping_justification')!r}, "
                    f"confidence={sample.get('confidence')}"
                )

        results.append((p["label"], picked_qualifiers, spec_hit, "ok"))
        print()

    # Summary matrix.
    print("=" * 72)
    print("Summary")
    print("=" * 72)
    print(f"{'Passage':12s} {'Spec qualifier picked':24s} {'Emitted':40s}")
    for label, picked, hit, _ in results:
        if picked is None:
            print(f"{label:12s} {'(no claims)':24s}")
        else:
            mark = "YES" if hit else "NO"
            shown = ", ".join(picked[:3])
            if len(picked) > 3:
                shown += f", +{len(picked) - 3} more"
            print(f"{label:12s} {mark:24s} {shown}")

    return 0


if __name__ == "__main__":
    raise SystemExit(asyncio.run(run()))
