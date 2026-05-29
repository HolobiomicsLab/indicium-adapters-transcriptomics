"""End-to-end integration test against Perspicacite's claim-extraction pipeline.

This script answers the question "does Perspicacite's extract_claims() flow
through BulkTranscriptomicsAdapter correctly?" without burning real LLM
calls. It mocks the LLM with a fixed JSON response shaped like what a real
extraction over Warner2018 would produce, then exercises the full contract:

* Perspicacite merges the adapter's qualifiers into its base SuperPattern set.
* Perspicacite appends extraction_context() to the prompt.
* Each coerced claim is passed through adapter.enrich_claim().
* Out-of-vocab qualifiers are dropped.
* Gene symbols in unexpected slots (e.g. relation) are silently ignored,
  since the adapter only inspects subject and object.

The synthetic claim list mixes:

* SuperPattern qualifiers (inhibits, activates) from Perspicacite's base set
* Bulk-adapter feature qualifiers (differentially_expressed_in)
* Bulk-adapter study-level qualifiers (aligned_with, consolidates_dataset)
* One gene-in-relation drift case
* One out-of-vocab qualifier that should be dropped

Run::

    $env:PYTHONPATH = "C:\\path\\to\\Perspicacite-AI\\src"
    py -3.13 examples/perspicacite_integration.py
"""
from __future__ import annotations

import asyncio
import json
import sys
from collections import Counter
from unittest.mock import AsyncMock

sys.stdout.reconfigure(encoding="utf-8")

# Guarded import block: Perspicacite lives in a sibling repo and may not be
# on PYTHONPATH. We catch ModuleNotFoundError to print a clear hint instead
# of a bare stack trace. The BulkTranscriptomicsAdapter import sits in the
# same try block for symmetry; it cannot fail in a working install of this
# package but if it did, the same message channel would surface the error.
try:
    from perspicacite.pipeline.claims import extract_claims

    from indicium_adapters_transcriptomics import BulkTranscriptomicsAdapter
except ModuleNotFoundError as exc:
    sys.stderr.write(
        f"[fatal] required module not importable: {exc}\n"
        "Set PYTHONPATH to your Perspicacite-AI/src directory before running, "
        "and install this package in editable mode.\n"
    )
    raise SystemExit(2) from None


# Synthetic LLM output that mimics what Perspicacite would return after
# running extract_claims over Warner et al. 2025 passages. The shape matches
# the real prompt contract (5-slot SuperPattern + evidence metadata).
WARNER_SYNTHETIC_CLAIMS = [
    # 1) SuperPattern claim: U0126 inhibits MAPK1.
    #    Bulk adapter's default priority should resolve U0126 (ChEBI) on
    #    the subject slot and MAPK1 (HGNC) on the object slot.
    {
        "context": "Nematostella oral regeneration",
        "subject": "U0126",
        "qualifier": "inhibits",
        "relation": "blocks MEK1/2 phosphorylation of",
        "object": "MAPK1",
        "claim_type": "explicit",
        "evidence_type": "data",
        "source_type": "text",
        "quote": "U0126 inhibited MAPK1/2 phosphorylation in regenerating polyps",
        "source_doi": "10.1242/dev.196782",
    },
    # 2) Bulk study-level qualifier with the STAR collision.
    #    Adapter's qualifier-aware priority should resolve STAR to
    #    biotools:star (aligner), not HGNC:11359 (the gene).
    {
        "context": "RNA-seq pipeline",
        "subject": "warner_2025_study",
        "qualifier": "aligned_with",
        "relation": "uses",
        "object": "STAR",
        "claim_type": "explicit",
        "evidence_type": "data",
        "source_type": "text",
        "quote": "Reads were aligned with STAR v2.7",
        "source_doi": "10.1242/dev.196782",
    },
    # 3) Bulk feature qualifier with a gene in subject.
    #    Default priority should resolve WNT4 to HGNC:12783.
    {
        "context": "regenerating oral pole",
        "subject": "WNT4",
        "qualifier": "differentially_expressed_in",
        "relation": "is upregulated in",
        "object": "regenerating_oral_pole",
        "claim_type": "explicit",
        "evidence_type": "data",
        "source_type": "text",
        "quote": "WNT4 was significantly upregulated in the regenerating oral pole",
        "source_doi": "10.1242/dev.196782",
    },
    # 4) Accession claim: GSE168938 on object slot should resolve via the
    #    accession regex to GEO:GSE168938.
    {
        "context": "data deposition",
        "subject": "warner_2025_study",
        "qualifier": "consolidates_dataset",
        "relation": "deposited as",
        "object": "GSE168938",
        "claim_type": "explicit",
        "evidence_type": "citation",
        "source_type": "text",
        "quote": "Raw sequencing data are available at GEO under accession GSE168938",
        "source_doi": "10.1242/dev.196782",
    },
    # 5) Drift case: gene symbol MAPK1 lives in the relation slot, not
    #    subject or object. The adapter inspects only subject and object,
    #    so MAPK1 here will NOT be resolved. This is a real risk when the
    #    LLM phrases a claim with the gene as a verb argument.
    {
        "context": "MEK pathway during regeneration",
        "subject": "FGF8 signaling",
        "qualifier": "activates",
        "relation": "MAPK1",
        "object": "regeneration",
        "claim_type": "implicit",
        "evidence_type": "inference",
        "source_type": "text",
        "quote": "FGF8 signaling drives MAPK activation during oral regeneration",
        "source_doi": "10.1242/dev.196782",
    },
    # 6) Out-of-vocab qualifier. Perspicacite's _coerce_claim should drop
    #    this entry entirely, since "frobnicates" is in neither the
    #    SuperPattern set nor the bulk-adapter qualifier set.
    {
        "context": "noise check",
        "subject": "U0126",
        "qualifier": "frobnicates",
        "relation": "does",
        "object": "stuff",
        "claim_type": "explicit",
        "evidence_type": "data",
        "source_type": "text",
        "quote": "noise",
        "source_doi": "10.1242/dev.196782",
    },
    # 7) Sequencer in object slot for a study-level qualifier. Should
    #    resolve to an EFO instrument CURIE via the OLS4 corpus harvest.
    {
        "context": "RNA-seq pipeline",
        "subject": "warner_2025_study",
        "qualifier": "sequenced_on",
        "relation": "uses",
        "object": "Illumina NovaSeq 6000",
        "claim_type": "explicit",
        "evidence_type": "data",
        "source_type": "text",
        "quote": "Libraries were sequenced on an Illumina NovaSeq 6000",
        "source_doi": "10.1242/dev.196782",
    },
]


def _make_mock_llm() -> AsyncMock:
    """Build an AsyncMock that mimics Perspicacite's AsyncLLMClient.complete()."""
    llm = AsyncMock()
    llm.complete = AsyncMock(
        return_value=json.dumps({"claims": WARNER_SYNTHETIC_CLAIMS})
    )
    return llm


async def run() -> int:
    print("=" * 72)
    print("Perspicacite + BulkTranscriptomicsAdapter end-to-end integration test")
    print("=" * 72)

    bulk = BulkTranscriptomicsAdapter()
    llm = _make_mock_llm()

    # Synthetic passages: a doi + a short text. Perspicacite's prompt formats
    # these as numbered chunks; the real LLM would use them as evidence, but
    # our mock ignores them and returns the canned WARNER_SYNTHETIC_CLAIMS.
    passages = [{
        "chunk_text": "Warner et al. 2025: U0126 inhibits MAPK during regeneration.",
        "source": {"doi": "10.1242/dev.196782", "title": "Warner 2025 abstract"},
    }]

    claims = await extract_claims(
        llm_client=llm,
        passages=passages,
        context="Nematostella regeneration transcriptomics",
        domain_adapter=bulk,
    )

    n_input = len(WARNER_SYNTHETIC_CLAIMS)
    n_coerced = len(claims)
    print(f"\nInput synthetic claims: {n_input}")
    print(f"After Perspicacite coercion + adapter enrichment: {n_coerced}")
    print(f"Dropped: {n_input - n_coerced} (expected: 1, the frobnicates entry)")

    # Verify the prompt actually contained the adapter's extraction_context.
    prompt_sent = llm.complete.call_args[1]["messages"][0]["content"]
    has_context = "Domain: bulk transcriptomics." in prompt_sent
    bulk_qualifier_in_prompt = "differentially_expressed_in" in prompt_sent
    superpattern_in_prompt = "causes" in prompt_sent and "inhibits" in prompt_sent
    print()
    print("Prompt composition checks:")
    print(f"  adapter extraction_context() appended : {has_context}")
    print(f"  bulk qualifier in valid set           : {bulk_qualifier_in_prompt}")
    print(f"  SuperPattern qualifiers still listed  : {superpattern_in_prompt}")

    # Per-claim enrichment report.
    print()
    print("Per-claim enrichment trace:")
    n_enriched = 0
    drift_log: list[tuple[int, str]] = []
    for i, c in enumerate(claims, 1):
        terms = c.get("ontology_terms") or {}
        marker = "[OK]" if terms else "[--]"
        slot_curies = ", ".join(f"{k}={v}" for k, v in terms.items()) or "(no resolution)"
        print(f"  [{i}] qualifier={c['qualifier']:30s} {marker} {slot_curies}")
        if terms:
            n_enriched += 1
        # Drift detection: scan the relation slot for tokens that look like
        # gene symbols or accessions. The adapter does not inspect this slot,
        # so any hit here is information lost.
        relation = (c.get("relation") or "").strip()
        looks_like_gene = (
            relation
            and len(relation) <= 12
            and relation.replace("-", "").isalnum()
            and (relation.isupper() or relation[0].isupper())
        )
        if looks_like_gene:
            drift_log.append((i, relation))

    print()
    print(f"Enriched: {n_enriched}/{n_coerced}")

    # Per-ontology coverage matrix.
    by_ontology: Counter = Counter()
    for c in claims:
        for curie in (c.get("ontology_terms") or {}).values():
            prefix = curie.split(":", 1)[0] if ":" in curie else "?"
            by_ontology[prefix] += 1
    if by_ontology:
        print()
        print("Per-ontology hit count:")
        for ont, n in by_ontology.most_common():
            print(f"  {ont:14s} {n}")

    # Drift findings.
    if drift_log:
        print()
        print("Drift signals (gene-like tokens in relation slot, not enriched):")
        for idx, rel in drift_log:
            print(f"  [{idx}] relation={rel!r}")
    else:
        print()
        print("No drift signals.")

    # Verify the STAR collision was routed correctly.
    print()
    print("Targeted verification (STAR collision):")
    star_claim = next(
        (c for c in claims if c["qualifier"] == "aligned_with"
         and c["object"] == "STAR"),
        None,
    )
    if star_claim is not None:
        star_curie = (star_claim.get("ontology_terms") or {}).get("object", "")
        expected = "biotools:star"
        ok = star_curie == expected
        print(f"  aligned_with STAR -> {star_curie!r} (expected {expected!r}): "
              f"{'OK' if ok else 'WRONG'}")

    # Try SHACL validation if indicium is importable.
    print()
    print("SHACL validation:")
    try:
        from perspicacite.pipeline.claims import validate_claims
        try:
            conforms, report = validate_claims(claims, domain_adapter=bulk)
            print(f"  conforms: {conforms}")
            if not conforms:
                print(f"  report excerpt: {report[:500]}")
        except (ImportError, ModuleNotFoundError) as exc:
            print(f"  skipped (indicium not installed): {exc}")
    except Exception as exc:
        print(f"  failed: {type(exc).__name__}: {exc}")

    return 0


def main() -> int:
    # asyncio.run because extract_claims is async.
    return asyncio.run(run())


if __name__ == "__main__":
    raise SystemExit(main())
