"""Maintainer script — harvest CURIEs for bioinformatics software via bio.tools.

Companion to ``build_corpus_concept_candidates.py``. The OLS4 harvest
catches assays / platforms / statistical concepts that live in
ontologies (EFO / OBI / NCIT / STATO), but it does **not** catch the
specific bioinformatics software referenced in the corpus
(``STAR``, ``Seurat``, ``DESeq2``, ``Cell Ranger`` …). Those tools are
registered at https://bio.tools/, the canonical registry for
bioinformatics software (each tool has a stable ``biotoolsID``).

This script reads ``corpus_concept_unresolved.txt`` (produced by the
OLS4 harvester), filters down to the **software axes** only, and tries
to resolve each remaining term via bio.tools::

  1. Direct API lookup ``GET https://bio.tools/api/<id>/?format=json``
     with a normalised ID (lowercased, spaces → hyphens). Cheap and
     usually authoritative for canonical tool names.
  2. If 404, fall back to the search endpoint
     ``GET https://bio.tools/api/tool/?name=<term>&format=json``,
     paginate through the response, and accept only entries whose
     ``name`` field matches the corpus term **exactly**, case-insensitive.

No CURIE is hand-picked. The accepted entries land in
``data/biotools_concept_index.tsv`` and are loaded into the same
``ConceptIndex`` singleton as the OLS4 indexes.

Usage (maintainer-only)::

    python -m indicium_adapters_transcriptomics.data.build_biotools_index
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path
from typing import Final, NamedTuple

import requests

BIOTOOLS_API = "https://bio.tools/api"

OUT_COLUMNS = ("term", "curie", "label", "description", "ontology", "source")

# Axes for which bio.tools is the right registry. Platforms, protocols and
# statistical concepts stay with OLS4 / EFO / NCIT / STATO.
SOFTWARE_AXES: Final[frozenset[str]] = frozenset({
    "aligner_quant",
    "de_tool",
    "sc_tool",
    "spatial_tool",
    "long_read_tool",
    "qc",
})


class Hit(NamedTuple):
    term: str         # original corpus term (preserves casing)
    axis: str         # software axis it came from (provenance)
    biotools_id: str  # canonical biotoolsID
    name: str         # bio.tools display name
    description: str  # short description (truncated to 500 chars)


def _retrying_get(
    url: str, *, params: dict[str, str] | None, timeout: int, retries: int
) -> requests.Response | None:
    """GET with backoff on 5xx and timeouts. ``None`` on persistent failure."""
    last_exc: Exception | None = None
    for attempt in range(1, retries + 2):
        try:
            r = requests.get(url, params=params, timeout=timeout)
            if 500 <= r.status_code < 600:
                last_exc = requests.HTTPError(f"HTTP {r.status_code}")
                print(
                    f"  ⚠ HTTP {r.status_code} on attempt "
                    f"{attempt}/{retries + 1}, retrying ...",
                    file=sys.stderr,
                )
                time.sleep(2 * attempt)
                continue
            return r
        except requests.exceptions.Timeout as exc:
            last_exc = exc
            print(
                f"  ⚠ timeout on attempt {attempt}/{retries + 1}, retrying ...",
                file=sys.stderr,
            )
            time.sleep(2 * attempt)
    print(f"  ✗ giving up after retries: {last_exc}", file=sys.stderr)
    return None


def _candidate_ids(term: str) -> list[str]:
    """Possible normalised biotoolsID candidates for a corpus term.

    bio.tools IDs are conventionally lower-case alnum + hyphen. We try
    a small set of common transforms; if none hit, the search fallback
    takes over.
    """
    t = term.strip()
    lowered = t.lower()
    seen: list[str] = []
    for c in (lowered, lowered.replace(" ", "-"), lowered.replace(" ", "")):
        if c and c not in seen:
            seen.append(c)
    return seen


def _normalize_for_match(s: str) -> str:
    """Collapse whitespace and hyphens so ``Bowtie2`` matches ``Bowtie 2`` and
    ``Cell Ranger`` matches ``cell-ranger``. Case-insensitive."""
    return "".join(ch for ch in s.strip().lower() if ch not in " -_")


def _accept_exact_name(
    term: str, payload: dict, axis: str
) -> Hit | None:
    """Convert a bio.tools tool payload into a ``Hit`` iff the tool ``name``
    matches the corpus term — either verbatim (case-insensitive) or after
    a whitespace / hyphen normalisation that catches stylistic variants
    (``Bowtie2`` vs ``Bowtie 2``)."""
    name = (payload.get("name") or "").strip()
    biotools_id = (payload.get("biotoolsID") or "").strip()
    if not name or not biotools_id:
        return None
    target = term.strip().lower()
    if name.lower() != target and _normalize_for_match(name) != _normalize_for_match(term):
        return None
    desc = (payload.get("description") or "").strip().replace("\n", " ")
    if len(desc) > 500:
        desc = desc[:497] + "..."
    return Hit(
        term=term,
        axis=axis,
        biotools_id=biotools_id,
        name=name,
        description=desc,
    )


def resolve_via_direct_id(
    term: str, axis: str, *, timeout: int = 30, retries: int = 2
) -> Hit | None:
    for cand in _candidate_ids(term):
        r = _retrying_get(
            f"{BIOTOOLS_API}/{cand}/",
            params={"format": "json"},
            timeout=timeout,
            retries=retries,
        )
        time.sleep(0.3)
        if r is None or r.status_code != 200:
            continue
        hit = _accept_exact_name(term, r.json(), axis)
        if hit is not None:
            return hit
    return None


def resolve_via_search(
    term: str, axis: str, *, timeout: int = 30, retries: int = 2,
    max_pages: int = 3,
) -> Hit | None:
    """Search-by-name, iterate up to ``max_pages`` of results, accept first
    entry whose ``name`` matches exactly (case-insensitive)."""
    page = 1
    while page <= max_pages:
        r = _retrying_get(
            f"{BIOTOOLS_API}/tool/",
            params={"name": term, "format": "json", "page": str(page)},
            timeout=timeout,
            retries=retries,
        )
        time.sleep(0.3)
        if r is None or r.status_code != 200:
            return None
        payload = r.json()
        tools = payload.get("list") or []
        for t in tools:
            hit = _accept_exact_name(term, t, axis)
            if hit is not None:
                return hit
        if not payload.get("next"):
            return None
        page += 1
    return None


def resolve_term(term: str, axis: str) -> Hit | None:
    hit = resolve_via_direct_id(term, axis)
    if hit is not None:
        return hit
    return resolve_via_search(term, axis)


def load_unresolved_terms(path: Path) -> list[tuple[str, str]]:
    """Parse the OLS4 harvester's unresolved log: ``<axis>\\t<term>`` per line,
    skipping blank lines and comments."""
    if not path.exists():
        return []
    out: list[tuple[str, str]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "\t" not in line:
            continue
        axis, term = line.split("\t", 1)
        out.append((term.strip(), axis.strip()))
    return out


def build(
    unresolved_path: Path,
    out_tsv: Path,
    still_unresolved_txt: Path,
) -> int:
    candidates = load_unresolved_terms(unresolved_path)
    software = [(t, a) for t, a in candidates if a in SOFTWARE_AXES]
    skipped = len(candidates) - len(software)
    print(
        f"[biotools] {len(candidates)} terms in {unresolved_path.name}; "
        f"{len(software)} are software (will query bio.tools); "
        f"{skipped} are non-software (will pass through unchanged).",
        file=sys.stderr,
    )

    hits: list[Hit] = []
    still_unresolved: list[tuple[str, str]] = []

    for term, axis in software:
        print(f"[biotools] resolving {term!r} ({axis}) ...", file=sys.stderr)
        hit = resolve_term(term, axis)
        if hit is None:
            print("  ✗ no exact-name match in bio.tools", file=sys.stderr)
            still_unresolved.append((term, axis))
            continue
        print(
            f"  ✓ biotools:{hit.biotools_id} {hit.name!r}",
            file=sys.stderr,
        )
        hits.append(hit)

    # Non-software terms pass through to the still-unresolved log so the
    # next iteration (alias tables, upstream contributions) can pick them up.
    for term, axis in candidates:
        if axis not in SOFTWARE_AXES:
            still_unresolved.append((term, axis))

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLUMNS, delimiter="\t")
        w.writeheader()
        for h in hits:
            w.writerow({
                "term": h.term.lower(),
                "curie": f"biotools:{h.biotools_id}",
                "label": h.name,
                "description": h.description,
                "ontology": "biotools",
                "source": f"bio.tools registry ({h.axis})",
            })

    with still_unresolved_txt.open("w", encoding="utf-8") as fh:
        fh.write(
            "# Terms still unresolved after OLS4 + bio.tools harvests.\n"
            "# Candidates for the v0.3.3 instrument-series alias table\n"
            "# (NovaSeq / HiSeq / PacBio / Oxford Nanopore...), upstream\n"
            "# ontology contributions, or hand-mapping.\n"
            "# Format: <axis>\\t<term>\n\n"
        )
        for term, axis in still_unresolved:
            fh.write(f"{axis}\t{term}\n")

    print(
        f"\n[out] wrote {out_tsv}: {len(hits)} resolved via bio.tools, "
        f"{len(still_unresolved)} still unresolved (of {len(candidates)} "
        f"input terms, {len(software)} software).",
        file=sys.stderr,
    )
    return len(hits)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--unresolved-in", type=Path,
        default=Path(__file__).parent / "corpus_concept_unresolved.txt",
        help="Unresolved-terms log from the OLS4 harvester.",
    )
    parser.add_argument(
        "--out", type=Path,
        default=Path(__file__).parent / "biotools_concept_index.tsv",
        help="Output TSV path (default: bundled location).",
    )
    parser.add_argument(
        "--unresolved-out", type=Path,
        default=Path(__file__).parent / "biotools_unresolved.txt",
        help="Output log for terms still unresolved after both harvests.",
    )
    args = parser.parse_args(argv)
    build(args.unresolved_in, args.out, args.unresolved_out)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
