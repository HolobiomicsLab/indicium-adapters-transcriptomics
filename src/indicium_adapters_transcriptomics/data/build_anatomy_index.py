"""Maintainer script: regenerate the bundled ``anatomy_index.tsv``.

Takes a curated list of anatomy terms (drawn from the Cheung et al. 2024
Tomo-seq paper on *Nematostella vectensis* whole-body remodelling and
adjacent general anatomy) and searches the EBI OLS4 UBERON index for each.
Only entries whose returned **UBERON label is the exact term (or starts
with it, picking the shortest such label)** are written to the TSV. Loose
matches that promote a hyponym (e.g. *"nasal tentacle"* instead of
*"tentacle"*) are rejected.

The output schema (kept narrow on purpose, anatomy lookups don't need
the full HGNC-style cross-reference matrix):

    term  uberon  label  description  source

``term`` is the lower-cased lookup key. ``source`` records which paper /
context the term was selected from so the bundle stays auditable.

Usage (maintainer-only, not run at install or import time)::

    python -m indicium_adapters_transcriptomics.data.build_anatomy_index

Re-run when new reference papers introduce additional anatomy terms.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import NamedTuple

import requests

OLS_SEARCH = "https://www.ebi.ac.uk/ols4/api/search"

OUT_COLUMNS = ("term", "uberon", "label", "description", "source")


class Query(NamedTuple):
    term: str           # lower-cased lookup key + the OLS4 query string
    source: str         # short string identifying provenance


# Anatomy terms drawn from Cheung et al. 2024 plus a couple of generic
# concepts useful as anatomical-region anchors. Each is searched against
# UBERON; entries that don't yield a label exactly matching (or starting
# with) the query are dropped; see TRIAGED_OUT at the bottom of this
# file for terms known not to have a UBERON entry.
CURATED_TERMS: tuple[Query, ...] = (
    Query("mesentery", "Cheung 2024"),
    Query("pharynx", "Cheung 2024"),
    Query("tentacle", "Cheung 2024"),
    Query("basement membrane", "Cheung 2024"),
    Query("body wall", "Cheung 2024"),
    Query("extracellular matrix", "Cheung 2024"),
    Query("mesoglea", "Cheung 2024"),
    Query("epidermis", "Cheung 2024"),
    Query("gastrodermis", "Cheung 2024"),
    Query("body cavity", "general"),
    Query("organism subdivision", "general"),
    Query("anatomical region", "general"),
)


TRIAGED_OUT_KNOWN_CNIDARIAN: tuple[tuple[str, str], ...] = (
    ("physa", "Cnidarian-specific (aboral end); no UBERON term as of 2026."),
    ("oral pole", "Cnidarian polar-axis term; no UBERON term."),
    ("aboral pole", "Cnidarian polar-axis term; no UBERON term."),
    ("body column", "Cnidarian-specific body-axis subdivision; no UBERON term."),
    ("polyp", "Cnidarian body form (life-cycle stage), not an anatomical region in UBERON."),
    ("planula", "Cnidarian larval stage; UBERON has no clean entry."),
)


def _extract(doc: dict) -> dict[str, str]:
    desc = doc.get("description") or ""
    if isinstance(desc, list):
        desc = desc[0] if desc else ""
    return {
        "uberon": doc.get("obo_id") or "",
        "label": doc.get("label") or "",
        "description": desc,
    }


def resolve_term(q: Query) -> dict[str, str] | None:
    """Search OLS4 UBERON for ``q.term`` and return only an exact label
    match (case-insensitive), else ``None``. The strict acceptance rule
    means hyponym promotion (e.g. *"nasal tentacle"* for query *"tentacle"*)
    cannot slip through.
    """
    params = {
        "q": q.term,
        "ontology": "uberon",
        "type": "class",
        "rows": 50,
    }
    r = requests.get(OLS_SEARCH, params=params, timeout=30)
    r.raise_for_status()
    docs = [
        d for d in r.json().get("response", {}).get("docs", [])
        if (d.get("obo_id") or "").startswith("UBERON:")
    ]
    needle = q.term.lower()
    for doc in docs:
        if (doc.get("label") or "").lower() == needle:
            return _extract(doc)
    return None


def build(out_path: Path) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_ok = 0
    n_unmatched: list[str] = []
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLUMNS, delimiter="\t")
        w.writeheader()
        for q in CURATED_TERMS:
            print(f"[ols4] resolving {q.term!r} ...", file=sys.stderr)
            hit = resolve_term(q)
            if hit is None:
                print("[ols4]   ✗ no exact UBERON label match", file=sys.stderr)
                n_unmatched.append(q.term)
                continue
            w.writerow({
                "term": q.term,
                "uberon": hit["uberon"],
                "label": hit["label"],
                "description": hit["description"],
                "source": q.source,
            })
            n_ok += 1
            print(f"[ols4]   ✓ {hit['uberon']} ({hit['label']!r})", file=sys.stderr)

    print(
        f"\n[out] wrote {out_path}: {n_ok} matched, "
        f"{len(n_unmatched)} unmatched (of {len(CURATED_TERMS)} curated)",
        file=sys.stderr,
    )
    if n_unmatched:
        print(
            f"[note] {len(n_unmatched)} curated terms had no exact UBERON "
            f"label match and were dropped from the bundle:",
            file=sys.stderr,
        )
        for t in n_unmatched:
            print(f"         - {t}", file=sys.stderr)
    if TRIAGED_OUT_KNOWN_CNIDARIAN:
        print(
            f"[note] {len(TRIAGED_OUT_KNOWN_CNIDARIAN)} cnidarian-specific "
            f"terms deliberately omitted from the curated set "
            f"(no UBERON entry exists):",
            file=sys.stderr,
        )
        for term, reason in TRIAGED_OUT_KNOWN_CNIDARIAN:
            print(f"         - {term}: {reason}", file=sys.stderr)
    return n_ok


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--out", type=Path,
        default=Path(__file__).parent / "anatomy_index.tsv",
        help="Output anatomy_index.tsv path (default: bundled location).",
    )
    args = parser.parse_args(argv)
    build(args.out)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
