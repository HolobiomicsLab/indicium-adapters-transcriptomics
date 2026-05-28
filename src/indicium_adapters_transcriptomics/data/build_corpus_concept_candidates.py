"""Maintainer script — harvest concept CURIEs from the mini-corpus extraction.

This script reads the per-paper / aggregate vocabulary JSON produced by
``articles/_extract_vocab.py`` over the 6-paper mini-corpus and uses the
**OLS4 search-by-label** endpoint to *propose* candidate CURIEs for each
extracted methodological term (STAR, Seurat, MARS-seq, Illumina, 10x …).

Every proposed CURIE is then re-verified by IRI lookup before being
accepted.

**Two-tier matching** (per axis): each axis declares both a preferred
ontology list *and* a match mode that tunes the strictness:

- ``exact`` — label or synonym must equal the corpus term verbatim
  (case-insensitive). Used as the default for ``de_tool`` / ``sc_tool``
  / ``aligner_quant`` axes where false positives are dangerous.
- ``prefix`` — corpus term must be the whole-word **prefix** of the
  OLS4 label (case-insensitive); allows descriptive suffixes
  ("Wilcoxon" → "Wilcoxon Signed Rank Test"; "bulk RNA-seq" → "bulk
  RNA-seq assay"; "Visium" → "Visium Spatial Gene Expression").
- ``substring`` — corpus term must appear as a whole-word **substring**
  of the OLS4 label (case-insensitive); used for instrument-series
  terms whose canonical EFO label prepends the vendor name
  ("NovaSeq" → "Illumina NovaSeq 6000"). Resolution picks the
  *first* OLS4 search result that matches — OLS4 search is ranked by
  relevance so the top hit is typically the canonical / oldest model.

The exact / prefix / substring modes are tried in cascade per axis;
the first hit wins. Hand-curation of CURIEs is still avoided — the
only curation is which mode applies to which axis.

The output schema matches ``concept_index.tsv``::

    term  curie  label  description  ontology  source

The ``source`` column tags the match mode used so the provenance is
auditable (e.g. ``corpus extraction (platform, substring-match)``).

Usage (maintainer-only — not run at install or import time)::

    python -m indicium_adapters_transcriptomics.data.build_corpus_concept_candidates

Re-run when the mini-corpus grows or when OLS4 indexes new tools / methods.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
from pathlib import Path
from typing import Literal, NamedTuple

import requests

MatchMode = Literal["exact", "prefix", "substring"]

OLS_SEARCH = "https://www.ebi.ac.uk/ols4/api/search"
OLS_TERM_BY_IRI = "https://www.ebi.ac.uk/ols4/api/ontologies"

OUT_COLUMNS = ("term", "curie", "label", "description", "ontology", "source")

# Default location of the vocab summary produced by articles/_extract_vocab.py.
DEFAULT_VOCAB_SUMMARY = (
    Path(__file__).resolve().parents[3]
    / "articles"
    / "_vocab_summary.json"
)


# Axis → ordered tuple of OLS4 ontology keys (lowercase) to search.
# The order encodes a *preference* — the first ontology that yields a
# match (under the axis's match modes) wins. This is one of the two
# places where domain knowledge is encoded (the other is the per-axis
# MATCH_MODES dict below). Specific CURIEs are NEVER hand-picked.
AXIS_ONTOLOGIES: dict[str, tuple[str, ...]] = {
    "platform":       ("efo", "obi", "ncit"),
    "protocol":       ("efo", "obi", "ncit"),
    "aligner_quant":  ("edam", "obi"),
    "de_tool":        ("edam", "obi"),
    "sc_tool":        ("edam", "obi"),
    "spatial_tool":   ("edam", "obi"),
    "long_read_tool": ("edam", "obi"),
    "qc":             ("edam", "obi"),
    "stat_test":      ("ncit", "stato", "obi"),
}


# Axis → cascade of match modes to try, in order. The first hit wins.
# - ``exact``: label or synonym equals the corpus term verbatim.
# - ``prefix``: label starts with the corpus term followed by a
#   whitespace / hyphen / end-of-string (whole-word prefix).
# - ``substring``: corpus term appears as a whole-word substring of
#   the label. Used for instrument families whose canonical EFO label
#   prepends a vendor name.
#
# Software axes (de_tool, sc_tool, aligner_quant, etc.) stay on
# ``exact`` only — they're handled by bio.tools downstream, and a
# loose match in EDAM/OBI would invent a CURIE for a different concept.
MATCH_MODES: dict[str, tuple[MatchMode, ...]] = {
    "platform":       ("exact", "prefix", "substring"),
    "protocol":       ("exact", "prefix"),
    "stat_test":      ("exact", "prefix"),
    "aligner_quant":  ("exact",),
    "de_tool":        ("exact",),
    "sc_tool":        ("exact",),
    "spatial_tool":   ("exact",),
    "long_read_tool": ("exact",),
    "qc":             ("exact",),
}


class Hit(NamedTuple):
    term: str            # original corpus term, exactly as extracted
    axis: str            # platform / protocol / ... (provenance)
    mode: MatchMode      # which match mode produced the hit
    curie: str           # OLS4-returned CURIE (e.g. "EFO:0008637")
    ontology: str        # uppercase ontology prefix
    label: str           # OLS4-canonical label
    description: str     # OLS4 definition (or "" if none)


def _retrying_get(
    url: str, *, params: dict[str, str], timeout: int, retries: int
) -> requests.Response | None:
    """GET with exponential-backoff retry on 5xx and timeouts.

    Returns the final ``Response`` on success (any 2xx / 4xx) or ``None``
    when the request kept failing past the retry budget."""
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
                f"  ⚠ timeout on attempt {attempt}/{retries + 1}, "
                f"retrying ...",
                file=sys.stderr,
            )
            time.sleep(2 * attempt)
    print(f"  ✗ giving up after retries: {last_exc}", file=sys.stderr)
    return None


def _search_ontology(
    term: str, ontology: str, *, exact: bool, rows: int = 10,
    timeout: int = 90, retries: int = 2,
) -> list[dict] | None:
    """OLS4 search restricted to ``ontology``. ``exact=True`` adds
    ``exact=true&queryFields=label,synonym``; ``False`` is the relaxed
    relevance search used by the prefix / substring modes.

    Returns the list of result docs (possibly empty) or ``None`` on
    network failure."""
    params: dict[str, str] = {
        "q": term,
        "ontology": ontology,
        "rows": str(rows),
    }
    if exact:
        params["exact"] = "true"
        params["queryFields"] = "label,synonym"
    r = _retrying_get(OLS_SEARCH, params=params, timeout=timeout, retries=retries)
    if r is None or r.status_code != 200:
        return None
    return ((r.json().get("response") or {}).get("docs") or [])


_WORD_BOUNDARY_AFTER_TERM = re.compile(r"[\s\-]")


def _prefix_match(corpus_term: str, label: str) -> bool:
    """``corpus_term`` is a whole-word prefix of ``label`` (case-insensitive)."""
    term_lc = corpus_term.strip().lower()
    label_lc = label.strip().lower()
    if not term_lc or not label_lc:
        return False
    if label_lc == term_lc:
        return True
    if not label_lc.startswith(term_lc):
        return False
    boundary = label_lc[len(term_lc):][:1]
    return bool(_WORD_BOUNDARY_AFTER_TERM.match(boundary))


def _substring_match(corpus_term: str, label: str) -> bool:
    """``corpus_term`` appears in ``label`` bounded by whitespace, hyphens or
    string ends (case-insensitive whole-word substring match)."""
    term_lc = corpus_term.strip().lower()
    label_lc = label.strip().lower()
    if not term_lc or not label_lc:
        return False
    pat = r"(?:^|[\s\-])" + re.escape(term_lc) + r"(?:[\s\-]|$)"
    return re.search(pat, label_lc) is not None


def verify_iri(
    iri: str, ontology: str, *, timeout: int = 90, retries: int = 2
) -> dict[str, str] | None:
    """Re-fetch a term by IRI from its declared ontology and return the
    authoritative label + description."""
    r = _retrying_get(
        f"{OLS_TERM_BY_IRI}/{ontology}/terms",
        params={"iri": iri},
        timeout=timeout,
        retries=retries,
    )
    if r is None or r.status_code != 200:
        return None
    terms = (r.json().get("_embedded") or {}).get("terms") or []
    if not terms:
        return None
    t = terms[0]
    desc = t.get("description") or []
    if isinstance(desc, list):
        desc = desc[0] if desc else ""
    return {"label": t.get("label") or "", "description": desc or ""}


def _try_match_in_docs(
    docs: list[dict], term: str, mode: MatchMode
) -> dict | None:
    """Walk the OLS4 result docs in returned order; return the first doc
    that satisfies the requested match mode.

    For ``exact`` mode the matcher does a two-pass over the docs:
    label-exact first, synonym-exact second. This protects against an
    OLS4 ranking quirk: e.g. searching ``macrophage`` against CL returns
    CL:0000394 (plasmatocyte, which lists macrophage as a synonym) ahead
    of CL:0000235 (whose label *is* ``macrophage``). The canonical entry
    must win even if OLS4's relevance score puts it second."""
    target = term.strip().lower()
    if mode == "exact":
        # Pass 1: label-exact.
        for d in docs:
            label = (d.get("label") or "").strip()
            if label.lower() == target:
                return d
        # Pass 2: synonym-exact. The OLS4 *search* endpoint exposes
        # synonyms under ``exact_synonyms`` (not the term endpoint's
        # ``synonym`` field); be defensive and check both, since some
        # documents may expose ``synonyms`` instead.
        for d in docs:
            synonyms: list[str] = []
            for key in ("exact_synonyms", "synonym", "synonyms"):
                val = d.get(key)
                if isinstance(val, list):
                    synonyms.extend(s.strip() for s in val if isinstance(s, str))
                elif isinstance(val, str):
                    synonyms.append(val.strip())
            if any(s.lower() == target for s in synonyms):
                return d
        return None
    if mode == "prefix":
        for d in docs:
            label = (d.get("label") or "").strip()
            if _prefix_match(term, label):
                return d
        return None
    if mode == "substring":
        for d in docs:
            label = (d.get("label") or "").strip()
            if _substring_match(term, label):
                return d
        return None
    return None


def resolve_term(term: str, axis: str) -> Hit | None:
    """Cascade the axis's match modes across the axis's preferred ontologies.

    The cascade is depth-first over modes: for each mode, every preferred
    ontology is tried before moving to the next mode. This guarantees
    an exact match in *any* preferred ontology wins over a prefix /
    substring match in the *first* preferred ontology — strictness beats
    preference."""
    ontologies = AXIS_ONTOLOGIES.get(axis, ())
    modes = MATCH_MODES.get(axis, ("exact",))
    if not ontologies or not modes:
        return None

    for mode in modes:
        for ont in ontologies:
            # ``exact`` uses the OLS4 exact endpoint; relaxed modes use a
            # ranked relevance search so the most relevant candidates
            # bubble up first.
            docs = _search_ontology(term, ont, exact=(mode == "exact"))
            time.sleep(0.3)
            if not docs:
                continue
            d = _try_match_in_docs(docs, term, mode)
            if d is None:
                continue
            iri = d.get("iri")
            obo_id = d.get("obo_id") or d.get("short_form")
            if not iri or not obo_id:
                continue
            info = verify_iri(iri, ont)
            time.sleep(0.3)
            if info is None:
                continue
            return Hit(
                term=term,
                axis=axis,
                mode=mode,
                curie=obo_id,
                ontology=ont.upper(),
                label=info["label"] or (d.get("label") or ""),
                description=info["description"],
            )
    return None


def load_corpus_terms(vocab_summary_path: Path) -> list[tuple[str, str]]:
    """Return ``[(term, axis), ...]`` from the aggregate vocab summary,
    deduped on ``(term.lower(), axis)``."""
    data = json.loads(vocab_summary_path.read_text(encoding="utf-8"))
    totals: dict[str, dict[str, int]] = data.get("totals") or {}
    seen: set[tuple[str, str]] = set()
    out: list[tuple[str, str]] = []
    for axis, by_term in totals.items():
        for term in by_term:
            key = (term.lower(), axis)
            if key in seen:
                continue
            seen.add(key)
            out.append((term, axis))
    return out


def already_in_existing_index(
    candidates: list[tuple[str, str]], existing_tsv: Path
) -> set[str]:
    """Return the set of lower-cased terms that are already present in the
    static ``concept_index.tsv`` (so we don't duplicate rows)."""
    if not existing_tsv.exists():
        return set()
    with existing_tsv.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return {(row.get("term") or "").strip().lower() for row in reader}


def build(
    vocab_summary_path: Path,
    out_tsv: Path,
    unresolved_txt: Path,
    existing_tsv: Path,
) -> int:
    candidates = load_corpus_terms(vocab_summary_path)
    skip_lower = already_in_existing_index(candidates, existing_tsv)
    print(
        f"[corpus] {len(candidates)} unique (term, axis) pairs from "
        f"{vocab_summary_path.name}; "
        f"{sum(1 for t, _ in candidates if t.lower() in skip_lower)} already "
        f"in concept_index.tsv (will skip).",
        file=sys.stderr,
    )

    hits: list[Hit] = []
    unresolved: list[tuple[str, str]] = []

    for term, axis in candidates:
        if term.lower() in skip_lower:
            print(f"[skip]   {term} ({axis}) — already in concept_index.tsv",
                  file=sys.stderr)
            continue
        print(f"[ols4] resolving {term!r} ({axis}) ...", file=sys.stderr)
        hit = resolve_term(term, axis)
        if hit is None:
            modes = MATCH_MODES.get(axis, ("exact",))
            print(
                f"  ✗ no {'/'.join(modes)} match across "
                f"{AXIS_ONTOLOGIES.get(axis, ())}",
                file=sys.stderr,
            )
            unresolved.append((term, axis))
            continue
        print(
            f"  ✓ [{hit.mode}] {hit.curie} {hit.label!r} ({hit.ontology})",
            file=sys.stderr,
        )
        hits.append(hit)

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLUMNS, delimiter="\t")
        w.writeheader()
        for h in hits:
            w.writerow({
                "term": h.term.lower(),
                "curie": h.curie,
                "label": h.label,
                "description": h.description,
                "ontology": h.ontology,
                "source": f"corpus extraction ({h.axis}, {h.mode}-match)",
            })

    with unresolved_txt.open("w", encoding="utf-8") as fh:
        fh.write(
            "# Corpus terms with NO exact label/synonym match in any axis "
            "ontology.\n"
            "# Candidates for bio.tools integration, upstream ontology "
            "contribution, or alias-table curation in a later release.\n"
            "# Format: <axis>\\t<term>\n\n"
        )
        for term, axis in unresolved:
            fh.write(f"{axis}\t{term}\n")

    print(
        f"\n[out] wrote {out_tsv}: {len(hits)} resolved, "
        f"{len(unresolved)} unresolved (of {len(candidates)} candidates, "
        f"{len(candidates) - len(hits) - len(unresolved)} skipped).",
        file=sys.stderr,
    )
    return len(hits)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--vocab-summary", type=Path,
        default=DEFAULT_VOCAB_SUMMARY,
        help=f"Vocab summary JSON (default: {DEFAULT_VOCAB_SUMMARY}).",
    )
    parser.add_argument(
        "--out", type=Path,
        default=Path(__file__).parent / "corpus_concept_index.tsv",
        help="Output TSV path (default: bundled location).",
    )
    parser.add_argument(
        "--unresolved", type=Path,
        default=Path(__file__).parent / "corpus_concept_unresolved.txt",
        help="Output unresolved-terms text path.",
    )
    parser.add_argument(
        "--existing", type=Path,
        default=Path(__file__).parent / "concept_index.tsv",
        help="Existing concept_index.tsv to skip duplicates against.",
    )
    args = parser.parse_args(argv)
    build(args.vocab_summary, args.out, args.unresolved, args.existing)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
