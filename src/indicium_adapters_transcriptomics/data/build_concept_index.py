"""Maintainer script, regenerate the bundled ``concept_index.tsv``.

Verifies a curated list of **analysis-, regulation- and protocol-concept
CURIEs** against the EBI OLS4 API and writes a TSV with the canonical
label and definition each ontology ships. The candidate set is drawn
directly from the ``TranscriptoKG`` mapping & vocabulary work (Boukhajou,
2026, every CURIE was already verified against OntoBee on 2026-05-19
during the KG-modelling sprint and is reused here for adapter-level
alignment).

The output schema:

    term  curie  label  description  ontology  source

``term`` is the lower-cased lookup key; ``ontology`` is the ontology
prefix (``NCIT`` / ``INO`` / ``OBI`` / ``STATO`` / ``GO`` / ``IAO`` /
``SIO``) so consumers can filter by axis.

Usage (maintainer-only, not run at install or import time)::

    python -m indicium_adapters_transcriptomics.data.build_concept_index

Re-run when the TranscriptoKG vocab grows or when an OLS4 label changes.
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path
from typing import NamedTuple

import requests

OLS_TERM_BY_IRI = "https://www.ebi.ac.uk/ols4/api/ontologies"

OUT_COLUMNS = ("term", "curie", "label", "description", "ontology", "source")


class Candidate(NamedTuple):
    curie: str        # e.g. "NCIT:C60693"
    ontology: str     # OLS4 ontology key, lowercase: "ncit", "ino", "obi", ...
    source: str       # short string identifying provenance


# CURIEs extracted from TranscriptoKG mapping + vocab.
# Every CURIE was verified against OntoBee during KG-modelling on 2026-05-19;
# the build script re-verifies via OLS4 to confirm the term still exists,
# and fills in canonical label + description at build time.
CANDIDATES: tuple[Candidate, ...] = (
    # --- INO, gene expression regulation processes ---
    Candidate("INO:0000094", "ino", "TranscriptoKG vocab (up-regulated)"),
    Candidate("INO:0000095", "ino", "TranscriptoKG vocab (down-regulated)"),
    # --- NCIT, analysis concepts and statistical outcomes ---
    Candidate("NCIT:C60693",  "ncit", "TranscriptoKG vocab (Fold Change)"),
    Candidate("NCIT:C64200",  "ncit", "TranscriptoKG vocab (log2)"),
    Candidate("NCIT:C68568",  "ncit", "TranscriptoKG vocab (Timepoint)"),
    Candidate("NCIT:C82576",  "ncit", "TranscriptoKG vocab (Reference Time Point)"),
    Candidate("NCIT:C201564", "ncit", "TranscriptoKG vocab (Not Significant)"),
    # --- NCIT, entity types used in TranscriptoKG mapping ---
    Candidate("NCIT:C14250",  "ncit", "TranscriptoKG mapping"),
    Candidate("NCIT:C16438",  "ncit", "TranscriptoKG mapping"),
    Candidate("NCIT:C61601",  "ncit", "TranscriptoKG mapping"),
    Candidate("NCIT:C70699",  "ncit", "TranscriptoKG mapping"),
    Candidate("NCIT:C73517",  "ncit", "TranscriptoKG mapping"),
    Candidate("NCIT:C93638",  "ncit", "TranscriptoKG mapping"),
    Candidate("NCIT:C142401", "ncit", "TranscriptoKG mapping"),
    Candidate("NCIT:C142485", "ncit", "TranscriptoKG mapping"),
    Candidate("NCIT:C164894", "ncit", "TranscriptoKG mapping (pathway id predicate)"),
    Candidate("NCIT:C164895", "ncit", "TranscriptoKG mapping (pathway name predicate)"),
    Candidate("NCIT:C175233", "ncit", "TranscriptoKG mapping"),
    Candidate("NCIT:C177610", "ncit", "TranscriptoKG mapping"),
    Candidate("NCIT:C188840", "ncit", "TranscriptoKG mapping"),
    Candidate("NCIT:C211600", "ncit", "TranscriptoKG mapping"),
    Candidate("NCIT:C213836", "ncit", "TranscriptoKG mapping"),
    # --- GO, three roots, useful as anchors ---
    Candidate("GO:0003674", "go", "TranscriptoKG mapping (molecular function root)"),
    Candidate("GO:0005575", "go", "TranscriptoKG mapping (cellular component root)"),
    Candidate("GO:0008150", "go", "TranscriptoKG mapping (biological process root)"),
    # --- OBI, investigation / role / protocol / DE transformation / assay ---
    # OBI:0000249 and OBI:0000650 were originally tagged "assay" and "data set"
    # in the TranscriptoKG vocab, but OLS4 verification (2026-05-27) shows
    # they actually denote "technical replicate role" and "differential
    # expression analysis data transformation", both still relevant to the
    # transcriptomics domain so we keep them, but the source notes reflect
    # the OLS4 truth.
    Candidate("OBI:0000066", "obi", "TranscriptoKG mapping (investigation)"),
    Candidate("OBI:0000249", "obi", "TranscriptoKG mapping (technical replicate role)"),
    Candidate("OBI:0000272", "obi", "TranscriptoKG mapping (protocol)"),
    Candidate("OBI:0000650", "obi", "TranscriptoKG mapping (DE analysis data transformation)"),
    Candidate("OBI:0001177", "obi", "TranscriptoKG mapping (RNA sequencing assay)"),
    # --- STATO, statistical method anchor ---
    Candidate("STATO:0000251", "stato", "TranscriptoKG mapping"),
)


def _curie_to_iri(curie: str) -> str:
    """Convert CURIE ``NS:LOCAL`` to its OBO PURL form."""
    ns, local = curie.split(":", 1)
    return f"http://purl.obolibrary.org/obo/{ns}_{local}"


def verify_curie(
    c: Candidate, *, timeout: int = 90, retries: int = 2
) -> dict[str, str] | None:
    """Return ``{"label": ..., "description": ...}`` for an OLS4-verified
    CURIE, or ``None`` when the term is not found in its declared ontology.

    The EBI endpoint can be slow on large ontologies (OBI, NCIT take 30-60s
    in cold cache), so the per-request timeout is generous and the call
    retries once on ``ReadTimeout``."""
    iri = _curie_to_iri(c.curie)
    url = f"{OLS_TERM_BY_IRI}/{c.ontology}/terms"
    last_exc: Exception | None = None
    for attempt in range(1, retries + 2):  # initial + retries
        try:
            r = requests.get(url, params={"iri": iri}, timeout=timeout)
            if r.status_code == 404:
                return None
            # Retry on 5xx (service unavailable / gateway timeout etc.).
            if 500 <= r.status_code < 600:
                last_exc = requests.HTTPError(f"HTTP {r.status_code}")
                print(
                    f"[ols4]   ⚠ HTTP {r.status_code} on attempt "
                    f"{attempt}/{retries + 1}, retrying ...",
                    file=sys.stderr,
                )
                time.sleep(2 * attempt)
                continue
            r.raise_for_status()
            terms = (r.json().get("_embedded") or {}).get("terms") or []
            if not terms:
                return None
            t = terms[0]
            desc = t.get("description") or []
            if isinstance(desc, list):
                desc = desc[0] if desc else ""
            return {"label": t.get("label") or "", "description": desc}
        except requests.exceptions.Timeout as exc:
            last_exc = exc
            print(
                f"[ols4]   ⚠ timeout on attempt {attempt}/{retries + 1}, "
                f"retrying ...",
                file=sys.stderr,
            )
            time.sleep(2 * attempt)
    # All retries exhausted.
    print(f"[ols4]   ✗ failed after retries: {last_exc}", file=sys.stderr)
    return None


def build(out_path: Path) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_ok = 0
    n_invalid: list[str] = []
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLUMNS, delimiter="\t")
        w.writeheader()
        for c in CANDIDATES:
            print(
                f"[ols4] verifying {c.curie} in {c.ontology.upper()} ...",
                file=sys.stderr,
            )
            info = verify_curie(c, retries=3)
            # Polite throttle between requests, EBI flapped 503s during
            # the v0.3.0 build.
            time.sleep(0.3)
            if info is None:
                print("[ols4]   ✗ not found", file=sys.stderr)
                n_invalid.append(c.curie)
                continue
            # The lookup term is the canonical label, lower-cased, that gives
            # us a stable, OLS4-sourced key (no risk of getting the label
            # wrong by hand).
            term = info["label"].lower()
            w.writerow({
                "term": term,
                "curie": c.curie,
                "label": info["label"],
                "description": info["description"],
                "ontology": c.ontology.upper(),
                "source": c.source,
            })
            n_ok += 1
            print(
                f"[ols4]   ✓ {info['label']!r}",
                file=sys.stderr,
            )

    print(
        f"\n[out] wrote {out_path}: {n_ok} verified, {len(n_invalid)} invalid "
        f"(of {len(CANDIDATES)} candidates)",
        file=sys.stderr,
    )
    if n_invalid:
        print(
            f"[note] {len(n_invalid)} CURIEs failed verification "
            f"(missing from declared ontology in OLS4):",
            file=sys.stderr,
        )
        for c in n_invalid:
            print(f"         - {c}", file=sys.stderr)
    return n_ok


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--out", type=Path,
        default=Path(__file__).parent / "concept_index.tsv",
        help="Output concept_index.tsv path (default: bundled location).",
    )
    args = parser.parse_args(argv)
    build(args.out)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
