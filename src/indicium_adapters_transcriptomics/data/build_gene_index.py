"""Maintainer script — regenerate the bundled ``gene_index.tsv``.

Joins three authoritative sources:

* The **G2V2 cross-reference table** (Nematostella vectensis).
* The **HGNC complete set** (HUGO Gene Nomenclature Committee, human).
* The **MGI marker list** (Mouse Genome Informatics, MRK_List2.rpt) —
  added in v0.6.0 to give Tabula Muris-style mouse claims a path to
  resolution.

Output schema (v0.6.0)::

    symbol  hgnc  mgi  ensembl  ncbigene  nemve1  nvertx  uniprot  species  description

The ``species`` column (``human`` / ``mouse`` / ``nematostella``)
drives the v0.6.0 species-aware lookup: callers may pass a species
hint on the claim, and a casing heuristic infers species from the
term itself (UPPER → human, Title-case → mouse) when no hint is
given. Disambiguation rule lives in ``GeneIndex.lookup``.

Usage (maintainer-only — not run at install or import time)::

    python -m indicium_adapters_transcriptomics.data.build_gene_index \\
        --g2v2-xlsx /path/to/G2V2_*.xlsx \\
        [--hgnc-tsv /path/to/hgnc_complete_set.txt] \\
        [--mgi-tsv /path/to/MRK_List2.rpt]

Re-run whenever G2V2, HGNC or MGI release a new version.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import pandas as pd
import requests

HGNC_URL = (
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/"
    "hgnc_complete_set.txt"
)
HGNC_URL_FALLBACK = (
    "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt"
)

OUT_COLUMNS = (
    "symbol", "hgnc", "mgi", "ensembl", "ncbigene",
    "nemve1", "nvertx", "uniprot", "species", "description",
)

MGI_URL = "https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt"


def _first_token(value: object, separators: str = ",|;") -> str:
    """Return the first non-empty token of a multi-value string, stripped."""
    if value is None or pd.isna(value):
        return ""
    s = str(value)
    for sep in separators:
        if sep in s:
            for tok in s.split(sep):
                t = tok.strip()
                if t:
                    return t
            return ""
    return s.strip()


def _download_hgnc(target: Path) -> None:
    """Download the HGNC complete set TSV. Tries the GCS bucket first, then
    falls back to the EBI mirror."""
    target.parent.mkdir(parents=True, exist_ok=True)
    for url in (HGNC_URL, HGNC_URL_FALLBACK):
        print(f"[hgnc] downloading {url} ...", file=sys.stderr)
        try:
            r = requests.get(url, stream=True, timeout=120)
            r.raise_for_status()
            with target.open("wb") as fh:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        fh.write(chunk)
            print(f"[hgnc] wrote {target} ({target.stat().st_size:,} bytes)", file=sys.stderr)
            return
        except Exception as exc:  # pragma: no cover - network failure path
            print(f"[hgnc] {url} failed: {exc}", file=sys.stderr)
    raise RuntimeError(
        "Could not download HGNC complete set from either mirror; "
        "pass --hgnc-tsv with a local copy."
    )


def _download_mgi(target: Path) -> None:
    """Download the MGI marker list (MRK_List2.rpt)."""
    target.parent.mkdir(parents=True, exist_ok=True)
    print(f"[mgi] downloading {MGI_URL} ...", file=sys.stderr)
    r = requests.get(MGI_URL, stream=True, timeout=180)
    r.raise_for_status()
    with target.open("wb") as fh:
        for chunk in r.iter_content(chunk_size=1024 * 1024):
            if chunk:
                fh.write(chunk)
    print(f"[mgi] wrote {target} ({target.stat().st_size:,} bytes)", file=sys.stderr)


def _load_mgi(path: Path) -> list[dict[str, str]]:
    """Parse MGI MRK_List2.rpt and return mouse-gene rows.

    Filters to ``Marker Type == "Gene"`` to drop BAC/YAC ends, syntenic
    markers, RIKEN cDNAs without a stable gene status, etc. Returned
    rows carry the minimum fields needed by the bundled TSV
    (``symbol``, ``mgi``, ``name``); ``ensembl`` / ``ncbigene`` /
    ``uniprot`` cross-refs are left empty in v0.6.0 — they live in
    different MGI report files (MRK_ENSEMBL / MGI_EntrezGene /
    MRK_SwissProt) and integrating those is on the v0.6.1 roadmap.
    """
    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        low_memory=False,
        na_filter=False,
        on_bad_lines="warn",
    )
    needed = {"MGI Accession ID", "Marker Symbol", "Marker Type", "Marker Name"}
    missing = needed - set(df.columns)
    if missing:
        raise RuntimeError(
            f"MGI rpt is missing expected columns: {sorted(missing)}. "
            f"Got: {sorted(df.columns)[:20]}..."
        )
    df = df[df["Marker Type"].str.strip() == "Gene"]
    rows: list[dict[str, str]] = []
    for _, r in df.iterrows():
        symbol = (r["Marker Symbol"] or "").strip()
        mgi_id = (r["MGI Accession ID"] or "").strip()
        if not symbol or not mgi_id:
            continue
        rows.append({
            "symbol": symbol,
            "mgi": mgi_id,
            "name": (r["Marker Name"] or "").strip(),
        })
    print(
        f"[mgi] indexed {len(rows):,} mouse Gene rows (filtered from {len(df):,} total)",
        file=sys.stderr,
    )
    return rows


def _load_hgnc(
    path: Path,
) -> tuple[
    dict[str, dict[str, str]],
    dict[str, dict[str, str]],
    list[dict[str, str]],
]:
    """Return three views of the HGNC complete set:

    * ``by_uniprot``  — UniProt accession → HGNC row
    * ``by_symbol``   — UPPER-cased HGNC symbol → HGNC row (case-insensitive
      fallback when no UniProt accession matched)
    * ``all_rows``    — every HGNC row, used to emit pure-human rows for genes
      that don't appear in G2V2

    The returned row dicts use the runtime-TSV column names so the writer can
    consume them uniformly.
    """
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False, na_filter=False)
    needed = {"hgnc_id", "symbol", "name", "uniprot_ids", "ensembl_gene_id", "entrez_id"}
    missing = needed - set(df.columns)
    if missing:
        raise RuntimeError(
            f"HGNC TSV is missing expected columns: {sorted(missing)}. "
            f"Got: {sorted(df.columns)[:20]}..."
        )
    by_uniprot: dict[str, dict[str, str]] = {}
    by_symbol: dict[str, dict[str, str]] = {}
    all_rows: list[dict[str, str]] = []
    for _, r in df.iterrows():
        if not r["symbol"]:
            continue
        row = {
            "hgnc": r["hgnc_id"],
            "symbol": r["symbol"],
            "ensembl": r["ensembl_gene_id"],
            "ncbigene": r["entrez_id"],
            "name": r["name"],
            "uniprot": (r["uniprot_ids"] or "").split("|")[0].strip(),
        }
        all_rows.append(row)
        by_symbol.setdefault(r["symbol"].upper(), row)
        for uniprot in (r["uniprot_ids"] or "").split("|"):
            uniprot = uniprot.strip()
            if uniprot:
                by_uniprot.setdefault(uniprot, row)
    print(
        f"[hgnc] indexed {len(all_rows):,} rows | "
        f"{len(by_uniprot):,} UniProt keys | "
        f"{len(by_symbol):,} symbol keys",
        file=sys.stderr,
    )
    return by_uniprot, by_symbol, all_rows


def build(
    g2v2_path: Path, hgnc_path: Path, mgi_path: Path, out_path: Path
) -> int:
    """Read G2V2 + HGNC + MGI, join, emit the bundled TSV. Return row count.

    Output composition (in this load order — affects ``GeneIndex.lookup``
    fallback ordering):

    1. One row per de-duplicated G2V2 Nematostella gene
       (``species=nematostella``).
    2. One row per HGNC gene not already paired with a Nematostella row
       (``species=human``).
    3. One row per MGI mouse marker of type ``Gene``
       (``species=mouse``).
    """
    print(f"[g2v2] reading {g2v2_path} ...", file=sys.stderr)
    g2v2 = pd.read_excel(g2v2_path, dtype=str, na_filter=False)
    print(f"[g2v2] {len(g2v2):,} rows, columns: {list(g2v2.columns)[:6]}...", file=sys.stderr)

    hgnc_by_uniprot, hgnc_by_symbol, hgnc_all_rows = _load_hgnc(hgnc_path)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_nem_total = 0
    n_via_uniprot = 0
    n_via_symbol = 0
    used_hgnc_ids: set[str] = set()
    seen_jgi: set[str] = set()
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLUMNS, delimiter="\t")
        w.writeheader()

        # Step 1 — Nematostella-anchored rows.
        for _, r in g2v2.iterrows():
            jgi = _first_token(r.get("JGI"))
            nvertx = _first_token(r.get("NvERTx_ID"))
            uniprot = _first_token(r.get("Uniprot_ID_NvERTx"))
            preferred = _first_token(r.get("Preferred_name"))
            g2v2_desc = _first_token(r.get("Uniprot_Description_NvERTx"))
            if not (jgi or nvertx or uniprot):
                continue
            if jgi and jgi in seen_jgi:
                continue
            if jgi:
                seen_jgi.add(jgi)

            # Primary: UniProt → HGNC. Secondary: Preferred_name → HGNC symbol.
            hgnc_row = hgnc_by_uniprot.get(uniprot) if uniprot else None
            if hgnc_row:
                n_via_uniprot += 1
            elif preferred and (hgnc_row := hgnc_by_symbol.get(preferred.upper())):
                n_via_symbol += 1

            if hgnc_row:
                used_hgnc_ids.add(hgnc_row["hgnc"])

            symbol = (hgnc_row or {}).get("symbol") or preferred or ""
            description = (hgnc_row or {}).get("name") or g2v2_desc or ""
            w.writerow({
                "symbol": symbol,
                "hgnc": (hgnc_row or {}).get("hgnc", ""),
                "mgi": "",
                "ensembl": (hgnc_row or {}).get("ensembl", ""),
                "ncbigene": (hgnc_row or {}).get("ncbigene", ""),
                "nemve1": jgi,
                "nvertx": nvertx,
                "uniprot": uniprot,
                "species": "nematostella",
                "description": description,
            })
            n_nem_total += 1

        # Step 2 — pure-HGNC rows for human genes not yet covered.
        n_pure_hgnc = 0
        for row in hgnc_all_rows:
            if row["hgnc"] in used_hgnc_ids:
                continue
            w.writerow({
                "symbol": row["symbol"],
                "hgnc": row["hgnc"],
                "mgi": "",
                "ensembl": row["ensembl"],
                "ncbigene": row["ncbigene"],
                "nemve1": "",
                "nvertx": "",
                "uniprot": row["uniprot"],
                "species": "human",
                "description": row["name"],
            })
            n_pure_hgnc += 1

        # Step 3 — MGI mouse rows (v0.6.0).
        mgi_rows = _load_mgi(mgi_path)
        n_mgi = 0
        for row in mgi_rows:
            w.writerow({
                "symbol": row["symbol"],
                "hgnc": "",
                "mgi": row["mgi"],
                "ensembl": "",
                "ncbigene": "",
                "nemve1": "",
                "nvertx": "",
                "uniprot": "",
                "species": "mouse",
                "description": row["name"],
            })
            n_mgi += 1

    n_total = n_nem_total + n_pure_hgnc + n_mgi
    print(
        f"[out] wrote {out_path}: {n_total:,} rows total\n"
        f"       Nematostella-anchored : {n_nem_total:,}\n"
        f"         via UniProt join    : {n_via_uniprot:,}\n"
        f"         via symbol fallback : {n_via_symbol:,}\n"
        f"         Nematostella-only   : {n_nem_total - n_via_uniprot - n_via_symbol:,}\n"
        f"       Pure-HGNC (no Nemve1) : {n_pure_hgnc:,}\n"
        f"       MGI (mouse)           : {n_mgi:,}",
        file=sys.stderr,
    )
    return n_total


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--g2v2-xlsx", type=Path, required=True,
        help="Path to the G2V2_*_NvERTx.4_JGI_ENSEMBL_NVE.xlsx cross-reference.",
    )
    parser.add_argument(
        "--hgnc-tsv", type=Path, default=None,
        help="Local HGNC complete set TSV. If omitted, downloaded to a cache.",
    )
    parser.add_argument(
        "--out", type=Path,
        default=Path(__file__).parent / "gene_index.tsv",
        help="Output gene_index.tsv path (default: bundled location).",
    )
    cache_dir = (
        Path.home() / ".cache" / "indicium-adapters-transcriptomics"
    )
    parser.add_argument(
        "--hgnc-cache", type=Path,
        default=cache_dir / "hgnc_complete_set.txt",
        help="Where to cache the downloaded HGNC TSV.",
    )
    parser.add_argument(
        "--mgi-tsv", type=Path, default=None,
        help="Local MGI MRK_List2.rpt. If omitted, downloaded to a cache.",
    )
    parser.add_argument(
        "--mgi-cache", type=Path,
        default=cache_dir / "MRK_List2.rpt",
        help="Where to cache the downloaded MGI report.",
    )
    args = parser.parse_args(argv)
    hgnc_path = args.hgnc_tsv
    if hgnc_path is None:
        if not args.hgnc_cache.is_file():
            _download_hgnc(args.hgnc_cache)
        hgnc_path = args.hgnc_cache
    mgi_path = args.mgi_tsv
    if mgi_path is None:
        if not args.mgi_cache.is_file():
            _download_mgi(args.mgi_cache)
        mgi_path = args.mgi_cache
    build(args.g2v2_xlsx, hgnc_path, mgi_path, args.out)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
