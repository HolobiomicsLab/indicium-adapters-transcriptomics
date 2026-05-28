"""Maintainer script, regenerate the bundled ``gene_index.tsv``.

Joins three authoritative sources:

* The **G2V2 cross-reference table** (Nematostella vectensis).
* The **HGNC complete set** (HUGO Gene Nomenclature Committee, human).
* The **MGI marker list** (Mouse Genome Informatics, MRK_List2.rpt),
  added in v0.6.0 to give Tabula Muris-style mouse claims a path to
  resolution.

Output schema (v0.6.0)::

    symbol  hgnc  mgi  ensembl  ncbigene  nemve1  nvertx  uniprot  species  description

The ``species`` column (``human`` / ``mouse`` / ``nematostella``)
drives the v0.6.0 species-aware lookup: callers may pass a species
hint on the claim, and a casing heuristic infers species from the
term itself (UPPER → human, Title-case → mouse) when no hint is
given. Disambiguation rule lives in ``GeneIndex.lookup``.

Usage (maintainer-only, not run at install or import time)::

    python -m indicium_adapters_transcriptomics.data.build_gene_index \\
        --g2v2-xlsx /path/to/G2V2_*.xlsx \\
        [--hgnc-tsv /path/to/hgnc_complete_set.txt] \\
        [--mgi-tsv /path/to/MRK_List2.rpt]

Re-run whenever G2V2, HGNC or MGI release a new version.
"""
from __future__ import annotations

import argparse
import csv
import gzip
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
    "symbol", "hgnc", "mgi", "flybase", "wormbase", "zfin",
    "ensembl", "ncbigene", "nemve1", "nvertx", "uniprot",
    "species", "description",
)

MGI_URL = "https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt"

# NCBI Gene catalogue per-organism gene_info files (v0.6.2). One TSV
# per species, gzip-compressed, with a consistent column layout:
#
#   #tax_id  GeneID  Symbol  LocusTag  Synonyms  dbXrefs  chromosome
#   map_location  description  type_of_gene  Symbol_from_nomenclature_authority
#   Full_name_from_nomenclature_authority  Nomenclature_status
#   Other_designations  Modification_date  Feature_type
#
# We read GeneID (Entrez), Symbol, Symbol_from_nomenclature_authority,
# description and dbXrefs. dbXrefs is pipe-separated and contains
# entries like ``FLYBASE:FBgn0040373`` / ``WormBase:WBGene00022277`` /
# ``ZFIN:ZDB-GENE-980526-104`` / ``Ensembl:ENSDARG00000021948``.
# UniProt per-organism id-mapping files (v0.6.3). Tab-separated
# ``<UniProtKB-AC>\t<TYPE>\t<ID>``; we filter to TYPE = FlyBase /
# WormBase / ZFIN to build a ``{organism id → UniProt accession}``
# dict and merge it into the gene_index rows.
UNIPROT_BY_ORGANISM_URLS: dict[str, str] = {
    "drosophila": (
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
        "knowledgebase/idmapping/by_organism/DROME_7227_idmapping.dat.gz"
    ),
    "celegans": (
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
        "knowledgebase/idmapping/by_organism/CAEEL_6239_idmapping.dat.gz"
    ),
    "zebrafish": (
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
        "knowledgebase/idmapping/by_organism/DANRE_7955_idmapping.dat.gz"
    ),
}

# Which idmapping TYPE column we read for each organism.
_SPECIES_UNIPROT_TYPE: dict[str, str] = {
    "drosophila": "FlyBase",
    "celegans":   "WormBase",
    "zebrafish":  "ZFIN",
}

NCBI_GENE_INFO_URLS: dict[str, str] = {
    "drosophila": (
        "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Invertebrates/"
        "Drosophila_melanogaster.gene_info.gz"
    ),
    "celegans": (
        "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Invertebrates/"
        "Caenorhabditis_elegans.gene_info.gz"
    ),
    "zebrafish": (
        "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Non-mammalian_vertebrates/"
        "Danio_rerio.gene_info.gz"
    ),
}

# Each species' organism-specific CURIE column. The parser writes the
# matching column when it finds the right dbXrefs prefix.
_SPECIES_ORG_COLUMN: dict[str, str] = {
    "drosophila": "flybase",
    "celegans":   "wormbase",
    "zebrafish":  "zfin",
}

# Prefixes in dbXrefs to extract organism-specific gene IDs. Keys are
# the species tag, values are the dbXrefs prefix (uppercase or
# Camel-cased as NCBI writes them, checked case-insensitively).
_SPECIES_ORG_DBXREF_PREFIXES: dict[str, tuple[str, ...]] = {
    "drosophila": ("FLYBASE",),
    "celegans":   ("WormBase",),
    "zebrafish":  ("ZFIN",),
}

# MGI cross-reference reports (v0.6.1). All three are header-less,
# tab-separated, positional. Schemas read from sampled rows on the
# upstream pages, re-verify if MGI ever changes their report layout.
MGI_ENSEMBL_URL = (
    "https://www.informatics.jax.org/downloads/reports/MRK_ENSEMBL.rpt"
)
MGI_ENTREZ_URL = (
    "https://www.informatics.jax.org/downloads/reports/MGI_EntrezGene.rpt"
)
MGI_UNIPROT_URL = (
    "https://www.informatics.jax.org/downloads/reports/MRK_SwissProt_TrEMBL.rpt"
)


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


def _download_url(url: str, target: Path, *, tag: str) -> None:
    """Generic resumable-ish file download to ``target``. ``tag`` is a short
    label used in log lines."""
    target.parent.mkdir(parents=True, exist_ok=True)
    print(f"[{tag}] downloading {url} ...", file=sys.stderr)
    r = requests.get(url, stream=True, timeout=180)
    r.raise_for_status()
    with target.open("wb") as fh:
        for chunk in r.iter_content(chunk_size=1024 * 1024):
            if chunk:
                fh.write(chunk)
    print(f"[{tag}] wrote {target} ({target.stat().st_size:,} bytes)", file=sys.stderr)


def _download_mgi(target: Path) -> None:
    """Download the MGI marker list (MRK_List2.rpt)."""
    _download_url(MGI_URL, target, tag="mgi")


def _load_uniprot_by_organism(
    path: Path, organism_type: str
) -> dict[str, str]:
    """Parse a UniProt per-organism idmapping.dat.gz file and return
    a ``{organism-specific id → UniProt accession}`` dict for the
    given ``organism_type`` (e.g. ``"FlyBase"`` / ``"WormBase"`` /
    ``"ZFIN"``).

    Each line of the file is ``<UniProtKB-AC>\\t<TYPE>\\t<ID>``;
    we filter to the requested TYPE. Multiple UniProt accessions per
    organism id are possible (paralog / fragment entries); the first
    one wins under load order, which is OK because UniProt sorts
    SwissProt-reviewed accessions before TrEMBL.
    """
    out: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            uniprot_ac, type_col, org_id = parts
            if type_col != organism_type:
                continue
            if org_id and uniprot_ac and org_id not in out:
                out[org_id] = uniprot_ac
    print(
        f"[uniprot] indexed {len(out):,} {organism_type}→UniProt rows",
        file=sys.stderr,
    )
    return out


def _load_ncbi_gene_info(
    path: Path, species: str
) -> list[dict[str, str]]:
    """Parse a gzip'd NCBI ``gene_info`` file for one organism and emit
    rows ready for ``build()`` step 4. Extracts:

    - ``symbol``: ``Symbol_from_nomenclature_authority`` if present
      (not ``-``), otherwise ``Symbol``.
    - ``ncbigene``: the GeneID column.
    - ``species``: the supplied tag (``drosophila`` / ``celegans`` /
      ``zebrafish``).
    - The organism-specific CURIE column (``flybase`` / ``wormbase``
      / ``zfin``) and ``ensembl``: extracted from ``dbXrefs``.
    - ``description``: from the description column.

    UniProt is not consistently surfaced in NCBI gene_info dbXrefs for
    these organisms; rows ship with ``uniprot=""`` in v0.6.2 and may be
    backfilled in a later release from UniProt's species-specific
    cross-reference files.
    """
    org_col = _SPECIES_ORG_COLUMN[species]
    org_prefixes = tuple(p.lower() for p in _SPECIES_ORG_DBXREF_PREFIXES[species])
    rows: list[dict[str, str]] = []
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 11:
                continue
            ncbigene = parts[1].strip()
            symbol_raw = parts[2].strip()
            dbxrefs_field = parts[5].strip()
            description = parts[8].strip()
            symbol_auth = parts[10].strip()
            symbol = symbol_auth if symbol_auth and symbol_auth != "-" else symbol_raw
            if not symbol or not ncbigene:
                continue

            org_id = ""
            ensembl = ""
            for xref in dbxrefs_field.split("|"):
                if not xref or ":" not in xref:
                    continue
                key, _, val = xref.partition(":")
                key_lc = key.strip().lower()
                if not org_id and key_lc in org_prefixes:
                    org_id = val.strip()
                elif not ensembl and key_lc in ("ensembl", "ensembl_metazoa"):
                    # v0.6.3: NCBI uses Ensembl_Metazoa for the invertebrate
                    # divisions (Drosophila / C. elegans); zebrafish stays on
                    # plain Ensembl. Accept both so the runtime sees a
                    # single ``ensembl`` column populated regardless of
                    # which division the organism comes from.
                    ensembl = val.strip()

            rows.append({
                "symbol": symbol,
                org_col: org_id,
                "ensembl": ensembl,
                "ncbigene": ncbigene,
                "description": description,
            })
    print(
        f"[{species}] indexed {len(rows):,} rows "
        f"({sum(1 for r in rows if r.get(org_col)):,} with {org_col} xref, "
        f"{sum(1 for r in rows if r.get('ensembl')):,} with ensembl xref)",
        file=sys.stderr,
    )
    return rows


def _load_mgi_xrefs(
    ensembl_path: Path | None,
    entrez_path: Path | None,
    uniprot_path: Path | None,
) -> dict[str, dict[str, str]]:
    """Build a ``{MGI accession → {ensembl, ncbigene, uniprot}}`` map by
    parsing the three MGI cross-reference reports.

    All three reports are header-less, tab-separated. Column positions
    are inlined here from upstream samples and noted at the URL
    constants; re-verify if MGI changes the layout. Empty / withdrawn
    rows are skipped.
    """
    out: dict[str, dict[str, str]] = {}

    def _ensure(mgi_id: str) -> dict[str, str]:
        if mgi_id not in out:
            out[mgi_id] = {"ensembl": "", "ncbigene": "", "uniprot": ""}
        return out[mgi_id]

    if ensembl_path is not None and ensembl_path.is_file():
        # MRK_ENSEMBL.rpt: col 0 = MGI id, col 5 = Ensembl Gene ID
        # (single id). Always one Ensembl gene per row.
        n = 0
        with ensembl_path.open("r", encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                mgi_id = parts[0].strip()
                ensembl = parts[5].strip()
                if mgi_id and ensembl:
                    _ensure(mgi_id)["ensembl"] = ensembl
                    n += 1
        print(f"[mgi] indexed {n:,} MGI→Ensembl rows", file=sys.stderr)

    if entrez_path is not None and entrez_path.is_file():
        # MGI_EntrezGene.rpt: col 0 = MGI id, col 2 = Status (O / W),
        # col 6 = Marker Type, col 8 = Entrez Gene ID. Filter to
        # Status="O" (official) AND Marker Type="Gene" (some lines
        # are withdrawn allele records with no Entrez id anyway).
        n = 0
        with entrez_path.open("r", encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                mgi_id = parts[0].strip()
                status = parts[2].strip()
                marker_type = parts[6].strip()
                entrez = parts[8].strip()
                if (
                    mgi_id
                    and entrez
                    and status == "O"
                    and marker_type == "Gene"
                ):
                    _ensure(mgi_id)["ncbigene"] = entrez
                    n += 1
        print(f"[mgi] indexed {n:,} MGI→Entrez rows", file=sys.stderr)

    if uniprot_path is not None and uniprot_path.is_file():
        # MRK_SwissProt_TrEMBL.rpt: col 0 = MGI id, col 6 = UniProt
        # accessions (space-separated). Take the first SwissProt-style
        # accession.
        n = 0
        with uniprot_path.open("r", encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 7:
                    continue
                mgi_id = parts[0].strip()
                uniprots = parts[6].strip().split()
                if mgi_id and uniprots:
                    _ensure(mgi_id)["uniprot"] = uniprots[0]
                    n += 1
        print(f"[mgi] indexed {n:,} MGI→UniProt rows", file=sys.stderr)

    print(
        f"[mgi] {len(out):,} MGI accessions carry at least one xref",
        file=sys.stderr,
    )
    return out


def _load_mgi(path: Path) -> list[dict[str, str]]:
    """Parse MGI MRK_List2.rpt and return mouse-gene rows.

    Filters to ``Marker Type == "Gene"`` to drop BAC/YAC ends, syntenic
    markers, RIKEN cDNAs without a stable gene status, etc. Returned
    rows carry the minimum fields needed by the bundled TSV
    (``symbol``, ``mgi``, ``name``); ``ensembl`` / ``ncbigene`` /
    ``uniprot`` cross-refs are left empty in v0.6.0, they live in
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

    * ``by_uniprot`` , UniProt accession → HGNC row
    * ``by_symbol``  , UPPER-cased HGNC symbol → HGNC row (case-insensitive
      fallback when no UniProt accession matched)
    * ``all_rows``   , every HGNC row, used to emit pure-human rows for genes
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
    g2v2_path: Path,
    hgnc_path: Path,
    mgi_path: Path,
    out_path: Path,
    *,
    mgi_xrefs: dict[str, dict[str, str]] | None = None,
    other_species_rows: dict[str, list[dict[str, str]]] | None = None,
    other_species_uniprot: dict[str, dict[str, str]] | None = None,
) -> int:
    """Read G2V2 + HGNC + MGI, join, emit the bundled TSV. Return row count.

    Output composition (in this load order, affects ``GeneIndex.lookup``
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

        # Step 1, Nematostella-anchored rows.
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
                "flybase": "",
                "wormbase": "",
                "zfin": "",
                "ensembl": (hgnc_row or {}).get("ensembl", ""),
                "ncbigene": (hgnc_row or {}).get("ncbigene", ""),
                "nemve1": jgi,
                "nvertx": nvertx,
                "uniprot": uniprot,
                "species": "nematostella",
                "description": description,
            })
            n_nem_total += 1

        # Step 2, pure-HGNC rows for human genes not yet covered.
        n_pure_hgnc = 0
        for row in hgnc_all_rows:
            if row["hgnc"] in used_hgnc_ids:
                continue
            w.writerow({
                "symbol": row["symbol"],
                "hgnc": row["hgnc"],
                "mgi": "",
                "flybase": "",
                "wormbase": "",
                "zfin": "",
                "ensembl": row["ensembl"],
                "ncbigene": row["ncbigene"],
                "nemve1": "",
                "nvertx": "",
                "uniprot": row["uniprot"],
                "species": "human",
                "description": row["name"],
            })
            n_pure_hgnc += 1

        # Step 3, MGI mouse rows (v0.6.0+).
        # v0.6.1 enriches each row with Ensembl / Entrez Gene /
        # UniProt cross-references when ``mgi_xrefs`` is provided.
        mgi_rows = _load_mgi(mgi_path)
        n_mgi = 0
        n_mgi_with_ensembl = 0
        n_mgi_with_entrez = 0
        n_mgi_with_uniprot = 0
        for row in mgi_rows:
            xrefs = (mgi_xrefs or {}).get(row["mgi"], {})
            ensembl = xrefs.get("ensembl", "")
            ncbigene = xrefs.get("ncbigene", "")
            uniprot = xrefs.get("uniprot", "")
            if ensembl:
                n_mgi_with_ensembl += 1
            if ncbigene:
                n_mgi_with_entrez += 1
            if uniprot:
                n_mgi_with_uniprot += 1
            w.writerow({
                "symbol": row["symbol"],
                "hgnc": "",
                "mgi": row["mgi"],
                "flybase": "",
                "wormbase": "",
                "zfin": "",
                "ensembl": ensembl,
                "ncbigene": ncbigene,
                "nemve1": "",
                "nvertx": "",
                "uniprot": uniprot,
                "species": "mouse",
                "description": row["name"],
            })
            n_mgi += 1

        # Step 4, other model organisms (v0.6.2): Drosophila / C. elegans
        # / zebrafish from NCBI gene_info files.
        n_other = {sp: 0 for sp in NCBI_GENE_INFO_URLS}
        n_other_uniprot = {sp: 0 for sp in NCBI_GENE_INFO_URLS}
        for species, rows in (other_species_rows or {}).items():
            org_col = _SPECIES_ORG_COLUMN.get(species)
            if org_col is None:
                continue
            uniprot_map = (other_species_uniprot or {}).get(species, {})
            for r in rows:
                org_id = r.get(org_col, "") if org_col in r else ""
                uniprot = uniprot_map.get(org_id, "") if org_id else ""
                if uniprot:
                    n_other_uniprot[species] += 1
                w.writerow({
                    "symbol": r["symbol"],
                    "hgnc": "",
                    "mgi": "",
                    "flybase": org_id if org_col == "flybase" else "",
                    "wormbase": org_id if org_col == "wormbase" else "",
                    "zfin": org_id if org_col == "zfin" else "",
                    "ensembl": r.get("ensembl", ""),
                    "ncbigene": r.get("ncbigene", ""),
                    "nemve1": "",
                    "nvertx": "",
                    "uniprot": uniprot,
                    "species": species,
                    "description": r.get("description", ""),
                })
                n_other[species] += 1

    n_total = n_nem_total + n_pure_hgnc + n_mgi + sum(n_other.values())
    print(
        f"[out] wrote {out_path}: {n_total:,} rows total\n"
        f"       Nematostella-anchored : {n_nem_total:,}\n"
        f"         via UniProt join    : {n_via_uniprot:,}\n"
        f"         via symbol fallback : {n_via_symbol:,}\n"
        f"         Nematostella-only   : {n_nem_total - n_via_uniprot - n_via_symbol:,}\n"
        f"       Pure-HGNC (no Nemve1) : {n_pure_hgnc:,}\n"
        f"       MGI (mouse)           : {n_mgi:,}\n"
        f"         with Ensembl xref   : {n_mgi_with_ensembl:,}\n"
        f"         with Entrez xref    : {n_mgi_with_entrez:,}\n"
        f"         with UniProt xref   : {n_mgi_with_uniprot:,}\n"
        f"       Drosophila            : {n_other.get('drosophila', 0):,}"
        f"  (+{n_other_uniprot.get('drosophila', 0):,} UniProt xref)\n"
        f"       C. elegans            : {n_other.get('celegans', 0):,}"
        f"  (+{n_other_uniprot.get('celegans', 0):,} UniProt xref)\n"
        f"       Zebrafish             : {n_other.get('zebrafish', 0):,}"
        f"  (+{n_other_uniprot.get('zebrafish', 0):,} UniProt xref)",
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
    # v0.6.1, MGI cross-reference reports (Ensembl / Entrez / UniProt).
    parser.add_argument(
        "--mgi-ensembl-cache", type=Path,
        default=cache_dir / "MRK_ENSEMBL.rpt",
        help="Where to cache MGI's MRK_ENSEMBL.rpt (mouse Ensembl xrefs).",
    )
    parser.add_argument(
        "--mgi-entrez-cache", type=Path,
        default=cache_dir / "MGI_EntrezGene.rpt",
        help="Where to cache MGI's MGI_EntrezGene.rpt (mouse NCBI Gene xrefs).",
    )
    parser.add_argument(
        "--mgi-uniprot-cache", type=Path,
        default=cache_dir / "MRK_SwissProt_TrEMBL.rpt",
        help="Where to cache MGI's MRK_SwissProt_TrEMBL.rpt (mouse UniProt xrefs).",
    )
    parser.add_argument(
        "--no-mgi-xrefs", action="store_true",
        help="Skip the v0.6.1 mouse Ensembl/Entrez/UniProt xref join "
             "(emits mouse rows with empty xref columns).",
    )
    # v0.6.2, other model organisms (Drosophila / C. elegans / Zebrafish).
    parser.add_argument(
        "--drosophila-cache", type=Path,
        default=cache_dir / "Drosophila_melanogaster.gene_info.gz",
        help="Where to cache the NCBI Drosophila gene_info.gz.",
    )
    parser.add_argument(
        "--celegans-cache", type=Path,
        default=cache_dir / "Caenorhabditis_elegans.gene_info.gz",
        help="Where to cache the NCBI C. elegans gene_info.gz.",
    )
    parser.add_argument(
        "--zebrafish-cache", type=Path,
        default=cache_dir / "Danio_rerio.gene_info.gz",
        help="Where to cache the NCBI zebrafish gene_info.gz.",
    )
    parser.add_argument(
        "--no-other-species", action="store_true",
        help="Skip the v0.6.2 Drosophila / C. elegans / zebrafish step "
             "(emits human + mouse + Nematostella only).",
    )
    # v0.6.3, UniProt per-organism idmapping for the three new species.
    parser.add_argument(
        "--drosophila-uniprot-cache", type=Path,
        default=cache_dir / "DROME_7227_idmapping.dat.gz",
        help="Where to cache the UniProt Drosophila idmapping file.",
    )
    parser.add_argument(
        "--celegans-uniprot-cache", type=Path,
        default=cache_dir / "CAEEL_6239_idmapping.dat.gz",
        help="Where to cache the UniProt C. elegans idmapping file.",
    )
    parser.add_argument(
        "--zebrafish-uniprot-cache", type=Path,
        default=cache_dir / "DANRE_7955_idmapping.dat.gz",
        help="Where to cache the UniProt zebrafish idmapping file.",
    )
    parser.add_argument(
        "--no-other-species-uniprot", action="store_true",
        help="Skip the v0.6.3 UniProt cross-ref pass for "
             "Drosophila / C. elegans / zebrafish.",
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

    mgi_xrefs: dict[str, dict[str, str]] | None = None
    if not args.no_mgi_xrefs:
        for url, cache in (
            (MGI_ENSEMBL_URL, args.mgi_ensembl_cache),
            (MGI_ENTREZ_URL, args.mgi_entrez_cache),
            (MGI_UNIPROT_URL, args.mgi_uniprot_cache),
        ):
            if not cache.is_file():
                _download_url(url, cache, tag="mgi-xref")
        mgi_xrefs = _load_mgi_xrefs(
            args.mgi_ensembl_cache,
            args.mgi_entrez_cache,
            args.mgi_uniprot_cache,
        )

    other_species_rows: dict[str, list[dict[str, str]]] | None = None
    other_species_uniprot: dict[str, dict[str, str]] | None = None
    if not args.no_other_species:
        other_species_rows = {}
        species_caches = {
            "drosophila": args.drosophila_cache,
            "celegans":   args.celegans_cache,
            "zebrafish":  args.zebrafish_cache,
        }
        for species, cache in species_caches.items():
            url = NCBI_GENE_INFO_URLS[species]
            if not cache.is_file():
                _download_url(url, cache, tag=species)
            other_species_rows[species] = _load_ncbi_gene_info(cache, species)

        if not args.no_other_species_uniprot:
            other_species_uniprot = {}
            uniprot_caches = {
                "drosophila": args.drosophila_uniprot_cache,
                "celegans":   args.celegans_uniprot_cache,
                "zebrafish":  args.zebrafish_uniprot_cache,
            }
            for species, cache in uniprot_caches.items():
                url = UNIPROT_BY_ORGANISM_URLS[species]
                if not cache.is_file():
                    _download_url(url, cache, tag=f"{species}-uniprot")
                other_species_uniprot[species] = _load_uniprot_by_organism(
                    cache, _SPECIES_UNIPROT_TYPE[species]
                )

    build(
        args.g2v2_xlsx, hgnc_path, mgi_path, args.out,
        mgi_xrefs=mgi_xrefs,
        other_species_rows=other_species_rows,
        other_species_uniprot=other_species_uniprot,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
