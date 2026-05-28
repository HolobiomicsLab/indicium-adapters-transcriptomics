"""Maintainer script, harvest Cell Ontology (CL) CURIEs for cell-type terms.

Companion to ``build_corpus_concept_candidates.py`` and
``build_biotools_index.py``. The OLS4 corpus harvester only resolves
methodological vocab; the bio.tools harvester only resolves software.
Neither catches **cell types**, yet single-cell claims like
``marker_of_cell_type: WNT4 cnidocyte`` need a CURIE for the object slot.

This script ships a hand-curated list of candidate cell-type names
(drawn from the three cell-atlas papers in the 6-paper mini-corpus,
Sebé-Pedrós 2018 *Nematostella* atlas, Tabula Muris 2018, Tabula
Sapiens 2022, plus the canonical mammalian immune / nervous /
epithelial cell types that show up in any transcriptomics study). Each
candidate is resolved against the **Cell Ontology** via OLS4 search
with the same cascade match modes as ``build_corpus_concept_candidates``
(exact label → exact synonym → prefix); ``substring`` is omitted
because cell-type labels in CL are typically already minimal and
substring matching would invite false positives.

Hand-curation of CURIEs is still avoided: only the *candidate name
list* is hand-picked. The CURIE that lands in the bundled TSV is
whatever OLS4 returns and re-verifies. If a candidate name resolves to
a non-canonical entry (cf. macrophage → plasmatocyte before the
label-first cascade refinement), the OLS4 cascade picks the canonical
label-match instead.

Pleasant surprise found while writing this: CL covers cnidarian cell
types like ``cnidocyte`` (CL:0011029), unlike UBERON which had no
entry for cnidarian-specific anatomy. So *Nematostella* single-cell
claims get real CL CURIEs, not deferred-to-future-release stubs.

Usage (maintainer-only)::

    python -m indicium_adapters_transcriptomics.data.build_cl_concept_index

Re-run when the cell-type candidate list grows or CL adds new entries.
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path
from typing import NamedTuple

from indicium_adapters_transcriptomics.data.build_corpus_concept_candidates import (
    _search_ontology,
    _try_match_in_docs,
    verify_iri,
)

OUT_COLUMNS = ("term", "curie", "label", "description", "ontology", "source")

# Cascade, same as build_corpus_concept_candidates for stat_test /
# protocol axes. We skip ``substring`` for cell types because the labels
# are already minimal, a substring match would invite false positives
# like ``cell`` matching ``B cell`` matching ``T cell`` matching anything.
_MATCH_MODES = ("exact", "prefix")


class Candidate(NamedTuple):
    name: str          # corpus-style mention (e.g. "T cell", "cnidocyte")
    source_hint: str   # short note on where the term comes from


# Curated candidate list, organised by source paper. Names are written
# as a human would type them, singular, lowercase except for proper
# nouns. The OLS4 cascade handles label/synonym matching with case-
# insensitive comparison.
CANDIDATES: tuple[Candidate, ...] = (
    # --- Nematostella vectensis (Sebé-Pedrós 2018 cell atlas) ---------------
    # cnidocyte is the iconic cnidarian cell type; nematocyte is a
    # synonym, both should resolve to CL:0011029.
    Candidate("cnidocyte",                  "Sebé-Pedrós 2018 (Nematostella atlas)"),
    Candidate("nematocyte",                 "Sebé-Pedrós 2018 (Nematostella atlas)"),
    Candidate("epidermal cell",             "Sebé-Pedrós 2018 + general"),
    Candidate("sensory neuron",             "Sebé-Pedrós 2018 + general"),
    Candidate("neuron",                     "general nervous-system"),
    Candidate("muscle cell",                "Sebé-Pedrós 2018 + general"),
    Candidate("gland cell",                 "Sebé-Pedrós 2018 (Nematostella secretory)"),
    Candidate("secretory cell",             "Sebé-Pedrós 2018"),
    Candidate("germ cell",                  "Sebé-Pedrós 2018"),
    Candidate("mucus secreting cell",       "Sebé-Pedrós 2018"),

    # --- Mammalian immune system (Tabula Muris / Tabula Sapiens) -----------
    Candidate("T cell",                     "Tabula Muris/Sapiens immune"),
    Candidate("B cell",                     "Tabula Muris/Sapiens immune"),
    Candidate("natural killer cell",        "Tabula Muris/Sapiens immune"),
    Candidate("dendritic cell",             "Tabula Muris/Sapiens immune"),
    Candidate("macrophage",                 "Tabula Muris/Sapiens immune"),
    Candidate("neutrophil",                 "Tabula Muris/Sapiens immune"),
    Candidate("monocyte",                   "Tabula Muris/Sapiens immune"),
    Candidate("eosinophil",                 "Tabula Muris/Sapiens immune"),
    Candidate("mast cell",                  "Tabula Muris/Sapiens immune"),
    Candidate("plasma cell",                "Tabula Muris/Sapiens immune"),
    Candidate("erythrocyte",                "Tabula Muris/Sapiens blood"),

    # --- Mammalian nervous system (Tabula Muris / Sapiens) -----------------
    Candidate("astrocyte",                  "Tabula Muris/Sapiens CNS"),
    Candidate("oligodendrocyte",            "Tabula Muris/Sapiens CNS"),
    Candidate("microglial cell",            "Tabula Muris/Sapiens CNS"),

    # --- Mammalian epithelial / connective (Tabula Muris/Sapiens) ----------
    Candidate("fibroblast",                 "Tabula Muris/Sapiens stromal"),
    Candidate("endothelial cell",           "Tabula Muris/Sapiens vasculature"),
    Candidate("hepatocyte",                 "Tabula Muris/Sapiens liver"),
    Candidate("cardiomyocyte",              "Tabula Muris/Sapiens cardiac"),
    Candidate("smooth muscle cell",         "Tabula Muris/Sapiens"),
    Candidate("skeletal muscle cell",       "Tabula Muris/Sapiens"),
    Candidate("enterocyte",                 "Tabula Muris/Sapiens GI"),
    Candidate("goblet cell",                "Tabula Muris/Sapiens GI"),
    Candidate("Paneth cell",                "Tabula Muris/Sapiens GI"),
    Candidate("adipocyte",                  "Tabula Muris/Sapiens"),
    Candidate("chondrocyte",                "Tabula Muris/Sapiens"),
    Candidate("osteocyte",                  "Tabula Muris/Sapiens"),
    Candidate("stem cell",                  "Tabula Muris/Sapiens stemness"),
)


def resolve_cell_type(cand: Candidate) -> dict | None:
    """Apply the exact → prefix cascade against CL for one candidate name.
    Returns ``{curie, label, description, ontology, mode}`` or ``None``."""
    for mode in _MATCH_MODES:
        docs = _search_ontology(cand.name, "cl", exact=(mode == "exact"))
        time.sleep(0.3)
        if not docs:
            continue
        d = _try_match_in_docs(docs, cand.name, mode)
        if d is None:
            continue
        iri = d.get("iri")
        obo_id = d.get("obo_id") or d.get("short_form")
        if not iri or not obo_id:
            continue
        info = verify_iri(iri, "cl")
        time.sleep(0.3)
        if info is None:
            continue
        return {
            "curie": obo_id,
            "label": info["label"] or (d.get("label") or ""),
            "description": info["description"],
            "ontology": "CL",
            "mode": mode,
        }
    return None


def build(out_path: Path, unresolved_path: Path) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_resolved = 0
    unresolved: list[Candidate] = []
    seen_curies: set[str] = set()

    with out_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLUMNS, delimiter="\t")
        w.writeheader()
        for cand in CANDIDATES:
            print(f"[ols4-cl] resolving {cand.name!r} ...", file=sys.stderr)
            info = resolve_cell_type(cand)
            if info is None:
                print("  ✗ no match in CL", file=sys.stderr)
                unresolved.append(cand)
                continue
            # Deduplicate on CURIE, cnidocyte + nematocyte both resolve to
            # CL:0011029. Keep the first (label-exact will land first).
            if info["curie"] in seen_curies:
                print(
                    f"  ⇣ duplicate CURIE {info['curie']} (already written for "
                    f"another candidate), skipping the row.",
                    file=sys.stderr,
                )
                continue
            seen_curies.add(info["curie"])
            print(
                f"  ✓ [{info['mode']}] {info['curie']} {info['label']!r}",
                file=sys.stderr,
            )
            w.writerow({
                "term": cand.name.lower(),
                "curie": info["curie"],
                "label": info["label"],
                "description": info["description"],
                "ontology": info["ontology"],
                "source": f"corpus extraction (cell_type, {info['mode']}-match)",
            })
            n_resolved += 1

    with unresolved_path.open("w", encoding="utf-8") as fh:
        fh.write(
            "# Cell-type candidates with no exact/prefix match in CL.\n"
            "# Each line: <term>\\t<source hint>.\n"
            "# Candidates for hand-mapping, upstream CL contribution,\n"
            "# or alternative ontology (e.g. BTO for cell lines).\n\n"
        )
        for cand in unresolved:
            fh.write(f"{cand.name}\t{cand.source_hint}\n")

    print(
        f"\n[out] wrote {out_path}: {n_resolved} resolved, "
        f"{len(unresolved)} unresolved (of {len(CANDIDATES)} candidates).",
        file=sys.stderr,
    )
    return n_resolved


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--out", type=Path,
        default=Path(__file__).parent / "cl_concept_index.tsv",
        help="Output TSV path.",
    )
    parser.add_argument(
        "--unresolved", type=Path,
        default=Path(__file__).parent / "cl_concept_unresolved.txt",
        help="Output log for cell-type candidates that didn't resolve.",
    )
    args = parser.parse_args(argv)
    build(args.out, args.unresolved)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
