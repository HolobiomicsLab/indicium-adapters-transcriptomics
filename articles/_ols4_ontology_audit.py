"""One-shot audit of the OLS4 ontology catalogue against the adapter's
existing coverage. Outputs a markdown report identifying which ontologies
are already covered, which are high-value adds for the ASB benchmark
phase, and which are niche or out-of-scope.

Run::

    py -3.13 articles/_ols4_ontology_audit.py > articles/_ols4_audit_report.md
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

CATALOG = Path.home() / ".cache" / "ols4_ontology_catalog.json"

# Ontology prefixes already in TRANSCRIPTOMICS_PREFIXES (or under active
# use via concept/CL/UBERON indexes).
ALREADY_COVERED = {
    "efo", "obi", "stato", "ncit", "cl", "uberon", "go", "so", "mondo",
    "pato", "edam", "ino", "sio", "ro", "iao", "ncbitaxon",
}

# Keywords that flag an ontology as high-value for transcriptomics
# benchmark exhaustiveness. Each tuple is (bucket_label, keywords).
# Iterated 2026-05-28 after a first pass left too many transcripto-
# adjacent ontologies in "Other".
HIGH_VALUE_BUCKETS: list[tuple[str, tuple[str, ...]]] = [
    # Disease / phenotype — clinical transcriptomics
    ("disease", ("disease", "disorder", "mondo", "doid", "omim",
                 "syndrome", "neoplasm", "cancer", "tumor", "tumour",
                 "pathology")),
    ("phenotype", ("phenotype", "phenotyp",  # phenotypic too
                   "trait", " hp ", "mammalian phenotype",
                   "zebrafish phenotype", "fly phenotype",
                   "worm phenotype", "behavior", "behaviour")),
    # Chemistry / drugs / lipids — perturbation experiments
    ("chemicals_drugs", ("chebi", "chemical entity", "drug",
                         "compound", "pharmacolog", "biomolecule",
                         "lipid", "lipidmaps", "metabolite")),
    # Protein / pathway / interaction — transcripto adjacent
    ("protein_pathway", ("protein ontology", "pathway", " pw ",
                         "reaction", "rxno", "kegg", "enzyme",
                         "molecular interaction", "interaction",
                         "cell cycle", "mhc", "complement system",
                         "signaling")),
    # Environment / sample collection — holobiont / ecology
    ("environment", ("environment", "envo", "habitat", "biome",
                     "ecosystem", "ecological")),
    # Sample / specimen / biobank metadata
    ("sample_biobank", ("biobank", "specimen", "sample collection",
                        "tissue bank", "cohort", "demographic",
                        "breed")),
    # Cell / anatomy adjacent (catches GeXO/ReXO/ReTO orgs and brain
    # atlases / model-organism anatomies / cell-atlas overlays)
    ("cell_anatomy_extra", ("cell line", "ccle", "cell ontology",
                            "anatomical", "anatomy", "tissue ontology",
                            "developmental stage", "life cycle",
                            "lifestage", "brain atlas",
                            "cell atlas", "embryonic", "developmental",
                            "gross anatomy")),
    # Methods / techniques (catches mass spec / proteomics / sequencing
    # methodology vocabs)
    ("methods", ("microarray", "sequencing", "ngs",
                 "single cell", "single-cell", "transcriptom",
                 "rna-seq", "rnaseq", "spatial transcriptom",
                 "experimental factor", "assay", "gene expression",
                 "regulation of transcription", "mass spectrometr",
                 "proteomics", "mirna target", "microrna",
                 "medical action")),
    # Statistics / measurement
    ("stat_measure", ("statistic", "measurement", "unit ontology",
                      "uo ", "outcome", "endpoint")),
    # Genomic features / variation
    ("genomic_features", ("sequence feature", "genomic feature",
                          "variation", " variant", "mutation",
                          "feature annotation", "splicing",
                          "transposable element", "conjugative element")),
    # Modifications / PTMs
    ("modifications", ("modification", "ptm", "post-translational",
                       "epigenetic", "methylation", "acetylation",
                       "glycan", "glycomics")),
    # Food / nutrition (diet studies)
    ("food_nutrition", ("food", "nutrition", "diet", "foodon")),
    # Plants (if multi-omics expands)
    ("plant", ("plant", "po ontology", "gramene", "agronomy",
               "agronomic", "agro", "flora")),
    # Upper / shared semantics
    ("upper_semantic", ("basic formal ontology", "bfo", "common core",
                        "relation", "information artifact",
                        "information resource")),
]

# Manually-curated short list — these are the ontologies most likely
# to fire on a diverse transcriptomics corpus, ordered by expected
# value. Each entry: (prefix, why-it-matters one-liner).
TOP_PICKS: tuple[tuple[str, str], ...] = (
    ("HP",       "Phenotype ontology — clinical transcripto papers (HPO is THE phenotype vocab for human disease)"),
    ("CHEBI",    "Small molecules — perturbation experiments (drug treatments, ligands, metabolites)"),
    ("DOID",     "Disease ontology — clinical / cancer transcriptomics (complements MONDO, more drug-perturbation-friendly)"),
    ("HPO",      "Same as HP if listed separately"),
    ("MP",       "Mammalian Phenotype — mouse studies (Tabula Muris, transgenic mice)"),
    ("OBA",      "Ontology of Biological Attributes — cross-kingdom traits, multi-species KGs"),
    ("PR",       "Protein Ontology — protein-level claims, isoforms, complexes"),
    ("ChEBI",    "Same as CHEBI if cased differently"),
    ("GeXO",     "Gene Expression Ontology — direct overlap with our domain, possibly redundant with our concept index"),
    ("ReXO",     "Regulation of Gene Expression Ontology — regulation-process vocabulary beyond INO"),
    ("ReTO",     "Regulation of Transcription Ontology — transcriptional control vocab"),
    ("EC",       "Enzyme Commission codes — enzymatic activity, multi-omics bridge"),
    ("MS",       "Mass spectrometry — proteomics/metabolomics overlap (cf. metabolomics adapter)"),
    ("OMIT",     "miRNA Target — for miRNA-focused transcriptomics studies"),
    ("HCAO",     "Human Cell Atlas Ontology — cell-type vocabulary harmonised across atlases"),
    ("EMAPA",    "Mouse Developmental Anatomy — Tabula Muris context, embryonic transcriptomics"),
    ("FBbt",     "Drosophila gross anatomy — paired with our v0.6.2 FlyBase gene index"),
    ("WBbt",     "C. elegans gross anatomy — paired with our v0.6.2 WormBase gene index"),
    ("ZFA",      "Zebrafish anatomy — paired with our v0.6.2 ZFIN gene index"),
    ("ZP",       "Zebrafish Phenotype — Zebrafish disease modelling"),
    ("XPO",      "Xenopus Phenotype — Xenopus regeneration studies"),
    ("ENVO",     "Environment Ontology — holobiont context, sample collection metadata"),
    ("NBO",      "Neuro Behavior Ontology — neural transcriptomics, behavioural phenotypes"),
    ("RXNO",     "Reaction Ontology — biochemical reaction context"),
    ("MOD",      "Protein modification (PTMs) — phospho/methyl/acetyl-proteomics overlap"),
    ("UO",       "Units of measurement — quantitative claim normalisation"),
    ("BTO",      "BRENDA Tissue / Cell Type — alternative to CL/UBERON, popular in enzyme literature"),
    ("OGSF",     "Genetic Susceptibility Factor — GWAS context"),
    ("CMPO",     "Cellular Microscopy Phenotype — imaging-derived single-cell phenotypes"),
    ("MEDDRA",   "Medical Dictionary for Regulatory Activities — drug-induced expression studies"),
)

# IDs we explicitly mark as "niche" — included for completeness but
# unlikely to fire on a transcriptomics paper unless the paper is in
# that exact sub-area.
NICHE_TAGS = {"covid", "rare", "specific", "regional"}


def categorise(o: dict) -> str:
    """Return a bucket tag: ``covered`` / ``high:<label>`` / ``niche`` /
    ``other``."""
    oid = (o.get("id") or "").lower()
    pref = (o.get("preferredPrefix") or "").lower()
    if oid in ALREADY_COVERED or pref in ALREADY_COVERED:
        return "covered"

    title = (o.get("title") or "").lower()
    desc = (o.get("description") or "").lower()
    blob = f"{oid} {pref} {title} {desc}"

    for bucket, keywords in HIGH_VALUE_BUCKETS:
        for kw in keywords:
            if kw in blob:
                return f"high:{bucket}"

    return "other"


def main() -> int:
    if not CATALOG.is_file():
        print(f"Catalog not found: {CATALOG}", file=sys.stderr)
        return 1
    cat = json.loads(CATALOG.read_text(encoding="utf-8"))

    buckets: dict[str, list[dict]] = {}
    for o in cat:
        b = categorise(o)
        buckets.setdefault(b, []).append(o)

    # Sort each bucket by term count descending (bigger ontology = more
    # likely to surface vocab).
    for b in buckets:
        buckets[b].sort(key=lambda o: -(o.get("numberOfTerms") or 0))

    # Render report
    def _row(o: dict) -> str:
        nt = o.get("numberOfTerms") or 0
        pref = o.get("preferredPrefix") or "—"
        title = (o.get("title") or "").strip() or "(untitled)"
        desc = (o.get("description") or "").strip()
        desc = desc.replace("\n", " ").replace("|", "\\|")
        if len(desc) > 200:
            desc = desc[:197] + "..."
        return f"| `{pref}` | {title} | {nt:,} | {desc} |"

    print("# OLS4 ontology audit for the ASB transcriptomics benchmark")
    print()
    print(f"Audit date: 2026-05-28 · Catalogue size: **{len(cat)} ontologies**")
    print()
    print("Two sections:")
    print()
    print("1. **Top picks** — curated short list of the ontologies most")
    print("   likely to fire on a diverse transcriptomics corpus, with")
    print("   one-liner justifications. **Start here.**")
    print("2. **Bucket dump** — heuristic keyword classification of every")
    print("   ontology in the OLS4 catalogue. Useful as raw material for")
    print("   building per-axis harvesters in future versions.")
    print()

    # --- TOP PICKS section ----------------------------------------------
    by_pref: dict[str, dict] = {}
    by_id: dict[str, dict] = {}
    for o in cat:
        pref = (o.get("preferredPrefix") or "").strip()
        oid = (o.get("id") or "").strip()
        if pref:
            by_pref[pref.lower()] = o
        if oid:
            by_id[oid.lower()] = o

    print("## 🎯 Top picks (curated)")
    print()
    print("| Prefix | Title | # terms | Why it matters |")
    print("|---|---|---:|---|")
    seen_picks: set[str] = set()
    for pref, why in TOP_PICKS:
        key = pref.lower()
        if key in seen_picks:
            continue
        seen_picks.add(key)
        o = by_pref.get(key) or by_id.get(key)
        if o is None:
            # Not found in OLS4 — note it but don't drop it.
            print(f"| `{pref}` | (not indexed in OLS4) | — | {why} |")
            continue
        title = (o.get("title") or "").strip() or "(untitled)"
        nt = o.get("numberOfTerms") or 0
        print(f"| `{o.get('preferredPrefix') or pref}` | {title} | {nt:,} | {why} |")
    print()
    print("---")
    print()
    print("## Bucket dump (heuristic keyword classification)")
    print()
    print("Each ontology is classified by keyword matching on its")
    print("`title + description + id`. The classification is rough — use")
    print("the **Top picks** above for short-term decisions, and dip into")
    print("the buckets when adding a new axis (e.g. cell-line or imaging")
    print("phenotype) to find candidates.")
    print()

    bucket_order = [
        ("covered", "Already covered by the adapter"),
        ("high:disease", "🎯 HIGH — Disease ontologies (clinical transcripto / cancer / disorders)"),
        ("high:phenotype", "🎯 HIGH — Phenotype ontologies (HP, MP, ZP, FBcv, …)"),
        ("high:chemicals_drugs", "🎯 HIGH — Chemicals / drugs (ChEBI, perturbation experiments)"),
        ("high:protein_pathway", "🎯 HIGH — Protein / pathway / interaction"),
        ("high:cell_anatomy_extra", "🎯 HIGH — Cell / anatomy extras (developmental stage, cell line, etc.)"),
        ("high:methods", "🎯 HIGH — Methods / techniques"),
        ("high:genomic_features", "Medium — Genomic features / variation"),
        ("high:modifications", "Medium — Modifications / PTMs"),
        ("high:stat_measure", "Medium — Statistics / measurement"),
        ("high:environment", "Medium — Environment / ecology (holobiont context)"),
        ("high:sample_biobank", "Medium — Sample / biobank metadata"),
        ("high:upper_semantic", "Medium — Upper-level / shared semantics"),
        ("high:food_nutrition", "Low — Food / nutrition (diet studies only)"),
        ("high:plant", "Low — Plant ontologies (only if multi-omics expands)"),
        ("other", "Other / Uncategorised"),
    ]

    for key, label in bucket_order:
        rows = buckets.get(key) or []
        if not rows:
            continue
        print(f"## {label} ({len(rows)} ontologies)")
        print()
        print("| Prefix | Title | # terms | Description |")
        print("|---|---|---:|---|")
        for o in rows:
            print(_row(o))
        print()

    # Summary line at the bottom
    total = sum(len(v) for v in buckets.values())
    cov_n = len(buckets.get("covered") or [])
    high_n = sum(len(v) for k, v in buckets.items() if k.startswith("high:") and "disease phenotype chemicals_drugs protein_pathway cell_anatomy_extra methods".find(k.split(":", 1)[1]) >= 0)
    print()
    print(f"**Total: {total}** · already covered: **{cov_n}** · "
          f"high-priority candidates (top 6 buckets): **{high_n}**")

    return 0


if __name__ == "__main__":
    sys.exit(main())
