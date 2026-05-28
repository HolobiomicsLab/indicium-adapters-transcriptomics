"""Maintainer script — harvest the v0.8.0 high-priority ontology indexes.

Following the 2026-05-28 OLS4 audit (see ``articles/_ols4_audit_report.md``),
seven ontologies were identified as the highest-yield additions for the
ASB transcriptomics benchmark phase:

- **CHEBI** — chemicals / drugs / ligands (perturbation experiments).
- **HP** — Human Phenotype Ontology (clinical transcriptomics).
- **MP** — Mammalian Phenotype Ontology (mouse studies).
- **DOID** — Human Disease Ontology (cancer / disorder vocab).
- **PR** — Protein Ontology (complex / isoform-level claims).
- **OBA** — Ontology of Biological Attributes (cross-species traits).
- **ENVO** — Environment Ontology (holobiont / ecology context).

Same cascade pattern as ``build_cl_concept_index.py``: candidate names
→ OLS4 exact-label / synonym / prefix search → re-verify by IRI lookup.
One bundled TSV per ontology so each can be re-harvested independently
when its candidate list grows.

Hand-curation of CURIEs is still avoided. The only curation here is
the candidate-NAME list per ontology — the OLS4 cascade picks the
matching CURIE, or the term lands in the per-ontology unresolved log
for the next iteration.

Usage (maintainer-only)::

    python -m indicium_adapters_transcriptomics.data.build_extra_ontology_indexes
    # or per-ontology:
    python -m indicium_adapters_transcriptomics.data.build_extra_ontology_indexes --only chebi
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
from dataclasses import dataclass
from pathlib import Path

from indicium_adapters_transcriptomics.data.build_corpus_concept_candidates import (
    _search_ontology,
    _try_match_in_docs,
    verify_iri,
)

OUT_COLUMNS = ("term", "curie", "label", "description", "ontology", "source")


@dataclass(frozen=True)
class OntologyConfig:
    """One harvest target: ontology id + axis label + candidate names."""
    ols_id: str             # OLS4 ontology key — lowercase: "chebi", "hp"
    ontology_label: str     # Canonical prefix written to TSV (e.g. "CHEBI")
    axis: str               # axis tag for the ``source`` column
    out_filename: str       # bundled TSV filename
    candidates: tuple[str, ...]


# --- Candidate lists -------------------------------------------------------
#
# Each list is small on purpose — ~20-40 terms per ontology covering the
# vocabulary most likely to fire on a diverse transcriptomics corpus.
# Expand the list (re-run) when benchmark misses surface new terms.


CHEBI_CANDIDATES: tuple[str, ...] = (
    # Perturbations / pathway inhibitors found in transcripto literature
    "U0126", "MEK inhibitor", "rapamycin", "wortmannin", "LY294002",
    # Hormones / signalling molecules
    "estradiol", "testosterone", "progesterone", "dexamethasone",
    "hydrocortisone", "retinoic acid", "all-trans-retinoic acid",
    # Cancer chemotherapy (common in clinical transcripto)
    "cisplatin", "paclitaxel", "doxorubicin", "5-fluorouracil",
    "tamoxifen", "imatinib",
    # Common metabolites / cofactors
    "ATP", "GTP", "NAD", "NADH", "glucose", "glutamine", "lactic acid",
    # Lipids / fatty acids
    "fatty acid", "cholesterol", "prostaglandin E2",
    # Neurotransmitters
    "dopamine", "serotonin", "GABA", "glutamic acid",
    # Toxicants / inducers
    "lipopolysaccharide", "caffeine", "ethanol", "nicotine",
    # Vitamins
    "vitamin D", "vitamin A", "ascorbic acid",
)


HP_CANDIDATES: tuple[str, ...] = (
    # Developmental / morphological
    "intellectual disability", "microcephaly", "macrocephaly",
    "short stature", "hypotonia", "polydactyly", "syndactyly",
    "developmental delay", "global developmental delay",
    # Neurological
    "seizure", "ataxia", "dystonia", "tremor", "neurodegeneration",
    # Cardiac / vascular
    "cardiomyopathy", "congenital heart defect", "hypertension",
    # Hepatic / renal
    "hepatomegaly", "splenomegaly", "elevated alanine aminotransferase",
    "renal failure",
    # Metabolic
    "hyperglycemia", "hypoglycemia", "obesity",
    # Hematologic
    "anemia", "thrombocytopenia", "leukopenia",
    # Ophthalmologic
    "macular degeneration", "retinitis pigmentosa",
    # General
    "weight loss", "fatigue", "fever",
)


MP_CANDIDATES: tuple[str, ...] = (
    # Body-level abnormalities
    "abnormal body weight", "decreased body weight",
    "increased body weight", "premature death",
    # Organ morphology
    "abnormal heart morphology", "abnormal liver morphology",
    "abnormal muscle morphology", "abnormal kidney morphology",
    "abnormal brain morphology",
    # Physiology / behavior
    "abnormal immune system physiology",
    "abnormal nervous system physiology",
    "abnormal locomotor activity", "abnormal feeding behavior",
    # Reproductive
    "decreased fertility", "infertility",
    # Cell-level
    "abnormal cell proliferation", "abnormal neuron differentiation",
    "abnormal T cell physiology", "abnormal B cell physiology",
    # Hematologic
    "abnormal blood circulation",
)


DOID_CANDIDATES: tuple[str, ...] = (
    # Cancer
    "cancer", "breast cancer", "lung cancer", "colorectal cancer",
    "pancreatic cancer", "melanoma", "leukemia", "lymphoma",
    "glioblastoma multiforme", "hepatocellular carcinoma",
    # Metabolic / endocrine
    "diabetes mellitus", "type 1 diabetes mellitus",
    "type 2 diabetes mellitus", "obesity",
    # Neurodegenerative
    "Alzheimer's disease", "Parkinson's disease",
    "amyotrophic lateral sclerosis",
    # Cardiovascular
    "cardiomyopathy", "atherosclerosis", "hypertension",
    # Infectious
    "viral infectious disease", "bacterial infectious disease",
    "COVID-19", "influenza",
    # Inflammatory / autoimmune
    "inflammatory bowel disease", "ulcerative colitis",
    "Crohn's disease", "rheumatoid arthritis", "lupus erythematosus",
)


PR_CANDIDATES: tuple[str, ...] = (
    # Functional groups (root-ish — high recall)
    "transcription factor", "enzyme", "kinase", "phosphatase",
    "G protein-coupled receptor", "ion channel", "ligand",
    # Common proteins under canonical names (the PR catalogue
    # uses both gene-symbol-style and full names)
    "p53", "BRCA1", "alpha-synuclein",
    # Protein complexes (PR's main differentiator from HGNC)
    "ribosomal protein", "histone H3", "actin",
    # Structural
    "collagen", "tubulin", "myosin",
    # Other
    "ubiquitin", "albumin", "insulin", "hemoglobin",
)


OBA_CANDIDATES: tuple[str, ...] = (
    # Body-level traits
    "body weight", "body size", "body length", "lifespan",
    # Metabolic / physiological
    "metabolic rate", "basal metabolic rate", "heart rate",
    "blood pressure",
    # Reproductive
    "reproductive output", "fertility",
    # Developmental / cellular
    "developmental rate", "cell proliferation rate",
    # Gene-expression-axis attributes (very on-topic for transcripto)
    "gene expression level", "mRNA abundance", "protein abundance",
    # Immunological
    "immune response",
    # Holobiont-relevant
    "regeneration capacity", "circadian rhythm",
)


ENVO_CANDIDATES: tuple[str, ...] = (
    # Major biomes
    "marine biome", "freshwater biome", "terrestrial biome",
    "ocean", "river", "lake", "estuary", "coral reef",
    "intertidal zone", "benthic zone",
    # Soil / sediment
    "soil", "sediment", "rhizosphere",
    # Abiotic factors
    "temperature", "salinity", "pH", "dissolved oxygen", "light",
    # Sample-context
    "host-associated environment", "anaerobic environment",
    "aerobic environment",
)


REGISTRY: tuple[OntologyConfig, ...] = (
    OntologyConfig(
        ols_id="chebi", ontology_label="CHEBI",  axis="chemical",
        out_filename="chebi_concept_index.tsv",
        candidates=CHEBI_CANDIDATES,
    ),
    OntologyConfig(
        ols_id="hp",    ontology_label="HP",     axis="phenotype",
        out_filename="hp_concept_index.tsv",
        candidates=HP_CANDIDATES,
    ),
    OntologyConfig(
        ols_id="mp",    ontology_label="MP",     axis="phenotype",
        out_filename="mp_concept_index.tsv",
        candidates=MP_CANDIDATES,
    ),
    OntologyConfig(
        ols_id="doid",  ontology_label="DOID",   axis="disease",
        out_filename="doid_concept_index.tsv",
        candidates=DOID_CANDIDATES,
    ),
    OntologyConfig(
        ols_id="pr",    ontology_label="PR",     axis="protein",
        out_filename="pr_concept_index.tsv",
        candidates=PR_CANDIDATES,
    ),
    OntologyConfig(
        ols_id="oba",   ontology_label="OBA",    axis="attribute",
        out_filename="oba_concept_index.tsv",
        candidates=OBA_CANDIDATES,
    ),
    OntologyConfig(
        ols_id="envo",  ontology_label="ENVO",   axis="environment",
        out_filename="envo_concept_index.tsv",
        candidates=ENVO_CANDIDATES,
    ),
)


# Same cascade as build_cl_concept_index — exact label/synonym then
# prefix; substring omitted because labels in these ontologies are
# typically minimal already.
_MATCH_MODES = ("exact", "prefix")


def resolve_term(term: str, cfg: OntologyConfig) -> dict | None:
    """Apply the exact → prefix cascade against ``cfg.ols_id``. Returns
    ``{curie, label, description, mode}`` or ``None``."""
    for mode in _MATCH_MODES:
        docs = _search_ontology(term, cfg.ols_id, exact=(mode == "exact"))
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
        info = verify_iri(iri, cfg.ols_id)
        time.sleep(0.3)
        if info is None:
            continue
        return {
            "curie": obo_id,
            "label": info["label"] or (d.get("label") or ""),
            "description": info["description"],
            "mode": mode,
        }
    return None


def harvest_ontology(cfg: OntologyConfig, data_dir: Path) -> tuple[int, int]:
    """Run the cascade against one ontology and write its TSV. Returns
    (resolved_count, unresolved_count)."""
    out_tsv = data_dir / cfg.out_filename
    unresolved_path = data_dir / cfg.out_filename.replace(
        "_concept_index.tsv", "_unresolved.txt"
    )

    n_resolved = 0
    unresolved: list[str] = []
    seen_curies: set[str] = set()

    print(
        f"\n[ {cfg.ontology_label} ] harvesting {len(cfg.candidates):,} "
        f"candidates ...",
        file=sys.stderr,
    )

    with out_tsv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLUMNS, delimiter="\t")
        w.writeheader()
        for term in cfg.candidates:
            print(f"  resolving {term!r} ...", file=sys.stderr)
            info = resolve_term(term, cfg)
            if info is None:
                print("    ✗ no match", file=sys.stderr)
                unresolved.append(term)
                continue
            if info["curie"] in seen_curies:
                print(
                    f"    ⇣ duplicate CURIE {info['curie']} — skipping row.",
                    file=sys.stderr,
                )
                continue
            seen_curies.add(info["curie"])
            print(
                f"    ✓ [{info['mode']}] {info['curie']} {info['label']!r}",
                file=sys.stderr,
            )
            w.writerow({
                "term": term.lower(),
                "curie": info["curie"],
                "label": info["label"],
                "description": info["description"],
                "ontology": cfg.ontology_label,
                "source": (
                    f"v0.8.0 ontology harvest ({cfg.axis}, "
                    f"{info['mode']}-match)"
                ),
            })
            n_resolved += 1

    with unresolved_path.open("w", encoding="utf-8") as fh:
        fh.write(
            f"# Candidates for {cfg.ontology_label} that did not match "
            "via exact or prefix mode.\n"
            "# Either the candidate name is too informal / non-canonical\n"
            "# (e.g. an acronym not declared as a synonym in OLS4) or\n"
            "# the ontology genuinely doesn't cover this term.\n"
            "# Format: one term per line.\n\n"
        )
        for term in unresolved:
            fh.write(f"{term}\n")

    print(
        f"[ {cfg.ontology_label} ] wrote {out_tsv.name}: "
        f"{n_resolved} resolved, {len(unresolved)} unresolved.",
        file=sys.stderr,
    )
    return n_resolved, len(unresolved)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--only", type=str, default=None,
        help="Run only one ontology by OLS4 id (chebi / hp / mp / doid / "
             "pr / oba / envo). Default: all 7.",
    )
    args = parser.parse_args(argv)

    data_dir = Path(__file__).parent
    targets = REGISTRY
    if args.only:
        targets = tuple(c for c in REGISTRY if c.ols_id == args.only)
        if not targets:
            print(f"[err] unknown ontology id: {args.only}", file=sys.stderr)
            return 2

    total_r, total_u = 0, 0
    for cfg in targets:
        r, u = harvest_ontology(cfg, data_dir)
        total_r += r
        total_u += u

    print(
        f"\n[overall] {total_r} resolved, {total_u} unresolved across "
        f"{len(targets)} ontology / ies.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
