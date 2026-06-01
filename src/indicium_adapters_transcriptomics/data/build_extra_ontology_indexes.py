"""Maintainer script, harvest the v0.8.0 high-priority ontology indexes.

Following the 2026-05-28 OLS4 audit (see ``articles/_ols4_audit_report.md``),
seven ontologies were identified as the highest-yield additions for the
ASB transcriptomics benchmark phase:

- **CHEBI**, chemicals / drugs / ligands (perturbation experiments).
- **HP**, Human Phenotype Ontology (clinical transcriptomics).
- **MP**, Mammalian Phenotype Ontology (mouse studies).
- **DOID**, Human Disease Ontology (cancer / disorder vocab).
- **PR**, Protein Ontology (complex / isoform-level claims).
- **OBA**, Ontology of Biological Attributes (cross-species traits).
- **ENVO**, Environment Ontology (holobiont / ecology context).

Same cascade pattern as ``build_cl_concept_index.py``: candidate names
→ OLS4 exact-label / synonym / prefix search → re-verify by IRI lookup.
One bundled TSV per ontology so each can be re-harvested independently
when its candidate list grows.

Hand-curation of CURIEs is still avoided. The only curation here is
the candidate-NAME list per ontology, the OLS4 cascade picks the
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
    ols_id: str             # OLS4 ontology key, lowercase: "chebi", "hp"
    ontology_label: str     # Canonical prefix written to TSV (e.g. "CHEBI")
    axis: str               # axis tag for the ``source`` column
    out_filename: str       # bundled TSV filename
    candidates: tuple[str, ...]


# --- Candidate lists -------------------------------------------------------
#
# Each list is small on purpose, ~20-40 terms per ontology covering the
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
    # Functional groups (root-ish, high recall)
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


# --- v0.9.0 Tier-1 universal: MS / EC / UO / MOD / RXNO / HCAO ------------


MS_CANDIDATES: tuple[str, ...] = (
    # Ionisation
    "electrospray ionization", "matrix-assisted laser desorption ionization",
    "chemical ionization", "atmospheric pressure chemical ionization",
    "fast atom bombardment",
    # Mass analysers
    "quadrupole", "time-of-flight", "ion trap", "orbitrap",
    "fourier transform ion cyclotron resonance",
    # Common combined methods (LC-MS, GC-MS spelt out)
    "liquid chromatography mass spectrometry",
    "gas chromatography mass spectrometry",
    "tandem mass spectrometry",
    # Quantification
    "label-free quantification", "isobaric tag for relative and absolute quantification",
    "stable isotope labeling by amino acids in cell culture",
    "tandem mass tag",
    # Acquisition modes
    "data-dependent acquisition", "data-independent acquisition",
    "selected reaction monitoring", "multiple reaction monitoring",
    # Concepts
    "mass-to-charge ratio", "precursor ion", "product ion",
    "fragment ion spectrum", "monoisotopic mass",
    # Sample prep
    "trypsin", "in-gel digestion", "in-solution digestion",
)


# EC_CANDIDATES kept here for reference but EC is deferred, see the
# REGISTRY comment. The OLS4 ``ec`` ontology re-exports GO molecular-
# function terms without human-readable labels, so cascade label-
# matching can't catch common enzyme class names. v0.9.x will land
# EC via a direct EC-number mapping file from ExPASy.
_EC_CANDIDATES_DEFERRED: tuple[str, ...] = (
    "oxidoreductase", "transferase", "hydrolase",
    "lyase", "isomerase", "ligase", "translocase",
    "kinase", "phosphatase", "protease", "peptidase",
    "nuclease", "polymerase", "helicase",
    "DNA polymerase", "RNA polymerase",
    "DNA helicase", "DNA ligase",
    "reverse transcriptase", "topoisomerase",
    "transposase", "integrase", "recombinase",
    "dehydrogenase", "oxidase", "reductase",
    "methyltransferase", "acetyltransferase",
)


UO_CANDIDATES: tuple[str, ...] = (
    # Mass units
    "gram", "kilogram", "milligram", "microgram", "nanogram", "picogram",
    # Volume units
    "liter", "milliliter", "microliter", "nanoliter",
    # Concentration units (molar), note: 'molar' might be ambiguous;
    # try the full canonical names
    "molar", "millimolar", "micromolar", "nanomolar", "picomolar",
    "mole", "millimole", "micromole",
    # Time
    "second", "minute", "hour", "day", "week", "month", "year",
    # Temperature
    "degree Celsius", "kelvin",
    # Counting / cellular units
    "molecule", "cell", "particle",
    # Length
    "meter", "centimeter", "millimeter", "micrometer", "nanometer",
    # Generic
    "percent", "fold change",
)


MOD_CANDIDATES: tuple[str, ...] = (
    # PSI-MOD uses "<verb-form> residue" / "<verb-form> L-<aa>" patterns
    # rather than the PTM-class names common in literature
    # ("phosphorylation"). Candidates aligned with that convention.
    # Class-level "residue" entries
    "phosphorylated residue", "acetylated residue",
    "methylated residue", "ubiquitinated residue",
    "sumoylated residue", "neddylated residue",
    "glycosylated residue", "N-glycosylated residue",
    "O-glycosylated residue",
    "palmitoylated residue", "myristoylated residue",
    "prenylated residue", "farnesylated residue",
    "geranylgeranylated residue",
    # Specific phospho / methyl / acetyl positions
    "phosphorylated L-serine", "phosphorylated L-threonine",
    "phosphorylated L-tyrosine", "phosphorylated L-histidine",
    "O-phosphorylated residue", "N-phosphorylated residue",
    "monomethylated L-lysine", "dimethylated L-lysine",
    "trimethylated L-lysine",
    "N6-acetyl-L-lysine",
    # Other common modifications
    "oxidized residue", "deamidated residue",
    "citrullinated residue", "nitrosylated residue",
    "protein modification",
)


RXNO_CANDIDATES: tuple[str, ...] = (
    # Generic reaction types
    "chemical reaction", "biochemical reaction",
    "oxidation reaction", "reduction reaction",
    "hydrolysis reaction", "condensation reaction",
    "phosphorylation reaction", "dephosphorylation reaction",
    "methylation reaction", "acetylation reaction",
    "ubiquitination reaction",
    # Biosynthesis / catabolism (might match as classes)
    "biosynthesis", "catabolism", "anabolism",
    # Common enzymatic reactions
    "transamination", "decarboxylation",
    # Pathway-adjacent
    "redox reaction", "esterification",
)


HCAO_CANDIDATES: tuple[str, ...] = (
    # HCA / Human Cell Atlas project-specific or harmonised terms
    "human cell", "peripheral blood mononuclear cell",
    "tumor infiltrating lymphocyte", "tissue-resident memory T cell",
    # Major organs as HCA references them
    "lung", "liver", "heart", "kidney", "spleen", "thymus",
    "pancreas", "brain", "bone marrow", "adipose tissue",
    # Specific HCA cell types
    "alveolar macrophage", "Kupffer cell", "dendritic cell",
    "myeloid cell", "lymphoid cell",
    # Developmental stage / project descriptors
    "adult human", "fetal human",
    # Common HCA experimental terms
    "single cell RNA sequencing",
)


# --- v0.9.0 Tier-2 organism-paired: anatomy + phenotype per species -------


FBBT_CANDIDATES: tuple[str, ...] = (
    # Life stages
    "adult", "embryo", "larva", "pupa", "first instar larva",
    "second instar larva", "third instar larva",
    # Imaginal discs (iconic Drosophila tissue)
    "imaginal disc", "wing disc", "eye disc", "leg disc",
    "antennal disc", "haltere disc", "genital disc",
    # Adult organ systems
    "ovary", "testis", "germline",
    "brain", "ventral nerve cord", "central nervous system",
    "thorax", "abdomen", "head",
    "gut", "midgut", "hindgut", "foregut",
    "Malpighian tubule", "fat body",
    "compound eye", "ocellus",
    # Tissues
    "muscle", "indirect flight muscle",
    "epidermis", "trachea",
)


WBBT_CANDIDATES: tuple[str, ...] = (
    # Life stages and sexes
    "adult hermaphrodite", "adult male",
    "L1 larva", "L2 larva", "L3 larva", "L4 larva",
    "dauer larva", "embryo",
    # Reproductive
    "germ line", "gonad", "oocyte", "spermatheca",
    # Digestive
    "intestine", "pharynx", "rectum",
    # Nervous system
    "neuron", "motor neuron", "sensory neuron", "interneuron",
    "nerve ring",
    # Muscle
    "body wall muscle", "pharyngeal muscle", "vulval muscle",
    # Other
    "vulva", "uterus", "hypodermis", "cuticle",
    "coelomocyte", "excretory cell",
)


ZFA_CANDIDATES: tuple[str, ...] = (
    # Life stages
    "embryo", "larva", "juvenile", "adult",
    # Early development structures
    "somite", "notochord", "neural tube", "neural crest",
    "yolk", "yolk syncytial layer",
    # Adult / late structures
    "heart", "atrium", "ventricle",
    "brain", "spinal cord", "forebrain", "midbrain", "hindbrain",
    "retina", "lens", "otic vesicle", "otolith",
    "liver", "intestine", "swim bladder", "pancreas",
    # Fins
    "fin", "pectoral fin", "pelvic fin", "dorsal fin",
    "anal fin", "caudal fin", "fin ray",
    # Pharyngeal apparatus
    "pharyngeal arch", "gill",
)


EMAPA_CANDIDATES: tuple[str, ...] = (
    # Life stages
    "embryo", "fetus", "neonate", "juvenile", "adult",
    # Early dev structures
    "somite", "neural tube", "neural crest",
    "neural plate", "primitive streak",
    # Limb development
    "limb bud", "forelimb bud", "hindlimb bud", "limb",
    # Cardio-vascular
    "heart", "heart tube", "atrium", "ventricle",
    # Brain regions
    "brain", "forebrain", "midbrain", "hindbrain",
    "spinal cord", "central nervous system",
    # Organ primordia
    "liver", "lung", "pancreas",
    "kidney", "gonad",
    # Maternal-fetal
    "placenta", "umbilical cord", "yolk sac",
    # Hematopoietic
    "fetal liver", "bone marrow",
)


ZP_CANDIDATES: tuple[str, ...] = (
    # Body morphology
    "abnormal somite morphology", "abnormal trunk morphology",
    "abnormal head morphology", "abnormal eye morphology",
    "abnormal heart morphology", "abnormal fin morphology",
    "abnormal pigmentation",
    # Developmental
    "decreased survival", "delayed development",
    "embryonic lethal", "larval lethal",
    # Behavior
    "abnormal swimming behavior", "abnormal touch response",
    # Tissue-level
    "abnormal vasculature morphology",
    "abnormal blood circulation",
    "abnormal jaw morphology", "abnormal swim bladder morphology",
    # Specific common ones
    "absent pigment cell", "decreased pigmentation",
)


XPO_CANDIDATES: tuple[str, ...] = (
    # Regeneration (Xenopus is the regen model)
    "abnormal limb regeneration", "abnormal tail regeneration",
    # Morphology
    "abnormal limb morphology", "abnormal tail morphology",
    "abnormal eye morphology", "abnormal heart morphology",
    "abnormal head morphology",
    # Developmental
    "delayed development", "abnormal metamorphosis",
    "embryonic lethal", "tadpole lethal",
    # Behavior
    "abnormal swimming behavior",
    # Tissue
    "abnormal pigmentation", "abnormal somite morphology",
    "abnormal blood vessel morphology",
)


# ZP (Zebrafish Phenotype Ontology) candidates: high-frequency abnormal-
# phenotype terms cited in zebrafish transcriptomics literature
# (Wagner 2018, knockdown / knockout phenotypes). ZP labels follow the
# pattern ``<entity> <quality> [, abnormal]`` so we list them in that
# canonical form. The cascade picks via exact label match.
ZP_CANDIDATES: tuple[str, ...] = (
    # Morphology of developmental structures
    "neural tube morphology, abnormal",
    "somite morphology, abnormal",
    "notochord morphology, abnormal",
    "fin morphology, abnormal",
    "pectoral fin morphology, abnormal",
    # Cardiac / circulation
    "heart morphology, abnormal",
    "blood circulation, abnormal",
    "blood vessel morphology, abnormal",
    # Eye / pigmentation
    "eye morphology, abnormal",
    "pigmentation, abnormal",
    # Whole-organism
    "body size, abnormal",
    "development, delayed",
    "embryonic lethal",
    # Swim-bladder + locomotion (zebrafish-specific)
    "swim bladder morphology, abnormal",
    "swimming behavior, abnormal",
)


# GO (Gene Ontology) candidates: high-frequency biological processes and
# molecular functions cited across the corpus. The three GO roots
# (molecular_function, cellular_component, biological_process) already
# ship in the static concept_index.tsv; we add the concrete terms here
# so claims like "WNT4 enriched_for Wnt signaling pathway" can resolve.
# Mix of process / function / pathway terms; ~50 entries to keep
# benchmark recall up on feature-level qualifiers.
GO_CANDIDATES: tuple[str, ...] = (
    # --- Cell cycle / proliferation / death --------------------------------
    "cell cycle", "mitotic cell cycle", "cell division", "DNA replication",
    "DNA repair", "apoptotic process", "programmed cell death",
    "autophagy", "cellular senescence",
    # --- Gene expression machinery -----------------------------------------
    "gene expression", "transcription, DNA-templated",
    "translation", "mRNA splicing, via spliceosome",
    "RNA processing", "rRNA processing",
    "regulation of transcription, DNA-templated",
    # --- Protein lifecycle -------------------------------------------------
    "protein folding", "protein phosphorylation",
    "ubiquitin-dependent protein catabolic process", "proteolysis",
    # --- Signal transduction (most common in transcriptomics) -------------
    "signal transduction", "Wnt signaling pathway",
    "MAPK cascade", "Notch signaling pathway",
    "transforming growth factor beta receptor signaling pathway",
    "BMP signaling pathway", "Hippo signaling pathway",
    "JAK-STAT cascade", "Hedgehog signaling pathway",
    # --- Metabolism --------------------------------------------------------
    "glycolytic process", "oxidative phosphorylation",
    "fatty acid biosynthetic process", "gluconeogenesis",
    "cellular respiration", "ATP biosynthetic process",
    # --- Immune / stress response -----------------------------------------
    "immune response", "innate immune response",
    "adaptive immune response", "inflammatory response",
    "response to oxidative stress", "cellular response to hypoxia",
    "cellular response to DNA damage stimulus",
    # --- Development / differentiation / morphogenesis --------------------
    "cell differentiation", "embryonic development",
    "neurogenesis", "gastrulation", "angiogenesis",
    "epithelial to mesenchymal transition",
    "wound healing", "regeneration",
    # --- Cell behaviour ---------------------------------------------------
    "cell migration", "cell adhesion",
    "cell proliferation",
    # --- Molecular function (top picks) -----------------------------------
    "DNA binding", "RNA binding", "protein binding",
    "kinase activity", "protein kinase activity",
    "ATP binding", "transcription factor activity",
)


REGISTRY: tuple[OntologyConfig, ...] = (
    # --- v0.8.0, Top-7 absolute priority -----------------------------------
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
    # --- v0.9.0 Tier-1 universal --------------------------------------------
    OntologyConfig(
        ols_id="ms",    ontology_label="MS",     axis="method",
        out_filename="ms_concept_index.tsv",
        candidates=MS_CANDIDATES,
    ),
    # EC removed from v0.9.0 REGISTRY, the OLS4 ``ec`` ontology re-
    # exports GO molecular-function terms without human-readable labels
    # so cascade label-matching can't catch the common enzyme class
    # names (kinase / phosphatase / etc.). 3 / 30 hits in the first
    # harvest run. Defer to v0.9.x with a direct EC-name → EC-number
    # mapping (e.g. from ExPASy enzyme.dat).
    # OntologyConfig(ols_id="ec", ontology_label="EC", axis="enzyme",
    #                out_filename="ec_concept_index.tsv",
    #                candidates=EC_CANDIDATES),
    OntologyConfig(
        ols_id="uo",    ontology_label="UO",     axis="unit",
        out_filename="uo_concept_index.tsv",
        candidates=UO_CANDIDATES,
    ),
    OntologyConfig(
        ols_id="mod",   ontology_label="MOD",    axis="modification",
        out_filename="mod_concept_index.tsv",
        candidates=MOD_CANDIDATES,
    ),
    # RXNO and HCAO deferred to v0.9.x, same issue as EC: in OLS4,
    # ``rxno`` is dominated by CHEBI / BFO imports without named
    # reaction labels, and ``hcao`` re-publishes CARO / CL / CLO / GO
    # terms (redundant with our existing CL index). 2 / 16 and 0 / 17
    # hits respectively in the first harvest. Defer until we have a
    # source that exposes reaction / atlas labels directly.
    # OntologyConfig(ols_id="rxno", ontology_label="RXNO", axis="reaction",
    #                out_filename="rxno_concept_index.tsv",
    #                candidates=RXNO_CANDIDATES),
    # OntologyConfig(ols_id="hcao", ontology_label="HCAO", axis="cell_atlas",
    #                out_filename="hcao_concept_index.tsv",
    #                candidates=HCAO_CANDIDATES),
    # --- v0.9.0 Tier-2 organism-paired (anatomy + phenotype) ---------------
    OntologyConfig(
        ols_id="fbbt",  ontology_label="FBbt",   axis="drosophila_anatomy",
        out_filename="fbbt_concept_index.tsv",
        candidates=FBBT_CANDIDATES,
    ),
    OntologyConfig(
        ols_id="wbbt",  ontology_label="WBbt",   axis="celegans_anatomy",
        out_filename="wbbt_concept_index.tsv",
        candidates=WBBT_CANDIDATES,
    ),
    OntologyConfig(
        ols_id="zfa",   ontology_label="ZFA",    axis="zebrafish_anatomy",
        out_filename="zfa_concept_index.tsv",
        candidates=ZFA_CANDIDATES,
    ),
    OntologyConfig(
        ols_id="emapa", ontology_label="EMAPA",  axis="mouse_anatomy",
        out_filename="emapa_concept_index.tsv",
        candidates=EMAPA_CANDIDATES,
    ),
    # --- v0.12 expansion: ZP (zebrafish phenotype) re-tried ---------------
    # The v0.9.0 attempt failed because the candidates were written as
    # natural English (``abnormal <entity>``) while ZP's canonical labels
    # follow ``<entity> <quality>, abnormal``. The v0.12 candidate list
    # is rewritten in the canonical form so the exact-label cascade
    # catches it. RXNO and HCAO remain deferred for the reasons
    # documented in v0.9.0 (CHEBI/CL re-publishers with low effective
    # transcripto yield).
    OntologyConfig(
        ols_id="zp",    ontology_label="ZP",    axis="zebrafish_phenotype",
        out_filename="zp_concept_index.tsv",
        candidates=ZP_CANDIDATES,
    ),
    # --- v0.12 expansion: GO (biological process + molecular function) ----
    # The static concept_index.tsv already ships the three GO roots; this
    # entry adds ~50 high-frequency concrete process / function terms so
    # feature-level claims (Wnt signaling pathway, apoptotic process,
    # ATP binding, ...) get a real GO CURIE.
    OntologyConfig(
        ols_id="go",    ontology_label="GO",    axis="biological_process",
        out_filename="go_concept_index.tsv",
        candidates=GO_CANDIDATES,
    ),
    OntologyConfig(
        ols_id="xpo",   ontology_label="XPO",    axis="xenopus_phenotype",
        out_filename="xpo_concept_index.tsv",
        candidates=XPO_CANDIDATES,
    ),
)


# Same cascade as build_cl_concept_index, exact label/synonym then
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
                    f"    ⇣ duplicate CURIE {info['curie']}, skipping row.",
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
