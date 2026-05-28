# OLS4 ontology audit for the ASB transcriptomics benchmark

Audit date: 2026-05-28 · Catalogue size: **280 ontologies**

Two sections:

1. **Top picks** — curated short list of the ontologies most
   likely to fire on a diverse transcriptomics corpus, with
   one-liner justifications. **Start here.**
2. **Bucket dump** — heuristic keyword classification of every
   ontology in the OLS4 catalogue. Useful as raw material for
   building per-axis harvesters in future versions.

## 🎯 Top picks (curated)

| Prefix | Title | # terms | Why it matters |
|---|---|---:|---|
| `HP` | Human Phenotype Ontology | 32,085 | Phenotype ontology — clinical transcripto papers (HPO is THE phenotype vocab for human disease) |
| `CHEBI` | Chemical Entities of Biological Interest | 224,691 | Small molecules — perturbation experiments (drug treatments, ligands, metabolites) |
| `DOID` | Human Disease Ontology | 19,470 | Disease ontology — clinical / cancer transcriptomics (complements MONDO, more drug-perturbation-friendly) |
| `HPO` | (not indexed in OLS4) | — | Same as HP if listed separately |
| `MP` | The Mammalian Phenotype Ontology | 35,181 | Mammalian Phenotype — mouse studies (Tabula Muris, transgenic mice) |
| `OBA` | Ontology of Biological Attributes (OBA) | 82,689 | Ontology of Biological Attributes — cross-kingdom traits, multi-species KGs |
| `PR` | PRotein Ontology (PRO) | 365,852 | Protein Ontology — protein-level claims, isoforms, complexes |
| `GeXO` | Gene Expression Ontology | 166,254 | Gene Expression Ontology — direct overlap with our domain, possibly redundant with our concept index |
| `ReXO` | Regulation of Gene Expression Ontology | 158,239 | Regulation of Gene Expression Ontology — regulation-process vocabulary beyond INO |
| `ReTO` | Regulation of Transcription Ontology | 147,738 | Regulation of Transcription Ontology — transcriptional control vocab |
| `EC` | Enzyme Nomenclature | 13,179 | Enzyme Commission codes — enzymatic activity, multi-omics bridge |
| `MS` | Mass spectrometry ontology | 4,097 | Mass spectrometry — proteomics/metabolomics overlap (cf. metabolomics adapter) |
| `OMIT` | Ontology for MIRNA Target | 90,928 | miRNA Target — for miRNA-focused transcriptomics studies |
| `HCAO` | Human Cell Atlas Ontology | 22,606 | Human Cell Atlas Ontology — cell-type vocabulary harmonised across atlases |
| `EMAPA` | Mouse Developmental Anatomy Ontology | 8,882 | Mouse Developmental Anatomy — Tabula Muris context, embryonic transcriptomics |
| `FBbt` | Drosophila gross anatomy | 29,751 | Drosophila gross anatomy — paired with our v0.6.2 FlyBase gene index |
| `WBbt` | C. elegans Gross Anatomy Ontology | 7,573 | C. elegans gross anatomy — paired with our v0.6.2 WormBase gene index |
| `ZFA` | Zebrafish anatomy and development ontology | 3,285 | Zebrafish anatomy — paired with our v0.6.2 ZFIN gene index |
| `ZP` | Zebrafish Phenotype Ontology (ZP) | 63,809 | Zebrafish Phenotype — Zebrafish disease modelling |
| `XPO` | Xenopus Phenotype Ontology | 30,225 | Xenopus Phenotype — Xenopus regeneration studies |
| `ENVO` | The Environment Ontology | 6,906 | Environment Ontology — holobiont context, sample collection metadata |
| `NBO` | Neuro Behavior Ontology | 4,546 | Neuro Behavior Ontology — neural transcriptomics, behavioural phenotypes |
| `RXNO` | RXNO | 1,019 | Reaction Ontology — biochemical reaction context |
| `MOD` | Protein modification | 2,116 | Protein modification (PTMs) — phospho/methyl/acetyl-proteomics overlap |
| `UO` | Units of measurement ontology | 574 | Units of measurement — quantitative claim normalisation |
| `BTO` | The BRENDA Tissue Ontology (BTO) | 6,569 | BRENDA Tissue / Cell Type — alternative to CL/UBERON, popular in enzyme literature |
| `OGSF` | Ontology of Genetic Susceptibility Factor | 406 | Genetic Susceptibility Factor — GWAS context |
| `CMPO` | Cellular Microscopy Phenotype Ontology | 1,058 | Cellular Microscopy Phenotype — imaging-derived single-cell phenotypes |
| `MEDDRA` | (not indexed in OLS4) | — | Medical Dictionary for Regulatory Activities — drug-induced expression studies |

---

## Bucket dump (heuristic keyword classification)

Each ontology is classified by keyword matching on its
`title + description + id`. The classification is rough — use
the **Top picks** above for short-term decisions, and dip into
the buckets when adding a new axis (e.g. cell-line or imaging
phenotype) to find candidates.

## Already covered by the adapter (16 ontologies)

| Prefix | Title | # terms | Description |
|---|---|---:|---|
| `NCBITaxon` | NCBI organismal classification | 2,739,624 | An ontology representation of the NCBI organismal taxonomy |
| `NCIT` | NCI Thesaurus OBO Edition | 209,910 | NCI Thesaurus (NCIt)is a reference terminology that includes broad coverage of the cancer domain, including cancer related diseases, findings and abnormalities. The NCIt OBO Edition aims to increas... |
| `EFO` | Experimental Factor Ontology | 93,358 | The Experimental Factor Ontology (EFO) provides a systematic description of many experimental variables available in EBI databases, and for projects such as the NHGRI-EBI GWAS catalog. It combines ... |
| `GO` | Gene Ontology | 84,523 | The Gene Ontology (GO) provides a framework and set of concepts for describing the functions of gene products from all organisms. |
| `MONDO` | Mondo Disease Ontology | 58,940 | A semi-automatically constructed ontology that merges in multiple disease resources to yield a coherent merged ontology. |
| `UBERON` | Uber-anatomy ontology | 27,295 | Uberon is an integrated cross-species anatomy ontology representing a variety of entities classified according to traditional anatomical criteria such as structure, function and developmental linea... |
| `CL` | Cell Ontology | 19,151 | An ontology of cell types. |
| `PATO` | PATO - the Phenotype And Trait Ontology | 8,625 | An ontology of phenotypic qualities (properties, attributes or characteristics). |
| `OBI` | Ontology for Biomedical Investigations | 5,219 | An ontology for representing biomedical investigations, including study designs, the collection and preparation of the targets of investigation, assays, instrumentation and reagents used, as well a... |
| `EDAM` | EDAM - The ontology of data analysis and management | 3,524 | EDAM is a domain ontology of data analysis and data management in bio- and other sciences, and science-based applications. It comprises concepts related to analysis, modelling, optimisation, and da... |
| `SO` | Sequence types and features ontology | 2,747 | A structured controlled vocabulary for sequence annotation, for the exchange of annotation data and for the description of sequence objects in databases. |
| `SIO` | Semanticscience Integrated Ontology | 1,585 | The Semanticscience Integrated Ontology (SIO) provides a simple, integrated ontology of types and relations for rich description of objects, processes and their attributes. |
| `STATO` | STATO: the statistical methods ontology | 1,053 | STATO is the statistical methods ontology. It contains concepts and properties related to statistical methods, probability distributions and other concepts related to statistical analysis, includin... |
| `INO` | INO: Interaction Network Ontology | 328 | he Interaction Network Ontology (INO) is an ontology in the domain of interactions and interaction networks. INO represents general and species-neutral types of interactions and interaction network... |
| `IAO` | Information Artifact Ontology (IAO) | 266 | The Information Artifact Ontology (IAO) is an ontology of information entities, originally driven by work by the OBI digital entity and realizable information entity branch. |
| `RO` | Relation Ontology | 58 | Relationship types shared across multiple ontologies |

## 🎯 HIGH — Disease ontologies (clinical transcripto / cancer / disorders) (22 ontologies)

| Prefix | Title | # terms | Description |
|---|---|---:|---|
| `HP` | Human Phenotype Ontology | 32,085 | The Human Phenotype Ontology (HPO) provides a standardized vocabulary of phenotypic abnormalities and clinical features encountered in human disease. |
| `CIDO` | CIDO: Ontology of Coronavirus Infectious Disease | 31,924 | The Ontology of Coronavirus Infectious Disease (CIDO) is a community-driven open-source biomedical ontology in the area of coronavirus infectious disease. The CIDO is developed to provide standardi... |
| `DOID` | Human Disease Ontology | 19,470 | The Disease Ontology has been developed as a standardized ontology for human disease with the purpose of providing the biomedical community with consistent, reusable and sustainable descriptions of... |
| `ORDO` | Orphanet Rare Disease Ontology | 16,036 | The Orphanet Rare Disease ontology (ORDO) is jointly developed by Orphanet and the EBI to provide a structured vocabulary for rare diseases capturing relationships between diseases, genes and other... |
| `CPONT` | Critical Path Ontology | 9,119 | An application ontology for integrating data from the Critical Path Institute's Rare Disease Cures Accelerator Data and Analytica Platform (RDCA-DAP). |
| `SCDO` | Sickle Cell Disease Ontology | 2,071 | An ontology for the standardization of terminology and integration of knowledge about Sickle Cell Disease. |
| `ADO` | Alzheimer's Disease Ontology (ADO) | 1,963 | Alzheimer's Disease Ontology is a knowledge-based ontology that encompasses varieties of concepts related to Alzheimer'S Disease, foundamentally structured by upper level Basic Formal Ontology(BFO)... |
| `EPIO` | Epilepsy Ontology (EPIO) | 1,916 | The Epilepsy Ontology (EPIO) is an assembly of structured knowledge on various aspects of epilepsy, developed according to basic formal ontology (BFO) and Open Biological and Biomedical Ontology (O... |
| `OHD` | Oral Health and Disease Ontology | 1,820 | The Oral Health and Disease Ontology is used for representing the diagnosis and treatment of dental maladies. |
| `APOLLO_SV` | Apollo Structured Vocabulary | 1,716 | An OWL2 ontology of phenomena in infectious disease epidemiology and population biology for use in epidemic simulation. |
| `MFOMD` | Mental Disease Ontology | 1,172 | The Mental Disease Ontology is developed to facilitate representation for all aspects of mental disease. It is an extension of the Ontology for General Medical Science (OGMS) and Mental Functioning... |
| `SYMP` | Symptom Ontology | 1,019 | The Symptom Ontology has been developed as a standardized ontology for symptoms of human diseases. |
| `MPATH` | Mouse pathology ontology | 891 | A structured controlled vocabulary of mutant and transgenic mouse pathology phenotypes |
| `CVDO` | Cardiovascular Disease Ontology | 720 | An ontology to describe entities related to cardiovascular diseases |
| `IDO` | Infectious Disease Ontology | 518 | A set of interoperable ontologies that will together provide coverage of the infectious disease domain. IDO core is the upper-level ontology that hosts terms of general relevance across the domain,... |
| `IDO-COVID-19` | The COVID-19 Infectious Disease Ontology | 485 | The COVID-19 Infectious Disease Ontology (IDO-COVID-19) is an extension of the Infectious Disease Ontology (IDO) and the Virus Infectious Disease Ontology (VIDO). IDO-COVID-19 follows OBO Foundry g... |
| `OGSF` | Ontology of Genetic Susceptibility Factor | 406 | An application ontology to represent genetic susceptibility to a specific disease, adverse event, or a pathological process. |
| `PCO` | Population and Community Ontology | 220 | The Population and Community Ontology (PCO) describes material entities, qualities, and processes related to collections of interacting organisms such as populations and communities. It is taxon ne... |
| `OGMS` | Ontology for General Medical Science | 186 | The Ontology for General Medical Science (OGMS) is an ontology of entities involved in a clinical encounter. OGMS includes very general terms that are used across medical disciplines, including: 'd... |
| `ExO` | Exposure ontology (ExO) | 168 | ExO is intended to bridge the gap between exposure science and diverse environmental health disciplines including toxicology, epidemiology, disease surveillance, and epigenetics. |
| `DISDRIV` | Disease Drivers | 84 | Drivers of human diseases including environmental, maternal and social exposures. |
| `TRANS` | Pathogen Transmission Ontology | 35 | The Pathogen Transmission Ontology describes the transmission methods of human disease pathogens describing how a pathogen is transmitted from one host, reservoir, or source to another host. The pa... |

## 🎯 HIGH — Phenotype ontologies (HP, MP, ZP, FBcv, …) (32 ontologies)

| Prefix | Title | # terms | Description |
|---|---|---:|---|
| `UPHENO` | Unified Phenotype Ontology (uPheno) | 193,296 | The uPheno ontology integrates multiple phenotype ontologies into a unified cross-species phenotype ontology. |
| `OBA` | Ontology of Biological Attributes (OBA) | 82,689 | A collection of biological attributes (traits) covering all kingdoms of life. Interoperable with VT (vertebrate trait ontology) and TO (plant trait ontology). Extends PATO. |
| `ZP` | Zebrafish Phenotype Ontology (ZP) | 63,809 | The Zebrafish Phenotype Ontology formally defines all phenotypes of the Zebrafish model organism. |
| `MP` | The Mammalian Phenotype Ontology | 35,181 | The Mammalian Phenotype Ontology is being developed by Cynthia L. Smith, Susan M. Bello, Anna V. Anagnostopoulos, Carroll W. Goldsmith and Janan T. Eppig, as part of the Mouse Genome Database (MGD)... |
| `FLOPO` | Flora Phenotype Ontology | 34,858 | Traits and phenotypes of flowering plants occurring in digitized Floras |
| `XPO` | Xenopus Phenotype Ontology | 30,225 | XPO represents anatomical, cellular, and gene function phenotypes occurring throughout the development of the African frogs Xenopus laevis and tropicalis. |
| `MICRO` | Ontology of Prokaryotic Phenotypic and Metabolic Characters | 19,246 | An ontology of prokaryotic phenotypic and metabolic characters |
| `GSSO` | GSSO - the Gender, Sex, and Sexual Orientation ontology | 13,192 | GSSO is the Gender, Sex, and Sex Orientation ontology, including terms related to gender identity and expression, sexual and romantic identity and orientation, and sexual and reproductive behavior. |
| `FOVT` | FuTRES Ontology of Vertebrate Traits | 11,980 | FuTRES Ontology of Vertebrate Traits is an application ontology used to convert vertebrate trait data in spreadsheet to triples. FOVT leverages the BioCollections Ontology (BCO) to link observation... |
| `PHIPO` | Pathogen Host Interactions Phenotype Ontology | 10,779 | Ontology of species-neutral phenotypes observed in pathogen-host interactions. |
| `FYPO` | Fission Yeast Phenotype Ontology (FYPO) | 10,541 | A formal ontology of phenotypes observed in fission yeast. |
| `WBPhenotype` | C elegans Phenotype Ontology | 6,546 | Ontology about C. elegans and other nematode phenotypes |
| `TO` | Plant Trait Ontology | 6,194 | A controlled vocabulary to describe phenotypic traits in plants. |
| `ECOCORE` | An ontology of core ecological entities | 5,586 | Ecocore is a community ontology for the concise and controlled description of ecological traits of organisms. |
| `NBO` | Neuro Behavior Ontology | 4,546 | An ontology of human and animal behaviours and behavioural phenotypes |
| `VT` | The Vertebrate Trait Ontology | 4,063 | The Vertebrate Trait Ontology is a controlled vocabulary for the description of traits (measurable or observable characteristics) pertaining to the morphology, physiology, or development of vertebr... |
| `PLANP` | Planarian Phenotype Ontology (PLANP) | 3,560 | Planarian Phenotype Ontology is an ontology of phenotypes observed in the planarian Schmidtea mediterranea. |
| `BCIO` | The Behaviour Change Intervention Ontology | 2,545 | The Behaviour Change Intervention Ontology is an ontology for all aspects of human behaviour change interventions and their evaluation. |
| `OMP` | Ontology of Microbial Phenotypes | 2,309 | An ontology of phenotypes covering microbes |
| `FBcv` | Drosophila Phenotype Ontology | 1,552 | An ontology of commonly encountered and/or high level Drosophila phenotypes. |
| `DDPHENO` | Dicty Phenotype Ontology (DDPHENO) | 1,373 | A structured controlled vocabulary of phenotypes of the slime-mould Dictyostelium discoideum. |
| `CMPO` | Cellular Microscopy Phenotype Ontology | 1,058 | CMPO is a species neutral ontology for describing general phenotypic observations relating to the whole cell, cellular components, cellular processes and cell populations. |
| `COLAO` | Coleoptera Anatomy Ontology | 894 | The Coleoptera Anatomy Ontology contains terms used for describing the anatomy and phenotype of beetles in biodiversity research. |
| `MMO` | Measurement method ontology | 869 | A representation of the variety of methods used to make clinical and phenotype measurements. |
| `SPD` | Spider Ontology | 832 | An ontology for spider comparative biology including anatomical parts (e.g. leg, claw), behavior (e.g. courtship, combing) and products (i.g. silk, web, borrow). |
| `LEPAO` | Lepidoptera Anatomy Ontology | 822 | The Lepidoptera Anatomy Ontology contains terms used for describing the anatomy and phenotype of moths and butterflies in biodiversity research. |
| `APO` | Ascomycete Phenotype Ontology (APO) | 619 | A structured controlled vocabulary for the phenotypes of Ascomycete fungi. |
| `GENO` | GENO ontology | 589 | GENO is an OWL model of genotypes, their more fundamental sequence components, and links to related biological and experimental entities.  At present many parts of the model are exploratory and set... |
| `SIBO` | Social Insect Behavior Ontology | 525 | Social Behavior in insects |
| `GALLONT` | Plant Gall Ontology | 394 | Ontology of plant gall phenotypes. Plant galls are novel plant structures, generated by plants in response to biotic stressors. This ontology is used to annotate gall phenotypes (e.g., their colors... |
| `MFMO` | Mammalian Feeding Muscle Ontology | 272 | The Mammalian Feeding Muscle Ontology is an antomy ontology for the muscles of the head and neck that participate in feeding, swallowing, and other oral-pharyngeal behaviors. |
| `TEDDY` | Terminology for Description of Dynamics | 167 | TEDDY is an ontology for dynamical behaviours, observable dynamical phenomena, and control elements of bio-models and biological systems in Systems and Synthetic Biology. |

## 🎯 HIGH — Chemicals / drugs (ChEBI, perturbation experiments) (12 ontologies)

| Prefix | Title | # terms | Description |
|---|---|---:|---|
| `SLM` | SwissLipids | 1,003,946 | SwissLipids is a curated resource that provides information about known lipids, including lipid structure, metabolism, interactions, and subcellular and tissue localization. Information is curated ... |
| `DRON` | The Drug Ontology | 778,648 | An ontology to support comparative effectiveness researchers studying claims data. |
| `CHEBI` | Chemical Entities of Biological Interest | 224,691 | A structured classification of molecular entities of biological interest focusing on 'small' chemical compounds. |
| `LIPIDMAPS` | LIPID MAPS | 50,413 |  |
| `ENM` | eNanoMapper ontology | 25,918 | The eNanoMapper project (https://www.enanomapper.net/), NanoCommons project (https://www.nanocommons.eu/) and ACEnano project (http://acenano-project.eu/) are creating a pan-European computational ... |
| `CHIRO` | CHEBI Integrated Role Ontology | 13,279 | CHEBI provides a distinct role hierarchy. Chemicals in the structural hierarchy are connected via a 'has role' relation. CHIRO provides links from these roles to useful other classes in other ontol... |
| `OAE` | OAE: Ontology of Adverse Events | 10,595 | The Ontology of Adverse Eventsy (OAE) is a biomedical ontology in the domain of adverse events. OAE aims to standardize adverse event annotation, integrate various adverse event data, and support c... |
| `TXPO` | TOXic Process Ontology (TXPO) | 9,406 | Elucidating the mechanism of toxicity is crucial in drug safety evaluations. TOXic Process Ontology (TXPO) systematizes a wide variety of terms involving toxicity courses and processes. The first v... |
| `FIDEO` | Food Interactions with Drugs Evidence Ontology | 3,543 | Food-Drug interactions automatically extracted from scientific literature |
| `DIDEO` | Drug-drug Interaction and Drug-drug Interaction Evidence Ontology | 684 | The Potential Drug-drug Interaction and Potential Drug-drug Interaction Evidence Ontology |
| `PDRO` | The Prescription of Drugs Ontology | 310 | An ontology to describe entities related to prescription of drugs |
| `MPIO` | Minimum PDDI Information Ontology | 51 | An ontology of minimum information regarding potential drug-drug interaction information. |

## 🎯 HIGH — Protein / pathway / interaction (15 ontologies)

| Prefix | Title | # terms | Description |
|---|---|---:|---|
| `PR` | PRotein Ontology (PRO) | 365,852 | An ontological representation of protein-related entities |
| `CCO` | Cell Cycle Ontology | 277,685 | The Cell Cycle Ontology extends existing ontologies for cell cycle knowledge building a resource that integrates and manages knowledge about the cell cycle components and regulatory aspects. |
| `MRO` | MHC Restriction Ontology | 51,436 | The MHC Restriction Ontology is an application ontology capturing how Major Histocompatibility Complex (MHC) restriction is defined in experiments, spanning exact protein complexes, individual prot... |
| `EC` | Enzyme Nomenclature | 13,179 | The Enzyme Nomenclature (also known as the Enzyme Commission Code) is a species-agnostic controlled vocabulary for specific enzymes and an associated hierarchical classification into 7 main categor... |
| `OHPI` | OHPI: Ontology of Host-Pathogen Interactions | 11,291 | OHPI is a biomedical ontology in the area of host-pathogen interactions. OHPI is developed by following the OBO Foundry Principles (e.g., openness and collaboration). |
| `BTO` | The BRENDA Tissue Ontology (BTO) | 6,569 | A structured controlled vocabulary for the source of an enzyme comprising tissues, cell lines, cell types and cell cultures. |
| `—` | Unipathway | 4,672 | A manually curated resource for the representation and annotation of metabolic pathways |
| `MOP` | MOP | 3,686 | MOP is the molecular process ontology. It contains the molecular processes that underlie the name reaction ontology RXNO, for example cyclization, methylation and demethylation. |
| `wikipathways` | WikiPathways | 2,899 | WikiPathways is a database of biological pathways maintained by and for the scientific community. |
| `PW` | Pathway ontology | 2,760 | A controlled vocabulary for annotating gene products to pathways. |
| `OHMI` | OHMI: Ontology of Host-Microbiome Interactions | 1,724 | OHMI is a biomedical ontology that represents the entities and relations in the domain of host-microbiome interactions. |
| `MI` | Molecular Interactions Controlled Vocabulary | 1,656 | A structured controlled vocabulary for the annotation of experiments concerned with protein-protein interactions. |
| `RXNO` | RXNO | 1,019 | RXNO is the name reaction ontology. It contains more than 500 classes representing organic reactions such as the Diels–Alder cyclization. |
| `OMRSE` | Ontology for Modeling and Representation of Social Entities | 474 | The Ontology for Modeling and Representation of Social Entities (OMRSE) is an OBO Foundry ontology that represents the various entities that arise from human social interactions, such as social act... |
| `PHI` | PHI-base Ontology | 42 | Pathogen-Host Interaction database Ontology used by Ensembl |

## 🎯 HIGH — Cell / anatomy extras (developmental stage, cell line, etc.) (43 ontologies)

| Prefix | Title | # terms | Description |
|---|---|---:|---|
| `FMA` | Foundational Model of Anatomy Ontology (subset) | 104,721 | This is currently a slimmed down version of FMA |
| `CLO` | CLO: Cell Line Ontology | 43,327 | The Cell Line Ontology (CLO) is a community-based ontology of cell lines. The CLO is developed to unify publicly available cell line entry data from multiple sources to a standardized logically def... |
| `PCL` | Provisional Cell Ontology | 37,210 | Cell types that are provisionally defined by experimental techniques such as single cell transcriptomics rather than a straightforward & coherent set of properties. |
| `FBbt` | Drosophila gross anatomy | 29,751 | An ontology representing the gross anatomy of Drosophila melanogaster. |
| `HCAO` | Human Cell Atlas Ontology | 22,606 | Application ontology for human cell types, anatomy and development stages for the Human Cell Atlas. |
| `—` | Mouse gross anatomy and development, timed | 19,444 | A structured controlled vocabulary of stage-specific anatomical structures of the mouse (Mus). |
| `EMAPA` | Mouse Developmental Anatomy Ontology | 8,882 | An ontology for mouse anatomy covering embryonic development and postnatal stages. |
| `AISM` | Ontology for the Anatomy of the Insect SkeletoMuscular system | 8,614 | The ontology for the Anatomy of the Insect SkeletoMuscular system (AISM) contains terms used to describe the cuticle - as a single anatomical structure - and the skeletal muscle system, to be used ... |
| `WBbt` | C. elegans Gross Anatomy Ontology | 7,573 | Ontology about the gross anatomy of the C. elegans |
| `ECAO` | Echinoderm Anatomy and Development Ontology | 7,309 | None |
| `DHBA` | Developing Human Brain Atlas Ontology | 6,843 | An application ontology built by combining ontologised versions of the Allen Institute Developing Human Brain Atlas (DHBA) Structure Graph mapped to Uberon. |
| `HBA` | Human Brain Atlas Ontology | 5,172 | An application ontology built by combining ontologised versions of the Allen Institute Human Brain Atlas (HBA) Structure Graph mapped to Uberon. |
| `DMBA` | Developing Mouse Brain Atlas Ontology | 4,238 | An application ontology built by combining ontologised versions of the Allen Institute Developing Mouse Brain Atlas (DMBA) StructureGraph mapped to Uberon. |
| `MA` | Mouse adult gross anatomy | 3,390 | A structured controlled vocabulary of the adult anatomy of the mouse (Mus). |
| `—` | Teleost Anatomy Ontology | 3,372 | Multispecies fish anatomy ontology. Originally seeded from ZFA, but intended to cover terms relevant to other taxa |
| `ZFA` | Zebrafish anatomy and development ontology | 3,285 | A structured controlled vocabulary of the anatomy and development of the Zebrafish |
| `ZFS` | Zebrafish developmental stages ontology | 3,285 | Developmental stages of the Zebrafish |
| `MBA` | Mouse Brain Atlas Ontology | 3,097 | An application ontology built by combining ontologised versions of the Allen Institute Mouse Brain Atlas (MBA) StructureGraph mapped to Uberon. |
| `—` | Human developmental anatomy, abstract | 2,743 | A structured controlled vocabulary of stage-specific anatomical structures of the developing human. |
| `HAO` | Hymenoptera Anatomy Ontology | 2,597 | A structured controlled vocabulary of the anatomy of the Hymenoptera (bees, wasps, and ants) |
| `PO` | Plant Ontology | 2,029 | The Plant Ontology is a structured vocabulary and database resource that links plant anatomy, morphology and growth and development to plant genomics data. |
| `—` | Mosquito gross anatomy ontology | 1,864 | A structured controlled vocabulary of the anatomy of mosquitoes. |
| `XAO` | Xenopus Anatomy Ontology | 1,830 | XAO represents the anatomy and development of the African frogs Xenopus laevis and tropicalis. |
| `—` | Cephalopod Ontology | 1,574 | An anatomical and developmental ontology for cephalopods |
| `CLAO` | Collembola Anatomy Ontology | 1,516 | CLAO is an ontology of anatomical terms employed in morphological descriptions for the Class Collembola (Arthropoda: Hexapoda). |
| `PORO` | Porifera (sponge) ontology | 932 | An ontology describing the anatomical structures and characteristics of Porifera (sponges) |
| `PLANA` | Planarian Anatomy Ontology (PLANA) | 875 | PLANA, the PLANarian Anatomy Ontology, encompasses the anatomy of developmental stages and adult biotypes of Schmidtea mediterranea. |
| `CARO` | Common Anatomy Reference Ontology | 720 | The Common Anatomy Reference Ontology (CARO) is being developed to facilitate interoperability between existing anatomy ontologies for different species, and will provide a template for building ne... |
| `CTENO` | Ctenophore Ontology | 709 | An anatomical and developmental ontology for ctenophores (Comb Jellies) |
| `—` | Tick Anatomy Ontology | 628 | The anatomy of the Tick, <i>Families: Ixodidae, Argassidae</i> |
| `OPL` | Ontology for Parasite Lifecycle | 560 | The Ontology for Parasite Lifecycle (OPL) models the life cycle stage details of various parasites, including Trypanosoma sp., Leishmania major, and Plasmodium sp., etc. In addition to life cycle s... |
| `AEO` | Anatomical Entity Ontology | 430 | AEO is an ontology of anatomical structures that expands CARO, the Common Anatomy Reference Ontology |
| `AMPHX` | Amphioxus Development and Anatomy Ontology (AMPHX) | 403 | An ontology for the development and anatomy of Amphioxus (Branchiostoma lanceolatum). |
| `—` | Vertebrate Skeletal Anatomy Ontology- | 314 | Vertebrate skeletal anatomy ontology. |
| `CLYH` | Clytia hemisphaerica Development and Anatomy Ontology (CLYH) | 265 | Anatomy, development and life cycle stages - planula, polyp, medusa/jellyfish - of the cnidarian hydrozoan species, Clytia hemiphaerica. |
| `HsapDv` | Human Developmental Stages | 260 | Life cycle stages for Human |
| `MmusDv` | Mouse Developmental Stages | 178 | Life cycle stages for Mus Musculus |
| `BSPO` | Biological Spatial Ontology | 169 | An ontology for respresenting spatial concepts, anatomical axes, gradients, regions, planes, sides and surfaces. These concepts can be used at multiple biological scales and in a diversity of taxa,... |
| `DDANAT` | Dicty Anatomy Ontology (DDANAT) | 135 | A structured controlled vocabulary of anatomies of the slime-mold Dictyostelium discoideum. |
| `FAO` | Fungal gross anatomy | 122 | A structured controlled vocabulary for the anatomy of fungi. |
| `OlatDv` | Medaka Developmental Stages | 47 | Life cycle stages for Medaka |
| `CCF` | Human Reference Atlas Common Coordinate Framework Ontology | 41 | The Common Coordinate Framework (CCF) Ontology is an application ontology built to support the development of the Human Reference Atlas (HRA).  It unifies vocabulary for HRA construction and usage—... |
| `PdumDv` | Platynereis Developmental Stages | 21 | Life cycle stages for Platynereis dumerilii |

## 🎯 HIGH — Methods / techniques (12 ontologies)

| Prefix | Title | # terms | Description |
|---|---|---:|---|
| `GeXO` | Gene Expression Ontology | 166,254 | Gene Expression Ontology |
| `ReXO` | Regulation of Gene Expression Ontology | 158,239 | Regulation of Gene Expression |
| `ReTO` | Regulation of Transcription Ontology | 147,738 | Regulation of Transcription |
| `PRIDE` | Proteomics Identification Database Ontology | 97,576 | Proteomics Identification Database Ontology, terms describing proteomics data and experimental metadata |
| `OMIT` | Ontology for MIRNA Target | 90,928 | Ontology to establish data exchange standards and common data elements in the microRNA (miR) domain |
| `BAO` | BioAssay Ontology | 8,936 | The BioAssay Ontology (BAO) describes biological screening assays and their results including high-throughput screening (HTS) data for the purpose of categorizing assays and data analysis. BAO is a... |
| `MAXO` | Medical Action Ontology | 7,120 | An ontology to represent medically relevant actions, procedures, therapies, interventions, and recommendations. |
| `MS` | Mass spectrometry ontology | 4,097 | A structured controlled vocabulary for the annotation of experiments concerned with proteomics mass spectrometry. |
| `UNIMOD` | Unimod protein modification database for mass spectrometry | 1,523 | Unimod is a community supported, comprehensive database of protein modifications for mass spectrometry applications. That is, accurate and verifiable values, derived from elemental compositions, fo... |
| `XLMOD` | HUPO-PSI cross-linking and derivatization reagents controlled vocabulary | 1,106 | A structured controlled vocabulary for cross-linking reagents used with proteomics mass spectrometry. |
| `MFOEM` | Emotion Ontology | 624 | An ontology of affective phenomena such as emotions, moods, appraisals and subjective feelings. |
| `ORNASEQ` | Ontology for RNA sequencing (ORNASEQ) | 163 | An application ontology designed to annotate next-generation sequencing experiments performed on RNA. |

## Medium — Genomic features / variation (3 ontologies)

| Prefix | Title | # terms | Description |
|---|---|---:|---|
| `ICEO` | ICEO: Ontology of Integrative and Conjugative Elements | 30,914 | A biological ontology to standardize and integrate Integrative and Conjugative Element (ICE) information and to support computer-assisted reasoning. |
| `ARO` | Antibiotic Resistance Ontology | 8,582 | Antibiotic resistance genes and mutations |
| `VariO` | Variation Ontology | 618 | Variation Ontology, VariO, is an ontology for standardized, systematic description of effects, consequences and mechanisms of variations. |

## Medium — Modifications / PTMs (2 ontologies)

| Prefix | Title | # terms | Description |
|---|---|---:|---|
| `GNO` | Glycan Naming Ontology | 198,340 | An ontology for glycans based on GlyTouCan, but organized by subsumption. |
| `MOD` | Protein modification | 2,116 | PSI-MOD is an ontology consisting of terms that describe protein chemical modifications |

## Medium — Statistics / measurement (4 ontologies)

| Prefix | Title | # terms | Description |
|---|---|---:|---|
| `CMO` | Clinical measurement ontology | 4,144 | Morphological and physiological measurement records generated from clinical and model organism research and health programs. |
| `XCO` | Experimental condition ontology | 1,834 | Conditions under which physiological and morphological measurements are made both in the clinic and in studies involving humans or model organisms. |
| `UO` | Units of measurement ontology | 574 | Metrical units for use in conjunction with PATO |
| `DUO` | Data Use Ontology | 311 | DUO is an ontology which represent data use conditions. |

## Medium — Environment / ecology (holobiont context) (16 ontologies)

| Prefix | Title | # terms | Description |
|---|---|---:|---|
| `mesh` | Medical Subject Headings | 355,028 | MeSH (Medical Subject Headings) is the National Library of Medicine's controlled vocabulary thesaurus. It consists of sets of terms naming descriptors in a hierarchical structure that permits searc... |
| `FOODON` | Food Ontology | 39,682 | A broadly scoped ontology representing entities which bear a “food role”. It encompasses materials in natural ecosystems and agriculture that are consumed by humans and domesticated animals. This i... |
| `VO` | Vaccine Ontology | 14,804 | The Vaccine Ontology (VO) is a community-based biomedical ontology in the domain of vaccines and vaccination. VO aims to standardize vaccine types and annotations, integrate various vaccine data, a... |
| `ECTO` | Environment Exposure Ontology | 11,770 | ECTO describes exposures to experimental treatments of plants and model organisms (e.g. exposures to modification of diet, lighting levels, temperature); exposures of humans or any other organisms ... |
| `RBO` | Radiation Biology Ontology | 9,209 | RBO is an ontology for the effects of radiation on biota in terrestrial and space environments. |
| `ENVO` | The Environment Ontology | 6,906 | ENVO is an ontology which represents knowledge about environments,environmental processes, ecosystems, habitats, and related entities |
| `EUPATH` | VEuPathDB Ontology | 5,109 | The VEuPathDB ontology is an application ontology developed to encode our understanding of what data is about in the public resources developed and maintained by the Eukaryotic Pathogen, Host and V... |
| `OPMI` | OPMI: Ontology of Precision Medicine and Investigation | 3,153 | OPMI is a biomedical ontology in the area of precision medicine and its related investigations. It is community-driven and developed by following the OBO Foundry ontology development principles. |
| `OBIB` | Ontology for BIoBanking (OBIB) | 1,803 | The Ontology for Biobanking (OBIB) is an ontology for the annotation and modeling of the activities, contents, and administration of a biobank. Biobanks are facilities that store specimens, such as... |
| `OCCO` | Occupation Ontology (OccO) | 1,801 | This Occupaton Ontology (OccO) is an ontology in the domain of human occupaitons. Current OccO version is developed based on the US Bureau of Labor Statistics Standard Occupation Classification (SO... |
| `OVAE` | OVAE: Ontology of Vaccine Adverse Events | 1,554 | OVAE is a biomedical ontology in the area of vaccine adverse events. OVAE is an extension of the community-based Ontology of Adverse Events (OAE). |
| `OBCS` | OBCS: Ontology of Biological and Clinical Statistics | 1,210 | OBCS stands for the Ontology of Biological and Clinical Statistics. OBCS is an ontology in the domain of biological and clinical statistics. It is aligned with the Basic Formal Ontology (BFO) and t... |
| `FBbi` | Biological Imaging Methods Ontology | 677 | An ontology for the description of sample preparation, visualization, and imaging methods used in biomedical research. |
| `ONTOAVIDA` | OntoAvida: ontology for Avida digital evolution platform. | 674 | The Ontology for Avida (OntoAvida) project aims to develop an integrated vocabulary for the description of the most widely used computational approach for performing experimental evolution using di... |
| `BCO` | Biological Collections Ontology | 253 | An ontology to support the interoperability of biodiversity data, including data on museum collections, environmental/metagenomic samples, and ecological surveys. |
| `COB` | Core Ontology for Biology and Biomedicine | 69 | COB brings together key terms from a wide range of OBO projects to improve interoperability. |

## Medium — Sample / biobank metadata (6 ontologies)

| Prefix | Title | # terms | Description |
|---|---|---:|---|
| `VBO` | Vertebrate Breed Ontology | 20,473 | Vertebrate Breed Ontology is an ontology created to serve as a single computable resource for vertebrate breed names. |
| `—` | Human Reference Atlas (HRA) | 10,479 | The Human Reference Atlas (HRA) is a comprehensive, high-resolution, three-dimensional atlas of all the cells in the healthy human body. The Human Reference Atlas provides standard terminologies an... |
| `NGBO` | Next generation biobanking ontology(NGBO). | 1,515 | Next Generation Biobanking Ontology (NGBO) is an open application ontology representing contextual data about omics digital assets in biobank. The ontology focuses on capturing the information abou... |
| `AfPO` | African Population Ontology | 473 | AfPO is an ontology that can be used in the study of diverse populations across Africa. It brings together publicly available demographic, anthropological and genetic data relating to African peopl... |
| `—` | Ontologized MIABIS | 427 | An ontological version of MIABIS (Minimum Information About BIobank data Sharing) |
| `GECKO` | Genomics Cohorts Knowledge Ontology | 154 | An ontology to represent genomics cohort attributes. |

## Medium — Upper-level / shared semantics (11 ontologies)

| Prefix | Title | # terms | Description |
|---|---|---:|---|
| `OGG` | OGG: Ontology of Genes and Genomes | 69,688 | OGG is a biological ontology in the area of genes and genomes. OGG uses the Basic Formal Ontology (BFO) as its upper level ontology. This OGG document contains the genes and genomes of a list of se... |
| `OM` | Ontology of units of Measure (OM) | 2,307 | The Ontology of units of Measure (OM) 2.0 models concepts and relations important to scientific research. It has a strong focus on units, quantities, measures, and dimensions. |
| `hgnc.genegroup` | HGNC Gene Group | 1,993 | The HGNC (HUGO Gene Nomenclature Committee) provides an approved gene name and symbol (short-form abbreviation) for each known human gene. All approved symbols are stored in the HGNC database, and ... |
| `PROCO` | PROCO: PROcess Chemistry Ontology | 922 | PROCO (PROcess Chemistry Ontology) is a formal ontology that aims to standardly represent entities and relations among entities in the domain of process chemistry. |
| `EXMO` | EXMO | 449 | A core reference ontology built upon BFO about exercise medicine. This ontology contains the related terms for healthy people, people with chronic conditions and people living with diability to exe... |
| `CDAO` | Comparative Data Analysis Ontology | 164 | a formalization of concepts and relations relevant to evolutionary comparative analysis |
| `MAMO` | Mathematical Modelling Ontology | 106 | The Mathematical Modelling Ontology (MAMO) is a classification of the types of mathematical models used mostly in the life sciences, their variables, relationships and other relevant features. |
| `semapv` | Semantic Mapping Vocabulary | 47 | The Semantic Mapping Vocabulary (SEMAPV) captures concepts and relations related to semantic matching, including types of matching (e.g. lexical, manually curated), preprocessing methods and specia... |
| `PROBONTO` | Probability Distribution Ontology | 39 | ProbOnto, is an ontology-based knowledge base of probability distributions, featuring more than eighty uni- and multivariate distributions with their defining functions, characteristics, relationsh... |
| `BFO` | Basic Formal Ontology | 35 | The upper level ontology upon which OBO Foundry ontologies are built. |
| `OMO` | OBO Metadata Ontology | 27 | An ontology specifies terms that are used to annotate ontology terms for all OBO ontologies. The ontology was developed as part of Information Artifact Ontology (IAO). |

## Low — Food / nutrition (diet studies only) (8 ontologies)

| Prefix | Title | # terms | Description |
|---|---|---:|---|
| `TFO` | TransformON | 10,288 | A domain ontology on food and bioproduct processing and characterization that is a specialization of the PO2 ontology (Process and Observation Ontology). |
| `GENEPIO` | Genomic Epidemiology Ontology | 9,159 | The Genomic Epidemiology Ontology (GenEpiO) covers vocabulary necessary to identify, document and research foodborne pathogens and associated outbreaks. |
| `LIFESTYLEFACTORS` | Lifestyle Factors Ontology | 5,598 | In this ontology, we defined nine distinct categories that encompass various aspects of lifestyle factors (LSFs), creating a unified Lifestyle Factors Ontology (LSFO) across these categories. Withi... |
| `ONS` | Ontology for Nutritional Studies | 4,834 | An ontology for description of concepts in the nutritional studies domain. |
| `ONE` | Ontology for Nutritional Epidemiology | 4,100 | An ontology to standardize research output of nutritional epidemiologic studies. |
| `CDNO` | Compositional Dietary Nutrition Ontology | 2,109 | The CDNO provides structured terminologies to describe nutritional attributes of material entities that contribute to human diet. These terms are intended primarily to be associated with datasets t... |
| `FOBI` | Food-Biomarker Ontology | 1,799 | FOBI (Food-Biomarker Ontology) is an ontology to represent food intake data and associate it with metabolomic data |
| `HSO` | Health Surveillance Ontology | 448 | The health Surveillance Ontology (HSO) focuses on "surveillance system level data", that is, data outputs from surveillance activities, such as number of samples collected, cases observed, etc. It ... |

## Low — Plant ontologies (only if multi-omics expands) (4 ontologies)

| Prefix | Title | # terms | Description |
|---|---|---:|---|
| `PSO` | Plant Stress Ontology | 4,899 | The Plant Stress Ontology describes biotic and abiotic stresses that a plant may encounter. |
| `AGRO` | Agronomy Ontology | 4,162 | AgrO is an ontlogy for representing agronomic practices, techniques, variables and related entities |
| `PECO` | Plant Experimental Conditions Ontology | 3,657 | A structured, controlled vocabulary which describes the treatments, growing conditions, and/or study types used in plant biology experiments. |
| `PPO` | Plant Phenology Ontology | 576 | An ontology for describing the phenology of individual plants and populations of plants, and for integrating plant phenological data across sources and scales. |

## Other / Uncategorised (74 ontologies)

| Prefix | Title | # terms | Description |
|---|---|---:|---|
| `GAZ` | Gazetteer | 668,838 | A gazetteer constructed on ontological principles. The countries are actively maintained. |
| `SNOMED` | SNOMED CT (International Edition) | 376,452 | SNOMED CT or SNOMED Clinical Terms is a systematically organized computer processable collection of medical terms providing codes, terms, synonyms and definitions used in clinical documentation and... |
| `—` | (untitled) | 195,343 |  |
| `hgnc` | HUGO Gene Nomenclature Committee | 163,229 | The HGNC (HUGO Gene Nomenclature Committee) provides an approved gene name and symbol (short-form abbreviation) for each known human gene.  All approved symbols are stored in the HGNC database, and... |
| `VTO` | Vertebrate Taxonomy Ontology | 107,254 | Comprehensive hierarchy of extinct and extant vertebrate taxa. |
| `TTO` | Teleost taxonomy ontology | 38,705 | An ontology covering the taxonomy of teleosts (bony fish) |
| `PREFER` | PREFER Ontology | 6,906 | An ontology for PREcision FERmentation. |
| `PBPKO` | Physiologically-Based Pharmacokinetic Ontology | 5,643 | An ontology for modeling and understanding of pharmacokinetics in biological systems. |
| `RS` | Rat Strain Ontology | 5,516 | Ontology of rat strains |
| `—` | Mosquito insecticide resistance | 4,409 | Application ontology for entities related to insecticide resistance in mosquitos |
| `SLSO` | Space Life Sciences Ontology | 3,978 | The Space Life Sciences Ontology is an application ontology and is intended to support the operation of NASA's Life Sciences Data Archive and other systems that contain space life science research ... |
| `AFO` | Allotrope Merged Ontology Suite | 3,516 | Allotrope Merged Ontology Suite |
| `MCO` | Microbial Conditions Ontology | 3,383 | Microbial Conditions Ontology is an ontology... |
| `ECO` | Evidence & Conclusion Ontology (ECO) | 3,358 | The Evidence & Conclusion Ontology (ECO) describes types of scientific evidence within the biological research domain that arise from laboratory experiments, computational methods, literature curat... |
| `CHMO` | Chemical Methods Ontology | 3,234 | CHMO, the chemical methods ontology, describes methods used to |
| `—` | Malaria Ontology | 3,158 | An application ontology to cover all aspects of malaria as well as the intervention attempts to control it. |
| `NCRO` | Non-Coding RNA Ontology | 3,078 | An ontology for non-coding RNA, both of biological origin, and engineered. |
| `FBcv` | FlyBase Controlled Vocabulary | 2,338 | A structured controlled vocabulary used for various aspects of annotation by FlyBase. |
| `—` | Beta Cell Genomics Ontology | 2,270 | An application ontology built for beta cell genomics studies. |
| `CHEMINF` | chemical information ontology (cheminf) - information entities about chemical entities | 2,141 | Includes terms for the descriptors commonly used in cheminformatics software applications and the algorithms which generate them. |
| `SWO` | Software ontology | 1,971 | The Software Ontology (SWO) is a resource for describing software tools, their types, tasks, versions, provenance and associated data. It contains detailed information on licensing and formats as w... |
| `—` | (untitled) | 1,824 |  |
| `ONTONEO` | Obstetric and Neonatal Ontology | 1,420 | The Obstetric and Neonatal Ontology is a structured controlled vocabulary to provide a representation of the data from electronic health records (EHRs) involved in the care of the pregnant woman, a... |
| `PSDO` | Performance Summary Display Ontology | 1,383 | Performance Summary Display Ontology (PSDO) is an application ontology about motivating information in performance summaries and the visual and textual entities that are used to communicate about p... |
| `ADDICTO` | Addiction Ontology | 1,363 | The Addiction Ontology (AddictO) is an ontology being developed all aspects of addiction research. |
| `gscmixs` | Minimum Information about any (x) Sequence (MIxS) standard | 1,335 | An OWL representation of the Minimum Information for any (x) Standard (MIxS), managed by the Genomic Standards Consortium. This file is released by the Genomic Standards Consortium (GSC; https://ww... |
| `HANCESTRO` | Human Ancestry Ontology | 1,310 | Human ancestry ontology for the NHGRI GWAS Catalog |
| `MSIO` | Metabolomics Standards Initiative Ontology (MSIO) | 1,294 | an application ontology for supporting description and annotation of mass-spectrometry and nmr-spectroscopy based metabolomics experiments and fluxomics studies. |
| `ICO` | Informed Consent Ontology (ICO) | 1,197 | The Informed Consent Ontology (ICO) is an ontology for the informed consent and informed consent process in the medical field. |
| `—` | (untitled) | 1,014 |  |
| `—` | (untitled) | 1,014 |  |
| `MCRO` | Model Card Report Ontology | 969 | An ontology representing the structure of model card reports - reports that describe basic characteristics of machine learning models for the public and consumers. |
| `OOSTT` | Ontology of Organizational Structures of Trauma centers and Trauma Systems | 844 | The Ontology of Organizational Structures of Trauma centers and Trauma systems (OOSTT) is a representation of the components of trauma centers and trauma systems coded in Web Ontology Language (OWL2). |
| `WBls` | C. elegans Development Ontology | 796 | Ontology about the development and life stages of the C. elegans |
| `SBO` | Systems Biology Ontology | 695 | Terms commonly used in Systems Biology, and in particular in computational modeling. |
| `—` | (untitled) | 663 |  |
| `BMONT` | Biomarker Ontology (BMONT) | 640 | An application ontology that represents comprehensive knowledge involving a variety of fields of medical and biological aspects. |
| `HTN` | Hypertension Ontology For Clinical Data | 603 | An ontology for representing clinical data about hypertension, intended to support classification of patients according to various diagnostic guidelines |
| `COVOC` | CoVoc Coronavirus Vocabulary | 563 | The COVID-19 Vocabulary (COVoc) is an ontology containing terms related to the research of the COVID-19 pandemic. This includes host organisms, pathogenicity, gene and gene products, barrier gestur... |
| `T4FS` | terms4FAIRskills | 563 | A terminology for the skills necessary to make data FAIR and to keep it FAIR. |
| `KISAO` | Kinetic Simulation Algorithm Ontology | 525 | A classification of algorithms for simulating biology, their parameters, and their outputs |
| `CAO` | Chemical Analysis Metadata Platform (ChAMP) Ontology version 0.2 | 445 | The Chemical Analysis Ontology (CAO) is an ontology developed as part of the Chemical Analaysis Metadata Project (ChAMP) as a resource to semantically annotate standards developed using the ChAMP p... |
| `SRAO` | FAIRsharing Subject Ontology (SRAO) | 445 | The FAIRsharing Subject Ontology (SRAO) is an application ontology for the categorization of research disciplines across all research domains, from the humanities to the natural sciences. It utiliz... |
| `MF` | Mental Functioning Ontology | 400 | The Mental Functioning Ontology is an overarching ontology for all aspects of mental functioning. |
| `REPR` | REPRODUCE-ME Ontology | 389 | The REPRODUCE-ME ontology is an extension of the PROV-O and the P-Plan ontology to describe a complete path of a scientific experiment. It expresses the REPRODUCE-ME Data Model using the OWL2 Web O... |
| `ENSGLOSS` | Ensembl Glossary | 377 | The Ensembl glossary lists the terms, data types and file types that are used in Ensembl and describes how they are used. |
| `OARCS` | Ontology of Arthropod Circulatory Systems | 308 | OArCS is an ontology describing the Arthropod ciruclatory system. |
| `NOMEN` | NOMEN - A nomenclatural ontology for biological names | 264 | NOMEN is a nomenclatural ontology for biological names (not concepts).  It encodes the goverened rules of nomenclature. |
| `FBdv` | Drosophila development | 250 | A structured controlled vocabulary of the development of Drosophila melanogaster. |
| `—` | (untitled) | 240 |  |
| `MIAPA` | Minimum Information for A Phylogenetic Analysis (MIAPA) Ontology | 238 | The MIAPA ontology is intended to be an application ontology for the purpose of semantic annotation of phylogenetic data according to the requirements and recommendations of the Minimum Information... |
| `gmho` | GALENOS Mental Health Ontology | 212 | The GALENOS Mental Health Ontology (GMHO) is being developed to represent the full range of concepts relevant to mental health. It forms a core component of the Global Alliance for Living Evidence ... |
| `LABO` | clinical LABoratory Ontology | 169 | LABO is an ontology of informational entities formalizing clinical laboratory tests prescriptions and reporting documents. |
| `EMI` | The Earth Metabolome Initiative (EMI) ontology | 161 | The Earth Metabolome Initiative (EMI) ontology is designed to capture and structure the metabolic content of all currently known species on our planet. |
| `ZECO` | Zebrafish Experimental Conditions Ontology | 161 | Ontology of Zebrafish Experimental Conditions |
| `CRO` | Contributor Role Ontology | 145 | A classification of the diverse roles performed in the work leading to a published research output in the sciences. Its purpose to provide transparency in contributions to scholarly published work,... |
| `SEPIO` | Ontology for Scientific Evidence and Provenance Information | 125 | An ontology for representing the provenance of scientific claims and the evidence that supports them. |
| `GEO` | Geographical Entity Ontology | 108 | An ontology of geographical entities |
| `SHARELOC` | ShareLoc | 95 | ShareLoc defines terms to annotate data sets from single molecule localization microscopy, including but not limited to: imaging technique, biological structures or molecules of interest, cell type... |
| `—` | (untitled) | 79 |  |
| `—` | (untitled) | 77 |  |
| `TAXRANK` | Taxonomic rank vocabulary | 77 | A vocabulary of taxonomic ranks (species, family, phylum, etc) |
| `HOM` | Homology Ontology | 66 | This ontology represents concepts related to homology, as well as other concepts used to describe similarity and non-homology. |
| `CIO` | Confidence Information Ontology | 44 | An ontology to capture confidence information about annotations. |
| `ORTH` | Orthology Ontology | 34 | The need of a common ontology for describing orthology information in biological research communities has led to the creation of the Orthology Ontology (ORTH). ORTH ontology is designed to describe... |
| `—` | The OWL 2 Schema vocabulary (OWL 2) | 32 |  |
| `PROV` | Provenance Ontology (PROV-O) | 31 |  |
| `—` | (untitled) | 22 |  |
| `credit` | Contributor Roles Taxonomy | 16 | CRediT (Contributor Roles Taxonomy) is high-level taxonomy, including 14 roles, that can be used to represent the roles typically played by contributors to scientific scholarly output. The roles de... |
| `ror` | Research Organization Registry | 10 | ROR (Research Organization Registry) is a global, community-led registry of open persistent identifiers for research organizations. ROR is jointly operated by California Digital Library, Crossref, ... |
| `—` | (untitled) | 6 |  |
| `—` | The RDF Schema vocabulary (RDFS) | 6 |  |
| `—` | (untitled) | 6 |  |
| `—` | (untitled) | 0 |  |


**Total: 280** · already covered: **16** · high-priority candidates (top 6 buckets): **136**
