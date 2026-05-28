"""Hardcoded namespace URI constants for transcriptomics ontologies.

We do not depend on any ontology Python package at runtime. The expansion
URIs are curated here and remain stable across releases. Update with care:
changing a URI here changes how every claim CURIE is expanded downstream.

The module exports one big dictionary ``TRANSCRIPTOMICS_PREFIXES`` mapping
CURIE prefixes to their URI expansions. The dict is what
``BulkTranscriptomicsAdapter.ontology_prefixes`` and
``SingleCellTranscriptomicsAdapter.ontology_prefixes`` expose.

Sections:

* Biomedical / experimental ontologies (EFO, OBI, STATO, ...).
* Concept / process / artifact ontologies (INO, NCIT, SIO, RO, IAO).
* Gene / transcript identifier authorities (HGNC, MGI, FlyBase, ...).
* Repository / accession authorities (GEO, SRA, ENA, ArrayExpress, ...).
"""
from __future__ import annotations

from typing import Final

# --- Biomedical / experimental ontologies ----------------------------------

EFO: Final[str] = "http://www.ebi.ac.uk/efo/EFO_"
OBI: Final[str] = "http://purl.obolibrary.org/obo/OBI_"
STATO: Final[str] = "http://purl.obolibrary.org/obo/STATO_"
NCBITAXON: Final[str] = "http://purl.obolibrary.org/obo/NCBITaxon_"
CL: Final[str] = "http://purl.obolibrary.org/obo/CL_"
UBERON: Final[str] = "http://purl.obolibrary.org/obo/UBERON_"
GO: Final[str] = "http://purl.obolibrary.org/obo/GO_"
SO: Final[str] = "http://purl.obolibrary.org/obo/SO_"
MONDO: Final[str] = "http://purl.obolibrary.org/obo/MONDO_"
PATO: Final[str] = "http://purl.obolibrary.org/obo/PATO_"
EDAM: Final[str] = "http://edamontology.org/"
# bio.tools registry: canonical home for bioinformatics software (STAR,
# Seurat, DESeq2, etc). Resolves to https://bio.tools/<id>.
BIOTOOLS: Final[str] = "https://bio.tools/"

# --- Concept / process / artifact ontologies -------------------------------
# Added in v0.3.0 to align with the TranscriptoKG modelling vocabulary,
# which uses NCIT for analysis concepts, INO for regulation processes, and
# SIO / IAO / RO for predicate-level metadata.

INO: Final[str] = "http://purl.obolibrary.org/obo/INO_"
NCIT: Final[str] = "http://purl.obolibrary.org/obo/NCIT_"
SIO: Final[str] = "http://semanticscience.org/resource/SIO_"
RO: Final[str] = "http://purl.obolibrary.org/obo/RO_"
IAO: Final[str] = "http://purl.obolibrary.org/obo/IAO_"

# --- Gene / transcript identifier authorities -----------------------------

HGNC: Final[str] = "https://identifiers.org/hgnc:"
MGI: Final[str] = "https://identifiers.org/MGI:"
FLYBASE: Final[str] = "https://identifiers.org/flybase:"
WORMBASE: Final[str] = "https://identifiers.org/wormbase:"
ZFIN: Final[str] = "https://identifiers.org/zfin:"
ENSEMBL: Final[str] = "https://identifiers.org/ensembl:"
ENTREZ_GENE: Final[str] = "https://identifiers.org/ncbigene:"
REFSEQ: Final[str] = "https://identifiers.org/refseq:"
UNIPROT: Final[str] = "https://identifiers.org/uniprot:"

# --- Repositories / accessions --------------------------------------------

GEO: Final[str] = "https://identifiers.org/geo:"
SRA: Final[str] = "https://identifiers.org/insdc.sra:"
BIOPROJECT: Final[str] = "https://identifiers.org/bioproject:"
ENA: Final[str] = "https://identifiers.org/ena.embl:"
ARRAYEXPRESS: Final[str] = "https://identifiers.org/arrayexpress:"
ENCODE: Final[str] = "https://identifiers.org/encode:"
DBGAP: Final[str] = "https://identifiers.org/dbgap:"
BIOSTUDIES: Final[str] = "https://identifiers.org/biostudies:"


# Combined CURIE-prefix dictionary exposed via the adapter classes'
# ``ontology_prefixes`` attribute. Consumers use it when rendering claims
# to RDF so a CURIE such as ``HGNC:12783`` expands to the full URI.
TRANSCRIPTOMICS_PREFIXES: Final[dict[str, str]] = {
    # Biomedical / experimental.
    "EFO": EFO,
    "OBI": OBI,
    "STATO": STATO,
    "NCBITaxon": NCBITAXON,
    "CL": CL,
    "UBERON": UBERON,
    "GO": GO,
    "SO": SO,
    "MONDO": MONDO,
    "PATO": PATO,
    "EDAM": EDAM,
    "biotools": BIOTOOLS,
    # Concept / process / artifact axes.
    "INO": INO,
    "NCIT": NCIT,
    "SIO": SIO,
    "RO": RO,
    "IAO": IAO,
    # Gene / transcript identifiers.
    "HGNC": HGNC,
    "MGI": MGI,
    "FlyBase": FLYBASE,
    "WormBase": WORMBASE,
    "ZFIN": ZFIN,
    "ENSEMBL": ENSEMBL,
    "NCBIGene": ENTREZ_GENE,
    "RefSeq": REFSEQ,
    "UniProt": UNIPROT,
    # Repositories.
    "GEO": GEO,
    "SRA": SRA,
    "BioProject": BIOPROJECT,
    "ENA": ENA,
    "ArrayExpress": ARRAYEXPRESS,
    "ENCODE": ENCODE,
    "dbGaP": DBGAP,
    "BioStudies": BIOSTUDIES,
}
