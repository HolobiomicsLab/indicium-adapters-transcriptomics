"""Tests for the v0.7.0 metadata parsers.

Each parser has a synthetic fixture that mirrors the real format
closely enough to exercise the field paths the parser cares about,
plus one end-to-end test that pipes the parser output through
``BulkTranscriptomicsAdapter.enrich_claim`` to prove the claim dicts
land in a valid enrichable shape.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from indicium_adapters_transcriptomics import BulkTranscriptomicsAdapter
from indicium_adapters_transcriptomics.parsers import (
    parse_biostudies_json,
    parse_ena_xml,
    parse_geo_soft,
    parse_isa_tab,
)

# --- GEO SOFT --------------------------------------------------------------


GSE168938_SOFT = """\
^SERIES = GSE168938
!Series_geo_accession = GSE168938
!Series_title = Warner et al. 2025 regeneration time course
!Series_pubmed_id = 12345678
!Series_sample_id = GSM4567890
!Series_sample_id = GSM4567891

^SAMPLE = GSM4567890
!Sample_title = oral_pole_0hpa
!Sample_organism_ch1 = Nematostella vectensis
!Sample_instrument_model = Illumina NovaSeq 6000

^SAMPLE = GSM4567891
!Sample_title = oral_pole_24hpa
!Sample_organism_ch1 = Nematostella vectensis
!Sample_instrument_model = Illumina NovaSeq 6000
"""


def test_parse_geo_soft_emits_consolidates_dataset_claims():
    claims = parse_geo_soft(GSE168938_SOFT)
    subjects = {c["subject"] for c in claims}
    assert subjects == {"GSE168938"}
    objs = {(c["qualifier"], c["object"]) for c in claims}
    assert ("consolidates_dataset", "GSM4567890") in objs
    assert ("consolidates_dataset", "GSM4567891") in objs


def test_parse_geo_soft_emits_sequenced_on_for_unique_instruments():
    claims = parse_geo_soft(GSE168938_SOFT)
    seq_objs = [c["object"] for c in claims if c["qualifier"] == "sequenced_on"]
    # Two samples list the same instrument — should dedupe.
    assert seq_objs == ["Illumina NovaSeq 6000"]


def test_parse_geo_soft_emits_uses_database_pubmed():
    claims = parse_geo_soft(GSE168938_SOFT)
    pm_objs = [c["object"] for c in claims if c["qualifier"] == "uses_database"]
    assert pm_objs == ["PubMed:12345678"]


def test_parse_geo_soft_attaches_species_hint_from_organism():
    claims = parse_geo_soft(GSE168938_SOFT)
    assert all(c.get("species") == "nematostella" for c in claims)


def test_parse_geo_soft_end_to_end_through_bulk_adapter():
    """The parser output must be directly consumable by enrich_claim."""
    a = BulkTranscriptomicsAdapter()
    claims = parse_geo_soft(GSE168938_SOFT)
    enriched_curies: dict[tuple[str, str], str] = {}
    for c in claims:
        out = a.enrich_claim(c)
        for slot in ("subject", "object"):
            curie = (out.get("ontology_terms") or {}).get(slot)
            if curie:
                enriched_curies[(c["qualifier"], slot)] = curie
    # Sample GSE/GSM accessions should resolve via the accession resolver.
    assert enriched_curies.get(("consolidates_dataset", "subject")) == "GEO:GSE168938"
    assert enriched_curies.get(("consolidates_dataset", "object")) == "GEO:GSM4567891"
    # Illumina NovaSeq 6000 resolves via the corpus index (substring match,
    # platform axis).
    assert enriched_curies.get(("sequenced_on", "object")) == "EFO:0008637"


# --- ENA / SRA XML ---------------------------------------------------------


ENA_XML = """<?xml version="1.0" encoding="UTF-8"?>
<EXPERIMENT_SET>
  <EXPERIMENT accession="ERX1234567">
    <STUDY_REF accession="PRJEB12345"/>
    <PLATFORM>
      <ILLUMINA>
        <INSTRUMENT_MODEL>Illumina NovaSeq 6000</INSTRUMENT_MODEL>
      </ILLUMINA>
    </PLATFORM>
  </EXPERIMENT>
  <EXPERIMENT accession="ERX1234568">
    <STUDY_REF accession="PRJEB12345"/>
    <PLATFORM>
      <OXFORD_NANOPORE>
        <INSTRUMENT_MODEL>PromethION</INSTRUMENT_MODEL>
      </OXFORD_NANOPORE>
    </PLATFORM>
  </EXPERIMENT>
</EXPERIMENT_SET>
"""


def test_parse_ena_xml_emits_consolidates_and_sequenced():
    claims = parse_ena_xml(ENA_XML)
    by_q = {(c["qualifier"], c["object"]) for c in claims}
    assert ("consolidates_dataset", "ERX1234567") in by_q
    assert ("consolidates_dataset", "ERX1234568") in by_q
    assert ("sequenced_on", "Illumina NovaSeq 6000") in by_q
    assert ("sequenced_on", "PromethION") in by_q


# --- BioStudies JSON -------------------------------------------------------


BIOSTUDIES_JSON = {
    "accno": "S-BIAD123",
    "section": {
        "type": "Study",
        "attributes": [
            {"name": "Title", "value": "Drosophila transcriptomics"},
            {"name": "Organism", "value": "Drosophila melanogaster"},
            {"name": "Pubmed", "value": "99887766"},
            {"name": "Instrument", "value": "Illumina HiSeq 2000"},
        ],
    },
}


def test_parse_biostudies_emits_pubmed_and_sequenced():
    claims = parse_biostudies_json(BIOSTUDIES_JSON)
    assert all(c["subject"] == "S-BIAD123" for c in claims)
    qual_obj = {(c["qualifier"], c["object"]) for c in claims}
    assert ("uses_database", "PubMed:99887766") in qual_obj
    assert ("sequenced_on", "Illumina HiSeq 2000") in qual_obj


def test_parse_biostudies_attaches_drosophila_species_hint():
    claims = parse_biostudies_json(BIOSTUDIES_JSON)
    assert all(c.get("species") == "drosophila" for c in claims)


def test_parse_biostudies_accepts_raw_json_string():
    claims = parse_biostudies_json(json.dumps(BIOSTUDIES_JSON))
    assert len(claims) >= 1


# --- ISA-Tab ---------------------------------------------------------------


ISATAB_INVESTIGATION = """\
ONTOLOGY SOURCE REFERENCE
Term Source Name\tEFO
Term Source File\thttp://www.ebi.ac.uk/efo
INVESTIGATION
Investigation Identifier\tINV001
STUDY
Study Identifier\tSTUDY-42
Study Title\tA mouse study
Study PubMed ID\t11223344
"""

ISATAB_STUDY = (
    "Source Name\tSample Name\tCharacteristics[organism]\n"
    "src1\tsamp1\tMus musculus\n"
    "src1\tsamp2\tMus musculus\n"
)

ISATAB_ASSAY = (
    "Sample Name\tProtocol REF\tParameter Value[instrument model]\n"
    "samp1\textraction\tIllumina NovaSeq 6000\n"
    "samp2\textraction\tIllumina NovaSeq 6000\n"
)


@pytest.fixture
def isa_tab_dir(tmp_path: Path) -> Path:
    (tmp_path / "i_investigation.txt").write_text(
        ISATAB_INVESTIGATION, encoding="utf-8"
    )
    (tmp_path / "s_study.txt").write_text(ISATAB_STUDY, encoding="utf-8")
    (tmp_path / "a_assay.txt").write_text(ISATAB_ASSAY, encoding="utf-8")
    return tmp_path


def test_parse_isa_tab_emits_pubmed_and_sequenced(isa_tab_dir: Path):
    claims = parse_isa_tab(isa_tab_dir)
    assert all(c["subject"] == "STUDY-42" for c in claims)
    qual_obj = {(c["qualifier"], c["object"]) for c in claims}
    assert ("uses_database", "PubMed:11223344") in qual_obj
    assert ("sequenced_on", "Illumina NovaSeq 6000") in qual_obj


def test_parse_isa_tab_detects_mouse_species(isa_tab_dir: Path):
    claims = parse_isa_tab(isa_tab_dir)
    assert claims, "expected at least one claim"
    assert all(c.get("species") == "mouse" for c in claims)


def test_parse_isa_tab_returns_empty_on_missing_investigation(tmp_path: Path):
    """No ``i_*.txt`` in the directory → no investigation to anchor on."""
    assert parse_isa_tab(tmp_path) == []
