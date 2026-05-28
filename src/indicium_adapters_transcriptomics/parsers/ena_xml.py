"""ENA / SRA XML parser (v0.7.0).

`ENA <https://www.ebi.ac.uk/ena>`__'s ``study`` and ``experiment``
endpoints return an XML envelope shaped like::

    <STUDY_SET>
      <STUDY accession="PRJEB...">
        <IDENTIFIERS>...</IDENTIFIERS>
        <DESCRIPTOR>
          <STUDY_TITLE>...</STUDY_TITLE>
          <STUDY_ABSTRACT>...</STUDY_ABSTRACT>
        </DESCRIPTOR>
      </STUDY>
    </STUDY_SET>

The ``experiment`` envelope adds ``LIBRARY_DESCRIPTOR`` (library
strategy / source / selection) and ``PLATFORM`` (Illumina / PacBio /
Nanopore vendor + instrument model).

This parser handles both envelopes via the same code path: it walks
all ``<STUDY>`` and ``<EXPERIMENT>`` elements and emits one claim per
useful field. ENA accessions (``PRJEB...``, ``SRP...``, ``ERR...``)
flow through the bulk adapter's accession resolver downstream.
"""
from __future__ import annotations

from xml.etree import ElementTree as ET


def parse_ena_xml(xml_text: str) -> list[dict]:
    """Parse ENA / SRA XML and return study-level claim dicts.

    Recognises ``<STUDY>`` and ``<EXPERIMENT>`` elements. For each
    study, emits ``consolidates_dataset`` claims linking the project
    accession to its experiments/samples; for each experiment, emits
    ``sequenced_on`` claims for the instrument model and
    ``aligned_with`` is **not** emitted (ENA descriptors don't
    typically include alignment-pipeline metadata).
    """
    root = ET.fromstring(xml_text)
    claims: list[dict] = []

    # Walk all STUDY elements (top-level or nested).
    for study in root.iter("STUDY"):
        study_acc = study.attrib.get("accession", "").strip()
        if not study_acc:
            # Fall back to the alias / PRIMARY_ID under IDENTIFIERS.
            for prim in study.iter("PRIMARY_ID"):
                if prim.text:
                    study_acc = prim.text.strip()
                    break
        if not study_acc:
            continue
        # Title — useful as descriptive data, but no qualifier maps
        # cleanly. Skipped in v0.7.0.

    # Walk all EXPERIMENT elements. Each carries a study reference and
    # a platform / instrument model.
    for exp in root.iter("EXPERIMENT"):
        exp_acc = exp.attrib.get("accession", "").strip()
        if not exp_acc:
            continue
        # Find the study reference.
        study_ref = ""
        for ref in exp.iter("STUDY_REF"):
            study_ref = ref.attrib.get("accession", "").strip()
            if study_ref:
                break
        if study_ref:
            claims.append({
                "subject": study_ref,
                "qualifier": "consolidates_dataset",
                "object": exp_acc,
            })
        # Platform / instrument.
        platform = exp.find("PLATFORM")
        if platform is not None:
            for vendor_elem in platform:
                # Children are vendor-named: <ILLUMINA>, <OXFORD_NANOPORE>,
                # <PACBIO_SMRT>, etc. Each contains an <INSTRUMENT_MODEL>.
                model_elem = vendor_elem.find("INSTRUMENT_MODEL")
                if model_elem is not None and model_elem.text:
                    claim_subject = study_ref or exp_acc
                    claims.append({
                        "subject": claim_subject,
                        "qualifier": "sequenced_on",
                        "object": model_elem.text.strip(),
                    })

    return claims
