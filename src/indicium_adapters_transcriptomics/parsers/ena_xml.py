"""ENA / SRA XML parser.

`ENA <https://www.ebi.ac.uk/ena>`__'s ``study`` and ``experiment``
endpoints return XML envelopes shaped like::

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
strategy / source / selection) and ``PLATFORM`` (vendor and instrument
model, e.g. Illumina / PacBio / Nanopore).

This parser handles both envelopes from the same code path: it walks
every ``<STUDY>`` and ``<EXPERIMENT>`` element and emits one claim per
useful field. ENA accessions (``PRJEB...``, ``SRP...``, ``ERR...``) are
recognised downstream by the bulk adapter's accession resolver.
"""
from __future__ import annotations

from xml.etree import ElementTree as ET


def parse_ena_xml(xml_text: str) -> list[dict]:
    """Parse ENA / SRA XML and return study-level claim dicts.

    Recognises ``<STUDY>`` and ``<EXPERIMENT>`` elements. For each
    experiment, emits a ``consolidates_dataset`` claim linking its
    project accession to the experiment, plus a ``sequenced_on`` claim
    for the instrument model. ``aligned_with`` is not emitted because
    ENA descriptors do not typically include alignment-pipeline metadata.
    """
    root = ET.fromstring(xml_text)
    claims: list[dict] = []

    # Walk every STUDY element (top-level or nested). The body is empty
    # for now: titles and abstracts do not map cleanly to a qualifier in
    # the current schema, so we record the accession only.
    for study in root.iter("STUDY"):
        study_acc = study.attrib.get("accession", "").strip()
        if not study_acc:
            # Fall back to the alias / PRIMARY_ID under IDENTIFIERS when
            # the accession attribute is missing.
            for prim in study.iter("PRIMARY_ID"):
                if prim.text:
                    study_acc = prim.text.strip()
                    break
        if not study_acc:
            continue

    # Walk every EXPERIMENT element. Each carries a study reference and
    # a platform / instrument model.
    for exp in root.iter("EXPERIMENT"):
        exp_acc = exp.attrib.get("accession", "").strip()
        if not exp_acc:
            continue
        # Find the parent study reference, if any.
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
        # Platform / instrument. PLATFORM children are vendor-named
        # (<ILLUMINA>, <OXFORD_NANOPORE>, <PACBIO_SMRT>, ...) and each
        # contains an <INSTRUMENT_MODEL> with the actual model string.
        platform = exp.find("PLATFORM")
        if platform is not None:
            for vendor_elem in platform:
                model_elem = vendor_elem.find("INSTRUMENT_MODEL")
                if model_elem is not None and model_elem.text:
                    claim_subject = study_ref or exp_acc
                    claims.append({
                        "subject": claim_subject,
                        "qualifier": "sequenced_on",
                        "object": model_elem.text.strip(),
                    })

    return claims
