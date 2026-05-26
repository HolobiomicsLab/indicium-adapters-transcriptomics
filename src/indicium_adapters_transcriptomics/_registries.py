"""Free-text → CURIE helpers for transcriptomics repository accessions.

Pure functions, no I/O. Used by adapter `enrich_claim()` implementations and by
upstream bridge layers that want to normalise study/run accessions before
emitting claims.
"""
from __future__ import annotations

import re
from typing import Final

# Accession patterns. Kept conservative: anchored, well-known prefixes only.
# To extend, add a (regex, prefix) row and keep the regex strict.
_ACCESSION_PATTERNS: Final[tuple[tuple[re.Pattern[str], str], ...]] = (
    (re.compile(r"^GSE\d{3,7}$"), "GEO"),
    (re.compile(r"^GSM\d{3,7}$"), "GEO"),
    (re.compile(r"^GPL\d{3,7}$"), "GEO"),
    (re.compile(r"^GDS\d{3,7}$"), "GEO"),
    (re.compile(r"^SRP\d{6,9}$"), "SRA"),
    (re.compile(r"^SRX\d{6,9}$"), "SRA"),
    (re.compile(r"^SRR\d{6,9}$"), "SRA"),
    (re.compile(r"^SRS\d{6,9}$"), "SRA"),
    (re.compile(r"^PRJ[A-Z]{1,2}\d{4,9}$"), "BioProject"),
    (re.compile(r"^ERP\d{6,9}$"), "ENA"),
    (re.compile(r"^ERX\d{6,9}$"), "ENA"),
    (re.compile(r"^ERR\d{6,9}$"), "ENA"),
    (re.compile(r"^ERS\d{6,9}$"), "ENA"),
    (re.compile(r"^E-[A-Z]{4}-\d+$"), "ArrayExpress"),
    (re.compile(r"^ENC[A-Z]{2}\d{3}[A-Z]{3}$"), "ENCODE"),
    (re.compile(r"^phs\d{6}(?:\.v\d+\.p\d+)?$"), "dbGaP"),
    (re.compile(r"^S-[A-Z]{4}\d+$"), "BioStudies"),
)


def accession_to_curie(accession: str) -> str | None:
    """Return ``"<prefix>:<accession>"`` if the accession matches a known
    transcriptomics repository pattern, else ``None``.

    >>> accession_to_curie("GSE168938")
    'GEO:GSE168938'
    >>> accession_to_curie("PRJEB104129")
    'BioProject:PRJEB104129'
    >>> accession_to_curie("not-a-real-accession") is None
    True
    """
    if not isinstance(accession, str):
        return None
    token = accession.strip()
    for pattern, prefix in _ACCESSION_PATTERNS:
        if pattern.match(token):
            return f"{prefix}:{token}"
    return None
