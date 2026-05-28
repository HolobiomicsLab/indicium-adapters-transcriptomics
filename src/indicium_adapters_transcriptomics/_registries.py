"""Free-text to CURIE helpers for transcriptomics repository accessions.

Pure functions, no I/O. Used by the bulk adapter's ``enrich_claim`` and by
upstream bridge layers that want to normalise study or run accessions
before emitting claims.

Adding a new repository: append a ``(regex, prefix)`` row to
``_ACCESSION_PATTERNS`` and keep the regex anchored and conservative. The
prefix must already be registered in ``_namespaces.TRANSCRIPTOMICS_PREFIXES``
so the produced CURIE expands correctly downstream.
"""
from __future__ import annotations

import re
from typing import Final

# Accession patterns. Anchored, well-known prefixes only. Each row is
# ``(compiled regex, CURIE prefix)``. The first matching regex wins.
_ACCESSION_PATTERNS: Final[tuple[tuple[re.Pattern[str], str], ...]] = (
    # NCBI GEO: series (GSE), sample (GSM), platform (GPL), data set (GDS).
    (re.compile(r"^GSE\d{3,7}$"), "GEO"),
    (re.compile(r"^GSM\d{3,7}$"), "GEO"),
    (re.compile(r"^GPL\d{3,7}$"), "GEO"),
    (re.compile(r"^GDS\d{3,7}$"), "GEO"),
    # NCBI SRA: project (SRP), experiment (SRX), run (SRR), sample (SRS).
    (re.compile(r"^SRP\d{6,9}$"), "SRA"),
    (re.compile(r"^SRX\d{6,9}$"), "SRA"),
    (re.compile(r"^SRR\d{6,9}$"), "SRA"),
    (re.compile(r"^SRS\d{6,9}$"), "SRA"),
    # BioProject: PRJNA / PRJEB / PRJDB.
    (re.compile(r"^PRJ[A-Z]{1,2}\d{4,9}$"), "BioProject"),
    # ENA: project (ERP), experiment (ERX), run (ERR), sample (ERS).
    (re.compile(r"^ERP\d{6,9}$"), "ENA"),
    (re.compile(r"^ERX\d{6,9}$"), "ENA"),
    (re.compile(r"^ERR\d{6,9}$"), "ENA"),
    (re.compile(r"^ERS\d{6,9}$"), "ENA"),
    # ArrayExpress: E-XXXX-N, e.g. E-MTAB-1234.
    (re.compile(r"^E-[A-Z]{4}-\d+$"), "ArrayExpress"),
    # ENCODE: ENCxxxxxxx ids, e.g. ENCSR000ABC.
    (re.compile(r"^ENC[A-Z]{2}\d{3}[A-Z]{3}$"), "ENCODE"),
    # dbGaP studies: phsNNNNNN with optional version / participant suffix.
    (re.compile(r"^phs\d{6}(?:\.v\d+\.p\d+)?$"), "dbGaP"),
    # BioStudies: S-XXXXN, e.g. S-BSST123.
    (re.compile(r"^S-[A-Z]{4}\d+$"), "BioStudies"),
)


def accession_to_curie(accession: str) -> str | None:
    """Return ``"<prefix>:<accession>"`` if the input matches a known
    transcriptomics repository pattern, else ``None``.

    Whitespace is stripped before matching. The accession is preserved
    verbatim in the returned CURIE (case included), since accession case
    is meaningful in every supported namespace.

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
