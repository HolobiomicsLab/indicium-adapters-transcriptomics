"""Free-text to CURIE helpers for transcriptomics repository accessions
and pathway-database identifiers.

Pure functions, no I/O. Used by the bulk adapter's ``enrich_claim`` and by
upstream bridge layers that want to normalise study or run accessions
before emitting claims.

The registry covers two families:

* Repository accessions (GEO / SRA / ENA / ArrayExpress / ENCODE / dbGaP /
  BioStudies / BioProject) — sample, run, study and experiment
  identifiers attached to claims via ``consolidates_dataset``,
  ``sequenced_on`` and similar study-level qualifiers.
* Pathway-database identifiers (REACTOME / KEGG / WikiPathways) — added
  to align with the indicium-adapters #1 spec so the
  ``enriched_in_pathway`` qualifier can resolve bare pathway IDs in its
  object slot.

Adding a new repository or DB: append a ``(regex, prefix)`` row to
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
    # --- Pathway databases (indicium-adapters #1 spec) --------------------
    # Reactome stable identifiers: ``R-{species 3-letter}-{6 to 7 digits}``,
    # e.g. ``R-HSA-9012999`` (canonical Wnt signaling), ``R-MMU-201681``.
    # The species code is uppercase 3-letter ISO-ish (HSA, MMU, DRE, DME, ...).
    (re.compile(r"^R-[A-Z]{3}-\d{4,9}$"), "REACTOME"),
    # KEGG pathway: 2-4 letter organism prefix followed by 5 digits, e.g.
    # ``hsa04310`` (human Wnt signaling), ``mmu04150``, ``map04310``
    # (reference map), ``ko04310`` (orthology map). Most KEGG taxon codes
    # are 3 letters; ``ko`` is the only widely-used 2-letter code, and a
    # handful of recent additions are 4 letters. The pattern accepts
    # the full 2-4 range. Resolution runs before the gene resolver, so a
    # 7-9 character ``xxxNNNNN`` mention will route to KEGG even if a
    # gene symbol happens to share that shape (none of HGNC/MGI/FlyBase/
    # WormBase/ZFIN have such a symbol today).
    (re.compile(r"^[a-z]{2,4}\d{5}$"), "KEGG"),
    # WikiPathways stable identifiers: ``WP`` followed by 1-7 digits,
    # optionally suffixed with a revision (e.g. ``WP554``, ``WP554_r118500``).
    # The ``_r<digits>`` revision tail is uncommon in claim text but kept
    # tolerant.
    (re.compile(r"^WP\d{1,7}(?:_r\d+)?$"), "WP"),
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
