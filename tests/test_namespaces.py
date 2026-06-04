from indicium_adapters_transcriptomics._namespaces import TRANSCRIPTOMICS_PREFIXES


def test_prefixes_present_and_iris():
    for key in ("ENSEMBL", "HGNC", "GO", "UBERON", "CL", "EFO"):
        assert key in TRANSCRIPTOMICS_PREFIXES
    assert all(v.startswith("http") for v in TRANSCRIPTOMICS_PREFIXES.values())


def test_pathway_database_prefixes_present():
    """indicium-adapters issue #1 alignment: REACTOME, KEGG and WP
    (WikiPathways) are the three pathway-DB CURIE prefixes the spec
    proposes for the ``enriched_in_pathway`` qualifier.

    Pin them so that resolving ``REACTOME:R-HSA-9012999`` or
    ``KEGG:hsa04310`` to its expanded IRI does not regress when a future
    edit to the prefix dict accidentally drops one.
    """
    for key in ("REACTOME", "KEGG", "WP"):
        assert key in TRANSCRIPTOMICS_PREFIXES, (
            f"pathway-DB prefix {key!r} missing from TRANSCRIPTOMICS_PREFIXES"
        )
        assert TRANSCRIPTOMICS_PREFIXES[key].startswith(
            "https://identifiers.org/"
        ), f"{key!r} should expand via the identifiers.org compact-URI scheme"
