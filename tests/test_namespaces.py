from indicium_adapters_transcriptomics._namespaces import TRANSCRIPTOMICS_PREFIXES


def test_prefixes_present_and_iris():
    for key in ("ENSEMBL", "HGNC", "GO", "UBERON", "CL", "EFO"):
        assert key in TRANSCRIPTOMICS_PREFIXES
    assert all(v.startswith("http") for v in TRANSCRIPTOMICS_PREFIXES.values())
