# Contributing

Thanks for your interest! The dev guide is `CLAUDE.md` (works for humans too).
The 4-quadrant filter for new features lives at the top of `ROADMAP.md`.

## Quickstart

```bash
git clone https://github.com/HolobiomicsLab/indicium-adapters-transcriptomics
cd indicium-adapters-transcriptomics
pip install -e .[dev]
pytest
```

`indicium` and `indicium-adapters` are private repos at the moment; pull them
from git first if you don't already have access via PyPI.

## Pull requests

- Branch off `main` with a name like `dina/short-description` (replace `dina`
  with your handle).
- Keep PRs focused — one bug or one feature per PR.
- Run `ruff check .` and `pytest --cov=indicium_adapters_transcriptomics`
  locally before pushing.
- Coverage threshold is 90%.
- New qualifiers / enrichment behaviour: add tests; update CHANGELOG.

## Code style

- `ruff` rules in `pyproject.toml`. Line length 100.
- Hardcoded namespace URIs in `_namespaces.py`; no runtime dep on ontology
  Python packages.
- `enrich_claim()` is a pure function — no I/O, no mutation of input.
