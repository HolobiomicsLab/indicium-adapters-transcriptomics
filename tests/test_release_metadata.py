"""The package version must agree across its release metadata.

``pyproject.toml`` is the single source of the version. The runtime
``__version__``, the ``version`` field of ``CITATION.cff``, the newest
release heading of ``CHANGELOG.md`` and the supported row of the
``SECURITY.md`` table are checked against it, so a bump that misses one
of them fails here instead of reaching a release.
"""
from __future__ import annotations

import re
import tomllib
from pathlib import Path

import indicium_adapters_transcriptomics

ROOT = Path(__file__).resolve().parent.parent


def _pyproject_version() -> str:
    with (ROOT / "pyproject.toml").open("rb") as fh:
        return tomllib.load(fh)["project"]["version"]


def test_runtime_version_matches_pyproject():
    assert indicium_adapters_transcriptomics.__version__ == _pyproject_version()


def test_citation_version_matches_pyproject():
    text = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    match = re.search(r'^version:\s*"?([^"\s]+)"?\s*$', text, re.MULTILINE)
    assert match, "CITATION.cff has no top-level version field"
    assert match.group(1) == _pyproject_version()


def test_newest_changelog_release_matches_pyproject():
    text = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    match = re.search(r"^## v?(\d+\.\d+\.\d+)\b", text, re.MULTILINE)
    assert match, "CHANGELOG.md has no release heading"
    assert match.group(1) == _pyproject_version()


def test_security_table_supports_current_minor():
    text = (ROOT / "SECURITY.md").read_text(encoding="utf-8")
    supported = set(re.findall(r"^\|\s*([^|]+?)\s*\|\s*yes\s*\|", text, re.MULTILINE))
    major, minor, _patch = _pyproject_version().split(".", 2)
    assert f"{major}.{minor}.x" in supported, (
        f"SECURITY.md supports {sorted(supported)}, not the current {major}.{minor}.x line"
    )
