"""SC5 (leafblower-l6h0 / phase 05-02): release-metadata version-drift gate.

Compares `DESCRIPTION`'s `Version:` field against `python/pyproject.toml`'s
`version = "..."` field. These are two independently-editable files with no
shared source of truth (CLAUDE.md: "Version sync: bump DESCRIPTION AND
python/pyproject.toml manually -- no automation"). A `DESCRIPTION`-only or
`pyproject.toml`-only version bump currently surfaces only as a silently
mismatched CRAN/PyPI release pair -- this test catches it as a plain pytest
assertion instead, on every local DoD gate run and in CI.

Verified baseline (2026-08-15): DESCRIPTION Version=0.1.0 == pyproject.toml
version=0.1.0 -- zero drift today.
"""

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_description_pyproject_version_match():
    desc_text = (REPO_ROOT / "DESCRIPTION").read_text()
    pyproject_text = (REPO_ROOT / "python" / "pyproject.toml").read_text()

    desc_version = re.search(r"^Version:\s*(\S+)", desc_text, re.M).group(1)
    py_version = re.search(r'^version\s*=\s*"([^"]+)"', pyproject_text, re.M).group(1)

    assert desc_version == py_version, (
        f"DESCRIPTION Version={desc_version!r} != pyproject.toml version={py_version!r}"
    )
