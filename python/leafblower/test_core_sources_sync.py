"""SC4 (leafblower-rywn / phase 02-02): build-list divergence gate.

Compares `src/*.cpp` (minus `r_bridge.cpp`, which is R-only — it includes
`Rinternals.h`, unavailable to the Python/pybind11 build) against
`python/CMakeLists.txt`'s `CORE_SOURCES` list. R auto-globs `src/*.cpp` via its
own build machinery; the Python build does NOT glob — `CORE_SOURCES` is an
explicit, hand-maintained list (CLAUDE.md: "Two build sites for src/*.cpp").
A new `src/*.cpp` added without a matching `CORE_SOURCES` entry currently
surfaces only as an undefined-symbol link error deep into the pybind11 build,
after R's tests have already gone green — this test catches it as a plain
pytest assertion instead.

Verified baseline (2026-08-15): src/*.cpp has 18 files; minus r_bridge.cpp,
17 remain, exactly matching CORE_SOURCES's 17 entries — zero drift today, so
this is pure regression prevention, not a fix for an existing gap.
"""

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

# Excluded from the src-side set: R-only, cannot build for Python (Rinternals.h).
_R_ONLY = {"r_bridge.cpp"}


def test_core_sources_matches_src_glob():
    src_files = {p.name for p in (REPO_ROOT / "src").glob("*.cpp")} - _R_ONLY
    cmake_text = (REPO_ROOT / "python" / "CMakeLists.txt").read_text()
    listed = set(re.findall(r"\.\./src/(\w+\.cpp)", cmake_text))
    assert src_files == listed, (
        f"src/*.cpp vs CORE_SOURCES drift: "
        f"missing from CORE_SOURCES={sorted(src_files - listed)}, "
        f"extra in CORE_SOURCES={sorted(listed - src_files)}"
    )
