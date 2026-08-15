"""SC3 (leafblower-rywn / phase 02-02): the unit-mode water-fill stays a single source.

Verifies the existing architecture (CONTEXT.md D-03): every bounded solver
delegates its `bounds_mode="unit"` per-cell finalization to the single shared
`calib_dispatch.hpp::finalize_weights`/`finalize_weights_buf` helper instead of
carrying a local copy. This is verification-only per D-03/D-04 -- no solver
code changes here. Two independent regressions this guards against:
  1. A solver stops calling the shared helper (re-copies the water-fill logic
     locally instead) -- Test 1.
  2. The helper itself gets duplicated into a second definition site -- Test 2.

Both are pure text scans: no compiled artifact, no BLAS threading dependence.
"""

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC = REPO_ROOT / "src"

# The 7 solvers verified (2026-08-15 direct read, CONTEXT.md D-03) to delegate
# per-cell bounds_mode="unit" finalization to the shared helper.
_SOLVERS_EXPECTED_TO_CALL_FINALIZE = {
    "oris_finalize.cpp",
    "raking.cpp",
    "chebyshev.cpp",
    "greenkhorn.cpp",
    "greg.cpp",
    "logit_calib.cpp",
    "sinkhorn.cpp",
}

# newton_calib.cpp is deliberately EXCLUDED, not silently omitted:
# newton_kl is a smooth-dual Newton method with no box-constrained inner step
# (unlike ORIS/raking's cell-table water-fill machinery) -- it only COUNTS
# bound violations (see newton_calib.cpp's `n_bounds_violated` accounting) and
# falls back to RK_ERR_NOCONV above a 5% violation fraction; it never
# redistributes per observation and never calls finalize_weights/
# finalize_weights_buf at all (RESEARCH.md Pitfall 2 / Assumption A1). This is
# an accepted pre-existing capability gap, NOT evidence that bounds_mode="unit"
# is enforced for newton_kl -- tracked on leafblower-og7d.5 (measured seed
# 7/11/23 out-of-bounds output on bounds_mode="unit"; report-not-clamp T4
# contract; no fix scheduled before Phase 2). No new duplicate ticket filed --
# leafblower-og7d.5 already records this exact gap; see its 2026-08-15 comment
# cross-referencing this test.
_EXCLUDED_NO_FINALIZE = "newton_calib.cpp"

# raking.cpp's OWN `water_fill_cat` (an inline per-margin IPF box projection
# integral to the raking algorithm itself) is different math from the shared
# post-hoc per-cell finalize step this test asserts about (CONTEXT.md D-03,
# RESEARCH.md "Water-fill consolidation (SC3)"). raking.cpp is asserted here
# purely for its (separate) call to the shared `finalize_weights` helper,
# not for water_fill_cat.


def test_every_bounded_solver_calls_shared_finalize():
    missing = []
    for name in sorted(_SOLVERS_EXPECTED_TO_CALL_FINALIZE):
        text = (SRC / name).read_text()
        if "finalize_weights" not in text:
            missing.append(name)
    assert not missing, (
        f"solver(s) stopped delegating to the shared finalize helper: {missing}"
    )


def test_finalize_weights_buf_defined_exactly_once():
    definition_re = re.compile(r"inline void finalize_weights_buf\(")
    defining_files = [
        p.name
        for p in list(SRC.glob("*.cpp")) + list(SRC.glob("*.hpp"))
        if definition_re.search(p.read_text())
    ]
    assert defining_files == ["calib_dispatch.hpp"], (
        f"finalize_weights_buf must be defined exactly once, in calib_dispatch.hpp; "
        f"found in: {defining_files}"
    )
