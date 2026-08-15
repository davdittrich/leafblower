"""SC1 (leafblower-rywn / phase 02-08): single-dispatch-site regression gate.

Phase 02 collapsed R's `strcmp(method_str, ...)` per-method dispatch chain
(`src/r_bridge.cpp`, plans 02-01 through 02-07) into calls through the single
shared `lbw::dispatch_solver` table (`src/calib_dispatch.hpp`). That
consolidation is correct today and undefended tomorrow unless a test fails
the moment someone re-introduces a second per-method chain -- this is that
test.

Pure text scan, same shape as its two SC3/SC4 siblings
(`test_finalize_weights_sync.py`, `test_core_sources_sync.py`): no import of
the compiled module, no BLAS dependence.

Verified baseline (2026-08-15, post plan 02-07): zero `strcmp(method_str, ...)`
occurrences remain in code (one mention survives inside a `//` comment
documenting the removal, at r_bridge.cpp:589); `lbw::dispatch_solver(` is
called from exactly 3 call sites -- the AUTO primary call, the AUTO
newton_kl fallback call, and the one unified explicit-method call (see
02-07-SUMMARY.md's D4 coverage entry, which enumerates these same 3 sites).
"""

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
R_BRIDGE = REPO_ROOT / "src" / "r_bridge.cpp"

# Bound on lbw::dispatch_solver(...) call sites: the AUTO branch dispatches
# up to twice (primary + newton_kl fallback), the non-AUTO branch dispatches
# once -- 3 total, never one call site per solver.
_MAX_DISPATCH_CALL_SITES = 3


def _strip_line_comments(text):
    """Drop every line whose first non-whitespace characters open a `//`
    comment, so documenting the retired pattern in a comment cannot make
    the guard self-invalidating."""
    return "\n".join(
        line for line in text.splitlines() if not line.lstrip().startswith("//")
    )


def test_r_bridge_has_no_method_dispatch_chain():
    code = _strip_line_comments(R_BRIDGE.read_text())

    strcmp_hits = code.count("strcmp(method_str")
    assert strcmp_hits == 0, (
        f"src/r_bridge.cpp still contains {strcmp_hits} strcmp(method_str, ...) "
        "comparison(s) -- that is the per-method dispatch chain phase 02 "
        "removed. A new solver belongs in the shared table in "
        "calib_dispatch.hpp, not a new branch here."
    )

    call_sites = len(re.findall(r"lbw::dispatch_solver\(", code))
    assert call_sites <= _MAX_DISPATCH_CALL_SITES, (
        f"src/r_bridge.cpp calls lbw::dispatch_solver() {call_sites} times; "
        f"expected at most {_MAX_DISPATCH_CALL_SITES} (AUTO primary, AUTO "
        "newton_kl fallback, and the one unified explicit-method call). A "
        "per-method dispatch chain would grow this toward one call per "
        "solver -- see calib_dispatch.hpp for where a new routing rule "
        "belongs instead."
    )
