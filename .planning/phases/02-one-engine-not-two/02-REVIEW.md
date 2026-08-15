---
phase: 02-one-engine-not-two
reviewed: 2026-08-15T00:00:00Z
depth: standard
files_reviewed: 9
files_reviewed_list:
  - src/calib_dispatch.hpp
  - src/c_api.cpp
  - src/r_bridge.cpp
  - python/CMakeLists.txt
  - python/leafblower/test_core_sources_sync.py
  - python/leafblower/test_finalize_weights_sync.py
  - python/leafblower/test_single_dispatch_site.py
  - python/leafblower/test_solver_parity.py
  - CLAUDE.md
findings:
  critical: 0
  warning: 4
  info: 0
  total: 4
status: issues_found
---

# Phase 02: Code Review Report

**Reviewed:** 2026-08-15T00:00:00Z
**Depth:** standard
**Files Reviewed:** 9
**Status:** issues_found

## Summary

This phase collapses R bridge / C API's per-solver `strcmp`/switch dispatch
chains into a single `lbw::dispatch_solver` table plus `lbw::route_auto` for
AUTO routing, both in `src/calib_dispatch.hpp`. I diffed every file against
`a8b2e57a9b692c6f10c7d725ce99efef2491aaea^` line-by-line (not just read the
final state) specifically hunting for the three failure modes called out in
the task: accidental solver-math changes, silently dropped/mismatched R↔Python
result fields, and the two named footguns (`CORE_SOURCES` desync,
premature post-water-fill renormalization).

**None of the three targeted failure modes are present:**

- `src/calib_dispatch.hpp`'s diff is purely additive (0 lines removed) — the
  pre-existing `finalize_weights_buf` (normalize-before-bounds order, no
  post-water-fill renormalize) and `compute_cell_metrics`/`apply_rule` bodies
  are byte-for-byte untouched by this phase.
- `python/CMakeLists.txt`'s `CORE_SOURCES` (17 entries) currently matches
  `src/*.cpp` minus `r_bridge.cpp` (17 files) exactly, and a new regression
  test (`test_core_sources_sync.py`) now guards this going forward.
- Two intentional behavior deltas were found (`st.oris_auto_selected` narrowed
  to `alg == RK_ALG_ORIS` in `c_api.cpp`, and `st.accelerate` now restored
  before the AUTO newton_kl fallback in `c_api.cpp`) — both are explicitly
  documented as adopting R's pre-existing, fixture-pinned behavior, and both
  are provably inert on solver *output* (`oris_auto_selected` only gates a
  verbose log prefix in `oris.cpp:353`; `newton_calib.cpp`/`.hpp` never reads
  `st.accelerate`, confirmed by grep). Neither is a solver-math change.
- The 49-field R result-list assembly in `r_bridge.cpp` (lines ~884-1030) is
  untouched by the diff; every `DispatchResult` field the new
  `pack_dispatch_result_c`/`pack_dispatch_oris_extras_c` (`c_api.cpp`) and the
  `dres.* -> res_*` copies (`r_bridge.cpp`) perform was traced against the old
  per-branch code and matches field-for-field, including the subtle
  pre-existing asymmetry where only the newton_kl branch resets
  `min_alpha_seen`/`final_alpha` to `1.0` (other non-ORIS solvers leave them
  at `rk_result_init`'s `memset`-0 default) — faithfully reproduced.

What remains are four maintainability/robustness findings below: the
migration replaced 9 *solver-selection* chains with 1, but the
*result-field-copy* code inside `dispatch_solver`'s 9 case arms is still
~30 lines duplicated 9 times, and two invariants ("always non-empty
`best_weights`" and "never dispatch an unmigrated/AUTO enum value") are
enforced only by comments, not by asserts or tests — the same gap class
this phase's own `test_core_sources_sync.py`/`test_finalize_weights_sync.py`
were written to close for other invariants.

## Warnings

### WR-01: `dispatch_solver`'s per-case field-copy blocks re-duplicate the pattern this phase set out to eliminate

**File:** `src/calib_dispatch.hpp:751-1122`
**Issue:** Each of the 9 `case RK_ALG_*` arms in `dispatch_solver` repeats an
near-identical ~25-30 line block copying `res.base.*` into `out.*`
(`status`, `iterations`, `max_error`, `mean_error`, `kl`, `chi2`,
`l1_weight_change`, `grake_norm`, `convergence_metric`, ... 20+ fields). This
is exactly the maintenance hazard the phase's stated goal ("eliminates
triplicated switch blocks", header comment at `calib_dispatch.hpp:1-3`)
targets — it has just moved from the *caller* (r_bridge.cpp/c_api.cpp) into
the *table* itself. Adding, renaming, or removing a `CalibResultBase` field
now requires 9 synchronized edits instead of 1, with no compiler or test
signal if one arm is missed (a field silently stays at `DispatchResult`'s
default in the missed arm).
**Fix:**
```cpp
// One helper, called at the top of every case arm instead of the 20-line block:
template <typename ResT>
inline void copy_base_fields(DispatchResult& out, const ResT& res) noexcept {
    out.status                       = res.base.status;
    out.iterations                   = res.base.iterations;
    out.max_error                    = res.base.max_error;
    out.mean_error                   = res.base.mean_error;
    out.kl                           = res.base.kl;
    out.chi2                         = res.base.chi2;
    out.l1_weight_change             = res.base.l1_weight_change;
    out.grake_norm                   = res.base.grake_norm;
    out.convergence_metric           = res.base.convergence_metric;
    out.convergence_rule             = res.base.convergence_rule;
    out.convergence_tol              = res.base.convergence_tol;
    out.convergence_iter             = res.base.convergence_iter;
    out.convergence_solver_objective = res.base.convergence_solver_objective;
    out.convergence_minimized_metric = res.base.convergence_minimized_metric;
    out.best_error                   = res.base.best_error;
    out.best_iter                    = res.base.best_iter;
    out.metric_first_check           = res.base.metric_first_check;
    out.metric_prev_check            = res.base.metric_prev_check;
    out.prev_check_iter              = res.base.prev_check_iter;
    out.stall_kind                   = res.base.stall_kind;
    out.n_bounds_violated            = res.n_bounds_violated;
    out.n_bounds_clamped             = res.n_bounds_clamped;
}
// case RK_ALG_SINKHORN: { auto res = lbw::sinkhorn_solve(st); copy_base_fields(out, res); ... }
```

### WR-02: `dispatch_solver`'s unreachable `default:` case fails silently, masquerading as ORIS

**File:** `src/calib_dispatch.hpp:1120-1122`, `DispatchResult` defaults at `src/calib_dispatch.hpp:667-671`
**Issue:** `default: break;` leaves `out` fully default-constructed:
`status=RK_ERR_NOCONV`, `alg_used=RK_ALG_ORIS`, `iterations=0`,
`max_error=1.0`. If `dispatch_solver` is ever called with `RK_ALG_AUTO`
(both current callers pre-resolve AUTO via `route_auto`, but nothing enforces
that at this call boundary) or a reserved/removed slot (2, 7), the caller
receives a result that is indistinguishable from "ORIS actually ran and
failed to converge" rather than "dispatch_solver was called with an invalid
algorithm" — a silent-wrong-result failure mode, not a crash, which is the
worst kind for a numerics library per this project's correctness philosophy.
Currently dead code (comment: "not yet migrated (D-01); caller's existing
branch handles it" — but there is no unmigrated case left after this phase),
so there is no live caller-side workaround to mask it either.
**Fix:**
```cpp
default:
    assert(false && "dispatch_solver: unhandled/unmigrated rk_algorithm_t");
    out.status = RK_ERR_BADARG;
    break;
```

### WR-03: "best_weights is never empty" is a comment-only invariant across 7 solver arms, unlike its 2 documented exceptions

**File:** `src/r_bridge.cpp:579-587` (comment), `src/calib_dispatch.hpp` (raking/oris/oris_soft/sinkhorn/greg/greenkhorn/logit case arms)
**Issue:** `chebyshev` and `newton_kl`'s dispatch arms explicitly guard against
an empty `res.base.best_weights` (violation-guard paths) and fall back to a
zero-filled sentinel of length `st.n` before assigning `out.best_weights`.
The other 7 arms (`sinkhorn`, `greg`, `greenkhorn`, `logit`, `raking`, `oris`,
`oris_soft`) do an unconditional `out.best_weights = std::move(res.base.best_weights);`
with no such guard, relying on the r_bridge.cpp comment's claim that these 7
solvers "never leave it empty on any path (verified by reading each solver's
best_weights assignment sites)". That verification is not machine-checked —
if a future edit to any of those 7 solvers introduces an empty-vector exit
path (plausible, since 2 of the other 9 solvers already have one), the R
result's `best_weights` SEXP element silently drops to length 0 while
`status` continues to report `RK_OK`/`RK_ERR_NOCONV` normally, and any R/harvest
code indexing `best_weights` by observation would break or misalign. This
project already added `test_core_sources_sync.py`, `test_finalize_weights_sync.py`,
and `test_single_dispatch_site.py` specifically to convert "verified by
reading the code" claims like this into machine-enforced regression gates.
**Fix:** Add the same guard uniformly across all 9 arms (cheapest, and
removes the asymmetry), or add a light assertion inside `dispatch_solver`
after populating `out.best_weights`:
```cpp
assert((out.best_weights.empty() || static_cast<int>(out.best_weights.size()) == st.n) &&
       "dispatch_solver: best_weights must be empty or length st.n");
```
combined with a length-`st.n` non-empty check before the `move` for the 7
currently-unguarded arms.

### WR-04: `DispatchResult::alg_used` defaults to a real algorithm (`RK_ALG_ORIS`) instead of a neutral sentinel

**File:** `src/calib_dispatch.hpp:671`
**Issue:** `rk_algorithm_t alg_used = RK_ALG_ORIS;` compounds WR-02: even if
the `default:` case above is hit, `out.alg_used` reads as a plausible,
real algorithm rather than an obviously-wrong sentinel. Every populated case
arm sets `alg_used` explicitly, so this default is never used unless WR-02's
dead path is (accidentally) reached — but a neutral default would make that
condition fail loudly (e.g., an "unknown" branch in the caller's algorithm-name
lookup) instead of quietly reporting a specific, wrong algorithm.
**Fix:**
```cpp
rk_algorithm_t alg_used = RK_ALG_AUTO;  // never a real dispatch outcome; a caller
                                          // that sees this back knows dispatch_solver
                                          // did not run a real case arm.
```

---

_Reviewed: 2026-08-15T00:00:00Z_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_
