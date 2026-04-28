# Convergence Status Redesign

**Date**: 2026-04-28
**Status**: Pending design review (rev 2 — gate round-1 fixes)
**Scope**: `src/raking.cpp` + `src/leafblower.h` + `R/harvest.R`
**Deferred**: ieppa.cpp stall changes (separate spec required — ieppa stall is PCT-based, not errRp-based; 5 regression fixes in last week make this high-risk without a dedicated plan)

---

## Problem

### Status 1 conflates two incompatible diagnostics

`RK_ERR_NOCONV = 1` fires when:
- **Budget exhausted** — solver was still improving; prescription: increase `max_iterations`
- **Stall** — loss function plateaued; prescription: done, optimum reached; weights are valid, accept them

Opposite prescriptions. A user checking `status == 1` cannot determine which case applies.

### Stall uses errRp, not the solver's loss function

Current stall detector (`min_errRp_window`, `kMaxNoImprove=5`): fires when errRp (max marginal error) fails to improve by ≥1% for 5 consecutive convergence-check intervals.

Problems:
1. **errRp is not monotone** near the constrained KL optimum. errRp oscillates as cyclic margin projections cycle through categories — it can stagnate while weight KL is still strictly decreasing. Result: false stall, early exit before the true optimum.
2. **errRp measures the wrong thing**. Raking minimizes weight KL. Water-filling bounded IPF is coordinate descent on weight KL, monotone by the KL projection theorem (Csiszar-Tusnady 1984). Each F_eval strictly decreases weight KL until the constrained minimum. **Weight KL plateau ↔ constrained KL minimum reached ↔ stall is correct signal.**

---

## Solution

### New status codes (appended — no renumbering of 0–3)

```c
// Existing (preserved for ABI/caller compatibility — unchanged):
#define RK_OK           0   // Converged: improvement criterion satisfied
#define RK_ERR_NOCONV   1   // Legacy: not emitted by new solvers; kept for ABI
#define RK_ERR_INFEAS   2   // Infeasible: empty cell with positive target
#define RK_ERR_BADARG   3   // Invalid argument

// New (appended):
#define RK_ERR_BUDGET   4   // Budget exhausted while loss still decreasing; increase max_iterations
#define RK_ERR_STALL    5   // Loss function plateau — at constrained optimum; weights are valid
```

**Breaking change**: callers using `status == 1` no longer receive it from `raking_solve`. They receive 4 or 5. Callers using `status > 0` for "not converged" continue to work. **This must be noted in NEWS/CHANGELOG.**

**ABI**: `rk_result_t.status` is int; `EXPECTED_RK_RESULT_BYTES` unchanged (no new struct fields).

### Stall criterion: weight KL for raking

Replace `min_errRp_window` + errRp-based stall in `raking_solve` with weight KL stall:

```cpp
// Ordering: convergence criterion check BEFORE stall detector on every convergence-check interval.
// This ensures perfect calibration (wkl=0) exits via criterion (status=OK), not stall.

// Primary convergence criterion (existing, unchanged):
lbw::CellMetrics m_conv; m_conv.errRp = errRp; ...
if (lbw::check_convergence(st.convergence_cfg, m_conv, prev_metric_for_rule, st.tol_abs)) {
    res.status = RK_OK; ... break;  // criterion fires first
}

// Stall detector (AFTER criterion check):
double wkl = compute_weight_kl();
// Guard: wkl ≤ tol_abs means effectively at optimum — exit as converged, not stalled.
if (wkl <= st.tol_abs) { res.status = RK_OK; res.convergence_iter = iter; break; }

if (!std::isfinite(min_loss_window)) {
    min_loss_window = wkl; n_no_improve = 0;
} else if (wkl < min_loss_window * (1.0 - st.convergence_cfg.pct_tol)) {
    min_loss_window = wkl; n_no_improve = 0;
} else {
    n_no_improve++;
}
if (n_no_improve >= kMaxNoImprove) { res.status = RK_ERR_STALL; break; }
```

**Ordering contract**: convergence criterion runs before stall check every iteration. This prevents false stall when wkl=0 (perfect calibration) or wkl is decreasing but below pct_tol threshold.

### Budget vs stall priority

When loop exits via condition (budget exhausted), status is set **after** the loop, only if not already set by stall or criterion:

```cpp
// After the loop (flat or SQUAREM):
if (res.status == RK_ERR_NOCONV)  // initial value; not set by stall or criterion
    res.status = RK_ERR_BUDGET;   // budget wins over stall only if stall never fired
```

**Priority**: criterion (0) > stall (fires during loop, breaks) > budget (set post-loop). Budget and stall cannot fire simultaneously: stall fires inside the loop with `break`, budget is set post-loop only for the initial NOCONV state.

### SQUAREM stall: errRp-based (pragmatic)

SQUAREM's CBB extrapolation is NOT monotone in weight KL. The accepted iterate's KL can increase after step-halving. Therefore: **SQUAREM stall uses errRp** (same as current), not weight KL.

The semantic difference from flat-loop stall=5:
- Flat loop status=5 (stall): **KL-guaranteed constrained optimum**
- SQUAREM status=5 (stall): **errRp plateau** — empirically near the optimum, not a KL-proof

To distinguish, `convergence_reason` includes `"stall_kl"` vs `"stall_errRp"`:
```r
convergence_reason =
    if (status == 0L) "criterion"
    else if (status == 4L) "budget"
    else if (status == 5L && accelerate_bool) "stall_errRp"  # SQUAREM: errRp plateau
    else if (status == 5L) "stall_kl"                         # flat: KL plateau (guaranteed optimum)
    else if (status == 2L) "infeasible"
    else if (status == 3L) "error"
    else "legacy"  # status=1: old solvers
```

### harvest.R: status handling

```r
# Hard stops (unchanged):
if (calib_result$status == 2L)
  stop("leafblower: infeasible problem — persistent empty cell with positive target.")
if (calib_result$status == 3L)
  stop("leafblower: invalid arguments — ", calib_result$message)

# Soft exits — weights are valid in all of these:
if (calib_result$status == 4L)
  warning("leafblower: budget exhausted — weights reflect best iterate; ",
          "increase max_iterations if further improvement is needed")
if (calib_result$status == 5L && isTRUE(accelerate_bool))
  warning("leafblower: SQUAREM errRp plateau — weights are valid; ",
          "try accelerate=FALSE for guaranteed KL-minimum convergence")
if (calib_result$status == 5L && !isTRUE(accelerate_bool))
  warning("leafblower: loss function plateau — at constrained optimum given bounds; ",
          "weights are valid; no further improvement is achievable")
if (calib_result$status == 1L)  # legacy, not emitted by new raking/ieppa
  warning("leafblower: did not converge (legacy status code from solver not yet updated)")
```

`convergence_reason` added to `result$convergence_used` (alongside existing `metric`, `rule`, `tol`, `fired_at_iter`). Also add to `@param` and `@return` documentation.

**@return addition**:
```r
#'         \item \code{convergence_used}: ... (existing) ...
#'         \item \code{convergence_used$convergence_reason}: Character string describing
#'           why the solver exited: \code{"criterion"} (improvement criterion satisfied),
#'           \code{"budget"} (max_iterations exhausted, still improving — increase budget),
#'           \code{"stall_kl"} (weight KL plateau — constrained KL minimum reached),
#'           \code{"stall_errRp"} (SQUAREM errRp plateau — empirically near optimum),
#'           \code{"infeasible"} (structural infeasibility detected),
#'           \code{"error"} (invalid arguments),
#'           \code{"legacy"} (old solver not yet emitting v2 codes).
```

---

## Acceptance Criteria

1. **AC1**: `status=0` fires only when primary convergence criterion (`pct_tol` improvement) satisfies. **Unchanged.**
2. **AC2**: `status=4` fires when loop exhausts budget while weight KL was still decreasing on the last check. Test: `max_iterations=5`, multi-margin problem far from converged → `status=4` AND `wkl_final < wkl_initial`.
3. **AC3**: `status=5` fires when weight KL fails to improve by `pct_tol` fraction for `kMaxNoImprove` consecutive check intervals (flat loop). Test: problem already at constrained minimum (all bounds active, errRp plateaued), `max_iterations=1000`, `pct_tol=0.01` → `status=5` before budget.
4. **AC4**: Perfect calibration (wkl=0) → `status=0`, not `status=5`. Guard: `if (wkl ≤ tol_abs) → RK_OK`.
5. **AC5**: `status=4` and `status=5` return valid weights (best iterate). `harvest()` emits `warning()`, does NOT `stop()`.
6. **AC6**: `status=2`, `status=3` still `stop()` in harvest.R. **Unchanged.**
7. **AC7**: `result$convergence_used$convergence_reason` is a non-NA character string for all exit paths including `status=1` (maps to `"legacy"`).
8. **AC8** *(local bench)*: On stepstone fulldata, flat raking with water-filling: `status=5` (KL stall) AND `max_err ≤ 2.97e-3` (at least as good as current errRp-stall exit of 2.98e-3).
9. **AC9**: `devtools::test()` FAIL ≤ 3 (pre-existing).

---

## Files Changed

| File | Change |
|------|--------|
| `src/leafblower.h` | Add `RK_ERR_BUDGET=4`, `RK_ERR_STALL=5` |
| `src/raking.cpp` | Replace errRp stall with weight KL stall; BUDGET at budget exit; ordering guarantee |
| `R/harvest.R` | Handle status 4/5 with warnings; add `convergence_reason` to result list; update `@return` |
| `tests/testthat/test-calibration-solvers.R` | Update 4 tests: `status %in% c(0L,1L)` → `c(0L,4L,5L)`, conditional `status==1L` → `%in% c(4L,5L)` |
| `tests/testthat/test-calib-linalg.R` | Same pattern update |
| `tests/testthat/test-convergence-criteria.R` | Same pattern update |
| `tests/testthat/test-ieppa-persistent-infeas.R` | Comment update only |

**Not in scope**: ieppa.cpp stall (requires separate spec — current ieppa stall is PCT-based at line 1188, not errRp-window based; separate high-risk change).

---

## Out of Scope

- ieppa.cpp stall criterion change (deferred — dedicated spec needed)
- sinkhorn stall criterion change
- SQUAREM stall on KL (CBB breaks monotonicity; errRp pragmatic and sufficient)
- Per-iteration loss trajectory exposed to R
