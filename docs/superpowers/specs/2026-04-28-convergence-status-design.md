# Convergence Status Redesign

**Date**: 2026-04-28
**Status**: Pending design review
**Files**: `src/leafblower.h`, `src/raking.cpp`, `src/ieppa.cpp`, `R/harvest.R`

---

## Problem

### Status 1 conflates two incompatible diagnostics

Current `RK_ERR_NOCONV = 1` fires when:
- Budget exhausted — solver was still improving; prescription: increase `max_iterations`
- Stall — loss function plateaued; prescription: done, optimum reached

These have **opposite prescriptions**. A user checking `status == 1` cannot know whether to add iterations or stop.

### Stall uses errRp, not the solver's loss function

Current stall detector (`min_errRp_window`, `kMaxNoImprove=5`): fires when errRp (max marginal error) fails to improve by ≥1% for 5 consecutive convergence-check intervals.

Problems:
1. **errRp is not monotone** near the constrained KL optimum. errRp oscillates as cyclic margin projections cycle through categories — it can stagnate or oscillate while weight KL is still strictly decreasing. Result: false stall detection, early exit before the true optimum.
2. **errRp measures the wrong thing**. Raking minimizes weight KL. IPF/water-filling is a coordinate descent on KL, monotone by the KL projection theorem (Csiszar-Tusnady 1984). Stalling on errRp misses this.
3. **Consequence on stepstone**: flat raking stalls at errRp≈2.98e-3 because errRp can't improve further in infeasible categories, while weight KL may still be descending in other categories.

**Mathematically correct stall**: `loss_function(iter k) - loss_function(iter k-1) < tol × loss_function(iter k-1)` where `loss_function` is:
- **raking**: weight KL = Σ_c X[c] log(X[c]/X_init[c]) / n
- **ieppa**: marginal KL = Σ_k Σ_j τ_kj log(τ_kj / achieved_kj) (already the default convergence metric)
- **sinkhorn**: weight KL

Weight KL is monotone for water-filling bounded IPF (each KL projection step decreases it; budget-exhausted iterates show it still descending). **Weight KL plateau = constrained KL minimum reached = stall is correct signal.**

---

## Solution

### New status codes (appended — no renumbering)

```c
// Existing (preserved for ABI/caller compatibility):
#define RK_OK           0   // Converged: improvement criterion satisfied
#define RK_ERR_NOCONV   1   // Legacy alias — mapped from old callers; not emitted by new solvers
#define RK_ERR_INFEAS   2   // Infeasible: empty cell with positive target (no solution exists)
#define RK_ERR_BADARG   3   // Invalid argument

// New (appended — no renumbering of existing codes):
#define RK_ERR_BUDGET   4   // Budget exhausted while loss was still decreasing; increase max_iterations
#define RK_ERR_STALL    5   // Loss function plateau — at constrained optimum; more iterations won't help
```

**ABI note**: `rk_result_t.status` is an int. Existing codes 0–3 are unchanged. New codes 4 and 5 are appended. External callers checking `status == 2` for INFEAS continue to work. Callers checking `status == 1` (NOCONV) will no longer receive it from new solvers — they receive 4 or 5 instead. Document this as a minor breaking change for `status == 1` users.

### Stall criterion: weight KL for raking, marginal KL for ieppa

Replace `min_errRp_window` + errRp-based stall with:

```cpp
// raking_solve flat loop and SQUAREM loop:
double wkl = compute_weight_kl();
if (!std::isfinite(min_loss_window)) {
    min_loss_window = wkl; n_no_improve = 0;
} else if (wkl < min_loss_window * (1.0 - st.convergence_cfg.pct_tol)) {
    min_loss_window = wkl; n_no_improve = 0;
} else {
    n_no_improve++;
}
if (n_no_improve >= kMaxNoImprove) { res.status = RK_ERR_STALL; break; }
```

Uses the same `pct_tol` as the primary improvement criterion — consistent tolerance.

**SQUAREM stall**: weight KL is NOT monotone for SQUAREM extrapolation steps (CBB can overshoot). Use errRp for SQUAREM stall (pragmatic — errRp of accepted iterates is approximately non-increasing with step-halving). Budget vs stall distinction still applies: loop-condition exit → BUDGET, n_no_improve → STALL.

**ieppa stall**: replace current errRp-based stall with marginal_kl-based stall (marginal_kl is already computed at every convergence check; use it for stall too).

### Budget vs stall at exit

```cpp
// Budget exhausted (while/for loop condition fails naturally):
res.status = RK_ERR_BUDGET;

// Stall detector fired:
res.status = RK_ERR_STALL;
```

Both are non-error exits: weights are valid (best iterate seen). Only INFEAS and BADARG are hard errors.

### harvest.R status handling

```r
# Hard stops (unchanged):
if (calib_result$status == 2L)  stop("leafblower: infeasible...")
if (calib_result$status == 3L)  stop("leafblower: invalid arguments...")

# Soft exits (warnings):
if (calib_result$status == 4L)
  warning("leafblower: budget exhausted — weights reflect best iterate; consider increasing max_iterations")
if (calib_result$status == 5L)
  warning("leafblower: loss function plateau — at constrained optimum given bounds and targets")

# Legacy status=1 (from old solvers not yet updated, or external C callers):
if (calib_result$status == 1L)
  warning("leafblower: did not converge (legacy status)")
```

New field in `result$convergence_used`:
```r
convergence_reason = c("criterion", "budget", "stall", "infeasible", "error")[
    match(calib_result$status, c(0L, 4L, 5L, 2L, 3L))]
```

---

## Mathematical justification of weight KL stall

For water-filling bounded IPF:
- Each category's KL projection satisfies: KL(X* || X_init) ≤ KL(X || X_init) − KL(X* || X) (Pythagorean theorem for KL)
- Weight KL strictly decreases with every F_eval call (until at the fixed point)
- Fixed point of bounded water-filling IPF = constrained KL minimum (Csiszar-Tusnady 1984)
- **Weight KL plateau ↔ fixed point reached ↔ constrained optimum**

For errRp (max marginal error):
- NOT monotone: as cyclic projections cycle through K margins, errRp on any single margin can increase temporarily while overall KL decreases
- errRp plateau ≠ optimum; may plateau mid-descent due to oscillation

**Practical difference on stepstone**: with errRp stall (current), flat raking exits at ≈2.98e-3 after 70 iterations because errRp oscillates. With weight KL stall, the solver should run longer (KL is still descending) and reach lower errRp before genuinely stalling.

---

## Acceptance Criteria

1. **AC1**: `status=0` fires only when primary convergence criterion (`pct_tol` improvement) satisfies. Unchanged.
2. **AC2**: `status=4` (BUDGET) fires when loop exhausts budget while loss function was still decreasing on the last check.
3. **AC3**: `status=5` (STALL) fires when weight KL (raking/sinkhorn) or marginal KL (ieppa) fails to improve by `pct_tol` fraction for `kMaxNoImprove` consecutive convergence-check intervals.
4. **AC4**: `status=4` and `status=5` both return valid weights (best iterate) and do NOT `stop()` in harvest.R.
5. **AC5**: `status=2` (INFEAS) and `status=3` (BADARG) still `stop()` in harvest.R. Unchanged.
6. **AC6**: `result$convergence_used$convergence_reason` exposed as character field: "criterion" | "budget" | "stall" | "infeasible" | "error".
7. **AC7**: On stepstone fulldata with flat raking+water-filling: `status=5` (stall at constrained optimum, better errRp than current errRp-stall exit).
8. **AC8**: `devtools::test()` FAIL ≤ 3 (pre-existing).

---

## Files Changed

| File | Change |
|------|--------|
| `src/leafblower.h` | Add `RK_ERR_BUDGET=4`, `RK_ERR_STALL=5`; update tripwire comment |
| `src/raking.cpp` | Replace errRp stall with weight KL stall; emit BUDGET vs STALL at exit |
| `src/ieppa.cpp` | Replace errRp stall with marginal KL stall; emit BUDGET vs STALL |
| `R/harvest.R` | Handle status 4/5 as warnings; add `convergence_reason` field |
| `tests/testthat/test-calibration-solvers.R` | Update tests checking status==1 to check status∈{4,5} |

No `rk_params_t` ABI change. `rk_result_t.status` int field semantics change (new values 4/5); `EXPECTED_RK_RESULT_BYTES` unchanged.

---

## Out of Scope

- sinkhorn stall (same pattern as raking but deferred)
- SQUAREM stall on marginal KL (errRp stall is pragmatic and sufficient for SQUAREM)
- Exposing per-iteration loss trajectory to R
