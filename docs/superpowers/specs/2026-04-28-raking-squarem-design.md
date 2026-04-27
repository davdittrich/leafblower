# Raking SQUAREM Acceleration Design

**Date**: 2026-04-28
**Status**: Pending design review (rev 3 — gate round-2 fixes)
**File**: `src/raking.cpp`, `R/harvest.R`

---

## Problem

Leafblower raking with `max_weight=5` on stepstone-fulldata (n=1.58M, K=9) converges in ~60 iterations (Greedy) to max_err=3.48e-3 — but improvement criterion fires before reaching the true KL minimum. Autumn achieves max_err=1.60e-3 using SQUAREM acceleration (Varadhan & Roland 2008). The gap is 2.2×.

---

## Solution: SQUAREM SqS3 outer loop

SQUAREM wraps the existing inner sweep as a fixed-point operator F and applies Cauchy-Barzilai-Borwein extrapolation to jump toward the fixed point. Same algorithm as `autumn::harvest(accelerate=TRUE)`.

### F definition

One call to F = one complete inner iteration:
- K-margin IPF sweep (Bregman Dykstra scaling, Greedy-sorted if `scheduler="greedy"`)
- Bregman box correction (p[c] multiplicative)
- Bregman hyperplane correction (q_hyp multiplicative)

### SQUAREM SqS3 super-step

```
w1 = F(X)
w2 = F(w1)
r  = w1 - X      # fixed-point residual at X
v  = w2 - w1     # change in residual direction

# Guard: already at fixed point → skip extrapolation
if ‖v‖₂ / (‖w2‖₂ + 1e-300) < 1e-10:
    X ← w2; continue to next super-step

α  = -‖r‖₂ / ‖v‖₂   (CBB step; capped to [-1000, -1])

# Snapshot BEFORE extrapolation (state at w2, not post-X*)
X_snap = copy(X); p_snap = copy(p); q_snap = q_hyp

X*     = X - 2α·r + α²·v    (quadratic extrapolation)
X*     = max(X*, 0)           (clamp structural zeros)
X_new  = F(X*)

# Step-halving: compare against reference residual ‖r‖₂ = ‖w1-X‖₂
# (already computed above; no extra F call needed)
for up to 16 halvings:
    if ‖F(X_new) - X_new‖₂ ≤ 1.01 × ‖r‖₂:
        break  # accept X_new
    α ← (α - 1) / 2        # midpoint toward -1; converges to -1 for any α ≤ 0
    X = copy(X_snap); p = copy(p_snap); q_hyp = q_snap  # restore
    X* = X - 2α·r + α²·v; X* = max(X*, 0)
    X_new = F(X*)
if all 16 halvings exhausted and still failing: accept w2 (α = -1 → X* = w2)

X ← X_new
```

**Key design decisions:**

- **Reference residual for step-halving**: `‖r‖₂ = ‖w1-X‖₂` (already computed, no extra F call). This is the fixed-point residual at the super-step start. Alternative `‖F(w2)-w2‖₂` would require a 4th F call per super-step.
- **Halving formula**: `α ← (α-1)/2` is the midpoint between α and -1. This is a contraction mapping with fixed point -1: converges monotonically to -1 for all α ≤ 0. No branching needed.
- **‖v‖₂ threshold**: relative `‖v‖₂ / (‖w2‖₂ + ε) < 1e-10` avoids false triggers at large weight magnitudes (n=1.58M cells, values O(1), ‖w2‖₂ ~ O(√n)).
- **α cap at -1**: if α from CBB formula is in (-1, 0), cap α = -1 (gentler than CBB, avoids extrapolating past w2).

**Cost per super-step**: 3 F-evaluations (F(X), F(w1), F(X*)) + 1 per halving iteration (F(X*)). Baseline (no halving) = 3 F-evals. Total = 3 + k where k = halvings triggered (typically 0).

### Snapshot timing

**Store snapshots of X[], p[], q_hyp AFTER computing w2 and AFTER the ‖v‖₂ guard, BEFORE computing X***. This is the correct restore point: on backtrack, we re-derive X* from the same X, r, v (which are fixed for this super-step) with a new α.

```cpp
// CORRECT order in pseudocode:
//   1. w1 = F(X); w2 = F(w1)
//   2. r = w1-X; v = w2-w1
//   3. if relative ‖v‖₂ guard → skip
//   4. compute α
//   5. SNAPSHOT HERE: X_snap=X; p_snap=p; q_snap=q_hyp
//   6. X* = extrapolate; X_new = F(X*)
//   7. step-halving loop (restores from snap, recomputes X* with new α)
```

### Integration

- **Greedy**: runs inside F. errRp_k from the final F(X*) used for next super-step's sort.
- **SOR**: not activated by default in accelerate mode (zero empirical benefit). Fires if `sor=list(auto=TRUE)` passed.
- **accelerate=TRUE with non-raking methods**: emits `warning("accelerate=TRUE is only supported for method='raking'; ignoring")` and falls back to standard solver.

---

## RED Test (TDD prerequisite)

`harvest()` uses `...` and passes unrecognized arguments to the C++ dispatch; they are currently silently ignored. Therefore `accelerate=TRUE` currently runs identically to `accelerate=FALSE`. The RED test exploits this:

```r
test_that("squarem-red: accelerate=TRUE produces different iterations than accelerate=FALSE", {
  # RED before SQUAREM: accelerate=TRUE silently ignored → same weights as FALSE.
  # GREEN after SQUAREM: accelerate=TRUE takes different (larger) steps → different weights.
  # Verified empirically: n=2000, K=2, max_weight=5 does NOT converge in 1 iteration,
  # so SQUAREM will produce different path than standard.
  set.seed(42L); n <- 2000L
  df  <- data.frame(v1 = factor(sample(3L, n, TRUE)), v2 = factor(sample(2L, n, TRUE)))
  tgt <- list(v1=c("1"=0.5,"2"=0.3,"3"=0.2), v2=c("1"=0.6,"2"=0.4))

  r_base <- leafblower::harvest(df, tgt, method="raking",
    accelerate=FALSE, max_weight=5, max_iterations=500L, attach_weights=FALSE)
  r_acc  <- leafblower::harvest(df, tgt, method="raking",
    accelerate=TRUE,  max_weight=5, max_iterations=500L, attach_weights=FALSE)

  iters_base <- attr(r_base, "result")$iterations
  iters_acc  <- attr(r_acc,  "result")$iterations

  # Before SQUAREM: same iterations (accelerate=TRUE ignored)
  # After  SQUAREM: fewer iterations (SQUAREM reaches convergence faster)
  expect_false(isTRUE(all.equal(iters_base, iters_acc)),
               label="accelerate=TRUE must use different (fewer) iterations than accelerate=FALSE")
})
```

Using `iterations` rather than weight vectors is more stable: even if SQUAREM reaches the same fixed point, it will do so in fewer super-steps, making `iters` the reliable discriminator.

---

## API

```r
harvest(df, tgt, method="raking", accelerate=FALSE, max_weight=5, ...)

#' @param accelerate Logical. If \code{TRUE}, applies SQUAREM SqS3 outer-loop
#'   acceleration (Varadhan & Roland 2008) to the raking fixed-point iteration.
#'   Only supported for \code{method="raking"}; silently ignored (with a warning)
#'   for all other methods. Default \code{FALSE} preserves pre-SQUAREM behavior.
```

---

## Acceptance Criteria

1. **AC1**: `accelerate=TRUE, method="raking"` runs without error
2. **AC2** *(local benchmark, manual gate)*: On stepstone-fulldata, `accelerate=TRUE` achieves max_err ≤ 2.50e-3. This must be verified locally before merge; it is not runnable in CI (dataset not committed). Log the result in the PR description.
3. **AC3**: `accelerate=FALSE` is bit-identical to pre-SQUAREM Approach C results. **Fixture**: generate `tests/testthat/fixtures/raking_squarem_baseline.rds` by running `accelerate=FALSE` BEFORE implementing SQUAREM. The test asserts `all.equal(w_new, w_baseline, tolerance=0)`.
4. **AC4**: Step-halving correctly restores X[], p[], q_hyp from snapshot taken after F(w2) (before extrapolation)
5. **AC5**: ‖v‖₂ / (‖w2‖₂ + ε) < 1e-10 triggers early return without NaN/Inf
6. **AC6**: `accelerate=TRUE` with `method="ieppa"` emits a warning and does not error
7. **AC7**: RED test passes after implementation; `devtools::test()` FAIL ≤ 2 (pre-existing)

---

## A/B Test Expectation

| Config | Expected max_err |
|--------|-----------------|
| Greedy only (baseline) | 3.48e-3 |
| SQUAREM + Greedy | ≤ 2.50e-3 |
| autumn (reference) | 1.60e-3 |

---

## Files Changed

| File | Change |
|------|--------|
| `src/types.hpp` | Add `bool accelerate = false` to CalibState (C++ only; not rk_params_t) |
| `src/r_bridge.cpp` | Add 31st SEXP arg (`accelerate_sexp`); update `R_registerRoutines` count 30→31 |
| `src/raking.cpp` | SQUAREM outer loop; read `st.accelerate` |
| `R/harvest.R` | Update `@param accelerate`; remove from ignored list; warn if non-raking; pass to `.Call` |
| `tests/testthat/test-calibration-solvers.R` | RED test + AC3 fixture test + AC4 step-halving test |
| `tests/testthat/fixtures/raking_squarem_baseline.rds` | Generated before implementation |

**ABI note**: `rk_params_t` (public C API in `leafblower.h`) is NOT modified. The `.Call` registration count changes from 30→31 — this is a package-internal R/C bridge change, not a public C ABI change. No `static_assert` tripwire needs updating.

---

## Out of Scope

- SQUAREM for ieppa (different solver structure)
- Nesterov acceleration (requires lf[] tracking not yet in raking)
- Python wrapper changes
