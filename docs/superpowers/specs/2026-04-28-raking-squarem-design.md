# Raking SQUAREM Acceleration Design

**Date**: 2026-04-28
**Status**: Pending design review (rev 2 — gate round-1 fixes)
**File**: `src/raking.cpp`, `R/harvest.R`

---

## Problem

Leafblower raking with `max_weight=5` on stepstone-fulldata (n=1.58M, K=9) converges in 60 iterations (Greedy) to max_err=3.48e-3 — but improvement criterion fires before reaching the true KL minimum. The trajectory stalls. Autumn achieves max_err=1.60e-3 using SQUAREM acceleration (Varadhan & Roland 2008). The gap is 2.2×.

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
w1  = F(X)
w2  = F(w1)

# Guard: if w1 == w2 (‖v‖₂ < ε), already at fixed point → return w2
r   = w1 - X
v   = w2 - w1
if ‖v‖₂ < 1e-15: return w2

α   = -‖r‖₂ / ‖v‖₂      (CBB step; capped at -1000)
X*  = X - 2α·r + α²·v    (quadratic extrapolation)
X*  = max(X*, 0)           (clamp: structural zeros cannot go negative)

# Snapshot X[], p[], q_hyp HERE (after w2, before extrapolation)
X_new = F(X*)

# Step-halving safety (up to 16 halvings):
# Residual = ‖F(X)-X‖₂ = ‖r‖₂  (fixed-point residual from first F call this super-step)
# Reference residual for w2 = ‖F(w2)-w2‖₂ (computed once before halving loop)
if ‖F(X_new)-X_new‖₂ > 1.01 × ‖F(w2)-w2‖₂:
    halve α toward -1 (α ← (α-1)/2 + (-1 if α < -1 else 1))
    restore X[], p[], q_hyp from snapshot
    recompute X* = X - 2α·r + α²·v; X* = max(X*, 0); X_new = F(X*)
    repeat up to 16 times; if still diverging after 16 halvings, accept w2
```

**Residual definition**: `‖F(X)-X‖₂` — the Euclidean norm of the fixed-point residual for the current iterate X. For the halving check, `‖F(X_new)-X_new‖₂` = one extra F evaluation. The reference for w2 is computed once as `‖F(w2)-w2‖₂` = `‖w1_next - w2‖₂` where w1_next is obtained by calling F(w2) once.

Cost per super-step: 3–4 F-evaluations + O(M_cell) arithmetic.

### Snapshot timing (B1 fix)

**Store snapshots of X[], p[], q_hyp AFTER computing w2 = F(w1)**, immediately before the extrapolation step. This is the correct restore point for step-halving: we restore to the state after w2 and re-try with a smaller α, without re-running F twice.

```cpp
// Correct snapshot point (AFTER w2 computed, BEFORE extrapolation):
auto X_snap    = X;
auto p_snap    = p;
auto q_snap    = q_hyp;

// ... extrapolation and step-halving use these snapshots ...
```

Snapshots at the start of the super-step (before w1) are WRONG — they would force re-running two F calls per halving step, tripling cost.

### Integration

- **Greedy**: runs inside F. Greedy sorts margins by per-margin errRp_k before each F call. The errRp_k from the final F(X*) is used for the next super-step's sort.
- **SOR**: kept wired but **not activated by default** in accelerate mode (zero empirical effect on stepstone). If `sor=list(auto=TRUE)` passed, SOR fires inside each F-call as usual.
- **Bregman Dykstra**: already in F; no change. p[] and q_hyp state saved/restored when stepping back during step-halving.
- **‖v‖₂ guard**: If ‖v‖₂ < 1e-15 (w1 ≈ w2 — already at fixed point), skip extrapolation and return w2. This avoids division-by-zero in α computation and terminates the outer loop immediately.

### Convergence criterion

Same as non-accelerated raking (improvement criterion on errRp). SQUAREM's larger steps mean the criterion fires much later = closer to the true KL minimum before stopping.

---

## RED Test (TDD prerequisite)

Write and commit BEFORE implementation. The test must be RED (error) before implementation and GREEN after.

```r
test_that("squarem-red: accelerate=TRUE signals not-yet-implemented before SQUAREM is added", {
  # This test is RED before SQUAREM is wired up:
  #   raking_solve() ignores accelerate=TRUE silently, so harvest() returns a result
  #   rather than erroring. The test below will pass once SQUAREM is implemented
  #   (no error) and fail (RED) before — confirming TDD prerequisite.
  #
  # The RED phase: calling accelerate=TRUE currently produces the same result as
  # accelerate=FALSE (parameter silently ignored). We assert the result is DIFFERENT
  # from accelerate=FALSE on a convergeable problem — which cannot pass before
  # SQUAREM is implemented.
  set.seed(42L); n <- 500L
  df  <- data.frame(v1 = factor(sample(3L, n, TRUE)), v2 = factor(sample(2L, n, TRUE)))
  tgt <- list(v1=c("1"=0.5,"2"=0.3,"3"=0.2), v2=c("1"=0.6,"2"=0.4))

  r_base <- leafblower::harvest(df, tgt, method="raking",
    accelerate=FALSE, max_weight=5, max_iterations=500L, attach_weights=FALSE)
  r_acc  <- leafblower::harvest(df, tgt, method="raking",
    accelerate=TRUE,  max_weight=5, max_iterations=500L, attach_weights=FALSE)

  w_base <- as.numeric(r_base)
  w_acc  <- as.numeric(r_acc)

  # Before SQUAREM: accelerate=TRUE is silently ignored → identical weights
  # After  SQUAREM: accelerate=TRUE takes different steps → different weights
  expect_false(isTRUE(all.equal(w_base, w_acc, tolerance=1e-8)),
               label="accelerate=TRUE must produce different weights than accelerate=FALSE (SQUAREM must fire)")
})
```

This is provably RED before implementation: with no SQUAREM wiring, `accelerate=TRUE` is silently ignored and the result equals `accelerate=FALSE`, causing `expect_false(all.equal(...))` to fail.

---

## API

```r
harvest(df, tgt, method="raking", accelerate=TRUE, max_weight=5, ...)
```

- `accelerate=FALSE` (default): current behavior unchanged
- `accelerate=TRUE`: SQUAREM outer loop; inner F uses same parameters (sor=, scheduler=, convergence=)
- Mirrors autumn's `accelerate=TRUE` API

---

## A/B Test Expectation

Benchmark on `stepstone-fulldata` (n=1.58M, K=9, max_weight=5):

| Config | Expected max_err |
|--------|-----------------|
| Greedy only (baseline) | 3.48e-3 |
| SQUAREM + Greedy | ≤ 2.50e-3 |
| autumn (reference) | 1.60e-3 |

---

## Acceptance Criteria

1. **AC1**: `accelerate=TRUE` for raking runs without error
2. **AC2**: On stepstone-fulldata, `accelerate=TRUE` achieves max_err ≤ 2.50e-3 (hard threshold; must beat Greedy-only baseline of 3.48e-3 by at least 28%)
3. **AC3**: `accelerate=FALSE` (default) is bit-identical to pre-SQUAREM Approach C results
4. **AC4**: Step-halving correctly restores X[], p[], q_hyp on backtrack (snapshots taken after F(w2))
5. **AC5**: ‖v‖₂ < 1e-15 guard triggers early return (no NaN/Inf in α)
6. **AC6**: RED test passes after implementation; `devtools::test()` FAIL ≤ 2 (pre-existing)

---

## Files Changed

| File | Change |
|------|--------|
| `src/raking.cpp` | SQUAREM outer loop wrapping existing inner sweep |
| `R/harvest.R` | Add `accelerate` parameter; wire to raking dispatch |
| `tests/testthat/test-calibration-solvers.R` | RED test (written before implementation) |

No new structs, no ABI change, no other files.

---

## Out of Scope

- SQUAREM for ieppa (different solver structure, different failure mode)
- Nesterov acceleration (requires lf[] tracking infrastructure not yet in raking)
- Python wrapper changes
