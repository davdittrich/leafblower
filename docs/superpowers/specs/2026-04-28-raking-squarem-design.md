# Raking SQUAREM Acceleration Design

**Date**: 2026-04-28
**Status**: Pending design review
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
r   = w1 - X
v   = w2 - w1
α   = -‖r‖₂ / ‖v‖₂      (CBB step; capped at -1000)
X*  = X - 2α·r + α²·v    (quadratic extrapolation)
X*  = max(X*, 0)           (clamp: structural zeros cannot go negative)
X_new = F(X*)

# Step-halving safety (up to 16 halvings):
if ‖residual(X_new)‖ > 1.01 × ‖residual(w2)‖:
    halve α toward -1; recompute X* and X_new
```

Cost per super-step: 3–4 F-evaluations + O(M_cell) arithmetic.

### Integration

- **Greedy**: runs inside F. Greedy sorts margins by per-margin errRp_k before each F call. The errRp_k from the final F(X*) is used for the next super-step's sort.
- **SOR**: kept wired but **not activated by default** in accelerate mode (zero empirical effect on stepstone). If `sor=list(auto=TRUE)` passed, SOR fires inside each F-call as usual.
- **Bregman Dykstra**: already in F; no change. p[] and q_hyp state saved/restored when stepping back during step-halving.
- **p[], q_hyp save/restore**: step-halving needs to restore X[], p[], q_hyp to the pre-extrapolation state before re-trying with α/2. Store snapshots at the start of each super-step.

### Convergence criterion

Same as non-accelerated raking (improvement criterion on errRp). SQUAREM's larger steps mean the criterion fires much later = closer to the true KL minimum before stopping.

---

## API

```r
harvest(df, tgt, method="raking", accelerate=TRUE, max_weight=5, ...)
```

- `accelerate=FALSE` (default): current behavior unchanged
- `accelerate=TRUE`: SQUAREM outer loop; inner F uses same parameters (sor=, scheduler=, convergence=)
- Mirrors autumn's `accelerate=TRUE` API

---

## Acceptance Criteria

1. **AC1**: `accelerate=TRUE` for raking runs without error
2. **AC2**: On stepstone-fulldata, `accelerate=TRUE` achieves max_err ≤ 3.48e-3 (at least as good as Greedy alone) and ideally closes gap toward autumn's 1.60e-3
3. **AC3**: `accelerate=FALSE` (default) is bit-identical to pre-SQUAREM Approach C results
4. **AC4**: Step-halving correctly restores X[], p[], q_hyp on backtrack
5. **AC5**: `devtools::test()` FAIL ≤ 2 (pre-existing)

---

## Files Changed

| File | Change |
|------|--------|
| `src/raking.cpp` | SQUAREM outer loop wrapping existing inner sweep |
| `R/harvest.R` | Add `accelerate` parameter; wire to raking dispatch |

No new structs, no ABI change, no other files.

---

## Out of Scope

- SQUAREM for ieppa (different solver structure, different failure mode)
- Nesterov acceleration (requires lf[] tracking infrastructure not yet in raking)
- Python wrapper changes
