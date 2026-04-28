# SQUAREM Geometry Fix: Obs-Level α + Weight-Change Stall

**Date**: 2026-04-28
**Status**: Pending design review
**File**: `src/raking.cpp`

---

## Problem

SQUAREM raking with water-filling underperforms the flat loop on stepstone-fulldata:

| Method | F-evals | max_err | DEFF |
|--------|---------|---------|------|
| raking flat (wf, KL stall) | 50 | **1.60e-3** | 1.958 |
| raking+SQUAREM (wf, errRp stall) | 113 | 1.71e-3 | 1.957 |
| autumn reference | ~100F | 1.60e-3 | 1.995 |

SQUAREM uses 2.3× more work and gets a worse answer — the opposite of its purpose.

---

## Root Cause Analysis

### Cause 1 (Primary): CBB step computed in cell space — wrong geometry

The Cauchy-Barzilai-Borwein step `α = -‖r‖/‖v‖` uses cell-level L2 norms:

```
‖r_cell‖² = Σ_c (w1[c] - X[c])²   (current implementation)
```

The geometrically correct norm for obs-level IPF is:

```
‖r_obs‖² = Σ_i (Δw_i)² = Σ_c n_c × (r_c/n_c)² = Σ_c r_c² / n_c
```

Since all obs in cell c get identical multipliers, `r_obs[i] = r_cell[c(i)] / n_{c(i)}`. The cell-level norm weights each cell equally; the obs-level norm weights each cell by `1/n_c`, downweighting large cells.

**Effect**: for stepstone (cells range from n_c=1 to n_c≈hundreds), large cells dominate `‖r_cell‖²` but correctly contribute little to `‖r_obs‖²`. Cell-level α is systematically larger in magnitude (more aggressive) than obs-level α → over-extrapolation → step-halving fires → α collapses to ≈-1 → no acceleration.

Autumn computes norms at obs level (n=1.58M), giving the correct geometry. We compute at cell level (M_cell≈29k), giving incorrect geometry.

### Cause 2 (Secondary): errRp stall fires during fine-convergence oscillation

Near the constrained KL minimum, SQUAREM's accepted iterate errRp oscillates (large CBB steps produce erratic trajectories). The errRp stall fires after 5 consecutive non-improvements at errRp≈1.71e-3, before the true KL minimum (errRp=1.60e-3) is reached.

The flat loop's KL stall is immune to errRp oscillation (KL is monotone). Autumn uses obs-level weight change `Σ|Δw_i|` which also bypasses errRp oscillation.

### Cause 3 (Structural): Fast convergence → α≈-1 regardless

Flat loop converges in 50 iters (ρ≈0.89/iter). CBB gives α≈-1/ρ≈-1.12 — barely any extrapolation. 3 F-evals/super-step for 12% boost = slower than flat (1 F-eval, 12% reduction). SQUAREM acceleration only helps when ρ is close to 1 (α is large). Fixing cause 1 may reveal α IS larger at obs level, providing real acceleration.

---

## Solution

### Fix 1: Obs-level weighted norms for α

Replace cell-level L2 norms with obs-level weighted norms:

```cpp
// Current (wrong — cell space):
for (int c = 0; c < ct.M_cell; c++) {
    double ri = w1[c] - X[c], vi = w2[c] - w1[c];
    norm_r += ri * ri;
    norm_v += vi * vi;
}

// Fixed (correct — obs space via 1/n_c weighting):
for (int c = 0; c < ct.M_cell; c++) {
    const double inv_nc = 1.0 / static_cast<double>(ct.n_per_cell[c]);
    double ri = w1[c] - X[c], vi = w2[c] - w1[c];
    norm_r += ri * ri * inv_nc;
    norm_v += vi * vi * inv_nc;
}
```

Also fix `norm_w2` used in the ‖v‖/‖w2‖ fixed-point guard:
```cpp
norm_w2 += w2[c] * w2[c] * inv_nc;  // obs-level norm
```

**Why this is correct**: within-cell obs all have identical per-obs changes `r_c/n_c`, so obs-level squared norm contribution is `n_c × (r_c/n_c)² = r_c²/n_c`. Summing over cells gives the correct obs-level SQUAREM geometry.

### Fix 2: Replace errRp stall with obs-level weight-change stall for SQUAREM

After each accepted super-step, compute obs-level L1 weight change:

```cpp
// Obs-level L1 weight change per super-step (autumn's convergence criterion)
double wchange_sq = 0.0;
for (int c = 0; c < ct.M_cell; c++)
    wchange_sq += std::fabs(X[c] - X_prev_sq[c]) / static_cast<double>(ct.n_per_cell[c]);
wchange_sq /= static_cast<double>(st.n);  // normalize to per-obs
```

Stall when `wchange_sq < pct_tol * wchange_sq_initial` (relative improvement below threshold for kMaxNoImprove consecutive super-steps). Or use the same absolute threshold as pct_tol.

Store `X_prev_sq` (previous accepted iterate snapshot) at the start of each super-step. This requires one extra O(M_cell) copy per super-step.

**Why weight change instead of errRp**: weight change goes to zero at the fixed point regardless of errRp trajectory. SQUAREM's errRp oscillates near the constrained minimum; weight change does not.

### Fix 3 (optional): Debug logging

Add verbose output per super-step: α, n_halvings, cand_resid/plain_resid ratio. Enables empirical verification that fix 1 changes α as predicted.

---

## Interaction with step-halving

The step-halving criterion `‖F(X*)-X*‖² ≤ 1.01 × ‖v‖²` remains **unchanged**. Both sides should use cell-level norms for consistency (they are compared to each other, so the geometry cancels). Only the α computation changes.

The `kVNormRel` guard `‖v‖/‖w2‖ < 1e-10` uses the fixed-point guard — both sides should use **the same norm** (either both cell-level or both obs-level). After fix 1, use obs-level norms throughout the norm computation block.

---

## Acceptance Criteria

1. **AC1**: `accelerate=TRUE` on stepstone-fulldata achieves max_err ≤ 1.60e-3 (flat loop result) — SQUAREM must reach the same constrained minimum
2. **AC2**: α on a typical super-step is measurably larger in magnitude with obs-level norms than cell-level norms (log at verbose=1)
3. **AC3**: `devtools::test()` FAIL ≤ 3 (pre-existing)
4. **AC4**: `accelerate=FALSE` results unchanged (no change to flat loop)

---

## Files Changed

| File | Change |
|------|--------|
| `src/raking.cpp` | SQUAREM norm computation: `ri²/n_c` and `vi²/n_c`; weight-change stall |

No other files. No API change. No ABI change.

---

## Out of Scope

- ieppa SQUAREM (different solver structure)
- Changing step-halving criterion geometry
- Convergence criterion for flat loop (already KL stall, working correctly)
