# SQUAREM Geometry Fix: Obs-Level α + Weight-Change Stall

**Date**: 2026-04-28
**Status**: Pending design review (rev 2 — gate round-1 fixes)
**File**: `src/raking.cpp` (SQUAREM branch only — flat loop unchanged)

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

**Effect**: for stepstone (cells range from n_c=1 to n_c≈hundreds), large cells dominate `‖r_cell‖²` but correctly contribute little to `‖r_obs‖²`. Cell-level α is systematically larger in magnitude than obs-level α → over-extrapolation → step-halving fires → α collapses to ≈-1 → no acceleration.

Autumn computes norms at obs level (n=1.58M), giving the correct geometry.

### Cause 2 (Secondary): errRp stall fires during fine-convergence oscillation

Near the constrained KL minimum, SQUAREM's accepted iterate errRp oscillates. The errRp stall fires after 5 consecutive non-improvements at errRp≈1.71e-3, before the true KL minimum (errRp=1.60e-3). Weight change `Σ_i|Δw_i|` monotonically goes to zero at the fixed point, bypassing this oscillation.

---

## Solution

### Fix 1: Obs-level weighted norms for α ONLY

**Key constraint**: two separate norm computations are needed:
- `norm_v_obs` (obs-level, `1/n_c` weighted) — used for α computation
- `norm_v_cell` (cell-level, unweighted) — used for step-halving (both sides cell-level = dimensionally consistent)

```cpp
// Compute both norms in a single pass:
double r_sq_obs = 0.0, v_sq_obs = 0.0;   // obs-level (for α)
double r_sq_cell = 0.0, v_sq_cell = 0.0; // cell-level (for halving)
double norm_w2_sq = 0.0;

for (int c = 0; c < ct.M_cell; c++) {
    // n_per_cell[c] >= 1 guaranteed by build_cell_table (cell allocated only if obs exist).
    const double inv_nc = 1.0 / static_cast<double>(ct.n_per_cell[c]);
    const double ri = w1[c] - X[c];
    const double vi = w2[c] - w1[c];
    r_sq_obs   += ri * ri * inv_nc;   // obs-level: Σ r_c²/n_c
    v_sq_obs   += vi * vi * inv_nc;   // obs-level: Σ v_c²/n_c
    r_sq_cell  += ri * ri;            // cell-level: Σ r_c²
    v_sq_cell  += vi * vi;            // cell-level: Σ v_c²
    norm_w2_sq += w2[c] * w2[c] * inv_nc;  // obs-level for fixed-point guard
}
const double norm_r_obs  = std::sqrt(r_sq_obs);
const double norm_v_obs  = std::sqrt(v_sq_obs);
const double norm_v_cell = std::sqrt(v_sq_cell);
const double norm_w2     = std::sqrt(norm_w2_sq);

// α uses obs-level geometry (correct CBB step):
double alpha = std::max(kAlphaMin, -norm_r_obs / (norm_v_obs + kVNormEps));

// Fixed-point guard uses obs-level:
if (norm_v_obs / (norm_w2 + kVNormEps) < kVNormRel) { ... }

// Step-halving criterion uses cell-level (both sides cell-level = consistent):
const double plain_resid = v_sq_cell;  // ‖v‖² cell-level
// ... cand_resid computed from cell-level diff as before ...
for (int h = 0; h < kMaxHalvings && cand_resid > kHalvingSlack * plain_resid; h++) { ... }
```

**Why separate norms for halving**: `cand_resid = ‖F(X*)-X*‖²` is computed as a cell-level norm (diff of cell-level vectors). `plain_resid = ‖v_cell‖²` must also be cell-level so both sides are in the same space. Obs-level `norm_v_obs` is ONLY used for α.

**n_per_cell[c] zero guard**: `build_cell_table` allocates cell c only when at least one obs belongs to it, so `n_per_cell[c] >= 1` for all c in `[0, M_cell)`. No runtime guard needed, but the code comment above documents this invariant.

### Fix 2: Weight-change stall replaces errRp stall for SQUAREM

After each accepted super-step, compute obs-level L1 weight change and use the same relative-improvement window as the KL stall in the flat loop:

```cpp
// X_prev_sq: snapshot of X at start of this super-step (O(M_cell) copy before w1=F(X))
// Initialize BEFORE the while loop: auto X_prev_sq = X;
// Update AFTER each accepted step: X_prev_sq = X;

double wchange = 0.0;
for (int c = 0; c < ct.M_cell; c++)
    wchange += std::fabs(X[c] - X_prev_sq[c]) / static_cast<double>(ct.n_per_cell[c]);
wchange /= static_cast<double>(st.n);  // per-obs L1 weight change

// Stall tracking: same relative improvement pattern as KL stall in flat loop
if (!std::isfinite(min_loss_window)) {
    min_loss_window = wchange; n_no_improve = 0;
} else if (wchange < min_loss_window * (1.0 - st.convergence_cfg.pct_tol)) {
    min_loss_window = wchange; n_no_improve = 0;
} else {
    n_no_improve++;
}
if (n_no_improve >= kMaxNoImprove) { res.status = RK_ERR_STALL; break; }

// Update snapshot for next super-step
X_prev_sq = X;
```

**`wchange_sq_initial` / `min_loss_window` initialization**: `min_loss_window` is initialized to `+∞` (same as the flat loop KL stall) before the while loop. On the first super-step, `wchange < ∞` is always true → sets `min_loss_window = wchange`, resets n_no_improve=0. This means the stall fires only after kMaxNoImprove consecutive super-steps where wchange fails to improve by `pct_tol` fraction. No separate "initial" baseline needed — the same sliding-window relative improvement pattern used throughout.

### Code path isolation (AC4 guarantee)

All changes in Fix 1 and Fix 2 are **inside `if (st.accelerate) { if (st.inner_max_iter >= 3) { while (...) {`}}}**. The flat loop is in `else { for (int iter ...) }`. The norm computation at lines 276–285 (SQUAREM while loop body) is SQUAREM-only. The flat loop uses its own `compute_weight_kl()` stall — unchanged.

| Code change | Path(s) touched |
|-------------|----------------|
| r_sq_obs, v_sq_obs, inv_nc computation | SQUAREM only |
| alpha = -norm_r_obs / norm_v_obs | SQUAREM only |
| plain_resid = v_sq_cell | SQUAREM only |
| wchange stall logic | SQUAREM only |
| X_prev_sq snapshot | SQUAREM only |
| Flat loop KL stall | **Flat loop only — UNCHANGED** |

---

## RED Test (TDD prerequisite)

On a synthetic problem with heterogeneous cell sizes (n_c=1 and n_c=100), cell-level and obs-level α differ by a factor of ≥10. Commit this test BEFORE implementation — it will fail (RED) because current code uses cell-level α:

```r
test_that("squarem-geo-red: obs-level α differs from cell-level α on heterogeneous cells", {
  # 2-margin problem where one cell has n_c=100 and another n_c=1.
  # Cell-level: ‖r_cell‖² = r_big² + r_small² ≈ r_big² (big cell dominates)
  # Obs-level:  ‖r_obs‖² = r_big²/100 + r_small²/1 ≈ r_small² (small cell dominates)
  # → α_obs << α_cell in magnitude → fewer halvings → real acceleration.
  #
  # RED: with cell-level α (current), SQUAREM stalls at max_err > 1.70e-3.
  # GREEN: with obs-level α (fixed), SQUAREM reaches max_err ≤ 1.65e-3.
  #
  # Construct: v1 with 2 categories of very different sizes.
  set.seed(42L)
  big   <- rep("A", 500L)   # 500 obs → big cell
  small <- rep("B", 5L)     # 5 obs → small cell
  df <- data.frame(v1 = factor(c(big, small)))
  tgt <- list(v1 = c("A" = 0.4, "B" = 0.6))  # B heavily overweight → forces large IPF step

  r <- leafblower::harvest(df, tgt, method = "raking", accelerate = TRUE,
    max_weight = 5, max_iterations = 200L, attach_weights = FALSE)
  res <- attr(r, "result")

  # With cell-level α: B gets huge (r_B dominates in cell space but not obs space),
  # SQUAREM over-extrapolates → step-halving brings α → -1 → status=stall_errRp early.
  # With obs-level α: correct geometry → fewer halvings → status=stall_kl later.
  # Binary discriminator: convergence_reason should be "stall_kl" not "stall_errRp".
  expect_equal(res$convergence_used$convergence_reason, "stall_kl",
               label = "geometry-fixed SQUAREM must converge to KL stall, not errRp stall")
})
```

---

## Acceptance Criteria

1. **AC1** *(local benchmark)*: `accelerate=TRUE` on stepstone-fulldata achieves max_err ≤ 1.60e-3. Flagged local-only (requires stepstone parquet, not in CI). Log result in commit message. Add `skip_if` guarded test:
   ```r
   test_that("squarem-geo-ac1: stepstone SQUAREM reaches flat-loop quality", {
     skip_if(!file.exists("benchmarks/stepstone_fulldata_bench_data.parquet"),
             "stepstone dataset not available")
     # ... run + assert max_err ≤ 1.60e-3 ...
   })
   ```

2. **AC2** *(synthetic unit test)*: On a problem with heterogeneous cell sizes (n_c=1 and n_c=100), `|α_obs| ≤ 0.5 × |α_cell|` for the first super-step — obs-level α is at least 2× smaller in magnitude than cell-level α. Test via RED test above (GREEN after fix).

3. **AC3**: `devtools::test()` FAIL ≤ 3 (pre-existing).

4. **AC4**: `accelerate=FALSE` results bit-identical to pre-fix baseline. All code changes are inside the SQUAREM-only branch (table above). Verified by the existing `squarem-ac3` test (which tests `accelerate=FALSE`).

5. **AC5**: `convergence_reason` for stepstone SQUAREM changes from `"stall_errRp"` to `"stall_kl"` — confirming the weight-change stall runs longer and reaches the KL minimum.

---

## Files Changed

| File | Change |
|------|--------|
| `src/raking.cpp` | SQUAREM branch only: obs-level norms for α; cell-level norm_v_cell for halving; weight-change stall; X_prev_sq snapshot |
| `tests/testthat/test-calibration-solvers.R` | RED test + AC1 skip_if-guarded test |

No other files. No API change. No ABI change.

---

## Out of Scope

- ieppa SQUAREM
- Changing step-halving criterion geometry (stays cell-level, dimensionally consistent)
- Convergence criterion for flat loop (KL stall, working correctly)
- Sinkhorn
