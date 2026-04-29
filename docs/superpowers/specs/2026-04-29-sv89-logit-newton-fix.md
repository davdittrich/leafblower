# sv89: Logit Newton Stabilization — Design Spec

**Date**: 2026-04-29
**Status**: Pending design review
**Ticket**: leafblower-sv89
**Files**: `src/logit_calib.cpp`, `tests/testthat/test-calibration-solvers.R`

---

## Problem

`method="logit"` diverges on large-K problems (stepstone: K=9, n=1.58M).
Benchmark: max_err=1.0, DEFF=527,577, ESS=3 — completely miscalibrated despite status=0.

**Root cause**: degenerate fixed-point convergence, not numerical error.

1. **Cold-start saturation**: λ=0 init gives w[c] = midpoint of [L,U] = 2.5·n_per_cell,
   far from calibrated solution. First Newton step is enormous.
2. **Saturation cascade**: Large Δλ → σ(z_c) → 0 or 1 for many cells →
   D_eff[c] = (U−L)·σ·(1−σ) → 0 → N nearly singular →
   LDLT Tikhonov regularization dominates → Δλ ≈ 0 → check_convergence fires.
3. **Convergence to wrong point**: Solver converges in 9 steps — to a degenerate
   solution where weights pile at bounds. The calibration error (max_err=1.0) is
   irrelevant to the convergence criterion being checked (||Δw|| < tol).

The identical root cause afflicts both the initialization quality AND the Newton step
size. Survey calibration literature (Deville-Sarndal, Kott 2006) prescribes both.

---

## Design: Three-Layer Fix

### Layer 1 — Armijo Line Search (primary fix)

After `ldlt_solve` gives Δλ in `b[]`, perform backtracking search for step size α ∈ (0,1]:

```
Acceptance condition: ||b_trial(λ + α·Δλ)||² < ||b_current||² × (1 - c·α)
  where c = 0.01 (Armijo sufficient-decrease constant)

Algorithm:
  alpha = 1.0
  for halv in 0..kMaxHalvings=10:
    compute w_trial[c] from lambda + alpha * delta_lambda (logit link)
    compute b_trial[k][j] = tau[k][j]*n - sum_{c in bucket} w_trial[c]
    if ||b_trial||² < ||b_current||² * (1 - 0.01 * alpha): break
    alpha *= 0.5
  lambda += alpha * delta_lambda
```

**Theoretical basis**: Newton with Armijo is globally convergent for strongly-convex
objectives. The D-S logit dual is strongly convex when D_eff[c] > 0. Armijo enforces
this: if a large step saturates cells (D_eff→0 next iter), the residual norm increases
and the step is halved until it doesn't. `survey::calib(method="logit")` uses step
halving and converges reliably. [Sufficient decrease guarantees from Nocedal-Wright §3.1]

**Cost per Newton step**: O(M_cell × K) per halving trial (recompute w_trial and
b_trial). Stepstone: ~1.8M ops per halving × max 10 halvings × ~9 Newton steps =
negligible vs total solver cost.

### Layer 2 — Design-Weight Initialization (cold-start prevention)

Instead of λ=0 (midpoint warm-start), initialize λ₀ so that σ(z_c) approximates
the normalized design weight position:

```
For each cell c:
  sigma_target[c] = clip((X_init[c] - L_cell[c]) / (U_cell[c] - L_cell[c]), ε, 1-ε)
  z_target[c] = log(sigma_target / (1 - sigma_target))  [logit transform]
  ε = 1e-4  (prevents log(0))

Find λ₀ minimizing ||A^T λ - z_target||²:
  Normal equations: (AA^T) λ₀ = A z_target
  where (Az_target)[k][j] = Σ_{c ∈ bucket(k,j)} z_target[c]
  and (AA^T) = compute_normal_equations(ct, D_eff=1, ...)  [reuse calib_linalg.hpp!]
```

This places the starting point where design weights match the logit link output —
Newton then operates within the convergence basin where D_eff[c] is well-conditioned.

**Degenerate cases**:
- L_cell[c] == U_cell[c]: set z_target[c] = 0 (contributes nothing to Az)
- X_init[c] < L_cell[c]: clip to L_cell[c] + ε·(U-L)
- X_init[c] > U_cell[c]: clip to U_cell[c] - ε·(U-L)

**Implementation**: Reuses existing `compute_normal_equations` and `ldlt_factor_inplace`
+ `ldlt_solve` from calib_linalg.hpp. One-time cost: same as one Newton step.

### Layer 3 — D_eff Floor (safety net)

```cpp
constexpr double kDeffEps = 1e-6;
double range = U_cell[c] - L_cell[c];
D_eff[c] = std::max(kDeffEps * range, range * sig * (1.0 - sig));
```

Prevents exact singularity of N if layers 1+2 still encounter ill-conditioning on
pathological inputs. Bias is O(1e-6) — negligible for all practical calibration.

**Why this is sound at 1e-6**: The minimum meaningful D_eff at non-saturated cells
is (U-L) × 0.5 × 0.5 = (U-L)/4. The floor at 1e-6×(U-L) is 4×10⁶ times smaller —
activated only when σ is within 1e-6 of 0 or 1, i.e., |z_c| > ~13.8. Cells this
saturated are genuinely at bounds and contribute nothing to calibration.

---

## Calibration of constants

| Constant | Value | Derivation |
|---|---|---|
| `c_armijo` | 0.01 | Standard Armijo constant; small → accepts many steps, large → forces improvement |
| `kMaxHalvings` | 10 | α_min = 2⁻¹⁰ ≈ 1e-3; below this, accept the step regardless (prevents infinite loop) |
| `kMaxDeltaZ` | 2.0 | Max z_c shift per Newton step; limits α_max before Armijo loop; logit(z+2) from 0.5 → 0.88 (safe range) |
| `kDeffEps` | 1e-6 | Bias < 1e-6 × (max_weight - min_weight) per cell — below double precision noise |
| `kInitSigmaEps` | 1e-4 | Clips σ_target away from 0/1; logit(1e-4) ≈ -9.2 (bounded z_target) |

## Fixes for design-review-gate findings

### Security Fix 1 — Step-norm guard before Armijo (prevents catastrophic fallthrough)

Even at α=2⁻¹⁰, a huge Newton step (large ||Δλ||) can produce a damaging update. Add a max-z-shift cap computed from the delta-z effect on each cell:

```cpp
// ORDERING: this block executes AFTER ldlt_solve (b[] now holds delta_lambda, not residuals)
// Cap alpha_max so no cell's z_c shifts by more than kMaxDeltaZ = 2.0 in one step
// delta_z[c] = sum_k delta_lambda[cat_offset[k] + g_per_cell[k][c]]
double max_delta_z = 0.0;
for (int c = 0; c < M; c++) {
    double dz = 0.0;
    for (int k = 0; k < K; k++) {
        int g = ct.g_per_cell[k][c];
        if (g >= 0 && g < st.cat_counts[k]) dz += std::abs(b[cat_offset[k]+g]);
    }
    max_delta_z = std::max(max_delta_z, dz);
}
double alpha_max = (max_delta_z > 0.0) ? std::min(1.0, kMaxDeltaZ / max_delta_z) : 1.0;
// Armijo searches in [0, alpha_max] instead of [0, 1]
double alpha = alpha_max;  // start from alpha_max
for (int halv = 0; halv < kMaxHalvings; halv++) { ... }
```

### Security Fix 2 — Pre-allocate Armijo scratch buffers

`w_trial` (size M_cell) and `b_trial` (size nct) must be allocated ONCE before the Newton loop, not per-halving. At M_cell=1.58M: per-halving alloc = 12MB × 10 halvings × 9 Newton steps = 1.08GB transient allocation. Pre-allocation cost: 12MB once.

```cpp
// Before Newton loop:
std::vector<double> w_trial(ct.M_cell, 0.0);   // Armijo trial weights
std::vector<double> b_trial(nct, 0.0);           // Armijo trial residuals
std::vector<double> lambda_trial(nct, 0.0);      // Trial lambda for Armijo eval
```

### Designer Fix 1 — Warning when Armijo stalls (all halvings exhausted)

When all 10 halvings fail (no sufficient decrease found), accept α (still makes a step) but populate the solver message:

```cpp
// After Armijo loop, if no improvement was found:
bool armijo_improved = (best_resid_sq < resid_sq_0);
if (!armijo_improved && st.verbose >= 1) {
    char msg[128];
    std::snprintf(msg, sizeof(msg),
        "[logit] Newton step: Armijo exhausted (alpha=%.2e), accepting best available",
        alpha);
    st.log(msg);
}
```

On STALL/BUDGET exit, the solver message includes this diagnostic.

### Designer Fix 2 — Convergence criterion is alpha-independent (clarification)

`check_convergence` from calib_dispatch.hpp evaluates `CellMetrics` (errRp, kl, chi2 etc.) computed from the CURRENT w[c] — not from ||Δw||. These calibration-quality metrics are alpha-independent: they measure how well margins are satisfied, regardless of step size. So:

- THRESHOLD rule: converges when max_err < tol. If alpha=2^-10 gave a tiny step but max_err is still 1.0, threshold doesn't fire. Correct.
- IMPROVEMENT rule: converges when metric barely improved. Small alpha → small improvement → fires. This IS correct behavior: IMPROVEMENT convergence on a tiny step means further Newton iterations won't help (solver has stalled). Status will be RK_ERR_STALL.

**No code change needed.** The convergence logic is alpha-independent by design. This finding is resolved by documentation, not code.

Note: The `prev_metric` update must happen AFTER the alpha-scaled step (current code already does this). If alpha=2^-10 gave a tiny improvement, the IMPROVEMENT rule correctly identifies a stall.

---

## Modified Newton Loop Structure

```cpp
// Layer 2: initialization (before main loop)
{
    std::vector<double> z_target(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) {
        if (U_cell[c] <= L_cell[c] + 1e-12) { z_target[c] = 0.0; continue; }
        double sig0 = (X_init[c] - L_cell[c]) / (U_cell[c] - L_cell[c]);
        sig0 = std::clamp(sig0, kInitSigmaEps, 1.0 - kInitSigmaEps);
        z_target[c] = std::log(sig0 / (1.0 - sig0));
    }
    // Az_target: sum z_target over each bucket
    std::fill(b_init.begin(), b_init.end(), 0.0);
    for (int k = 0; k < K; k++)
        for (int j = 0; j < st.cat_counts[k]; j++)
            for (int c : cells_per_cat[k][j]) b_init[cat_offset[k]+j] += z_target[c];
    // (AA^T) lambda_0 = b_init, using D_eff=1 everywhere
    std::fill(D_ones.begin(), D_ones.end(), 1.0);
    compute_normal_equations(ct, D_ones.data(), N.data(), cat_offset.data(), K, (size_t)nct);
    ldlt_factor_inplace(N.data(), (size_t)nct, 1e-10);
    ldlt_solve(N.data(), (size_t)nct, b_init.data());  // b_init = lambda_0
    lambda = b_init;
}

// Main Newton loop
for (int iter = 0; iter < kMaxNewtonIters; iter++) {

    // Layer 3: D_eff floor applied in sigma computation
    for (int c = 0; c < ct.M_cell; c++) {
        // ... z computation ... sigma computation ...
        D_eff[c] = std::max(kDeffEps * range, range * sig * (1.0 - sig));
    }

    // Residuals, build N, factor, solve -> b = delta_lambda

    // Layer 1: Armijo line search
    double resid_sq_0 = 0.0;
    for (double bj : b_resid) resid_sq_0 += bj * bj;  // ||b_current||²

    double alpha = 1.0;
    for (int halv = 0; halv < kMaxHalvings; halv++) {
        // Compute w_trial from lambda + alpha * b (delta_lambda)
        // Compute b_trial residuals
        double resid_sq_trial = compute_resid_sq(b_trial);
        if (resid_sq_trial < resid_sq_0 * (1.0 - kArmijoC * alpha)) break;
        alpha *= 0.5;
    }
    for (int j = 0; j < nct; j++) lambda[j] += alpha * b[j];  // b = delta_lambda

    // check_convergence (unchanged)
}
```

---

## TDD Requirements

### T_logit_armijo — Armijo prevents divergence on large-K synthetic problem

```r
test_that("T_logit_armijo: logit converges on K=5 tight-bound problem", {
  # Problem designed to stress Newton: 5 margins, many conflicting constraints
  set.seed(42); n <- 20000L
  df <- data.frame(
    a = factor(sample(letters[1:4], n, TRUE)),
    b = factor(sample(LETTERS[1:5], n, TRUE)),
    c = factor(sample(c("x","y","z"), n, TRUE)),
    d = factor(sample(c("M","F"), n, TRUE)),
    e = factor(sample(c("Y","O"), n, TRUE))
  )
  tgt <- list(
    a = setNames(c(0.3,0.2,0.3,0.2), letters[1:4]),
    b = setNames(rep(0.2,5), LETTERS[1:5]),
    c = c(x=0.4,y=0.35,z=0.25),
    d = c(M=0.48,F=0.52),
    e = c(Y=0.55,O=0.45)
  )
  r <- harvest(df, tgt, method="logit", max_weight=4.0, min_weight=0.1)
  me <- attr(r,"result")$max_error
  expect_lt(me, 1e-3, label=sprintf("logit K=5 must converge: got max_err=%.2e", me))
  # DEFF must be in reasonable range (not 527k like pre-fix)
  w <- r$weights
  expect_true(max(w) <= 4.0 + 1e-9)
  expect_true(min(w) >= 0.1 - 1e-9)
})
```

### T_logit_init — Design-weight initialization places lambda in convergence basin

This is tested implicitly by the above (both A+B are active). A diagnostic test:

```r
test_that("T_logit_init: logit converges in few Newton steps from design-weight init", {
  set.seed(7); n <- 5000L
  df <- data.frame(
    x = factor(sample(letters[1:3],n,TRUE)),
    y = factor(sample(c("M","F"),n,TRUE))
  )
  tgt <- list(x=c(a=0.3,b=0.4,c=0.3), y=c(M=0.5,F=0.5))
  r <- harvest(df, tgt, method="logit", max_weight=3.0, min_weight=0.1,
               convergence=list(absolute=1e-8))
  n_iters <- attr(r,"result")$iterations
  me <- attr(r,"result")$max_error
  expect_lt(me, 1e-8)
  # Good initialization → should converge in < 20 Newton steps even with tight tol
  expect_lt(n_iters, 20L,
    label=sprintf("design-weight init should converge fast: got %d steps", n_iters))
})
```

### T5-T8 regression guard
Existing T5-T8 tests must remain GREEN (bounds, quality, backward compat).

---

## Acceptance Criteria

| # | Criterion | How to verify |
|---|-----------|---------------|
| AC1 | T_logit_armijo GREEN: K=5 tight-bound problem, max_err < 1e-3 | `devtools::test(filter="calibration-solvers")` |
| AC2 | T_logit_init GREEN: converges in < 20 Newton steps with absolute=1e-8 | Same |
| AC3 | T5-T8 still GREEN | Same |
| AC4 | FAIL count unchanged (3) | `devtools::test()` |
| AC5 (benchmark-only) | Stepstone logit: max_err < 0.01, DEFF < 5 | Manual `Rscript benchmarks/stepstone_all_methods.R` |

---

## Out of Scope

- `method="greenkhorn"` SQUAREM fix (separate: leafblower-i0am)
- Logit with per-obs bounds mode (cell-aggregate only; same limitation as all solvers)
- Exposing Armijo constants to R user (internal tuning; hardcoded is simpler)
- Continuation from GREG (higher complexity, not needed given A+B suffice)
