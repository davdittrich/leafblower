# Chebyshev IPM Rewrite + GREG Quality Warning

**Date**: 2026-04-29
**Status**: Pending design review
**Tickets**: leafblower-yned (chebyshev), leafblower-scii adjacent
**Files**: `src/chebyshev.hpp`, `src/chebyshev.cpp`, `src/r_bridge.cpp`, `R/harvest.R`,
           `tests/testthat/`

---

## Problem

Three methods perform poorly on stepstone K=9 (n=1.58M):

| Method | max_err | status | issue |
|---|---|---|---|
| greg | 1.10e-01 (11%) | converged | converges to wrong solution for K≥5/tight bounds |
| chebyshev | NaN | RK_ERR_NOCONV | 500-iter hard cap hit; cold-start too far from optimum |
| grake | NaN | RK_ERR_NOCONV | same as chebyshev (one-line wrapper) |

**GREG root cause**: GREG's regression model assumes independent margins and chi-squared distance minimization. For K≥5 overlapping margins with tight bounds, the design matrix Gram inverse is near-singular and Newton converges to the wrong interior point.

**Chebyshev root cause**: Two compounding failures:
1. Cold uniform initialization → barrier path too long for 500 iterations
2. Ill-conditioned normal equations (`A*D*A^T`) for K=9 overlapping margins → Schur complement degenerates, dual explosion

**Value of fixing chebyshev**: Chebyshev minimizes L∞ error **directly** — it is the only method that certifiably produces the globally optimal max_err solution. If it works, every benchmark metric we report becomes an upper bound on the true optimum.

---

## Design

### Change 1: GREG quality-check warning

After greg dispatch, if `max_err > 0.05` (5%), emit a diagnostic warning. User keeps their result — warn-and-continue, not fallback.

**`R/harvest.R`** (after greg result assembled, before return):
```r
if (method == "greg" &&
    !is.null(calib_result$max_error) &&
    is.finite(calib_result$max_error) &&
    calib_result$max_error > 0.05) {
  warning(sprintf(
    paste0("greg converged but max_err=%.4g (>5%%). ",
           "greg may be unreliable for K=%d margins or tight bounds ",
           "(max_weight=%.4g). Consider method='raking' or 'ieppa'."),
    calib_result$max_error, length(target), max_weight),
    call. = FALSE)
}
```

No algorithmic change. One testthat test: `expect_warning(harvest(..., method="greg", ...), "greg.*unreliable")`.

---

### Change 2: Chebyshev warm-start from ieppa (Option A)

**Step A — r_bridge.cpp**: Before dispatching chebyshev/grake, run a fast ieppa pre-solve:

```cpp
// For chebyshev or grake: warm-start from ieppa result
std::vector<double> X_warm;        // cell masses (empty = cold start)
double delta_warm = -1.0;          // -1 = use default delta_0 = 1.0
if (strcmp(method_str, "chebyshev") == 0 || strcmp(method_str, "grake") == 0) {
    CalibState st_warm = st;
    st_warm.inner_max_iter = std::min(100, st.inner_max_iter / 10);
    auto ieppa_res = lbw::ieppa_solve(st_warm);
    if (!ieppa_res.best_weights.empty() &&
        (int)ieppa_res.best_weights.size() == st.n) {
        // Build temporary cell table to aggregate obs-level → cell masses
        CellTable ct_tmp;
        if (build_cell_table(st.n, st.K, st.group_ids, st.cat_counts,
                             st.weights, ct_tmp) == RK_OK) {
            X_warm.assign(ct_tmp.M_cell, 0.0);
            for (int i = 0; i < st.n; i++)
                X_warm[ct_tmp.cell_of[i]] += ieppa_res.best_weights[i];
            delta_warm = ieppa_res.max_error * 1.5;
        }
    }
}
// Pass X_warm and delta_warm to chebyshev_ipm (see Step B)
```

**Step B — chebyshev.hpp signature change**:
```cpp
ChebyshevResult chebyshev_ipm(
    CalibState& st,
    LpVariant variant,
    const std::vector<double>& X_warm = {},  // cell-level masses from ieppa
    double delta_warm = -1.0);               // ieppa_max_err * 1.5; -1 → default 1.0
```

**Step C — chebyshev_ipm initialization**:
```cpp
// Primal initialization
std::vector<double> X_primal = X_warm.empty() ? X_init_uniform : X_warm;
// Bound X_warm to cell bounds (ieppa may have clamped differently)
for (int c = 0; c < M; c++)
    X_primal[c] = std::clamp(X_primal[c], L_cell[c], U_cell[c]);

// Delta initialization: start slightly above ieppa quality (or at 1.0 cold)
double delta_0 = (delta_warm > 0.0) ? delta_warm : 1.0;
```

Memory overhead: one extra CellTable build + X_warm vector (M_cell doubles ≈ negligible).
Latency overhead: ieppa pre-solve at max_iterations=100 ≈ 0.5s on stepstone.

---

### Change 3: Mehrotra predictor-corrector IPM

Replace the current single-step barrier update with Mehrotra's two-phase predictor-corrector. This achieves near-quadratic convergence compared to the current O(n^0.5) path.

**Current iteration (pseudocode)**:
```
Solve: (A*D*A^T) Δλ = rhs  [one solve per iter]
Update: x += α*Δx, λ += α*Δλ, s += α*Δs
Reduce: μ *= 0.1  [fixed schedule]
```

**Mehrotra iteration**:
```
Phase A — Affine predictor (barrier terms = 0):
  Solve normal equations once (r_aff, no centering)
  Compute: α_aff = max α s.t. (x + α*Δx_aff) ≥ 0 AND (s + α*Δs_aff) ≥ 0
  Compute: μ_aff = (x + α_aff*Δx_aff)·(s + α_aff*Δs_aff) / (2*m)

Compute centering parameter:
  σ = (μ_aff / μ)³    [adaptive: large σ → more centering when affine step short]

Phase B — Corrector (centering + second-order):
  RHS correction: += σ*μ*e - Δx_aff·Δs_aff   [element-wise]
  Solve normal equations again (same factorization, different RHS)
  Compute corrected step lengths α_p, α_d with 0.99 damping

Update: x, λ, s with corrected step
Recompute: μ = x·s / (2*m)  [adaptive — no explicit schedule]
```

The second-order term `Δx_aff·Δs_aff` removes the linearization error that causes slow convergence near the solution. In practice, Mehrotra reduces iteration count from O(100-500) to O(20-50) for well-conditioned problems.

**Implementation detail**: Phase A and Phase B use the same LDLT factorization of `(A*D*A^T)` — the matrix is factored once per iteration, not twice. Only the RHS differs. This is O(nct²) factorization + 2×O(nct²) solves per iteration instead of the current O(nct²) factorization + 1 solve.

---

### Change 4: Jacobi diagonal preconditioning

For K=9 overlapping margins, the normal equations matrix `N = A*D*A^T` is ill-conditioned. Diagonal (Jacobi) preconditioning scales each row/column by `1/sqrt(N[j][j])`:

```cpp
// Before LDLT factorization:
std::vector<double> D_jac(nct);
for (int j = 0; j < nct; j++)
    D_jac[j] = 1.0 / std::sqrt(std::max(N[j*nct+j], 1e-12));

// Scale N: N_scaled[i][j] = D_jac[i] * N[i][j] * D_jac[j]
for (int i = 0; i < nct; i++)
    for (int j = 0; j < nct; j++)
        N[i*nct+j] *= D_jac[i] * D_jac[j];

// Scale RHS: rhs_scaled[j] = D_jac[j] * rhs[j]
// After solve: Δλ_unscaled[j] = D_jac[j] * Δλ_scaled[j]
```

This costs 2×O(nct²) extra per iteration but reduces condition number dramatically for overlapping margins. For K=9 stepstone with nct≈80: negligible relative to the O(nct³) LDLT factorization.

---

## Architecture: Files changed

| File | Change |
|---|---|
| `R/harvest.R` | Add greg quality-check warning (3 lines) |
| `src/chebyshev.hpp` | Add X_warm + delta_warm parameters to chebyshev_ipm |
| `src/chebyshev.cpp` | (1) Warm-start initialization; (2) Mehrotra loop; (3) Jacobi preconditioning |
| `src/r_bridge.cpp` | For chebyshev/grake: run ieppa pre-solve, extract X_warm, pass to chebyshev_ipm |
| `tests/testthat/test-chebyshev.R` | New test file with T_greg_warn, T_cheby_warm, T_cheby_convergence |

---

## TDD Requirements

### T_greg_warn
```r
test_that("T_greg_warn: greg warns when max_err > 5% on tight-bounds K=5 problem", {
  set.seed(99); n <- 3000L
  df <- data.frame(
    a=factor(sample(letters[1:4],n,TRUE)), b=factor(sample(LETTERS[1:3],n,TRUE)),
    c=factor(sample(c("x","y"),n,TRUE)),   d=factor(sample(c("M","F"),n,TRUE)),
    e=factor(sample(c("Y","O"),n,TRUE))
  )
  tgt <- list(a=setNames(c(0.4,0.3,0.2,0.1),letters[1:4]),
              b=setNames(c(0.5,0.3,0.2),LETTERS[1:3]),
              c=c(x=0.6,y=0.4), d=c(M=0.45,F=0.55), e=c(Y=0.55,O=0.45))
  expect_warning(
    harvest(df, tgt, method="greg", max_weight=1.8, min_weight=0,
            max_iterations=50, attach_weights=FALSE),
    regexp="greg.*unreliable|greg.*max_err")
})
```

### T_cheby_warm (K=3, verifies warm-start reduces iterations)
```r
test_that("T_cheby_warm: chebyshev with ieppa warm-start converges on K=3 problem", {
  set.seed(7); n <- 5000L
  df <- data.frame(
    a=factor(sample(letters[1:3],n,TRUE)),
    b=factor(sample(LETTERS[1:4],n,TRUE)),
    c=factor(sample(c("M","F"),n,TRUE))
  )
  tgt <- list(a=setNames(c(0.4,0.35,0.25),letters[1:3]),
              b=setNames(c(0.3,0.3,0.2,0.2),LETTERS[1:4]),
              c=c(M=0.48,F=0.52))
  r_cheby <- suppressWarnings(harvest(df,tgt,method="chebyshev",
                                       max_iterations=500,attach_weights=FALSE))
  r_raking <- suppressWarnings(harvest(df,tgt,method="raking",
                                        max_iterations=500,attach_weights=FALSE))
  me_c <- attr(r_cheby,"result")$max_error
  me_r <- attr(r_raking,"result")$max_error
  # Chebyshev (L∞ optimal) must be at least as good as raking
  expect_lte(me_c, me_r * 1.001 + 1e-10,
    label=sprintf("chebyshev (%.4e) must not exceed raking (%.4e)", me_c, me_r))
  expect_lt(me_c, 1e-3, label="chebyshev must converge to <1e-3 on K=3")
})
```

### T_cheby_K9 (stepstone — skip in CI, manual verify)
```r
test_that("T_cheby_K9: chebyshev K=9 stepstone max_err <= greenkhorn", {
  skip_if(!file.exists("benchmarks/stepstone_fulldata_bench_data.parquet"),
          "stepstone not available")
  # ... same pattern as T_sraa_adaptive_K9
})
```

---

## Acceptance Criteria

| # | Criterion | Verify |
|---|-----------|--------|
| AC1 | T_greg_warn GREEN | `devtools::test(filter="chebyshev")` |
| AC2 | T_cheby_warm GREEN (K=3 converges, ≤ raking) | Same |
| AC3 | T_cheby_K9 GREEN (stepstone chebyshev ≤ greenkhorn) | Manual benchmark |
| AC4 | T_sraa_grk, T_sraa_rk, T_sraa_outer_revert still GREEN | `devtools::test()` |
| AC5 | FAIL count = 3 (unchanged) | Same |
| AC6 (benchmark) | chebyshev max_err ≤ 1.57e-3 on K=9 | `Rscript benchmarks/stepstone_all_methods.R` |

---

## Calibration of constants

| Constant | Value | Derivation |
|---|---|---|
| Greg warning threshold | `0.05` | 5× raking quality (1e-3 → flag at 5e-2); catches systematic failure |
| ieppa pre-solve iterations | `min(100, max_iterations/10)` | 100 iters ≈ 0.5s; enough for ieppa to get to ~1e-3 quality |
| delta_warm multiplier | `1.5` | 50% above ieppa max_err; gives IPM room to improve |
| Mehrotra damping | `0.99` | Standard value from Mehrotra (1992) preventing boundary-touching |
| Jacobi eps | `1e-12` | Guards against zero diagonal in unscaled N |

---

## Out of Scope

- grake separate fix (inherits chebyshev improvements automatically)
- Exposing Mehrotra parameters to R user
- ADMM reformulation of L∞ LP (future research if Mehrotra insufficient)
- GREG algorithmic improvements (correct routing is sufficient)
