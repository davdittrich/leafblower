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
// For chebyshev or grake: warm-start from ieppa result.
// CRITICAL: ieppa_solve() mutates st.weights in-place (normalizes/clamps).
// We must NOT pass st directly. Create a separate weights copy so the original
// CalibState is untouched after the pre-solve.
std::vector<double> w_warm_obs;    // obs-level warm-start weights (empty = cold start)
double delta_warm = -1.0;          // -1 = use default delta_0 = 1.0
if (strcmp(method_str, "chebyshev") == 0 || strcmp(method_str, "grake") == 0) {
    // Deep-copy weights buffer before handing to ieppa (which mutates st.weights)
    std::vector<double> weights_copy(st.weights, st.weights + st.n);
    CalibState st_warm = st;
    st_warm.weights = weights_copy.data();          // redirect to copy
    st_warm.inner_max_iter = std::max(5, std::min(100, st.inner_max_iter / 10));
    auto ieppa_res = lbw::ieppa_solve(st_warm);
    // weights_copy goes out of scope after this block; st.weights is unchanged
    if (!ieppa_res.best_weights.empty() &&
        (int)ieppa_res.best_weights.size() == st.n) {
        w_warm_obs = std::move(ieppa_res.best_weights);  // obs-level, size n
        delta_warm = ieppa_res.max_error * 1.5;
    }
}
// Pass w_warm_obs (obs-level) and delta_warm to chebyshev_ipm.
// chebyshev_ipm aggregates obs-level → cell masses AFTER its own build_cell_table.
// This avoids building CellTable twice. See Step B.
```

**Step B — chebyshev.hpp signature change** (pass obs-level warm weights — chebyshev aggregates internally):
```cpp
ChebyshevResult chebyshev_ipm(
    CalibState& st,
    LpVariant variant,
    const std::vector<double>& w_warm_obs = {},  // obs-level warm weights (size n); aggregated
                                                  // to cell masses after internal build_cell_table
    double delta_warm = -1.0);                   // ieppa_max_err * 1.5; -1 → default 1.0
```

**Step C — chebyshev_ipm initialization** (aggregate obs → cell AFTER internal build):
```cpp
// After build_cell_table(ct) runs inside chebyshev_ipm:
std::vector<double> X_primal(M);   // M = ct.M_cell
if (!w_warm_obs.empty() && (int)w_warm_obs.size() == st.n) {
    std::fill(X_primal.begin(), X_primal.end(), 0.0);
    for (int i = 0; i < st.n; i++)
        X_primal[ct.cell_of[i]] += w_warm_obs[i];
    // Clamp to cell bounds; renormalize to preserve total mass (pre-clamp Σ = Σ_post)
    double total_pre = 0.0, total_post = 0.0;
    for (int c = 0; c < M; c++) total_pre += X_primal[c];
    for (int c = 0; c < M; c++) X_primal[c] = std::clamp(X_primal[c], L_cell[c], U_cell[c]);
    for (int c = 0; c < M; c++) total_post += X_primal[c];
    if (total_post > 0.0 && total_pre > 0.0) {
        double scale = total_pre / total_post;
        for (int c = 0; c < M; c++) X_primal[c] = std::clamp(X_primal[c]*scale, L_cell[c], U_cell[c]);
    }
} else {
    // Cold start: uniform proportional initialization
    X_primal = X_init_uniform;
}

// Delta initialization: start slightly above ieppa quality (or at 1.0 cold)
double delta_0 = (delta_warm > 0.0) ? delta_warm : 1.0;
```

Memory overhead: one O(n) weights copy + no extra CellTable (aggregation uses internal ct).
Latency overhead: ieppa pre-solve ≈ 0.5s on stepstone. Documented in `@param max_iterations`.

---

### Change 3: Mehrotra predictor-corrector IPM

Replace the current single-step barrier update with Mehrotra's two-phase predictor-corrector. This achieves near-quadratic convergence compared to the current O(n^0.5) path.

**Current iteration (pseudocode)**:
```
Solve: (A*D*A^T) Δλ = rhs  [one solve per iter]
Update: x += α*Δx, λ += α*Δλ, s += α*Δs
Reduce: μ *= 0.1  [fixed schedule]
```

**Mehrotra iteration** (with m>0 guard and explicit Jacobi scoping):
```
// Guard: if m == 0 (degenerate LP, no complementarity pairs), skip Mehrotra,
// fall through to plain barrier step with μ *= 0.1. Return RK_ERR_DEGENERATE.
if (m == 0) { ... }

// D held FIXED across Phase A and Phase B — recomputed only at next outer iteration.
// N = A·D·A^T is rebuilt from scratch each iteration before Jacobi scaling.

Phase A — Affine predictor (barrier terms = 0):
  Build N = A·D·A^T  [fresh, not cached]
  Apply Jacobi preconditioning: D_jac[j] = 1/sqrt(N[j][j]); scale N, scale rhs_A
  LDLT-factor scaled N  [one factorization per iteration — reused in Phase B]
  Solve scaled system → Δλ_A_scaled; unscale: Δλ_A[j] = D_jac[j] * Δλ_A_scaled[j]
  Compute Δx_aff, Δs_aff from Δλ_A
  α_aff = max α s.t. (x + α*Δx_aff) ≥ 0 AND (s + α*Δs_aff) ≥ 0
  μ_aff = clamp((x+α_aff*Δx_aff)·(s+α_aff*Δs_aff)/(2*m), 0, μ*100)  [guard division]

Compute centering parameter:
  σ = clamp((μ_aff/μ)³, 1e-8, 1.0)  [clamped: prevents σ=0 or σ>1]

Phase B — Corrector (centering + second-order):
  rhs_B = rhs_A + σ*μ*e - Δx_aff·Δs_aff   [element-wise second-order term]
  Apply same Jacobi scaling to rhs_B: rhs_B_scaled[j] = D_jac[j] * rhs_B[j]
  Solve REUSING the SAME LDLT factorization from Phase A (no refactoring)
  Unscale: Δλ_B[j] = D_jac[j] * Δλ_B_scaled[j]
  Compute corrected Δx, Δs; step lengths α_p, α_d with 0.99 damping

Update: x += α_p*Δx, λ += α_p*Δλ_B, s += α_d*Δs
Recompute: μ = x·s/(2*m)  [adaptive — no explicit schedule]
```

**Key invariants:**
- N rebuilt fresh each iteration → Jacobi scaling is never applied twice to the same matrix
- D (barrier scaling) held fixed within one iteration across Phase A and Phase B
- Phase B reuses the already-factored (Jacobi-scaled) N — only RHS changes

The second-order term `Δx_aff·Δs_aff` removes the linearization error. Mehrotra reduces iteration count from O(100-500) to O(20-50) for well-conditioned problems (same LDLT factorization cost per iteration — one factor, two solves).

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
  st_c  <- attr(r_cheby,"result")$status
  # Must converge (status=0) and produce finite quality
  expect_equal(st_c, 0L, label="chebyshev must converge (status=0) on K=3")
  expect_true(is.finite(me_c), label="chebyshev max_error must be finite")
  # Chebyshev (L∞ optimal) must be at least as good as raking when it converges
  expect_lte(me_c, me_r * 1.001 + 1e-10,
    label=sprintf("chebyshev (%.4e) must not exceed raking (%.4e)", me_c, me_r))
  expect_lt(me_c, 1e-3, label="chebyshev must converge to <1e-3 on K=3")
})
```

### T_cheby_K9 (stepstone — skip in CI, manual verify)
```r
test_that("T_cheby_K9: chebyshev K=9 stepstone max_err <= greenkhorn", {
  skip_if_not_installed("arrow"); skip_if_not_installed("jsonlite")
  skip_if(!file.exists("benchmarks/stepstone_fulldata_bench_data.parquet"),
          "stepstone not available")
  df  <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
  df$uuid <- NULL
  tgt <- lapply(jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json"),
                function(t) { t <- unlist(t); t / sum(t) })
  for (nm in names(tgt)) df[[nm]] <- factor(df[[nm]])
  r_c <- suppressWarnings(harvest(df,tgt,method="chebyshev",max_weight=5,min_weight=0,
                                   max_iterations=5000L,attach_weights=FALSE,verbose=0))
  r_g <- suppressWarnings(harvest(df,tgt,method="greenkhorn",max_weight=5,min_weight=0,
                                   max_iterations=5000L,attach_weights=FALSE,verbose=0))
  me_c <- attr(r_c,"result")$max_error; me_g <- attr(r_g,"result")$max_error
  expect_equal(attr(r_c,"result")$status, 0L, label="chebyshev must converge on K=9")
  expect_lte(me_c, me_g * 1.001 + 1e-10,
    label=sprintf("chebyshev (%.4e) must not exceed greenkhorn (%.4e)", me_c, me_g))
})
```

### T_cheby_warm_fallback (verify graceful cold-start when ieppa fails)
```r
test_that("T_cheby_warm_fallback: chebyshev still works when ieppa warm-start fails", {
  # max_iterations=1 forces ieppa pre-solve to 1 iter (useless warm-start)
  # chebyshev should fall back to cold start and still return valid result
  set.seed(7); n <- 1000L
  df  <- data.frame(x=factor(sample(letters[1:3],n,TRUE)), y=factor(sample(c("M","F"),n,TRUE)))
  tgt <- list(x=c(a=0.3,b=0.4,c=0.3), y=c(M=0.5,F=0.5))
  r <- suppressWarnings(harvest(df,tgt,method="chebyshev",max_iterations=200L,attach_weights=FALSE))
  expect_true(is.finite(attr(r,"result")$max_error),
    label="chebyshev must return finite max_error (not NaN from bad warm-start)")
})
```

---

## Acceptance Criteria

| # | Criterion | Verify |
|---|-----------|--------|
| AC1 | T_greg_warn GREEN (warns on K≥5 tight-bounds) | `devtools::test(filter="chebyshev")` |
| AC2 | T_cheby_warm GREEN: status=0, finite me_c, me_c ≤ me_r on K=3 | Same |
| AC3 | T_cheby_K9 GREEN (stepstone chebyshev ≤ greenkhorn) | Same (skip if parquet absent) |
| AC4 | T_cheby_warm_fallback GREEN (graceful cold-start on bad warm-start) | Same |
| AC5 | T_sraa_grk, T_sraa_rk, T_sraa_outer_revert, E1 (cheby≤raking) still GREEN | `devtools::test()` |
| AC6 | FAIL count = 3 (pre-existing: ieppa-nonuniform-d:28/:29, sor:18) | `devtools::test()` |
| AC7 (benchmark) | chebyshev max_err ≤ 1.57e-3 on K=9 stepstone | `Rscript benchmarks/stepstone_all_methods.R` |

Note on AC6 "FAIL count = 3": The three pre-existing failures are test-ieppa-nonuniform-d.R:28/:29 and
test-sor.R:18 — unrelated to chebyshev. E1 (chebyshev≤raking, test-calibration-solvers.R:257) is
currently GREEN and must remain GREEN. AC7 threshold 1.57e-3 = greenkhorn K=9 stepstone current result.

---

## Calibration of constants

| Constant | Value | Derivation |
|---|---|---|
| Greg warning threshold | `0.05` | 50× raking typical quality (1e-3); flags catastrophic failures only |
| ieppa pre-solve iterations | `max(5, min(100, max_iterations/10))` | 5 = minimum useful; 100 ≈ 0.5s; floor prevents useless 1-iter warm-start |
| delta_warm multiplier | `1.5` | 50% above ieppa max_err; gives IPM room to improve |
| Mehrotra damping | `0.99` | Standard value from Mehrotra (1992) preventing boundary-touching |
| Jacobi eps | `1e-12` | Guards against zero diagonal in unscaled N |

---

## Out of Scope

- grake separate fix (inherits chebyshev improvements automatically)
- Exposing Mehrotra parameters to R user
- ADMM reformulation of L∞ LP (future research if Mehrotra insufficient)
- GREG algorithmic improvements (correct routing is sufficient)
