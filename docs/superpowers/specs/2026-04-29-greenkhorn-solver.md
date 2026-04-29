# Greenkhorn + Logit Calibration Solvers — Design Spec (rev 5 — gate approved)

**Date**: 2026-04-29
**Status**: Pending design review (rev 3 adds Part 2: Logit calibration)

## Two new methods

| Method | Algorithm | Autumn analogue | Bounds |
|--------|-----------|-----------------|--------|
| `method="greenkhorn"` | Greedy coordinate-descent IPF (Part 1) | `autumn::harvest()` greedy | clamp |
| `method="logit"` | Deville-Särndal logit Newton-Raphson (Part 2) | `autumn::calibrate()` | **by construction** |

Greenkhorn ≡ autumn's greedy raking. Logit ≡ autumn's calibrate with `"logit"` distance
(Deville & Särndal 1992), which "enforces the cap by construction at every Newton step."

---

## Problem / Motivation

The benchmark shows autumn's greedy IPF achieves `max_err=1.60e-3` in 33.6s on stepstone
while leafblower's raking achieves the same quality in 2.5s via round-robin. The goal is NOT
to replicate autumn's speed profile (leafblower already wins by 13×). The goals are:

1. **Machine precision on simple margins**: on a trivial 1-margin problem, Greenkhorn
   converges in 1 step to floating-point precision. Raking also does this in 1 round (K=1),
   so this benefit is only meaningful when K≥2 and margins have very different difficulty.
   More precisely: Greenkhorn's greedy selection ensures the hardest margin is always
   addressed first — simple margins drift only slightly before being revisited.

2. **Research baseline**: Greenkhorn + bounded IPF is the foundation for more advanced
   variants (per-margin step sizes, multi-margin selection). Having a clean C++
   implementation enables future research.

3. **Explicit method name**: `method="greenkhorn"` for users who explicitly want greedy
   coordinate descent (distinct from `scheduler="greedy"` on raking — see §API Notes).

**Use case**: K≥2 margins with highly unequal difficulty (e.g., one 2-cat margin and one
K=100-cat margin) where greedy scheduling polishes the easy margin to machine precision while
concentrating budget on the hard one. Raking and Greenkhorn should converge to the same
fixed point in infinite budget; the difference is the path.

---

## Algorithm

**Greenkhorn IPF** (Altschuler-Weed-Rigollet 2017, NIPS, arxiv:1705.09634):

```
Initialize:
  X[c] = Σ_{i: cell_of[i]=c} st.weights[i]   ∀c          (cell initial masses)
  X_init[c] = X[c]                                          (saved for reconstruction)
  W = Σ_c X[c]                                              (total mass, maintained)
  S[k][j] = Σ_{c ∈ cells[k][j]} X[c]   ∀k,j               (bucket sums)
  errRp[k] = max_j |S[k][j]/W − τ[k][j]|   ∀k

Repeat until max_k errRp[k] < tol OR iter = max_iter:
  (1) k* = argmax_k errRp[k]

  (2) Guard against W ≤ 0 (all cells clamped to zero — structurally infeasible):
        if (W <= 0.0) {
            // All cells clamped to zero — structurally infeasible
            res.status = RK_ERR_INFEAS;
            std::snprintf(res.message, sizeof(res.message),
                "greenkhorn: total mass W<=0 (all cells at zero bound)");
            break;
        }

      For each j ∈ {0…cat_counts[k*]−1}:
        if S[k*][j] < ε·W: continue          (structural empty)
        f = (W × τ[k*][j]) / S[k*][j]       (scale factor)

        For each c ∈ cells_per_cat[k*][j]:   (local per-bucket cell list)
          X_old = X[c]
          X_new = clamp(X_old × f, L_cell[c], U_cell[c])
          δ = X_new − X_old
          X[c] = X_new
          W   += δ                            (maintain total mass)

          For each k2 ≠ k*:
            g2 = g_per_cell[k2][c]
            if g2 ≥ 0: S[k2][g2] += δ

        S[k*][j] recomputed: Σ_{c∈cells_per_cat[k*][j]} X[c]

  (3) errRp[k] = max_j |S[k][j]/W − τ[k][j]|  for all k  (O(Σ cat_counts))
  (4) Best-iterate tracking (guards against returning worse-than-best on budget exhaustion):
        double curr_max_errRp = max_k errRp[k]
        if curr_max_errRp < best_errRp:
            best_errRp = curr_max_errRp
            X_best = X  // snapshot
            res.best_iter = iter
  (5) iter++

Return X_best (not X) as final weights — X_best is the iterate with lowest max errRp seen,
protecting against regression when the budget runs out mid-oscillation.
```

---

## Data structure notes (verified against codebase)

### What CellTable actually provides (from src/cell_table.hpp)

```cpp
struct CellTable {
    int    M_cell;
    std::vector<int>    cell_of;          // cell_of[i] ∈ [0, M_cell)
    std::vector<int>    n_per_cell;       // n_per_cell[c] = obs count in cell c
    std::vector<std::vector<int>> g_per_cell;  // g_per_cell[k][c] ∈ [-1, cat_counts[k])
    double W_input;
    double capacity_mu_auto;
};
```

**NOT available**: `cell_members[c]` (obs per cell) and `cells[k][j]` (cells per bucket).
These must be built locally inside `greenkhorn_solve`, mirroring raking.cpp's pattern.

### Local construction (mirrors raking.cpp lines 82-97)

```cpp
// X init from design weights:
std::vector<double> X(ct.M_cell, 0.0);
for (int i = 0; i < st.n; i++) X[ct.cell_of[i]] += st.weights[i];
const std::vector<double> X_init = X;   // save for reconstruction

// Per-bucket cell lists — built once, used every iteration:
// cells_per_cat[k][j] = list of cell indices c where g_per_cell[k][c] == j
std::vector<std::vector<std::vector<int>>> cells_per_cat(st.K);
for (int k = 0; k < st.K; k++) {
    cells_per_cat[k].assign(st.cat_counts[k], {});
    for (int c = 0; c < ct.M_cell; c++) {
        int g = ct.g_per_cell[k][c];
        if (g >= 0 && g < st.cat_counts[k]) cells_per_cat[k][g].push_back(c);
    }
}
```

---

## W maintenance

W = Σ X[c] changes when cells are clamped. Maintain incrementally:
```cpp
double W = 0.0;
for (int c = 0; c < ct.M_cell; c++) W += X[c];

// Inside bucket update, per cell:
double delta = X_new - X_old;
W += delta;
```

Re-scale `errRp` and `S/W` using the updated W at every check. This ensures the
convergence criterion reflects the actual constrained optimum, not a biased approximation.

---

## Convergence API (from src/calib_dispatch.hpp)

Use `check_convergence(cfg, m, prev_metric, st.tol_abs)` from calib_dispatch.hpp.
`prev_metric` is passed by reference and updated inside `check_convergence`.

```cpp
#include "calib_dispatch.hpp"
double prev_metric = std::numeric_limits<double>::infinity();
CellMetrics m = compute_cell_metrics(st, ct, X, W, bucket_scratch);
// ...
bool converged = check_convergence(st.convergence_cfg, m, prev_metric, st.tol_abs);
```

`compute_cell_metrics` (calib_dispatch.hpp) fills all 6 metrics (errRp, mean_err, kl,
chi2, grake_norm, l1) in one pass. Use this for the convergence check AND for final
result population (to fill all `rk_result_t` fields correctly).

Note: `compute_cell_metrics` needs a `bucket` scratch vector of size ≥ max(cat_counts[k]).

---

## Input validation (pre-Greenkhorn, in r_bridge.cpp dispatch)

```cpp
// Validate before calling greenkhorn_solve():
if (p.min_weight >= p.max_weight && p.max_weight > 0.0) {
    Rf_error("leafblower: min_weight (%.4g) >= max_weight (%.4g)",
             p.min_weight, p.max_weight);
}
```

This prevents UB in `std::clamp(x, L, U)` when `L > U` (undefined behavior in C++17 per
[alg.clamp] precondition).

---

## GreenkornResult / result mapping

`GreenkornResult` maps 1:1 to the fields that `pack_solver_result` expects in raking. Use
the SAME pattern as `RakingResult` — check `src/raking.hpp` for the full struct definition
and copy it. Then r_bridge.cpp dispatch mirrors the raking dispatch block:

```cpp
} else if (strcmp(method_str, "greenkhorn") == 0) {
    auto res = lbw::greenkhorn_solve(st);
    pack_solver_result(res);
    res_status     = res.status;
    res_iterations = res.iterations;
    res_max_error  = res.max_error;
    res_alg_used   = (int)RK_ALG_GREENKHORN;   // = 9
    if (!res.best_weights.empty())
        res_best_weights = std::move(res.best_weights);
    else
        res_best_weights.assign(st.n, 0.0);
```

**harvest.R `alg_names` extension** (harvest.R line ~437):
```r
alg_names <- c("", "ieppa", "lbfgsb", "raking", "sinkhorn",
                "chebyshev", "greg", "grake", "ieppa_soft", "greenkhorn")
#                                                             ^ index 9 (RK_ALG_GREENKHORN=9)
```

---

## API Notes

**`method="greenkhorn"` vs `scheduler="greedy"`**:
- `scheduler="greedy"` on `method="raking"`: sweeps ALL K margins per iteration in
  residual-sorted order (greedy-sorted round-robin). Documented in harvest.R as
  "Greenkhorn priority" (an approximation).
- `method="greenkhorn"`: pure single-margin coordinate descent. ONE margin per step.
  These are different algorithms with different convergence paths. Add to `@param method`
  and `@param scheduler` @details clarifying this distinction.

**Unsupported params** — emit named warning if any of the following are supplied non-NULL:
`accelerate`, `scheduler`, `homotopy_levels`, `homotopy_start_factor`,
`homotopy_end_factor`, `eta_schedule`, `capacity_penalty`, `sor_enabled`.
Exact warning: `warning("method='greenkhorn' does not support '", nm, "'; ignored", call.=FALSE)`

---

## New enum value

In `src/leafblower.h`, after `RK_ALG_IEPPA_SOFT = 8`:
```c
RK_ALG_GREENKHORN = 9
```

`EXPECTED_RK_PARAMS_BYTES` and `EXPECTED_RK_RESULT_BYTES` are NOT affected
(no new fields in `rk_params_t` or `rk_result_t`).

---

## R/harvest.R critical changes

Three changes are required in `R/harvest.R` — without them ALL tests fail:

1. **`map_method()` must include both `"greenkhorn"` and `"logit"`** — if either is missing, `harvest()` rejects the method string before dispatch reaches C++ and every T1/T5 test fails with a wrong-method error.

2. **`alg_names` must extend to 11 elements (indices 0–10)**, adding `"greenkhorn"` at index 9 and `"logit"` at index 10:
   ```r
   alg_names <- c("", "ieppa", "lbfgsb", "raking", "sinkhorn",
                   "chebyshev", "greg", "grake", "ieppa_soft",
                   "greenkhorn", "logit")
   #               ^ 9           ^ 10
   ```
   Without this, `algorithm_used` returns `""` and T1/T5 `algorithm_used` assertions fail.

3. **`harvest.R` `stop()` for `status==2` must use `res$message`** (the C++ solver's actual message) instead of a hardcoded empty-cell string. The logit solver returns `RK_ERR_INFEAS` with a specific message (e.g., "logit: singular normal equations") — a hardcoded override discards that diagnostic and confuses T5/T8 failures.

---

## Files to create / modify (both methods combined)

| File | Change |
|------|--------|
| `src/greenkhorn.cpp` | New solver (greenkhorn_solve) |
| `src/greenkhorn.hpp` | GreenkornResult struct (mirrors RakingResult) |
| `src/logit_calib.cpp` | New solver (logit_calibrate) |
| `src/logit_calib.hpp` | LogitCalibResult struct (mirrors GregResult) |
| `src/leafblower.h` | Add RK_ALG_GREENKHORN=9, RK_ALG_LOGIT=10 to enum |
| `src/Makevars.in` | Add both `.cpp` files to PKG_SOURCES |
| `src/r_bridge.cpp` | Add "greenkhorn" and "logit" dispatch blocks; alg_names to 11 elements |
| `src/c_api.cpp` | Add GREENKHORN and LOGIT cases |
| `R/harvest.R` | Add both to map_method(); update @param method; alg_names to 11 |
| `man/harvest.Rd` | Auto-regenerate via devtools::document() |
| `tests/testthat/test-calibration-solvers.R` | Add T1-T4 (Greenkhorn) + T5-T8 (Logit) |

---

## TDD (RED before implementation)

### T1 — Greenkhorn available, calibrates, returns correct algorithm_used

```r
test_that("T1: greenkhorn available and calibrates", {
  set.seed(1); n <- 1000L
  df  <- data.frame(sex=factor(sample(c("M","F"),n,TRUE)),
                    age=factor(sample(c("Y","O"),n,TRUE)))
  tgt <- list(sex=c(M=0.5,F=0.5), age=c(Y=0.6,O=0.4))
  r   <- harvest(df, tgt, method="greenkhorn", max_iterations=500L)
  expect_lt(attr(r,"result")$max_error, 1e-3)
  expect_equal(attr(r,"result")$algorithm_used, "greenkhorn",
    label="algorithm_used must be 'greenkhorn' (requires alg_names extended to index 9)")
})
```

### T2 — Machine precision on trivial 1-margin problem

```r
test_that("T2: greenkhorn reaches near-machine precision on 1-margin 2-cat", {
  # Greenkhorn on a single 2-cat margin with no bounds converges in 1 step
  # (one exact scale factor). Use absolute convergence to not stop early.
  set.seed(2); n <- 5000L
  df  <- data.frame(g=factor(sample(c("A","B"), n, TRUE, prob=c(0.3,0.7))))
  tgt <- list(g=c(A=0.5, B=0.5))
  r   <- harvest(df, tgt, method="greenkhorn",
                 max_iterations=100L,
                 convergence=list(absolute=1e-12))  # REQUIRED: default pct=1e-4 stops at 1e-4
  me  <- attr(r,"result")$max_error
  # 1 step → exact scale → machine precision
  expect_lt(me, 1e-10,
    label=sprintf("greenkhorn should reach ~machine precision on 1-margin (got %.2e)", me))
})
```

### T3 — Quality within 2× of raking on 3-margin problem

```r
test_that("T3: greenkhorn max_err within 2x of raking", {
  set.seed(42); n <- 10000L
  df <- data.frame(
    a=factor(sample(letters[1:3],n,TRUE)),
    b=factor(sample(LETTERS[1:4],n,TRUE)),
    c=factor(sample(c("x","y"),n,TRUE))
  )
  tgt <- list(a=c(a=0.3,b=0.4,c=0.3),
              b=c(A=0.25,B=0.25,C=0.25,D=0.25),
              c=c(x=0.6,y=0.4))
  r_rk  <- harvest(df, tgt, method="raking",
                   convergence=list(absolute=1e-6))
  r_grk <- harvest(df, tgt, method="greenkhorn",
                   convergence=list(absolute=1e-6))
  me_rk  <- attr(r_rk,  "result")$max_error
  me_grk <- attr(r_grk, "result")$max_error
  expect_lt(me_grk, 2.0 * me_rk + 1e-6)
})
```

### T4 — Bounds respected

```r
test_that("T4: greenkhorn respects max_weight and min_weight exactly", {
  set.seed(5); n <- 2000L
  df  <- data.frame(v=factor(sample(5, n, TRUE)))
  tgt <- list(v=setNames(c(0.4,0.3,0.15,0.1,0.05), as.character(1:5)))
  r   <- harvest(df, tgt, method="greenkhorn", max_weight=2.0, min_weight=0.1)
  w   <- r$weights
  expect_true(max(w) <= 2.0 + 1e-9)
  expect_true(min(w) >= 0.1 - 1e-9)
})
```

---

## Acceptance Criteria

1. T1 GREEN: max_error < 1e-3; `algorithm_used == "greenkhorn"` (alg_names extended).
2. T2 GREEN: max_error < 1e-10 with `convergence=list(absolute=1e-12)` on 1-margin.
3. T3 GREEN: max_err within 2× of raking on 3-margin problem.
4. T4 GREEN: bounds respected exactly.
5. `min_weight > max_weight` → `Rf_error` with informative message (validation test).
6. Status codes: RK_OK=0, RK_ERR_BUDGET=4, RK_ERR_STALL=5. Never RK_ERR_NOCONV=1.
7. `devtools::test()` FAIL count unchanged (currently 3).
8. `R CMD INSTALL --preclean .` compiles clean.
9. `method="greenkhorn"` emits named warning for unsupported params (accelerate, scheduler, etc.).
10. *(Benchmark-only, not CI-blocking)*: Stepstone stepstone benchmark: greenkhorn
    max_err ≤ 2× raking max_err AND wall time < 30s. Run manually: `OMP_NUM_THREADS=1
    Rscript benchmarks/stepstone_all_methods.R`.

---

## Out of Scope (Part 1 — Greenkhorn)

- SQUAREM (non-stationary F)
- Top-k multi-margin selection per step
- Log-space Greenkhorn (KL-transport entropic regularization)
- Homotopy / ALM / SOR
- AUTO routing (Greenkhorn is explicit opt-in only)
- Autumn API compatibility (convergence_used field differences — out of scope)

---

# Part 2: Logit Calibration (Deville-Särndal 1992)

`method="logit"` implements bounded calibration via logit distance — `autumn::calibrate()`.

## Mathematical formulation

**Objective**: minimize generalized distance
```
Σ_c d_c × F(w_c / d_c)
```
where `d_c = X_init[c]` (cell initial mass), `w_c` is the calibrated cell mass, and F is the
logit distance with bounds A = min_weight, B = max_weight (absolute, cell-aggregate):

```
F(u; A·n_c, B·n_c) = convex separation function whose gradient gives the logit link
```

The key property: **the optimizer writes the solution as a closed-form logistic function**:

```
w_c = L_cell[c] + (U_cell[c] − L_cell[c]) × σ(z_c)
```
where `σ(z) = 1/(1+exp(−z))` is the logistic sigmoid and
`z_c = Σ_k λ[cat_offset[k] + g_per_cell[k][c]]` is the linear predictor for cell c.

This guarantees `w_c ∈ [L_cell[c], U_cell[c]]` for ALL λ — no clamping, no water-fill.

**Calibration constraints**: `Σ_{c ∈ bucket(k,j)} w_c = τ[k][j] × n` for all (k,j).

## Newton-Raphson (dual problem)

Starting from λ=0 (→ w_c = midpoint of [L,U]), iterate:

```
Δλ = N^{-1} × b

where:
  N[j1][j2] = (A·diag(D_eff)·Aᵀ)[j1][j2]     (normal equations, via calib_linalg.cpp)
  b[k][j]   = τ[k][j]×n − Σ_{c ∈ bucket(k,j)} w_c   (calibration residual)
  D_eff[c]  = (U_cell[c] − L_cell[c]) × σ(z_c) × (1−σ(z_c))   (Newton weight)
              // D_eff[c] = ∂w_c/∂z_c = (U-L)·σ(z)·(1-σ(z)) — this IS the correct Newton
              // weight because w_c = L+(U-L)·σ(z_c) and ∂²/∂z² of the logit distance equals
              // ∂w/∂z (Deville-Särndal 1992 eq. 8).

λ += Δλ
```

Then recompute `w_c` and `D_eff[c]` from updated λ. Repeat until `max|b| < tol`.

## Comparison with existing greg.cpp

| Property | greg.cpp (chi-square) | logit_calib.cpp (logit) |
|---|---|---|
| Distance | chi-square | logit |
| Bounds | active-set clamp | analytical (by construction) |
| D_eff[c] | `X_init[c]` (constant) | `(U−L)·σ·(1−σ)` (varies per step) |
| Newton iters | 5-50 (may cycle on tight bounds) | 5-20 (smooth, provably convergent) |
| Degeneracy | ill-conditioned N at bounds | N well-conditioned (σ(1-σ)>0 always) |

## Existing infrastructure reuse (verified)

`src/calib_linalg.cpp` already provides exactly what's needed:
- `compute_normal_equations(ct, D_eff, N, cat_offset, K, nct)` — builds N = A·diag(D)·Aᵀ
- `ldlt_factor_inplace(N, nct, eps)` — factorizes N
- `ldlt_solve(N, nct, b)` — solves NΔλ = b in-place

These are the SAME functions used by greg.cpp. Logit calibration adds only the logit link
function on top of the same normal-equation skeleton.

## Core implementation (greenkhorn_solve → logit_calibrate)

```cpp
// src/logit_calib.cpp

LogitCalibResult logit_calibrate(CalibState& st) {
    // ... cell table build, X_init, L_cell/U_cell as in greg.cpp ...

    // nct = Σ cat_counts[k]; cat_offset[k] = Σ_{k'<k} cat_counts[k']
    int nct = 0;
    std::vector<int> cat_offset(st.K);
    for (int k = 0; k < st.K; k++) { cat_offset[k] = nct; nct += st.cat_counts[k]; }

    std::vector<double> lambda(nct, 0.0);  // dual variables; init = 0
    std::vector<double> w(ct.M_cell);      // cell calibrated masses
    std::vector<double> D_eff(ct.M_cell);  // Newton weights
    std::vector<double> N(nct * nct);      // normal equations matrix (symmetric)
    std::vector<double> b(nct);            // calibration residuals / Newton RHS

    constexpr int    kMaxNewtonIters = 50;
    constexpr double kRegularization = 1e-10;

    // Track initial and best residuals for post-loop status classification:
    double initial_resid = std::numeric_limits<double>::infinity();
    double best_resid    = std::numeric_limits<double>::infinity();

    for (int iter = 0; iter < kMaxNewtonIters; iter++) {

        // (1) Compute w[c] and D_eff[c] from current lambda:
        for (int c = 0; c < ct.M_cell; c++) {
            double z = 0.0;
            for (int k = 0; k < st.K; k++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k])
                    z += lambda[cat_offset[k] + g];
            }
            // Clamp z to avoid exp overflow (z = ±700 → sigma = 0 or 1)
            z = std::clamp(z, -700.0, 700.0);
            double sig = 1.0 / (1.0 + std::exp(-z));
            double range = U_cell[c] - L_cell[c];
            w[c]     = L_cell[c] + range * sig;
            D_eff[c] = range * sig * (1.0 - sig);  // ≥ 0 always
        }

        // (2) Compute residuals b[k][j] = τ[k][j]×n − Σ_{c∈bucket(k,j)} w[c]
        std::fill(b.begin(), b.end(), 0.0);
        for (int k = 0; k < st.K; k++) {
            for (int j = 0; j < st.cat_counts[k]; j++) {
                double target = st.targets[k][j] * static_cast<double>(st.n);
                double S_kj   = 0.0;
                for (int c : cells_per_cat[k][j]) S_kj += w[c];
                b[cat_offset[k] + j] = target - S_kj;
            }
        }

        // (3) Check convergence via shared infrastructure:
        std::fill(bucket_scratch.begin(), bucket_scratch.end(), 0.0);
        lbw::CellMetrics m = lbw::compute_cell_metrics(st, ct, w_cell, W, bucket_scratch);
        bool converged = lbw::check_convergence(st.convergence_cfg, m, prev_metric, st.tol_abs);
        // NOTE: This correctly handles absolute, pct, and improvement rules —
        // including convergence=list(absolute=1e-12) used in T2/T3.
        // Do NOT hand-roll `if (max_resid < pct_tol * n)` — that only implements
        // one convergence rule and silently ignores the others.
        if (converged) {
            res.status = RK_OK;
            break;
        }

        // (4) Build N = A·diag(D_eff)·Aᵀ and solve NΔλ = b
        if (compute_normal_equations(ct, D_eff.data(), N.data(),
                                     cat_offset.data(), st.K,
                                     static_cast<size_t>(nct)) != RK_OK) {
            res.status = RK_ERR_INFEAS;  // singular system
            break;
        }
        if (ldlt_factor_inplace(N.data(), static_cast<size_t>(nct), kRegularization)
                != RK_OK) {
            res.status = RK_ERR_INFEAS;
            break;
        }
        ldlt_solve(N.data(), static_cast<size_t>(nct), b.data());  // b = Δλ

        // (5) Update lambda
        for (int j = 0; j < nct; j++) lambda[j] += b[j];
    }

    // After Newton loop exits without convergence:
    // Classify: if residual improved vs initial → BUDGET; if plateau → STALL
    if (res.status == RK_ERR_NOCONV) {
        res.status = (best_resid < initial_resid * 0.999) ? RK_ERR_BUDGET : RK_ERR_STALL;
        std::snprintf(res.message, sizeof(res.message),
            "logit: %s after %d Newton steps; best max_err=%.4e",
            res.status == RK_ERR_BUDGET ? "budget exhausted" : "stall",
            kMaxNewtonIters, res.best_error);
    }
    // NOTE: Status RK_ERR_NOCONV=1 is NEVER returned; caller sees either
    // RK_OK=0, RK_ERR_INFEAS=2, RK_ERR_BUDGET=4, or RK_ERR_STALL=5.

    // ... result population, weight reconstruction ...
}
```

## Weight reconstruction (cell → obs)

```cpp
res.best_weights.resize(st.n);
for (int i = 0; i < st.n; i++) {
    int c = ct.cell_of[i];
    double X_init_c = 0.0;
    for (int c2 = 0; c2 < ct.M_cell; c2++) { /* aggregate */ }
    // Simpler: use precomputed X_init[c]:
    res.best_weights[i] = (X_init[c] > 0.0)
        ? st.weights[i] * w[c] / X_init[c]
        : st.weights[i];
    // Result w[i] ∈ [min_weight, max_weight] cell-aggregate (guaranteed by logit)
}
```

## LogitCalibResult struct (logit_calib.hpp)

`LogitCalibResult` has IDENTICAL fields to `GregResult` (src/greg.hpp) — copy the struct
verbatim and rename. No logit-specific extra fields are added. λ (dual variables) is an
internal solver quantity — NOT exposed in the result (consistent with leafblower's
existing solvers that also discard duals). This is explicitly out of scope.

```cpp
// logit_calib.hpp
struct LogitCalibResult {
    int    status           = RK_ERR_NOCONV;
    int    iterations       = 0;
    double max_error        = 1.0;
    double best_error       = 1.0;
    int    best_iter        = 0;
    char   message[256]     = {};
    double mean_error       = 0.0;
    double kl               = 0.0;
    double chi2             = 0.0;
    double l1_weight_change = 0.0;
    double grake_norm       = 0.0;
    int    convergence_metric = 0;
    int    convergence_rule   = 0;
    double convergence_tol    = 0.0;
    int    convergence_iter   = -1;
    double convergence_solver_objective = 0.0;
    int    convergence_minimized_metric = 0;
    double sor_min_omega    = 1.0;
    int    sor_n_damped     = 0;
    int    n_bounds_violated = 0;
    int    n_bounds_clamped  = 0;
    int    n_xcur_writes_per_iter_linear = 0;
    double min_alpha_seen   = 1.0;
    double final_alpha      = 1.0;
    double alm_capacity_mu_final  = 0.0;
    int    alm_n_growth_events    = 0;
    double alm_max_dual_norm      = 0.0;
    double alm_sum_drift          = 0.0;
    std::vector<double> best_weights;
};
```

## Degenerate cell handling (L_cell[c] = U_cell[c])

When `min_weight × n_per_cell[c] = max_weight × n_per_cell[c]` (i.e., min=max, or
n_per_cell[c]=0), `D_eff[c] = (U-L)·σ·(1-σ) = 0`. The cell contributes **zero** to
the normal equations matrix N. If ALL cells in a particular bucket have D_eff=0, that
column of A·diag(D_eff)·Aᵀ is zero — N is singular.

**Handling**: `ldlt_factor_inplace` uses Tikhonov regularization (`kRegularization=1e-10`)
which handles near-singular N gracefully. After factorization, verify the solution quality:
if `max|b| > tol` after the solve (meaning the "zero" constraint could not be matched),
return `RK_ERR_INFEAS` with message "logit: singular normal equations (degenerate bounds)".

Additionally, if `n_per_cell[c] == 0` for any cell: skip (contributes nothing to any
margin, D_eff=0 naturally).

## Files (logit calibration additions)

| File | Change |
|------|--------|
| `src/logit_calib.cpp` | New solver (logit_calibrate) |
| `src/logit_calib.hpp` | LogitCalibResult struct (full field list above) |
| `src/leafblower.h` | Add RK_ALG_LOGIT = 10; extend alg_names comment |
| `src/Makevars.in` | Add `logit_calib.cpp` to PKG_SOURCES |
| `src/r_bridge.cpp` | Add "logit" dispatch block; alg_names to 11 elements |
| `src/c_api.cpp` | Add LOGIT case |
| `R/harvest.R` | Add "logit" to map_method(); update @param method; alg_names vector |
| `man/harvest.Rd` | Auto-regenerate |
| `tests/testthat/test-calibration-solvers.R` | Add T5-T8 tests |

`alg_names` in harvest.R extended to 11:
```r
alg_names <- c("", "ieppa", "lbfgsb", "raking", "sinkhorn",
                "chebyshev", "greg", "grake", "ieppa_soft",
                "greenkhorn", "logit")
#               ^ index 9            ^ index 10
```

## TDD — Logit calibration

### T5 — Logit available and calibrates

```r
test_that("T5: logit available and calibrates", {
  set.seed(5); n <- 1000L
  df  <- data.frame(sex=factor(sample(c("M","F"),n,TRUE)),
                    age=factor(sample(c("Y","O"),n,TRUE)))
  tgt <- list(sex=c(M=0.5,F=0.5), age=c(Y=0.6,O=0.4))
  r   <- harvest(df, tgt, method="logit", max_iterations=50L)
  expect_lt(attr(r,"result")$max_error, 1e-3)
  expect_equal(attr(r,"result")$algorithm_used, "logit")
})
```

### T6 — Logit bounds by construction (Newton steps < raking rounds on K=2 tight problem)

```r
test_that("T6: logit bounds by construction (Newton steps < raking rounds on K=2 tight problem)", {
  set.seed(6); n <- 5000L
  df  <- data.frame(
    v=factor(sample(5, n, TRUE)),
    g=factor(sample(c("M","F"), n, TRUE))  # second margin
  )
  tgt <- list(v=setNames(c(0.4,0.3,0.15,0.1,0.05),as.character(1:5)),
              g=c(M=0.55, F=0.45))
  r   <- harvest(df, tgt, method="logit", max_weight=1.5, min_weight=0.1)
  w   <- r$weights
  expect_true(max(w) <= 1.5 + 1e-9)
  expect_true(min(w) >= 0.1 - 1e-9)
  n_iters <- attr(r,"result")$iterations
  expect_lt(n_iters, 50L)
  # K=2 means raking needs ≥2 rounds/sweeps; logit Newton should converge faster:
  r_rk <- harvest(df, tgt, method="raking", max_weight=1.5, min_weight=0.1,
                  convergence=list(absolute=1e-6))
  n_rk <- attr(r_rk,"result")$iterations
  expect_lt(n_iters, n_rk,
    label=sprintf("logit Newton (%d steps) < raking (%d rounds) on K=2 tight problem", n_iters, n_rk))
})
```

### T7 — Logit matches GREG on unconstrained problem (same chi-square minimizer)

```r
test_that("T7: logit and raking reach same calibration target (max_err < 1e-4)", {
  # Without tight bounds, logit and raking should both converge to same weights
  set.seed(7); n <- 5000L
  df  <- data.frame(
    a=factor(sample(letters[1:3],n,TRUE)),
    b=factor(sample(LETTERS[1:4],n,TRUE))
  )
  tgt <- list(a=c(a=0.3,b=0.4,c=0.3), b=c(A=0.25,B=0.25,C=0.25,D=0.25))
  r_logit <- harvest(df, tgt, method="logit",
                     convergence=list(absolute=1e-6))
  me_logit <- attr(r_logit,"result")$max_error
  expect_lt(me_logit, 1e-4)
})
```

### T8 — Logit vs raking on tight-bounds problem (logit must be competitive)

```r
test_that("T8: logit max_err within 2x of raking on tight-bounds problem", {
  set.seed(8); n <- 5000L
  df  <- data.frame(v=factor(sample(5, n, TRUE)))
  tgt <- list(v=setNames(c(0.4,0.3,0.15,0.1,0.05), as.character(1:5)))
  r_rk    <- harvest(df, tgt, method="raking", max_weight=1.8, min_weight=0,
                     convergence=list(absolute=1e-6))
  r_logit <- harvest(df, tgt, method="logit", max_weight=1.8, min_weight=0,
                     convergence=list(absolute=1e-6))
  me_rk    <- attr(r_rk,    "result")$max_error
  me_logit <- attr(r_logit, "result")$max_error
  expect_lt(me_logit, 2.0 * me_rk + 1e-6)
})
```

## Acceptance Criteria (logit additions — Part 2)

AC-L1: T5 GREEN — logit calibrates, `algorithm_used == "logit"`.
AC-L2: T6 GREEN — bounds respected; Newton iters < 50.
AC-L3: T7 GREEN — max_err < 1e-4 on 2-margin unconstrained problem.
AC-L4: T8 GREEN — logit max_err within 2× of raking on tight-bounds problem.
AC-L5: Status codes: RK_OK=0 on convergence; RK_ERR_INFEAS=2 on singular N (infeasible). Never RK_ERR_NOCONV=1.
AC-L6: `devtools::test()` FAIL count unchanged (currently 3).
AC-L7: `R CMD INSTALL --preclean .` compiles clean.
AC-L8: *(Benchmark-only)* Stepstone logit: max_err ≤ raking max_err AND convergence in <30 Newton steps.

## Key difference from greg.cpp (no active-set, no clamping)

Greg clamps weights and iterates active-set steps. Logit never clamps — the logit link
guarantees bounds analytically. When the bounds constraint is structurally infeasible
(sum(L_cell) > n or sum(U_cell) < n), the normal equations N become singular and the
solver returns RK_ERR_INFEAS with a message. Otherwise, Newton converges smoothly.

## λ warm-start note

`λ=0` initializes `w_c = L_cell[c] + (U_cell[c]-L_cell[c]) × 0.5` — the midpoint of
[L, U]. When `min_weight=0` (common), midpoint = `max_weight × n_per_cell[c] / 2`.
This is a valid interior point for Newton to start from. Large `max_weight` means the
initial weights are far from the constrained optimum — Newton still converges quickly
(typically 10-20 steps) because the logit Hessian is well-conditioned at the midpoint.

For future warm-starting from a prior solution: out of scope.

## Out of Scope (Part 2 — Logit)

- **λ (dual variables) returned to R**: internal solver quantity, not exposed in result.
  Leafblower's existing solvers (greg, ieppa, raking) also discard their duals.
  Users needing λ can use the `survey` package's `calib()` function.
- Per-obs logit calibration (cell-aggregate bounds only, same as bounds_mode=CELL)
- Other Deville-Särndal distances (raking=KL, truncated logit)
- Homotopy warmstart for initial λ
- SQUAREM / acceleration (not needed — Newton converges in O(10-20) steps)
- AUTO routing (explicit opt-in only)

---

## Amendment: SQUAREM acceleration for Greenkhorn (added 2026-04-29)

`method="greenkhorn"` with `accelerate=TRUE` adds SQUAREM at the **round level**:

- One round = K greedy Greenkhorn steps, margins sorted **once at round entry** (stationary F_eval)
- SQUAREM CBB step accelerates the sequence of rounds
- This is autumn's `harvest(accelerate=TRUE, scheduler="greedy")` equivalent
- Valid SQUAREM usage: F_eval is stationary because sort order is fixed at entry, not re-computed per step

Differs from leafblower's R8 fix (greedy disabled under raking SQUAREM) because:
- R8 issue: greedy re-sorts inside F_eval → different calls get different orders → non-stationary
- Greenkhorn fix: sort once at F_eval entry → same input → same sort → same K updates → stationary

Default: `accelerate=FALSE` (pure Greenkhorn, one margin per step). When `accelerate=TRUE`: round-level SQUAREM with greedy-sorted F_eval.
