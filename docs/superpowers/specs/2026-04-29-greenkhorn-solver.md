# Greenkhorn Calibration Solver — Design Spec (rev 2)

**Date**: 2026-04-29
**Status**: Pending design review (rev 2 addresses all round-1 gate findings)

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

  (2) For each j ∈ {0…cat_counts[k*]−1}:
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
  (4) iter++

Return best-iterate weights (per-obs reconstruction from X_best)
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

## Files to create / modify

| File | Change |
|------|--------|
| `src/greenkhorn.cpp` | New solver (greenkhorn_solve) |
| `src/greenkhorn.hpp` | GreenkornResult struct (mirrors RakingResult) |
| `src/leafblower.h` | Add RK_ALG_GREENKHORN = 9; add "greenkhorn" to rk_algorithm_t comment |
| `src/Makevars.in` | Add `greenkhorn.cpp` to PKG_SOURCES |
| `src/r_bridge.cpp` | Add "greenkhorn" dispatch block; extend alg_names to 10 elements in R result |
| `src/c_api.cpp` | Add GREENKHORN case |
| `R/harvest.R` | Add "greenkhorn" to map_method(); update @param method; update @param scheduler; extend alg_names vector |
| `man/harvest.Rd` | Auto-regenerate via devtools::document() |
| `tests/testthat/test-calibration-solvers.R` | Add T1-T4 tests |

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

## Out of Scope

- SQUAREM (non-stationary F)
- Top-k multi-margin selection per step
- Log-space Greenkhorn (KL-transport entropic regularization)
- Homotopy / ALM / SOR
- AUTO routing (Greenkhorn is explicit opt-in only)
- Autumn API compatibility (convergence_used field differences — out of scope)
