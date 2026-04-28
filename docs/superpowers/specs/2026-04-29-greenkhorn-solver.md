# Greenkhorn Calibration Solver — Design Spec

**Date**: 2026-04-29
**Status**: Pending design review
**Motivation**: Implement autumn's greedy-coordinate-descent calibration (Greenkhorn algorithm)
  in leafblower's C++ / cell-table architecture, matching leafblower's performance patterns.

---

## Problem

Autumn achieves `max_err=1.60e-3` on stepstone (K=9, n=1.58M) in 33.6s by using
greedy margin selection — picking the ONE margin with highest residual per step
(Greenkhorn). Leafblower's raking achieves the SAME `max_err=1.60e-3` in 2.5s but via
round-robin: it sweeps ALL K margins per iteration. The architectures differ:

| Property | Leafblower raking | Autumn (Greenkhorn) |
|---|---|---|
| Per-step cost | O(K × M_cell) | O(K + bucket_size × K) |
| Steps to converge | Fewer (full rounds) | More (single-margin steps) |
| Margin focus | Uniform across K | Greedy: hardest first |
| Simple margin polish | Not guaranteed | Machine precision possible |
| SQUAREM compatible | Yes | No (non-stationary F) |

Leafblower already beats autumn on wall time (13×). The value of implementing
Greenkhorn is NOT speed (we already win) but:
1. **Correctness parity**: simple margins polished to machine precision
2. **Alternative convergence path**: useful when one or two margins dominate residual
3. **Research baseline**: enables comparison with Greenkhorn variants (P-B scheduler
   from design docs; per-margin adaptive step sizes)
4. **Autumn drop-in**: `method="greenkhorn"` as a documented algorithm

---

## Algorithm

**Greenkhorn IPF** (Altschuler-Weed-Rigollet 2017, arxiv:1705.09634):

```
Initialize:
  X[c] = X_init[c]  ∀ c ∈ {1…M_cell}
  S[k][j] = Σ_{c ∈ bucket(k,j)} X[c]   ∀ k,j   (maintained incrementally)
  W = Σ_c X[c] = n                               (total mass = n)

Repeat until max_k errRp_k < tol OR iter = max_iter:
  (1) k* = argmax_k  max_j |S[k][j]/W − τ[k][j]|
                       errRp_k

  (2) For each j ∈ {0…cat_counts[k*]−1}:
        f_kj = (W × τ[k*][j]) / S[k*][j]         (scale factor; skip if S ≈ 0)
        For each c ∈ bucket(k*, j):
          δ_c  = X[c] × (f_kj − 1)               (weight change)
          X_new = clamp(X[c] + δ_c, L_cell[c], U_cell[c])   (bounded clamp)
          δ_actual = X_new − X[c]
          X[c] = X_new

          // Incremental bucket-sum update for ALL margins affected by cell c:
          For each k' ≠ k*:
            g = g_per_cell[k'][c]
            if g ≥ 0: S[k'][g] += δ_actual

        S[k*][j] = Σ_{c ∈ bucket(k*,j)} X[c]     (recompute exactly for margin k*)

  (3) iter++

Return X[c], reconstruct per-obs weights: w[i] = X[cell_of[i]] × w_init[i] / X_init[cell_of[i]]
```

**Convergence guarantee** (A-W-R 2017 Thm 1): For the unconstrained case,
Greenkhorn achieves KL(π | π*) ≤ ε in O(K/ε²) steps, same complexity as alternating
Sinkhorn in O(1/ε²) rounds (each round = K steps). Bounded case: same order.

---

## Why incremental bucket-sum update is critical

After updating cells in bucket (k*, j), each cell c in that bucket belongs to ONE
bucket in each of the other K−1 margins. So the incremental update is:

- Step cost: O(|bucket(k*,j)| × K)
- Compare with full recompute: O(M_cell × K)
- Ratio: |bucket(k*,j)| / M_cell ≪ 1 for sparse problems

For stepstone: M_cell ≈ 100k+, but |bucket(k*,j)| for a simple margin like
`rk_gender` (j=M or F) = M_cell / 2 ≈ 50k. For complex margins like
`rk_i_loc_time_age10_gender`, individual buckets ≈ 100-1000 cells.

The argmax over K margins after each step requires recomputing errRp only for the K−1
margins whose bucket sums changed (all do when the update δ is nonzero). Full errRp
recompute is O(K + Σ_k n_j^k) = O(M_cell) per step — still fast.

---

## Bounds integration

Bounds enforced per-cell inside the update loop:
```
X_new = clamp(X[c] × f_kj, L_cell[c], U_cell[c])
δ_actual = X_new − X[c]  (may differ from X[c]×(f_kj−1) due to clamp)
```

After clamping, the bucket sum S[k*][j] will be ≤ target (if caps hit), which
prevents exact margin matching — this is the same constraint as raking's water-fill.

When all cells in a bucket are at their upper bound (`Σ U_cell < target_mass`): the
update leaves S[k*][j] at max achievable mass. errRp_k* remains nonzero; Greenkhorn
will revisit k* at the next step, but no progress can be made. Detect and break.

**Structural infeasibility detection**: identical to raking — track is_infeasible flag
in the bucket update loop. Post-loop: if infeasible AND stalled → RK_ERR_STALL (the
capacity constraint is binding but the problem is still feasible in the raking sense).

---

## Data structures

### GreenkornState (solver-local)
```cpp
struct GreenkornState {
    // Bucket sums: S[k][j] = Σ_{c ∈ bucket(k,j)} X[c]
    // Layout: flat vector, indexed as S[k * max_cats + j]
    std::vector<double> S;
    int S_stride;  // max(cat_counts[k]) for flat indexing

    // Per-margin residuals (maintained incrementally)
    std::vector<double> errRp;  // size K

    // Cell masses
    std::vector<double> X;  // size M_cell

    // Best iterate
    double best_errRp;
    int    best_iter;
    std::vector<double> X_best;  // size M_cell

    // Total weight (= n when no clipping or perfect normalization)
    double W;
};
```

### Pre-computed per-cell group lookup
`ct.g_per_cell[k][c]` — already exists in CellTable. Provides O(1) lookup of which
bucket cell c belongs to in each margin k.

---

## Implementation

### New files
- `src/greenkhorn.cpp` — solver implementation
- `src/greenkhorn.hpp` — GreenkornResult struct + `greenkhorn_solve(st)` declaration

### Modified files
- `src/leafblower.h` — add `RK_ALG_GREENKHORN = 9`
- `src/r_bridge.cpp` — add `"greenkhorn"` dispatch (same pattern as raking)
- `src/c_api.cpp` — add GREENKHORN case
- `src/Makevars.in` — add `greenkhorn.cpp` to `PKG_SOURCES`
- `R/harvest.R` — add `"greenkhorn"` to method enum and dispatch
- `tests/testthat/test-calibration-solvers.R` — add Greenkhorn tests

### greenkhorn.hpp
```cpp
#pragma once
#include "lbw_config.h"
#include "types.hpp"
#include "cell_table.hpp"
#include "calib_dispatch.hpp"

namespace lbw {

struct GreenkornResult {
    int    status         = RK_ERR_NOCONV;
    int    iterations     = 0;
    double max_error      = 1.0;
    double best_error     = 1.0;
    int    best_iter      = 0;
    char   message[256]   = {};
    std::vector<double> best_weights;  // per-obs weights at best iterate
    // Standard fields mirrored by r_bridge:
    double sor_min_omega  = 1.0;
    int    sor_n_damped   = 0;
    double convergence_solver_objective = 0.0;
    int    convergence_minimized_metric = 0;
    double convergence_tol = 0.0;
    int    convergence_iter = -1;
    double mean_error      = 0.0;
    double kl              = 0.0;
    double chi2            = 0.0;
    double l1_weight_change = 0.0;
    double grake_norm       = 0.0;
    int    convergence_metric = 0;
    int    convergence_rule   = 0;
    double min_alpha_seen     = 1.0;
    double final_alpha        = 1.0;
    int    n_bounds_violated  = 0;
    int    n_bounds_clamped   = 0;
    int    n_xcur_writes_per_iter_linear = 0;
    // ALM fields (zero for Greenkhorn)
    double alm_capacity_mu_final  = 0.0;
    int    alm_n_growth_events    = 0;
    double alm_max_dual_norm      = 0.0;
    double alm_sum_drift          = 0.0;
};

GreenkornResult greenkhorn_solve(CalibState& st);

}  // namespace lbw
```

### greenkhorn.cpp — core loop
```cpp
#include "greenkhorn.hpp"
#include "cell_table.hpp"
#include "calib_dispatch.hpp"
#include "lbw_math.hpp"
#include <limits>
#include <algorithm>
#include <cstring>
#include <cmath>

namespace lbw {

GreenkornResult greenkhorn_solve(CalibState& st) {
    GreenkornResult res;

    // ── Build cell table ──────────────────────────────────────────────────────
    CellTable ct;
    {
        std::vector<const int32_t*> gids(st.K);
        for (int k = 0; k < st.K; k++) gids[k] = st.group_ids[k];
        if (build_cell_table(st.n, st.K, gids.data(), st.cat_counts,
                             st.weights, ct) != RK_OK) {
            res.status = RK_ERR_INFEAS;
            std::snprintf(res.message, sizeof(res.message),
                "greenkhorn: cell table build failed (structural infeasibility)");
            return res;
        }
    }

    // ── Initialize cell masses from design weights ────────────────────────────
    const int  M = ct.M_cell;
    const double W = static_cast<double>(st.n);  // sum of weights = n (normalized)

    std::vector<double> X(M);
    for (int c = 0; c < M; c++) {
        double cell_w = 0.0;
        for (int i : ct.cell_members[c]) cell_w += st.weights[i];  // or use n_per_cell × mean
        X[c] = cell_w;
    }

    // Capacity bounds per cell
    std::vector<double> L_cell(M), U_cell(M);
    for (int c = 0; c < M; c++) {
        L_cell[c] = st.min_weight * ct.n_per_cell[c];
        U_cell[c] = st.max_weight * ct.n_per_cell[c];
        // Clamp to [L, U] initially
        X[c] = std::clamp(X[c], L_cell[c], U_cell[c]);
    }

    // ── Initialize bucket sums ────────────────────────────────────────────────
    // S_flat[k][j] = S[k * S_stride + j]
    int max_cats = *std::max_element(st.cat_counts, st.cat_counts + st.K);
    const int S_stride = max_cats;
    std::vector<double> S_flat(st.K * S_stride, 0.0);
    auto S = [&](int k, int j) -> double& { return S_flat[k * S_stride + j]; };

    for (int k = 0; k < st.K; k++) {
        for (int j = 0; j < st.cat_counts[k]; j++) {
            for (int c : ct.cells[k][j]) S(k,j) += X[c];
        }
    }

    // ── Initialize errRp ─────────────────────────────────────────────────────
    std::vector<double> errRp(st.K, 0.0);
    auto compute_errRp_k = [&](int k) {
        double e = 0.0;
        for (int j = 0; j < st.cat_counts[k]; j++) {
            if (j < 0) continue;
            double achieved = S(k,j) / W;
            e = std::max(e, std::abs(achieved - st.targets[k][j]));
        }
        return e;
    };
    for (int k = 0; k < st.K; k++) errRp[k] = compute_errRp_k(k);

    // ── Best-iterate tracking ─────────────────────────────────────────────────
    double best_errRp = *std::max_element(errRp.begin(), errRp.end());
    res.best_error = best_errRp;
    res.best_iter  = 0;
    std::vector<double> X_best = X;

    constexpr double kEmptyThreshold = 1e-15;
    constexpr int    kErrCheckInterval = 10;
    const CalibConvergenceCfg& cfg = st.convergence_cfg;
    double prev_metric = std::numeric_limits<double>::infinity();

    // ── Main Greenkhorn loop ──────────────────────────────────────────────────
    for (int iter = 0; iter < st.inner_max_iter; iter++) {
        // (1) Pick k* = argmax_k errRp_k
        int k_star = std::max_element(errRp.begin(), errRp.end()) - errRp.begin();

        // Early exit if already converged
        if (errRp[k_star] < cfg.pct_tol && cfg.rule == CalibRule::THRESHOLD) {
            res.status = RK_OK;
            res.convergence_iter = iter;
            break;
        }

        // (2) Apply single-margin IPF update for margin k*
        bool any_capped = false;
        for (int j = 0; j < st.cat_counts[k_star]; j++) {
            double S_kj = S(k_star, j);
            double target_mass = st.targets[k_star][j] * W;

            if (S_kj < kEmptyThreshold * W) continue;  // structural empty

            double f = target_mass / S_kj;  // scale factor

            for (int c : ct.cells[k_star][j]) {
                double X_old = X[c];
                double X_new = std::clamp(X[c] * f, L_cell[c], U_cell[c]);
                double delta = X_new - X_old;
                if (X_new >= U_cell[c] * (1.0 - 1e-12)) any_capped = true;
                X[c] = X_new;

                if (std::abs(delta) < 1e-300) continue;  // no change

                // Incremental bucket-sum update for all other margins
                for (int k2 = 0; k2 < st.K; k2++) {
                    if (k2 == k_star) continue;
                    int g2 = ct.g_per_cell[k2][c];
                    if (g2 >= 0 && g2 < st.cat_counts[k2]) S(k2, g2) += delta;
                }
            }
            // Recompute S(k*, j) exactly (may differ from target due to clamp)
            S(k_star, j) = 0.0;
            for (int c : ct.cells[k_star][j]) S(k_star, j) += X[c];
        }

        // (3) Update errRp for k* and ALL other margins (incremental: all changed)
        for (int k = 0; k < st.K; k++) errRp[k] = compute_errRp_k(k);

        // (4) Convergence check (every kErrCheckInterval steps)
        res.iterations = iter + 1;
        double curr_max_errRp = *std::max_element(errRp.begin(), errRp.end());

        if (curr_max_errRp < best_errRp) {
            best_errRp = curr_max_errRp;
            res.best_error = best_errRp;
            res.best_iter  = iter + 1;
            X_best = X;
        }

        if ((iter + 1) % kErrCheckInterval == 0 || iter == st.inner_max_iter - 1) {
            bool converged = apply_convergence_rule(cfg, curr_max_errRp, prev_metric);
            if (converged) {
                res.status = RK_OK;
                res.convergence_iter = iter + 1;
                break;
            }
            prev_metric = curr_max_errRp;

            if (st.verbose >= 1) {
                char msg[128];
                std::snprintf(msg, sizeof(msg),
                    "[greenkhorn] iter=%d k*=%d errRp=%.4e",
                    iter + 1, k_star, curr_max_errRp);
                st.log(msg);
            }
        }
    }

    // ── Status classification ─────────────────────────────────────────────────
    if (res.status == RK_ERR_NOCONV) {
        // Distinguish BUDGET (still improving) vs STALL (plateau)
        double final_errRp = *std::max_element(errRp.begin(), errRp.end());
        res.status = (final_errRp < prev_metric) ? RK_ERR_BUDGET : RK_ERR_STALL;
        std::snprintf(res.message, sizeof(res.message),
            "greenkhorn: %s after %d steps; best max_err=%.4e",
            (res.status == RK_ERR_BUDGET) ? "budget exhausted" : "stall",
            res.iterations, res.best_error);
    }

    // ── Reconstruct per-obs weights from best iterate ─────────────────────────
    res.best_weights.resize(st.n);
    for (int i = 0; i < st.n; i++) {
        int c = ct.cell_of[i];
        double X_init_c = 0.0;
        for (int i2 : ct.cell_members[c]) X_init_c += st.weights[i2];
        res.best_weights[i] = (X_init_c > 0.0)
            ? st.weights[i] * X_best[c] / X_init_c
            : st.weights[i];
    }

    // ── Final metrics ─────────────────────────────────────────────────────────
    res.max_error = *std::max_element(errRp.begin(), errRp.end());
    res.mean_error = std::accumulate(errRp.begin(), errRp.end(), 0.0) / st.K;

    return res;
}

}  // namespace lbw
```

---

## Performance note: cell_members vs cells[k][j]

The current `CellTable` stores `cells[k][j]` — the list of cell indices in each bucket.
For Greenkhorn's inner loop this is all we need. We do NOT need `cell_members[c]` (the
list of observation indices within each cell) at the inner loop — it's only needed at
the end for weight reconstruction.

The weight reconstruction O(n) is done once at the end. The inner loop is O(M_cell×K)
per step for errRp + O(|bucket| × K) for the incremental update — much less than O(n×K).

---

## R API

```r
# harvest.R: add to method dispatch
harvest(df, target,
        method = "greenkhorn",   # new method
        max_weight = 5,
        min_weight = 0,
        max_iterations = 5000L,
        convergence = list(pct = 1e-4),
        ...)
```

- **Unsupported params** (silently ignored with warning if supplied):
  `accelerate`, `scheduler`, `homotopy_*`, `eta_*`, `sor_*`, `capacity_penalty`
- **Supported params**: `max_iterations`, `max_weight`, `min_weight`,
  `convergence`, `verbose`, `bounds_mode`, `start_weights`, `design_weights`

---

## TDD (RED before implementation)

### T1 — Greenkhorn available and calibrates
```r
test_that("greenkhorn: method available and produces calibrated weights", {
  set.seed(1); n <- 1000L
  df <- data.frame(
    sex = factor(sample(c("M","F"), n, TRUE)),
    age = factor(sample(c("Y","O"), n, TRUE))
  )
  tgt <- list(sex=c(M=0.5,F=0.5), age=c(Y=0.6,O=0.4))
  r <- harvest(df, tgt, method="greenkhorn", max_iterations=500L)
  expect_lt(attr(r,"result")$max_error, 1e-3)
  expect_equal(attr(r,"result")$algorithm_used, 9L)  # RK_ALG_GREENKHORN
})
```

### T2 — Greenkhorn polishes simple margins to machine precision
```r
test_that("greenkhorn: simple 2-cat margin reaches machine precision", {
  set.seed(2); n <- 5000L
  df <- data.frame(g = factor(sample(c("A","B"), n, TRUE, prob=c(0.3,0.7))))
  tgt <- list(g = c(A=0.5, B=0.5))
  r <- harvest(df, tgt, method="greenkhorn", max_iterations=1000L)
  me <- attr(r,"result")$max_error
  # Greenkhorn should reach near-machine precision on a single 2-cat margin
  expect_lt(me, 1e-10,
    label=sprintf("greenkhorn max_err=%.2e, expected < 1e-10 on trivial 2-cat margin", me))
})
```

### T3 — Greenkhorn vs raking: same quality on stepstone-like problem
```r
test_that("greenkhorn: max_err within 2x of raking on multi-margin problem", {
  set.seed(42); n <- 10000L
  df <- data.frame(
    a = factor(sample(letters[1:3], n, TRUE)),
    b = factor(sample(LETTERS[1:4], n, TRUE)),
    c = factor(sample(c("x","y"), n, TRUE))
  )
  tgt <- list(a=c(a=0.3,b=0.4,c=0.3), b=c(A=0.25,B=0.25,C=0.25,D=0.25), c=c(x=0.6,y=0.4))
  r_raking <- harvest(df, tgt, method="raking")
  r_grk    <- harvest(df, tgt, method="greenkhorn")
  me_raking <- attr(r_raking,"result")$max_error
  me_grk    <- attr(r_grk,"result")$max_error
  expect_lt(me_grk, 2 * me_raking + 1e-6,
    label=sprintf("greenkhorn=%.2e, raking=%.2e; must be within 2x", me_grk, me_raking))
})
```

### T4 — Greenkhorn respects bounds
```r
test_that("greenkhorn: final weights respect max_weight and min_weight", {
  set.seed(5); n <- 2000L
  df <- data.frame(v=factor(sample(5, n, TRUE)))
  tgt <- list(v=setNames(c(0.4,0.3,0.15,0.1,0.05), as.character(1:5)))
  r <- harvest(df, tgt, method="greenkhorn", max_weight=2.0, min_weight=0.1)
  w <- r$weights
  expect_true(max(w) <= 2.0 + 1e-9)
  expect_true(min(w) >= 0.1 - 1e-9)
})
```

---

## Acceptance Criteria

1. `method="greenkhorn"` available in harvest() — no error on valid input.
2. **T1 GREEN**: max_error < 1e-3 on 2-margin 2-cat problem.
3. **T2 GREEN**: max_error < 1e-10 on trivial 1-margin 2-cat problem (machine precision).
4. **T3 GREEN**: max_err within 2× of raking on 3-margin mixed problem.
5. **T4 GREEN**: bounds respected exactly.
6. Status codes: RK_OK=0 on convergence; RK_ERR_BUDGET=4 on iter exhaustion with improving trend; RK_ERR_STALL=5 on plateau.
7. `devtools::test()` FAIL count unchanged (currently 3).
8. `R CMD INSTALL --preclean .` compiles clean.
9. Stepstone benchmark: greenkhorn max_err ≤ 2 × raking max_err (same quality class).
10. Wall time on stepstone: greenkhorn < 10s (faster than Python's ipfn=17.2s; autumn=33.6s).

---

## Out of Scope

- SQUAREM acceleration (incompatible with non-stationary F)
- Top-k greedy (select multiple margins per step) — simplicity first
- Log-space Greenkhorn (KL-transport variant with entropic regularization)
- Homotopy / eta scheduling (no capacity relaxation needed for IPF)
- ALM/SOR stabilization (Greenkhorn has no oscillation; stable by construction)
- auto-routing: AUTO shall NOT select Greenkhorn (explicit opt-in only)
