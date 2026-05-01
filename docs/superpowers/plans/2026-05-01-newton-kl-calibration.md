# Newton-KL Calibration Solver — Implementation Plan

> For agentic workers: REQUIRED SUB-SKILL: superpowers:subagent-driven-development

**Goal:** Implement method="newton_kl" in harvest() — standalone Newton solver on the KL dual for the zero-compression regime (K≥5, M_cell/n≥0.9). Target < 2s on kk1204 (n=1M, K=20).

**Architecture:** Single obs-level pass accumulates n_λ×n_λ Fisher Hessian; LDLT solve; Armijo line search. Reuses ldlt_factor_inplace/ldlt_solve from calib_linalg.hpp and mark_converged from calib_dispatch.hpp.

**Tech Stack:** C++17, R package build, testthat 3.

**Spec:** `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md`

**Beads epics + tasks:**
- Epic A — Core solver: `leafblower-3wp9`
  - N1: `leafblower-hvbw` — infrastructure (enum + Makevars + c_api dispatch arm)
  - N2: `leafblower-197c` — implementation (newton_calib.hpp/.cpp algorithm)
- Epic B — Integration + tests: `leafblower-zpkd`
  - N3: `leafblower-0qp0` — r_bridge.cpp dispatch + AUTO + result fields
  - N4: `leafblower-2qhw` — R/harvest.R match.arg + docstring
  - N5: `leafblower-58ce` — test-newton-kl.R 6-test suite
  - N6: `leafblower-cn5s` — kk1204 benchmark gate

**Dependency order (linear, no parallelism):** N1 → N2 → N3 → N4 → N5 → N6.

---

## N1 — Infrastructure scaffolding (`leafblower-hvbw`)

**Files:** `src/leafblower.h`, `src/Makevars.in`, `src/c_api.cpp`, NEW `src/newton_calib.hpp`, NEW `src/newton_calib.cpp` (stub)

### Edit 1.1 — `src/leafblower.h:50` enum

**Before** (lines 40–51):
```c
typedef enum {
    RK_ALG_AUTO   = 0,
    RK_ALG_IEPPA  = 1,
    RK_ALG_LBFGSB = 2,
    RK_ALG_RAKING    = 3,
    RK_ALG_SINKHORN  = 4,
    RK_ALG_CHEBYSHEV = 5,
    RK_ALG_GREG      = 6,
    RK_ALG_IEPPA_SOFT = 8,   /* ieppa + ADMM soft capacity enforcement */
    RK_ALG_GREENKHORN = 9,   /* greedy coordinate-descent IPF (autumn::harvest style) */
    RK_ALG_LOGIT      = 10   /* Deville-Sarndal 1992 logit Newton calibration (autumn::calibrate style) */
} rk_algorithm_t;
```

**After:**
```c
typedef enum {
    RK_ALG_AUTO   = 0,
    RK_ALG_IEPPA  = 1,
    RK_ALG_LBFGSB = 2,
    RK_ALG_RAKING    = 3,
    RK_ALG_SINKHORN  = 4,
    RK_ALG_CHEBYSHEV = 5,
    RK_ALG_GREG      = 6,
    RK_ALG_IEPPA_SOFT = 8,   /* ieppa + ADMM soft capacity enforcement */
    RK_ALG_GREENKHORN = 9,   /* greedy coordinate-descent IPF (autumn::harvest style) */
    RK_ALG_LOGIT      = 10,  /* Deville-Sarndal 1992 logit Newton calibration (autumn::calibrate style) */
    RK_ALG_NEWTON_KL  = 11   /* Newton on KL dual; obs-level single-pass Hessian */
} rk_algorithm_t;
```

### Edit 1.2 — `src/Makevars.in:3` PKG_SOURCES

**Before:**
```
PKG_SOURCES = c_api.cpp logit.cpp lbfgsb_solver.cpp ieppa.cpp cell_table.cpp r_bridge.cpp calib_validate.cpp sinkhorn.cpp calib_linalg.cpp greg.cpp chebyshev.cpp greenkhorn.cpp logit_calib.cpp
```

**After:**
```
PKG_SOURCES = c_api.cpp logit.cpp lbfgsb_solver.cpp ieppa.cpp cell_table.cpp r_bridge.cpp calib_validate.cpp sinkhorn.cpp calib_linalg.cpp greg.cpp chebyshev.cpp greenkhorn.cpp logit_calib.cpp newton_calib.cpp
```

### Edit 1.3 — `src/c_api.cpp:170` enum dispatch

**Before** (line 170):
```cpp
case RK_ALG_LOGIT:       alg = RK_ALG_LOGIT;       break;
```

**After:**
```cpp
case RK_ALG_LOGIT:       alg = RK_ALG_LOGIT;       break;
case RK_ALG_NEWTON_KL:   alg = RK_ALG_NEWTON_KL;   break;
```

### Edit 1.4 — `src/c_api.cpp:171–181` AUTO branch

**Before:**
```cpp
            case RK_ALG_AUTO:
            default: {
                // Route to raking when cell table is nearly incompressible (M_cell/n > 0.9).
                // At high ratios iEPPA has no compression benefit; raking is equivalent and simpler.
                int M_cell_est = lbw::estimate_M_cell(n, K, group_ids, cat_counts);
                // Exact integer comparison: M_cell_est / n > 0.9  ↔  M_cell_est * 10 > n * 9
                alg = (static_cast<int64_t>(M_cell_est) * 10 > static_cast<int64_t>(n) * 9)
                      ? RK_ALG_RAKING : RK_ALG_IEPPA;
                auto_selected = true;
                break;
            }
```

**After:**
```cpp
            case RK_ALG_AUTO:
            default: {
                // Newton-KL: K>=5 AND M_cell/n>=0.9 (zero-compression, dense margins)
                // Raking: M_cell/n>=0.9 AND K<5 (zero-compression, narrow margins)
                // iEPPA: M_cell/n<0.9 (compressible cell table)
                int M_cell_est = lbw::estimate_M_cell(n, K, group_ids, cat_counts);
                bool incompressible = (static_cast<int64_t>(M_cell_est) * 10 > static_cast<int64_t>(n) * 9);
                if (incompressible && K >= 5)      alg = RK_ALG_NEWTON_KL;
                else if (incompressible)           alg = RK_ALG_RAKING;
                else                               alg = RK_ALG_IEPPA;
                auto_selected = true;
                break;
            }
```

### Create 1.5 — `src/newton_calib.hpp`

```cpp
#pragma once
#include "calib_state.hpp"
#include "calib_dispatch.hpp"   // for CalibResult

namespace lbw {

struct NewtonCalibResult {
    CalibResult base;
    int    n_lambda  = 0;
    double dual_gap  = 0.0;
    double step_norm = 0.0;
    double line_alpha = 1.0;
};

NewtonCalibResult newton_calibrate(CalibState& st);

} // namespace lbw
```

### Create 1.6 — `src/newton_calib.cpp` (stub)

```cpp
#include "newton_calib.hpp"
#include "leafblower.h"
#include <cstring>

namespace lbw {

NewtonCalibResult newton_calibrate(CalibState& /*st*/) {
    NewtonCalibResult res;
    res.base.status = RK_ERR_BADARG;
    std::strncpy(res.message, "newton_kl unimplemented (N2)", sizeof(res.message) - 1);
    return res;
}

} // namespace lbw
```

### Edit 1.7 — `src/c_api.cpp` dispatch switch

Add `case RK_ALG_NEWTON_KL:` arm next to other solver dispatches in the main switch (calls `lbw::newton_calibrate(st)` and packs `result->status`, `result->iterations`, `result->max_error`, `result->algorithm_used = RK_ALG_NEWTON_KL`). Mirror the existing `case RK_ALG_LOGIT:` block exactly.

### Compile gate

```sh
R CMD INSTALL --preclean . 2>&1 | tail -2
```
Expected: `* DONE (leafblower)`.

### Commit

`feat(newton_kl): scaffold enum, Makevars, c_api dispatch arm`

---

## N2 — Newton solver implementation (`leafblower-197c`)

**Files:** `src/newton_calib.hpp` (extend), `src/newton_calib.cpp` (full implementation, replacing N1 stub)

### Algorithm reference (verbatim from spec §"Single-Pass Per-Step Algorithm")

```cpp
// lam_off[k] = Σ_{k'<k} (cat_counts[k']-1); n_lam = lam_off[K]
// For each obs i: u_i = Σ_k (j=group_ids[k][i]>0 ? lam[lam_off[k]+j-1] : 0)
// f_i = st.weights[i] * exp(u_i)
// Gradient: G[lam_off[k]+j-1] += f_i for j>0
// Hessian upper triangle: H[lam_off[k1]+j1-1][lam_off[k2]+j2-1] += f_i (j1>0,j2>0)
// Normalize: G/=Z; H/=Z; H -= G⊗G; G -= T
// Solve: LDLT(H) δ = G → Newton step
// Armijo: find α such that g(λ+αδ) > g(λ)
```

### Edit 2.1 — `src/newton_calib.cpp` (replace stub)

```cpp
#include "newton_calib.hpp"
#include "leafblower.h"
#include "calib_linalg.hpp"
#include "calib_dispatch.hpp"
#include <vector>
#include <cmath>
#include <cstring>
#include <algorithm>

namespace lbw {

NewtonCalibResult newton_calibrate(CalibState& st) {
    NewtonCalibResult res;
    res.base.status = RK_ERR_BUDGET;  // assume budget unless converged

    // 1. Setup — cell-table validation preamble (do NOT use ct in the Newton loop).
    CellTable ct;
    std::vector<double> X_init, L_cell, U_cell;
    std::vector<int> cat_offset;
    double hi_eff;
    int n_cats_total = 0;
    if (solver_setup_ct(st, ct, X_init, hi_eff, L_cell, U_cell,
                        cat_offset, n_cats_total, res) != RK_OK) {
        return res;
    }

    const int n = st.n;
    const int K = st.K;

    // 2. Reference elimination — fix lam_{k,0} = 0; free vars: n_lam = Σ (cat_counts[k]-1).
    std::vector<int> lam_off(K + 1, 0);
    for (int k = 0; k < K; ++k)
        lam_off[k + 1] = lam_off[k] + std::max(0, st.cat_counts[k] - 1);
    const int n_lam = lam_off[K];
    res.n_lambda = n_lam;
    if (n_lam <= 0 || n_lam > 256) {
        res.base.status = RK_ERR_BADARG;
        std::strncpy(res.message, "newton_kl: n_lambda out of [1,256]", sizeof(res.message) - 1);
        return res;
    }

    // 3. State allocation.
    std::vector<double> lam(n_lam, 0.0);
    std::vector<double> H(static_cast<size_t>(n_lam) * n_lam, 0.0);
    std::vector<double> G(n_lam, 0.0);
    std::vector<double> delta(n_lam, 0.0);
    std::vector<double> T(n_lam, 0.0);
    for (int k = 0; k < K; ++k) {
        for (int j = 1; j < st.cat_counts[k]; ++j)
            T[lam_off[k] + j - 1] = st.targets[k][j];
    }

    const int max_iter = std::min(st.outer_max_iter, 50);
    const double tol_abs = st.tol_abs;
    int iter = 0;
    double Z_curr = 0.0;
    double g_curr = 0.0;  // dual objective value

    auto eval_dual = [&](const std::vector<double>& lam_eval, double& Z_out) -> double {
        double Z = 0.0;
        for (int i = 0; i < n; ++i) {
            double u = 0.0;
            for (int k = 0; k < K; ++k) {
                int j = st.group_ids[k][i];
                if (j > 0) u += lam_eval[lam_off[k] + j - 1];
            }
            Z += st.weights[i] * std::exp(u);
        }
        Z_out = Z;
        // g(lam) = log Z - sum T*lam
        double Tlam = 0.0;
        for (int a = 0; a < n_lam; ++a) Tlam += T[a] * lam_eval[a];
        return std::log(Z) - Tlam;
    };

    g_curr = eval_dual(lam, Z_curr);

    for (iter = 0; iter < max_iter; ++iter) {
        // 4a. Zero accumulators.
        std::fill(G.begin(), G.end(), 0.0);
        std::fill(H.begin(), H.end(), 0.0);
        double Z = 0.0;

        // 4b. Single obs-level pass.
        for (int i = 0; i < n; ++i) {
            double u = 0.0;
            for (int k = 0; k < K; ++k) {
                int j = st.group_ids[k][i];
                if (j > 0) u += lam[lam_off[k] + j - 1];
            }
            const double f_i = st.weights[i] * std::exp(u);
            Z += f_i;
            for (int k = 0; k < K; ++k) {
                int j = st.group_ids[k][i];
                if (j > 0) G[lam_off[k] + j - 1] += f_i;
            }
            for (int k1 = 0; k1 < K; ++k1) {
                int j1 = st.group_ids[k1][i];
                if (j1 <= 0) continue;
                for (int k2 = k1; k2 < K; ++k2) {
                    int j2 = st.group_ids[k2][i];
                    if (j2 <= 0) continue;
                    H[(lam_off[k1] + j1 - 1) * n_lam + (lam_off[k2] + j2 - 1)] += f_i;
                }
            }
        }
        Z_curr = Z;

        // 4c. Normalize.
        for (int a = 0; a < n_lam; ++a) G[a] /= Z;
        for (size_t idx = 0; idx < H.size(); ++idx) H[idx] /= Z;
        // 4d. Subtract outer product (upper triangle).
        for (int a = 0; a < n_lam; ++a)
            for (int b = a; b < n_lam; ++b)
                H[a * n_lam + b] -= G[a] * G[b];
        // 4e. Mirror upper to lower.
        for (int a = 0; a < n_lam; ++a)
            for (int b = a + 1; b < n_lam; ++b)
                H[b * n_lam + a] = H[a * n_lam + b];
        // 4f. Subtract targets.
        for (int a = 0; a < n_lam; ++a) G[a] -= T[a];

        // 4g. Convergence: ||G||_inf < tol_abs.
        double dual_gap = 0.0;
        for (int a = 0; a < n_lam; ++a) dual_gap = std::max(dual_gap, std::fabs(G[a]));
        res.dual_gap = dual_gap;
        if (dual_gap < tol_abs) {
            mark_converged(res, st.conv_cfg, iter);
            break;
        }

        // 4h. Newton solve: H δ = G.
        std::copy(G.begin(), G.end(), delta.begin());
        if (ldlt_factor_inplace(H.data(), n_lam, 1e-12) != RK_OK) {
            res.base.status = RK_ERR_BADARG;
            std::strncpy(res.message, "newton_kl: ldlt_factor failed", sizeof(res.message) - 1);
            return res;
        }
        ldlt_solve(H.data(), n_lam, delta.data());

        // 4i. Armijo line search (gradient-ascent direction = -delta).
        double g_dot_d = 0.0;
        for (int a = 0; a < n_lam; ++a) g_dot_d += G[a] * delta[a];
        double alpha = 1.0;
        std::vector<double> lam_trial(n_lam);
        double Z_trial = 0.0, g_trial = 0.0;
        const double c1 = 1e-4;
        bool accepted = false;
        for (int ls = 0; ls < 20; ++ls) {
            for (int a = 0; a < n_lam; ++a) lam_trial[a] = lam[a] - alpha * delta[a];
            g_trial = eval_dual(lam_trial, Z_trial);
            if (g_trial >= g_curr - c1 * alpha * g_dot_d) { accepted = true; break; }
            alpha *= 0.5;
        }
        if (!accepted) alpha = 0.0;  // no descent; will trigger step_norm break next iter
        res.line_alpha = alpha;

        // 4j. Update.
        double step2 = 0.0;
        for (int a = 0; a < n_lam; ++a) {
            double d = alpha * delta[a];
            lam[a] -= d;
            step2 += d * d;
        }
        res.step_norm = std::sqrt(step2);
        g_curr = g_trial;
        Z_curr = Z_trial;

        // 4k. Step-norm secondary stop.
        if (res.step_norm < 1e-12) {
            mark_converged(res, st.conv_cfg, iter);
            break;
        }
    }
    res.base.iterations = iter + 1;

    // 5. Compute weights w_i = f_i / Z * n  and bounds-fallback check.
    std::vector<double> w(n, 0.0);
    int n_violated = 0;
    for (int i = 0; i < n; ++i) {
        double u = 0.0;
        for (int k = 0; k < K; ++k) {
            int j = st.group_ids[k][i];
            if (j > 0) u += lam[lam_off[k] + j - 1];
        }
        double w_i = (st.weights[i] * std::exp(u)) / Z_curr * n;
        w[i] = w_i;
        if (w_i > st.max_weight || w_i < st.min_weight) ++n_violated;
    }
    const double frac_violated = static_cast<double>(n_violated) / n;
    if (frac_violated > 0.05) {
        res.base.status = RK_ERR_NOCONV;
        std::snprintf(res.message, sizeof(res.message),
            "newton_kl: %.1f%% obs violate bounds; Newton solution may be inadmissible",
            frac_violated * 100.0);
    }

    // 6. Write weights and compute max_error (per-margin achieved vs target).
    for (int i = 0; i < n; ++i) st.weights[i] = w[i];
    double max_err = 0.0;
    for (int k = 0; k < K; ++k) {
        std::vector<double> achieved(st.cat_counts[k], 0.0);
        double total = 0.0;
        for (int i = 0; i < n; ++i) {
            int j = st.group_ids[k][i];
            if (j >= 0) { achieved[j] += w[i]; total += w[i]; }
        }
        for (int j = 0; j < st.cat_counts[k]; ++j) {
            double err = std::fabs(achieved[j] / total - st.targets[k][j]);
            if (err > max_err) max_err = err;
        }
    }
    res.base.max_error = max_err;
    res.base.best_weights = std::move(w);
    return res;
}

} // namespace lbw
```

### Compile gate

```sh
R CMD INSTALL --preclean . 2>&1 | tail -2
```

### Commit

`feat(newton_kl): obs-level Newton solver on KL dual`

---

## N3 — r_bridge.cpp dispatch + AUTO + result fields (`leafblower-0qp0`)

**File:** `src/r_bridge.cpp`

### Edit 3.1 — `src/r_bridge.cpp:24-34` kAlgMap

**Before:**
```cpp
const std::unordered_map<std::string_view, rk_algorithm_t> kAlgMap = {
    {"ieppa",      RK_ALG_IEPPA},
    {"ieppa_soft", RK_ALG_IEPPA_SOFT},
    {"lbfgsb",     RK_ALG_LBFGSB},
    {"raking",     RK_ALG_RAKING},
    {"greg",       RK_ALG_GREG},
    {"chebyshev",  RK_ALG_CHEBYSHEV},
    {"sinkhorn",   RK_ALG_SINKHORN},
    {"auto",       RK_ALG_AUTO},
    {"greenkhorn", RK_ALG_GREENKHORN},
    {"logit",      RK_ALG_LOGIT},
};
```

**After:**
```cpp
const std::unordered_map<std::string_view, rk_algorithm_t> kAlgMap = {
    {"ieppa",      RK_ALG_IEPPA},
    {"ieppa_soft", RK_ALG_IEPPA_SOFT},
    {"lbfgsb",     RK_ALG_LBFGSB},
    {"raking",     RK_ALG_RAKING},
    {"greg",       RK_ALG_GREG},
    {"chebyshev",  RK_ALG_CHEBYSHEV},
    {"sinkhorn",   RK_ALG_SINKHORN},
    {"auto",       RK_ALG_AUTO},
    {"greenkhorn", RK_ALG_GREENKHORN},
    {"logit",      RK_ALG_LOGIT},
    {"newton_kl",  RK_ALG_NEWTON_KL},
};
```

Add at the top of r_bridge.cpp near other solver headers: `#include "newton_calib.hpp"`.

### Edit 3.2 — `src/r_bridge.cpp:425-480` AUTO arm three-way

**Before** (key lines 426–448):
```cpp
    } else if (strcmp(method_str, "auto") == 0) {
        // AUTO routing: select raking when M_cell/n > 0.9, else iEPPA.
        int M_cell_est = lbw::estimate_M_cell(n, K,
            const_cast<const int32_t**>(group_ids.data()),
            cat_counts.data());
        // Exact integer comparison: M_cell_est / n > 0.9  ↔  M_cell_est * 10 > n * 9
        bool use_raking = (static_cast<int64_t>(M_cell_est) * 10 > static_cast<int64_t>(n) * 9);
        // Save for auto-fallback: only st.weights is mutated by solvers in-place
        const std::vector<double> weights_backup(weights);
        if (use_raking) {
            auto res = lbw::raking_solve(st);
            ...
        } else {
            st.ieppa_auto_selected = true;
            auto res = lbw::ieppa_solve(st);
            ...
        }
```

**After** (insert newton branch ahead of raking; condition reordered):
```cpp
    } else if (strcmp(method_str, "auto") == 0) {
        // AUTO routing: newton_kl when M_cell/n>=0.9 AND K>=5; raking when
        // M_cell/n>=0.9 AND K<5; iEPPA otherwise.
        int M_cell_est = lbw::estimate_M_cell(n, K,
            const_cast<const int32_t**>(group_ids.data()),
            cat_counts.data());
        bool incompressible = (static_cast<int64_t>(M_cell_est) * 10 > static_cast<int64_t>(n) * 9);
        bool use_newton = incompressible && K >= 5;
        bool use_raking = incompressible && !use_newton;
        const std::vector<double> weights_backup(weights);
        if (use_newton) {
            auto res = lbw::newton_calibrate(st);
            res_status     = res.base.status;
            res_iterations = res.base.iterations;
            res_max_error  = res.base.max_error;
            res_alg_used   = (int)RK_ALG_NEWTON_KL;
            pack_solver_result(res);
            res_n_lambda   = res.n_lambda;
            res_dual_gap   = res.dual_gap;
            res_step_norm  = res.step_norm;
            res_line_alpha = res.line_alpha;
            res_best_weights = std::move(res.base.best_weights);
        } else if (use_raking) {
            auto res = lbw::raking_solve(st);
            // ... unchanged raking branch
        } else {
            st.ieppa_auto_selected = true;
            auto res = lbw::ieppa_solve(st);
            // ... unchanged ieppa branch
        }
```

### Edit 3.3 — `src/r_bridge.cpp:527` Newton dispatch arm

Insert after the `logit` arm (line 527, before the `else` at line 528):

```cpp
    } else if (strcmp(method_str, "newton_kl") == 0) {
        auto res = lbw::newton_calibrate(st);
        pack_solver_result(res);
        res_status     = res.base.status;
        res_iterations = res.base.iterations;
        res_max_error  = res.base.max_error;
        res_alg_used   = (int)RK_ALG_NEWTON_KL;
        res_n_lambda   = res.n_lambda;
        res_dual_gap   = res.dual_gap;
        res_step_norm  = res.step_norm;
        res_line_alpha = res.line_alpha;
        if (!res.base.best_weights.empty())
            res_best_weights = std::move(res.base.best_weights);
        else
            res_best_weights.assign(st.n, 0.0);
```

### Edit 3.4 — `src/r_bridge.cpp:625-632` alg_name lookup

**Before:**
```cpp
const char* alg_name = (res_alg_used == (int)RK_ALG_LBFGSB)                       ? "L-BFGS-B"
                     : (res_alg_used == (int)RK_ALG_RAKING)                       ? "raking"
                     : (res_alg_used == static_cast<int>(RK_ALG_SINKHORN))        ? "sinkhorn"
                     : (res_alg_used == static_cast<int>(RK_ALG_GREG))            ? "greg"
                     : (res_alg_used == static_cast<int>(RK_ALG_CHEBYSHEV))       ? "chebyshev"
                     : (res_alg_used == static_cast<int>(RK_ALG_GREENKHORN))      ? "greenkhorn"
                     : (res_alg_used == static_cast<int>(RK_ALG_LOGIT))           ? "logit"
                     : "iEPPA";
```

**After:** same chain, insert `: (res_alg_used == (int)RK_ALG_NEWTON_KL)                  ? "newton_kl"` before `: "iEPPA"`.

### Edit 3.5 — Result list packing

Declare four new locals near the existing `res_*` declarations:
```cpp
int    res_n_lambda   = 0;
double res_dual_gap   = 0.0;
double res_step_norm  = 0.0;
double res_line_alpha = 1.0;
```
Append them to the result list (mirror the SET_VECTOR_ELT idiom used for `alm_capacity_mu_final`); add the four name strings to the `setNames()` call.

### Compile gate

```sh
R CMD INSTALL --preclean . 2>&1 | tail -2
```

### Commit

`feat(newton_kl): r_bridge dispatch + AUTO + result fields`

---

## N4 — R/harvest.R `method="newton_kl"` (`leafblower-2qhw`)

**File:** `R/harvest.R`

### Edit 4.1 — `R/harvest.R:577` `match.arg`

**Before:**
```r
match.arg(method, c("auto", "ieppa", "ieppa_soft", "lbfgsb", "raking", "sinkhorn", "chebyshev", "greg", "greenkhorn", "logit"))
```

**After:**
```r
match.arg(method, c("auto", "ieppa", "ieppa_soft", "lbfgsb", "raking", "sinkhorn", "chebyshev", "greg", "greenkhorn", "logit", "newton_kl"))
```

### Edit 4.2 — `R/harvest.R:505` `alg_names`

**Before:**
```r
alg_names <- c("", "ieppa", "lbfgsb", "raking", "sinkhorn", "chebyshev", "greg", "", "ieppa_soft", "greenkhorn", "logit")
```

**After:**
```r
alg_names <- c("", "ieppa", "lbfgsb", "raking", "sinkhorn", "chebyshev", "greg", "", "ieppa_soft", "greenkhorn", "logit", "newton_kl")
```

### Edit 4.3 — `@param method` docstring

After the `\code{"logit"}` entry (around line 30), append:
```
#'   \code{"newton_kl"} (Newton solver on the KL dual; obs-level single-pass
#'   Hessian; targets the zero-compression regime K >= 5, M_cell/n >= 0.9;
#'   typically converges in 5-10 quadratic steps; falls back to RK_ERR_NOCONV
#'   when more than 5% of observations violate weight bounds).
```

### Edit 4.4 — `@param method` AUTO comment

**Before (line 16-17):**
```
#' @param method Calibration method. One of \code{"auto"} (default: iEPPA or
#'   raking based on M_cell/n ratio), \code{"ieppa"} ...
```

**After:**
```
#' @param method Calibration method. One of \code{"auto"} (default:
#'   \code{"newton_kl"} when M_cell/n >= 0.9 AND K >= 5; \code{"raking"} when
#'   M_cell/n >= 0.9 AND K < 5; \code{"ieppa"} otherwise), \code{"ieppa"} ...
```

### Compile / smoke gate

```sh
R CMD INSTALL --preclean . 2>&1 | tail -2
Rscript -e 'library(leafblower); df <- data.frame(x=factor(c("a","b","a","b"))); harvest(df, list(x=c(a=0.5,b=0.5)), method="newton_kl")' 2>&1 | tail -5
```

### Commit

`docs+feat(harvest): expose method="newton_kl"`

---

## N5 — `tests/testthat/test-newton-kl.R` (`leafblower-58ce`)

**File (new):** `tests/testthat/test-newton-kl.R`

### Content

```r
# Newton-KL Calibration Tests
# Spec: docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md

test_that("newton_kl K=3 small problem converges", {
  set.seed(1); n <- 5000
  df <- data.frame(
    a = factor(sample(letters[1:3], n, TRUE)),
    b = factor(sample(letters[1:3], n, TRUE)),
    c = factor(sample(letters[1:3], n, TRUE))
  )
  tgt <- list(
    a = c(a = 0.4, b = 0.3, c = 0.3),
    b = c(a = 0.5, b = 0.3, c = 0.2),
    c = c(a = 0.33, b = 0.33, c = 0.34)
  )
  r <- harvest(df, tgt, method = "newton_kl")
  expect_equal(attr(r, "result")$status, 0L)
  expect_lt(attr(r, "result")$max_error, 1e-6)
})

test_that("newton_kl K=9 stepstone-like fixture", {
  skip_if_not(file.exists(test_path("fixtures", "stepstone_small.rds")))
  fix <- readRDS(test_path("fixtures", "stepstone_small.rds"))
  r <- harvest(fix$df, fix$tgt, method = "newton_kl")
  expect_equal(attr(r, "result")$status, 0L)
  expect_lt(attr(r, "result")$max_error, 1e-4)
})

test_that("newton_kl bounds fallback returns RK_ERR_NOCONV when >5% violated", {
  set.seed(2); n <- 5000
  df <- data.frame(x = factor(sample(letters[1:3], n, TRUE,
                                     prob = c(0.85, 0.10, 0.05))))
  tgt <- list(x = c(a = 0.05, b = 0.05, c = 0.90))
  r <- harvest(df, tgt, method = "newton_kl", max_weight = 2)
  expect_equal(attr(r, "result")$status, 1L)  # RK_ERR_NOCONV
})

test_that("newton_kl produces lower max_error than greg (KL vs chi2)", {
  set.seed(1); n <- 5000
  df <- data.frame(
    a = factor(sample(letters[1:3], n, TRUE)),
    b = factor(sample(letters[1:3], n, TRUE))
  )
  tgt <- list(
    a = c(a = 0.4, b = 0.3, c = 0.3),
    b = c(a = 0.5, b = 0.3, c = 0.2)
  )
  r_nk <- harvest(df, tgt, method = "newton_kl")
  r_gr <- harvest(df, tgt, method = "greg")
  expect_lt(attr(r_nk, "result")$max_error,
            attr(r_gr, "result")$max_error)
})

test_that("auto routes to newton_kl when K>=5 and M_cell/n>=0.9", {
  set.seed(3); n <- 2000
  # K=5 dense factors -> M_cell/n approaches 1
  df <- as.data.frame(lapply(1:5, function(k)
    factor(sample(letters[1:5], n, TRUE))))
  names(df) <- paste0("v", 1:5)
  tgt <- lapply(df, function(f) {
    p <- runif(5); p <- p / sum(p)
    setNames(p, levels(f))
  })
  r <- harvest(df, tgt, method = "auto")
  expect_equal(attr(r, "algorithm"), "newton_kl")
})

test_that("kk1204 Newton-KL converges in budget (manual; skipped on CRAN)", {
  skip_on_cran()
  skip_if_not(Sys.getenv("LBW_HEAVY_TESTS") == "1")
  # Heavy fixture: defer to N6 benchmark; this is a smoke check only.
  set.seed(4); n <- 1e5  # smaller than kk1204; kk1204 lives in benchmarks/
  K <- 20; nj <- 5
  df <- as.data.frame(lapply(1:K, function(k)
    factor(sample(letters[1:nj], n, TRUE))))
  names(df) <- paste0("m", 1:K)
  tgt <- lapply(df, function(f) {
    p <- runif(nj); p <- p / sum(p)
    setNames(p, levels(f))
  })
  t0 <- Sys.time()
  r <- harvest(df, tgt, method = "newton_kl", max_weight = 3,
               max_iterations = 50)
  dt <- as.numeric(Sys.time() - t0, units = "secs")
  expect_lt(dt, 1.0)  # n=1e5 should be ≤ 1s; n=1M wall budget is N6
  expect_lt(attr(r, "result")$max_error, 1e-4)
})
```

### Test gate

```sh
Rscript -e 'library(leafblower); testthat::test_dir("tests/testthat")' 2>&1 | grep -E "FAIL|PASS" | tail -1
```
Expected: `[ FAIL 0 |` line.

### Commit

`test(newton_kl): K=3, K=9, kk1204, bounds-fallback, KL-vs-chi2`

---

## N6 — kk1204 benchmark (`leafblower-cn5s`)

**File (new):** `benchmarks/newton_kl_bench.R`

### Content

```r
# Newton-KL kk1204 Benchmark
# Spec: docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md
# Target: wall < 2s, max_error < 1e-4 on n=1M, K=20.

Sys.setenv(OMP_NUM_THREADS = "1")
suppressPackageStartupMessages({
  library(leafblower)
  library(bench)
})

set.seed(1)

# kk1204 fixture: n=1M, K=20, nj=5, max_weight=3, skewed targets.
n  <- 1e6L; K <- 20L; nj <- 5L
df <- as.data.frame(lapply(seq_len(K), function(k)
  factor(sample(letters[seq_len(nj)], n, TRUE))))
names(df) <- paste0("m", seq_len(K))

# Skewed targets: each margin has one dominant cell.
tgt <- lapply(df, function(f) {
  p <- c(0.6, 0.2, 0.1, 0.07, 0.03)  # skewed
  p <- p / sum(p)
  setNames(p, levels(f))
})

cat(sprintf("[bench] kk1204 newton_kl: n=%d K=%d nj=%d\n", n, K, nj))

res <- bench::mark(
  newton_kl = harvest(df, tgt, method = "newton_kl",
                      max_weight = 3, max_iterations = 50),
  iterations = 3, check = FALSE, memory = FALSE, filter_gc = FALSE
)

# Reference run for max_error
r <- harvest(df, tgt, method = "newton_kl", max_weight = 3, max_iterations = 50)
wall <- as.numeric(res$median)
me   <- attr(r, "result")$max_error
it   <- attr(r, "result")$iterations
cat(sprintf("[bench] wall=%.3fs max_err=%.3e iters=%d\n", wall, me, it))

dir.create("benchmarks/results", showWarnings = FALSE, recursive = TRUE)
write.csv(
  data.frame(method = "newton_kl", wall_s = wall, max_error = me,
             iterations = it, n = n, K = K, nj = nj),
  "benchmarks/results/newton_kl_kk1204.csv",
  row.names = FALSE
)

stopifnot(wall < 2.0)
stopifnot(me   < 1e-4)
cat("[bench] PASS: wall<2s and max_err<1e-4\n")
```

### Test gate

```sh
OMP_NUM_THREADS=1 Rscript benchmarks/newton_kl_bench.R 2>&1 | tail -5
```
Expected: `[bench] PASS: wall<2s and max_err<1e-4`.

### Commit

`bench(newton_kl): kk1204 < 2s gate`

---

## Final verification (after N6)

```sh
R CMD INSTALL --preclean . 2>&1 | tail -2
Rscript -e 'library(leafblower); testthat::test_dir("tests/testthat")' 2>&1 | grep -E "FAIL|PASS" | tail -1
OMP_NUM_THREADS=1 Rscript benchmarks/newton_kl_bench.R 2>&1 | tail -5
```

All three must pass before closing epic `leafblower-zpkd`.
