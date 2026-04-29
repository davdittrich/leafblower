# Epic C — Logit Calibration Solver (logit_calib.{hpp,cpp})

**Date**: 2026-04-29
**Branch**: `fix/correctness-performance-2026-04-28`
**Spec**: `docs/superpowers/specs/2026-04-29-greenkhorn-solver.md` (Part 2)
**Scope**: 2 tasks — header + full solver. NO bridges, NO Makevars, NO R changes (those land
in Epic D).

## Plan Header

- **Mechanism**: Deville–Särndal 1992 logit calibration via dual Newton–Raphson on the
  reduced normal-equations system `N·Δλ = b`, with `N = A·diag(D_eff)·Aᵀ`,
  `D_eff[c] = (U−L)·σ(z_c)·(1−σ(z_c))`, and weights reconstructed by
  `w_c = L + (U−L)·σ(z_c)`. Reuses `compute_normal_equations`, `ldlt_factor_inplace`,
  `ldlt_solve` (from `calib_linalg.hpp`) and `compute_cell_metrics`/`check_convergence`
  (from `calib_dispatch.hpp`).
- **Forbidden**:
  - NO clamping of `w[c]` (bounds are by construction).
  - NO active-set logic (logit has none — that is greg.cpp's pattern, not logit's).
  - NO water-fill / per-obs clamping.
  - NO new fields in `LogitCalibResult` beyond the verbatim copy of `GregResult` plus the
    `alm_*` zero-init fields the spec mandates.
  - NO RK_ERR_NOCONV path — post-loop must classify as RK_ERR_BUDGET or RK_ERR_STALL only.
  - NO bridge wiring (r_bridge.cpp / harvest.R / Makevars / leafblower.h enum) in this epic.
  - NO hand-rolled `if (max_resid < tol)` convergence — must use `check_convergence`.
- **Audit**: Compile-only verification via `R CMD INSTALL --preclean .` after each task.
  Tests T5–T8 from the spec are RED until Epic D lands the bridge dispatch (`map_method`,
  `alg_names`, r_bridge.cpp dispatch block). After C2, the test file should still compile
  and tests should fail with a recognizable "method not found" / dispatch error — NOT with
  link errors from `lbw::logit_calibrate` being unresolved (because Epic C does not extend
  `Makevars.in`; `logit_calib.cpp` is unreferenced TU and is not yet linked).

## Pre-conditions verified against codebase

Read confirmed:
- `src/greg.hpp` — `GregResult` field list (status, iterations, max_error, mean_error, kl,
  chi2, grake_norm, l1_weight_change, convergence_metric, convergence_rule,
  convergence_tol, convergence_iter, best_objective_seen, convergence_solver_objective,
  convergence_minimized_metric, best_error, best_iter, best_weights, M_cell, message[256]).
- `src/greg.cpp` lines 25–164 — pattern: build_cell_table → X_init → L_cell/U_cell →
  cat_offset/n_cats_total → max_cats scratch → Newton loop with compute_normal_equations →
  ldlt_factor_inplace → ldlt_solve → metrics via `compute_cell_metrics` → obs expansion via
  `apply_obs_expansion`.
- `src/calib_linalg.hpp` — exact signatures verified: `compute_normal_equations(const
  CellTable&, const double* D, double* N, const int* cat_offset, int K, size_t nct)`,
  `ldlt_factor_inplace(double*, size_t, double eps)`, `ldlt_solve(const double*, size_t,
  double*)`. All return `int` status; `ldlt_solve` is `void`.
- `src/calib_dispatch.hpp` — `compute_cell_metrics(st, ct, X, W, bucket)` returns
  `CellMetrics`; `check_convergence(cfg, m, prev_metric, tol_abs)` updates `prev_metric` in
  place. `apply_obs_expansion(ct, X, X_init, n, lo, hi, weights)` for obs reconstruction.

Spec quote (rev 5, lines 580–617) defines `LogitCalibResult` with the GregResult fields
PLUS `alm_capacity_mu_final`, `alm_n_growth_events`, `alm_max_dual_norm`, `alm_sum_drift`
zero-init. Plus solver-bookkeeping fields (`sor_min_omega`, `sor_n_damped`,
`n_bounds_violated`, `n_bounds_clamped`, `n_xcur_writes_per_iter_linear`, `min_alpha_seen`,
`final_alpha`). The header below is a verbatim transcription of the spec struct.

Confidence: 95 (exact code/signatures read from source, spec literal struct).

---

## Task C1 — `src/logit_calib.hpp`

### Goal

Create new header declaring `LogitCalibResult` struct (verbatim from spec rev 5 §580–617)
and the `logit_calibrate(CalibState&)` function.

### Exact code to write

```cpp
// src/logit_calib.hpp
#pragma once
#include "types.hpp"
#include "cell_table.hpp"
#include <vector>
#include <limits>

namespace lbw {

// Deville-Särndal 1992 bounded logit calibration.
// Bounds-by-construction via w_c = L + (U-L)·σ(z_c) — no clamping, no active set.
// Newton-Raphson on dual variables λ; uses compute_normal_equations + LDLT
// (same skeleton as greg.cpp).
struct LogitCalibResult {
    int    status              = RK_ERR_NOCONV;
    int    iterations          = 0;
    double max_error           = 1.0;
    double best_error          = 1.0;
    int    best_iter           = 0;
    char   message[256]        = {};
    double mean_error          = 0.0;
    double kl                  = 0.0;
    double chi2                = 0.0;
    double l1_weight_change    = 0.0;
    double grake_norm          = 0.0;
    int    convergence_metric  = static_cast<int>(CalibMetric::CHI2);
    int    convergence_rule    = 0;
    double convergence_tol     = 0.0;
    int    convergence_iter    = -1;
    double convergence_solver_objective = std::numeric_limits<double>::infinity();
    int    convergence_minimized_metric = static_cast<int>(CalibMetric::CHI2);
    double best_objective_seen = std::numeric_limits<double>::infinity();
    double sor_min_omega       = 1.0;
    int    sor_n_damped        = 0;
    int    n_bounds_violated   = 0;
    int    n_bounds_clamped    = 0;
    int    n_xcur_writes_per_iter_linear = 0;
    double min_alpha_seen      = 1.0;
    double final_alpha         = 1.0;
    double alm_capacity_mu_final = 0.0;
    int    alm_n_growth_events   = 0;
    double alm_max_dual_norm     = 0.0;
    double alm_sum_drift         = 0.0;
    int    M_cell                = 0;
    std::vector<double> best_weights;
};

// Newton dual ascent for logit calibration (Deville-Särndal 1992).
// Maximizes the logit-distance dual subject to margin constraints.
// w_c = L_cell[c] + (U_cell[c]-L_cell[c])·σ(z_c) guarantees bounds by construction.
// Returns RK_OK | RK_ERR_INFEAS | RK_ERR_BUDGET | RK_ERR_STALL.  Never RK_ERR_NOCONV.
LogitCalibResult logit_calibrate(CalibState& st);

} // namespace lbw
```

### Verification gate (after C1)

```bash
R CMD INSTALL --preclean .
```

Header is unreferenced (no .cpp includes it yet), so this MUST compile clean. If it fails,
fix root cause before proceeding to C2 — DO NOT pivot to "skip header" or modify Makevars.

### Commit (C1)

```
feat(logit): add LogitCalibResult struct and logit_calibrate declaration

Header for upcoming Deville-Särndal 1992 bounded logit calibration solver.
Field layout mirrors GregResult plus alm_* zeros (per spec rev 5 §580-617).
No .cpp / Makevars / bridge wiring yet — those land in Epic D.
```

---

## Task C2 — `src/logit_calib.cpp`

### Goal

Implement `logit_calibrate(CalibState&)` — full Newton solver mirroring greg.cpp's
prologue/epilogue but with logit-specific Newton inner loop (sigmoid weights, no
active set, no clamping).

### Exact code to write

```cpp
// src/logit_calib.cpp
#include "logit_calib.hpp"
#include "calib_linalg.hpp"
#include "calib_dispatch.hpp"
#include "calib_validate.hpp"     // kNCatsTotalMax
#include "leafblower.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

namespace lbw {

LogitCalibResult logit_calibrate(CalibState& st) {
    LogitCalibResult res;
    res.status = RK_ERR_NOCONV;

    // --- Build cell table (same prologue as greg.cpp) ---
    CellTable ct;
    if (build_cell_table(st.n, st.K, st.group_ids, st.cat_counts, st.weights, ct) != 0) {
        res.status = RK_ERR_BADARG;
        std::strncpy(res.message, "build_cell_table failed", sizeof(res.message) - 1);
        return res;
    }
    res.M_cell = ct.M_cell;

    // --- Initial cell masses + bounds ---
    std::vector<double> X_init(ct.M_cell, 0.0);
    for (int i = 0; i < st.n; i++) X_init[ct.cell_of[i]] += st.weights[i];

    const double lo = st.min_weight;
    const double hi = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
    std::vector<double> L_cell(ct.M_cell), U_cell(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) {
        L_cell[c] = lo * static_cast<double>(ct.n_per_cell[c]);
        U_cell[c] = hi * static_cast<double>(ct.n_per_cell[c]);
    }

    // --- Per-bucket cell lists (built once, used every iteration) ---
    std::vector<std::vector<std::vector<int>>> cells_per_cat(st.K);
    for (int k = 0; k < st.K; k++) {
        cells_per_cat[k].assign(st.cat_counts[k], {});
        for (int c = 0; c < ct.M_cell; c++) {
            int g = ct.g_per_cell[k][c];
            if (g >= 0 && g < st.cat_counts[k]) cells_per_cat[k][g].push_back(c);
        }
    }

    // --- Constraint indexing ---
    int nct = 0;
    std::vector<int> cat_offset(st.K);
    for (int k = 0; k < st.K; k++) { cat_offset[k] = nct; nct += st.cat_counts[k]; }

    if (nct > kNCatsTotalMax) {
        res.status = RK_ERR_BADARG;
        std::snprintf(res.message, sizeof(res.message),
                      "n_cats_total=%d exceeds limit %d; use method='ieppa'",
                      nct, kNCatsTotalMax);
        return res;
    }

    // --- State vectors ---
    std::vector<double> lambda(nct, 0.0);                    // dual variables
    std::vector<double> w(ct.M_cell, 0.0);                   // cell calibrated masses
    std::vector<double> D_eff(ct.M_cell, 0.0);               // Newton weights
    std::vector<double> N(static_cast<size_t>(nct) * static_cast<size_t>(nct), 0.0);
    std::vector<double> b(nct, 0.0);                         // residuals → Δλ
    const int max_cats =
        *std::max_element(st.cat_counts, st.cat_counts + st.K);
    std::vector<double> bucket_scratch(max_cats, 0.0);

    // --- Best-iterate tracking + convergence state ---
    std::vector<double> w_best(ct.M_cell, 0.0);
    double best_resid = std::numeric_limits<double>::infinity();
    double initial_resid = std::numeric_limits<double>::infinity();
    double prev_metric = std::numeric_limits<double>::infinity();

    // --- Newton iteration budget: capped at 50, possibly tighter via inner_max_iter ---
    constexpr int kMaxNewtonIters = 50;
    const int newton_budget =
        (st.inner_max_iter > 0 && st.inner_max_iter < kMaxNewtonIters)
        ? st.inner_max_iter : kMaxNewtonIters;

    const auto& cfg = st.convergence_cfg;

    bool converged = false;
    for (int iter = 0; iter < newton_budget; iter++) {
        res.iterations = iter + 1;

        // (1) Compute w[c] and D_eff[c] from current λ (logit link).
        for (int c = 0; c < ct.M_cell; c++) {
            double z = 0.0;
            for (int k = 0; k < st.K; k++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k]) z += lambda[cat_offset[k] + g];
            }
            z = std::clamp(z, -700.0, 700.0);             // exp overflow guard
            double sig = 1.0 / (1.0 + std::exp(-z));
            double range = U_cell[c] - L_cell[c];
            w[c]     = L_cell[c] + range * sig;            // bounds by construction
            D_eff[c] = range * sig * (1.0 - sig);          // Deville-Särndal 1992 eq 8
        }

        // (2) Residuals b[k][j] = τ[k][j]·n − Σ_{c∈bucket(k,j)} w[c]
        const double n_total = static_cast<double>(st.n);
        std::fill(b.begin(), b.end(), 0.0);
        for (int k = 0; k < st.K; k++) {
            for (int j = 0; j < st.cat_counts[k]; j++) {
                double target = st.targets[k][j] * n_total;
                double S_kj = 0.0;
                for (int c : cells_per_cat[k][j]) S_kj += w[c];
                b[cat_offset[k] + j] = target - S_kj;
            }
        }

        // (3) Best-iterate + convergence via shared infrastructure.
        double W_total = 0.0;
        for (int c = 0; c < ct.M_cell; c++) W_total += w[c];
        lbw::CellMetrics m =
            lbw::compute_cell_metrics(st, ct, w, W_total, bucket_scratch);

        if (iter == 0) initial_resid = m.errRp;
        if (m.errRp < best_resid) {
            best_resid = m.errRp;
            w_best     = w;
            res.best_error = m.errRp;
            res.best_iter  = iter + 1;
        }

        if (lbw::check_convergence(cfg, m, prev_metric, st.tol_abs)) {
            res.status     = RK_OK;
            res.max_error  = m.errRp;
            res.mean_error = m.mean_err;
            res.kl         = m.kl;
            res.chi2       = m.chi2;
            res.grake_norm = m.grake_norm;
            res.convergence_iter = iter + 1;
            converged = true;
            break;
        }

        // (4) Build N = A·diag(D_eff)·Aᵀ and factor.
        std::fill(N.begin(), N.end(), 0.0);
        if (compute_normal_equations(ct, D_eff.data(), N.data(),
                                     cat_offset.data(), st.K,
                                     static_cast<size_t>(nct)) != RK_OK) {
            res.status = RK_ERR_INFEAS;
            std::snprintf(res.message, sizeof(res.message),
                "logit: singular normal equations (degenerate bounds — L=U cells)");
            break;
        }
        if (ldlt_factor_inplace(N.data(), static_cast<size_t>(nct), 1e-10) != RK_OK) {
            res.status = RK_ERR_INFEAS;
            std::snprintf(res.message, sizeof(res.message),
                "logit: LDLT factorization failed (degenerate bounds)");
            break;
        }

        // (5) Solve NΔλ = b in-place; b ← Δλ.
        ldlt_solve(N.data(), static_cast<size_t>(nct), b.data());

        // (6) Update λ (full Newton step — no line search needed; logit Hessian SPD).
        for (int j = 0; j < nct; j++) lambda[j] += b[j];
    }

    // --- Post-loop status classification ---
    if (!converged && res.status == RK_ERR_NOCONV) {
        // Spec line 549: BUDGET if residual improved >0.1% vs initial; else STALL.
        res.status = (best_resid < initial_resid * 0.999)
                     ? RK_ERR_BUDGET : RK_ERR_STALL;
        std::snprintf(res.message, sizeof(res.message),
            "logit: %s after %d Newton steps; best max_err=%.4e",
            res.status == RK_ERR_BUDGET ? "budget exhausted" : "stall",
            res.iterations, best_resid);
    }

    // --- Final metrics from best iterate ---
    {
        double W_best_total = 0.0;
        for (int c = 0; c < ct.M_cell; c++) W_best_total += w_best[c];
        lbw::CellMetrics m =
            lbw::compute_cell_metrics(st, ct, w_best, W_best_total, bucket_scratch);
        res.max_error  = m.errRp;
        res.mean_error = m.mean_err;
        res.kl         = m.kl;
        res.chi2       = m.chi2;
        res.grake_norm = m.grake_norm;
        res.convergence_solver_objective = m.errRp;
    }

    // --- Weight reconstruction (cell → obs) ---
    // w_obs[i] = st.weights[i] · w_best[cell_of[i]] / X_init[cell_of[i]]
    // Bounds preserved automatically because w_best ∈ [L_cell, U_cell] by construction.
    res.best_weights.resize(st.n);
    for (int i = 0; i < st.n; i++) {
        int c = ct.cell_of[i];
        double mult = (X_init[c] > 1e-10) ? w_best[c] / X_init[c] : 1.0;
        st.weights[i] = st.weights[i] * mult;
        res.best_weights[i] = st.weights[i];
    }

    return res;
}

} // namespace lbw
```

### Algorithm trace (one Newton step on T5 fixture)

T5 input: n=1000, K=2, cat_counts=[2,2], targets sex=(0.5,0.5), age=(0.6,0.4),
default bounds (`min_weight=0`, `max_weight=Inf` → `hi=1e300`).

- Initialization: λ=[0,0,0,0]; for every cell c, z=0 → σ=0.5 → w=L+0.5·(U−L). With
  `min_weight=0`, `max_weight=Inf` → `range=1e300·n_per_cell[c]` → `w[c] ≈ 0.5e300·n_per_cell[c]`.
- Note: This blows up. **This is the case that requires Epic D's input validation** —
  `max_weight=Inf` MUST be rejected at the bridge for `method="logit"` because logit needs
  finite bounds. T5 must be revised in Epic D to set `max_weight=2.0` (or similar finite
  value) — or the bridge must default `max_weight` to e.g. `5.0` for logit. Plan flag, not
  this epic's problem.

For T6/T7/T8 (which have finite `max_weight`), the trace is well-defined:
- Iter 0: w[c] = midpoint, residuals b reflect distance from target, `D_eff[c] = (U−L)·0.25`.
- Solve N·Δλ = b → Δλ in `b`. λ += Δλ.
- Iter 1: z[c] now nonzero, σ moves toward target side, residuals shrink.
- Convergence typically in 5–20 iters per spec line 446.

### Verification gate (after C2)

```bash
R CMD INSTALL --preclean .
```

The TU compiles standalone (it is NOT yet in `Makevars.in` PKG_SOURCES — that is Epic D).
We rely on direct compilation by manually invoking the compiler, OR we accept that until
Epic D adds the file to Makevars, this TU is not yet linked into `leafblower.so`. To
exercise the compile gate without Epic D, run:

```bash
cd /home/dd/Gemini/leafblower
PKG_CXXFLAGS="$(R CMD config --cppflags) -I./src -std=c++17 -Wall -Wextra -O2" \
R CMD COMPILE src/logit_calib.cpp
```

Must produce `src/logit_calib.o` with zero warnings. If it fails, fix root cause — DO NOT
pivot to suppression flags.

### RED test verification (after C2)

Tests T5–T8 (per spec lines 661–732) live in
`tests/testthat/test-calibration-solvers.R`. After C2 alone:

- `R/harvest.R::map_method` does NOT yet recognize `"logit"`. Therefore the tests are
  expected to fail at the R-level dispatch with an "unknown method" error. This is RED
  but for the **wrong reason** (bridge missing, not solver missing).
- This is acceptable because Epic C is intentionally scoped to solver-only. Epic D adds
  `map_method` + `alg_names` + r_bridge.cpp dispatch + Makevars, after which T5–T8 should
  exercise the new code path.
- Verification command: `Rscript -e 'devtools::test(filter="calibration-solvers")'` — must
  show T5–T8 FAIL with messages indicating dispatch (not link) error.

### Commit (C2)

```
feat(logit): implement logit_calibrate Newton solver

Deville-Särndal 1992 bounded logit calibration via dual Newton-Raphson.
- w_c = L + (U-L)·σ(z_c): bounds by construction, no clamp, no active set.
- D_eff[c] = (U-L)·σ(1-σ): Newton weight from D-S 1992 eq 8.
- Reuses compute_normal_equations + ldlt_factor_inplace + ldlt_solve.
- Reuses compute_cell_metrics + check_convergence (no hand-rolled tol logic).
- Status: RK_OK | RK_ERR_INFEAS | RK_ERR_BUDGET | RK_ERR_STALL.
- Best-iterate tracking guards against worse-than-best on budget exhaust.

Bridge wiring (Makevars, r_bridge, harvest.R, leafblower.h enum) lands in
Epic D. T5-T8 expected RED with dispatch error until then.
```

---

## Compliance Bead Lock

```
Task C1 [logit_calib.hpp = verbatim spec struct + decl] !
        [no extra fields, no clamping, no premature bridge wiring]
Task C2 [logit_calib.cpp = Newton solver reusing calib_linalg + calib_dispatch] !
        [no active-set, no manual tol, no RK_ERR_NOCONV exit, no Makevars edit]
```

If `R CMD COMPILE` fails on C2 → **SPEC_FAILURE**. Halt, log finding, do not pivot.

## Confidence: 92

(Exact code traced from spec rev 5 + greg.cpp; struct field list verbatim from spec
§580-617; algorithm steps quoted from spec §460-559. Deduction at 92 because Epic D's
input validation contract for `max_weight=Inf` under logit is unspecified — flagged inline
above for Epic D to resolve, not blocking for this epic.)
