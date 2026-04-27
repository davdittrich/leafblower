# Calibration Solvers — Plan C: method="sinkhorn" (KL Bregman Dykstra)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `method="sinkhorn"` — true KL minimum subject to hard capacity constraints. Removes the `RK_ERR_BADARG` stub and makes the A1 test GREEN.

**Algorithm:** Alternating KL Bregman projections (spec §4). K Sinkhorn sweeps per outer iteration + bisection-based KL projection onto capacity box with log-domain Dykstra corrections. Monotone KL decrease guaranteed each outer iteration.

**Architecture:** New `src/sinkhorn.cpp` + `src/sinkhorn.hpp`. Wire into c_api.cpp (replace SINKHORN stub) and r_bridge.cpp (remove "sinkhorn" early-exit; add dispatch branch). T1a test updated to remove "sinkhorn" from stub loop.

**Spec:** `docs/superpowers/specs/2026-04-25-calibration-solvers-design.md` §4

**Acceptance:** A1 test GREEN (sinkhorn KL ≤ iEPPA KL at best_iter on stepstone-fulldata).

**Baseline:** FAIL 0 | PASS 355 | SKIP 5

**KL metric definition (critical):** All six metrics use the SAME formulas as raking.cpp and ieppa.cpp. The `kl` metric is `max_k Σ_j T[k][j] * log(T[k][j] / S_p[k][j])` (KL of target proportions from achieved proportions, max over margins) — NOT the cell-mass KL. This matches the iEPPA reference fixture.

---

## File Structure

| File | Action |
|---|---|
| `src/sinkhorn.hpp` | New — SinkhornResult struct + sinkhorn_solve declaration |
| `src/sinkhorn.cpp` | New — full implementation |
| `src/c_api.cpp` | Replace SINKHORN stub with real dispatch |
| `src/r_bridge.cpp` | Remove "sinkhorn" early-exit; add sinkhorn dispatch branch |
| `src/Makevars` + `src/Makevars.in` | Add sinkhorn.cpp to PKG_SOURCES |
| `tests/testthat/test-calibration-solvers.R` | Update T1a: remove "sinkhorn" from stub loop |

---

## Task 1 — Update T1a test (RED step for Plan C)

**Before implementing sinkhorn, update T1a so that when sinkhorn succeeds, T1a doesn't fail.**

```bash
grep -n "T1a\|sinkhorn.*chebyshev\|for.*m in" tests/testthat/test-calibration-solvers.R | head -5
```

Find the T1a loop `for (m in c("sinkhorn", "chebyshev", "greg", "grake"))`. Change to:
```r
for (m in c("chebyshev", "greg", "grake")) {  # sinkhorn removed — implemented in Plan C
```

Run to confirm T1a still passes with remaining 3 stubs:
```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | grep -E "T1a|FAIL|PASS" | head -5
```

Commit:
```bash
git add tests/testthat/test-calibration-solvers.R
git commit -m "test(T1a): remove sinkhorn from stub loop — Plan C implements it"
```

---

## Task 2 — src/sinkhorn.hpp + src/sinkhorn.cpp

**Ticket:**
```bash
bd create --title="feat(sinkhorn): KL Bregman Dykstra solver" --type=feature --priority=1 2>&1 | tail -2
```
Note the ticket ID from output.

### Step 2.1: Create src/sinkhorn.hpp

```cpp
#pragma once
#include "types.hpp"
#include "cell_table.hpp"
#include <vector>
#include <limits>

namespace lbw {

struct SinkhornResult {
    int    status        = RK_ERR_NOCONV;
    int    iterations    = 0;
    double max_error     = 1.0;
    double mean_error    = 0.0;
    double kl            = 0.0;
    double chi2          = 0.0;
    double grake_norm    = 0.0;
    double l1_weight_change = 0.0;
    int    convergence_metric = 0;
    int    convergence_rule   = 0;
    double convergence_tol    = 0.001;
    int    convergence_iter   = -1;
    double convergence_objective        = std::numeric_limits<double>::infinity();
    int    convergence_minimized_metric = 0;
    double best_error   = std::numeric_limits<double>::infinity();
    int    best_iter    = 0;
    std::vector<double> best_weights;
    int    M_cell       = 0;
    char   message[256] = {};
};

// KL Bregman Dykstra calibration solver.
// Minimizes KL(target || achieved) subject to margin constraints and hard capacity bounds.
// Monotone KL decrease each outer iteration (spec §4).
SinkhornResult sinkhorn_solve(CalibState& st);

} // namespace lbw
```

### Step 2.2: Create src/sinkhorn.cpp

```cpp
#include "sinkhorn.hpp"
#include "calib_dispatch.hpp"
#include "leafblower.h"
#include <cmath>
#include <cstring>
#include <algorithm>
#include <limits>

namespace lbw {

// Bisect on log-scale multiplier μ to project X onto capacity box.
// Finds μ s.t. Σ_c clamp(X[c]*exp(a[c]+μ), L[c], U[c]) = target_mass.
// Returns false if infeasible.
static bool bisect_capacity(const std::vector<double>& X,
                             const std::vector<double>& a,
                             const std::vector<double>& L,
                             const std::vector<double>& U,
                             int M_cell,
                             double target_mass,
                             double& mu_out,
                             std::vector<double>& X_proj)
{
    double sum_L = 0.0, sum_U = 0.0;
    for (int c = 0; c < M_cell; c++) { sum_L += L[c]; sum_U += U[c]; }
    if (sum_L > target_mass + 1e-9 || sum_U < target_mass - 1e-9) return false;

    auto f = [&](double mu) -> double {
        double s = 0.0;
        for (int c = 0; c < M_cell; c++)
            s += std::clamp(X[c] * std::exp(a[c] + mu), L[c], U[c]);
        return s - target_mass;
    };

    double lo = -50.0, hi = 50.0;
    while (lo > -500.0 && f(lo) > 0.0) lo *= 2.0;
    while (hi < 500.0  && f(hi) < 0.0) hi *= 2.0;
    for (int i = 0; i < 80; i++) {
        double mid = 0.5 * (lo + hi);
        if (f(mid) < 0.0) lo = mid; else hi = mid;
        if (hi - lo < 1e-12) break;
    }
    mu_out = 0.5 * (lo + hi);
    for (int c = 0; c < M_cell; c++)
        X_proj[c] = std::clamp(X[c] * std::exp(a[c] + mu_out), L[c], U[c]);
    return true;
}

SinkhornResult sinkhorn_solve(CalibState& st) {
    static constexpr int    kErrCheckInterval = 10;
    static constexpr double kMetricEps        = 1e-10;
    static constexpr double kChi2Floor        = 1.0;

    SinkhornResult res;
    res.status = RK_ERR_NOCONV;

    CellTable ct;
    if (build_cell_table(st.n, st.K, st.group_ids, st.cat_counts, st.weights, ct) != 0) {
        res.status = RK_ERR_BADARG;
        std::strncpy(res.message, "build_cell_table failed", sizeof(res.message) - 1);
        return res;
    }
    res.M_cell = ct.M_cell;

    std::vector<double> X(ct.M_cell, 0.0);
    for (int i = 0; i < st.n; i++) X[ct.cell_of[i]] += st.weights[i];
    const std::vector<double> X_init(X);

    const double lo = st.min_weight;
    const double hi = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
    std::vector<double> L_cell(ct.M_cell), U_cell(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) {
        L_cell[c] = lo * static_cast<double>(ct.n_per_cell[c]);
        U_cell[c] = hi * static_cast<double>(ct.n_per_cell[c]);
    }

    std::vector<double> a(ct.M_cell, 0.0);
    int max_cats = *std::max_element(st.cat_counts, st.cat_counts + st.K);
    // Note: no lf (log Sinkhorn factors) — state is fully captured by in-place X[c] updates.

    double best_metric_seen = std::numeric_limits<double>::infinity();
    int    best_iter_val    = 0;
    std::vector<double> W_best(ct.M_cell, 0.0);

    std::vector<double> bucket(max_cats);
    std::vector<double> scale(max_cats);
    std::vector<double> X_proj(ct.M_cell);
    std::vector<double> X_prev(X);
    double prev_metric_for_rule = std::numeric_limits<double>::infinity();

    for (int iter = 1; iter <= st.inner_max_iter; iter++) {
        res.iterations = iter;

        // ① K Sinkhorn sweeps
        double W_total = 0.0;
        for (int c = 0; c < ct.M_cell; c++) W_total += X[c];

        // O(K × M_cell): one bucket-accumulation pass + one scale-apply pass per margin.
        for (int k = 0; k < st.K; k++) {
            std::fill(bucket.begin(), bucket.begin() + st.cat_counts[k], 0.0);
            for (int c = 0; c < ct.M_cell; c++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k]) bucket[g] += X[c];
            }
            std::fill(scale.begin(), scale.begin() + st.cat_counts[k], 1.0);
            for (int j = 0; j < st.cat_counts[k]; j++) {
                if (bucket[j] < 1e-300) continue;
                double ratio = st.targets[k][j] * W_total / bucket[j];
                if (ratio <= 0.0) continue;
                scale[j] = ratio;
            }
            for (int c = 0; c < ct.M_cell; c++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k]) X[c] *= scale[g];
            }
        }

        // ② KL capacity projection via bisection
        double target_mass = 0.0;
        for (int c = 0; c < ct.M_cell; c++) target_mass += X[c];

        double mu;
        if (!bisect_capacity(X, a, L_cell, U_cell, ct.M_cell, target_mass, mu, X_proj)) {
            res.status = RK_ERR_INFEAS;
            break;
        }
        for (int c = 0; c < ct.M_cell; c++) {
            if (X[c] > 1e-300 && X_proj[c] > 1e-300)
                a[c] += std::log(X[c]) - std::log(X_proj[c]);
            X[c] = X_proj[c];
        }

        // Convergence check
        if (iter == 1 || iter % kErrCheckInterval == 0 || iter == st.inner_max_iter) {
            double W = 0.0;
            for (int c = 0; c < ct.M_cell; c++) W += X[c];

            // Reuse bucket for all 6 metrics (mirrors raking.cpp exactly)
            double errRp = 0.0;
            double mean_err_sum = 0.0, kl_max = 0.0, chi2_total = 0.0, grake_norm = 0.0;
            for (int k = 0; k < st.K; k++) {
                const int nj = st.cat_counts[k];
                std::fill(bucket.begin(), bucket.begin() + nj, 0.0);
                for (int c = 0; c < ct.M_cell; c++) {
                    int g = ct.g_per_cell[k][c];
                    if (g >= 0 && g < nj) bucket[g] += X[c];
                }
                double max_k = 0.0, kl_k = 0.0;
                for (int j = 0; j < nj; j++) {
                    double S_p = bucket[j] / W;
                    double T   = st.targets[k][j];
                    double err = std::fabs(S_p - T);
                    if (err > max_k) max_k = err;
                    if (err > errRp) errRp = err;
                    if (T > 0.0)
                        kl_k += T * std::log((T + kMetricEps) / (S_p + kMetricEps));
                    double obs    = bucket[j];
                    double pop_kj = T * W;
                    chi2_total += (obs - pop_kj) * (obs - pop_kj) / (pop_kj + kChi2Floor);
                    double nm = std::fabs(obs - pop_kj) / (1.0 + std::fabs(pop_kj));
                    if (nm > grake_norm) grake_norm = nm;
                }
                mean_err_sum += max_k;
                if (kl_k > kl_max) kl_max = kl_k;
            }
            double mean_err = (st.K > 0) ? mean_err_sum / static_cast<double>(st.K) : 0.0;

            // l1_weight (cell-level approximation: Σ|ΔX[c]|/n)
            double l1_sum = 0.0;
            for (int c = 0; c < ct.M_cell; c++)
                l1_sum += std::fabs(X[c] - X_prev[c]);
            double l1_weight = l1_sum / static_cast<double>(st.n);

            res.max_error   = errRp;
            res.kl          = kl_max;
            res.mean_error  = mean_err;
            res.chi2        = chi2_total;
            res.grake_norm  = grake_norm;
            res.l1_weight_change = l1_weight;
            for (int c = 0; c < ct.M_cell; c++) X_prev[c] = X[c];

            // Best-iterate tracking (sinkhorn is monotone: best_iter == last_iter)
            const double curr_best = lbw::select_metric(
                st.convergence_cfg.metric,
                errRp, mean_err, kl_max, chi2_total, grake_norm, l1_weight);
            if (std::isfinite(curr_best) && curr_best < best_metric_seen) {
                best_metric_seen = curr_best;
                best_iter_val    = iter;
                W_best           = X;
            }

            // Convergence dispatch (identical to raking.cpp)
            const auto& cfg = st.convergence_cfg;
            const double curr_metric = lbw::select_metric(
                cfg.metric, errRp, mean_err, kl_max, chi2_total, grake_norm, l1_weight);
            bool converged_abs = (cfg.absolute_tol > 0.0) && (curr_metric < cfg.absolute_tol);
            const bool converged_pct = lbw::apply_rule(
                cfg.rule, curr_metric, prev_metric_for_rule, cfg.pct_tol);
            bool have_pct = (cfg.pct_tol > 0.0), have_abs = (cfg.absolute_tol > 0.0);
            bool converged = false;
            if (have_pct && have_abs)
                converged = (cfg.stop_when == lbw::CalibStopWhen::ALL)
                            ? (converged_pct && converged_abs)
                            : (converged_pct || converged_abs);
            else if (have_pct)  converged = converged_pct;
            else if (have_abs)  converged = converged_abs;
            else                converged = (errRp < st.tol_abs);

            if (converged) {
                res.status             = RK_OK;
                res.convergence_metric = static_cast<int>(cfg.metric);
                res.convergence_rule   = static_cast<int>(cfg.rule);
                res.convergence_tol    = cfg.pct_tol;
                res.convergence_iter   = iter;
                break;
            }
        }
    }

    // Exit: populate best-iterate
    res.convergence_objective        = best_metric_seen;
    res.convergence_minimized_metric = static_cast<int>(st.convergence_cfg.metric);
    res.best_error = best_metric_seen;
    res.best_iter  = best_iter_val;

    if (std::isfinite(best_metric_seen) && !W_best.empty()) {
        double s = 0.0;
        for (int c = 0; c < ct.M_cell; c++) s += W_best[c];
        if (s > 0.0) {
            const double sc = static_cast<double>(st.n) / s;
            for (int c = 0; c < ct.M_cell; c++) W_best[c] *= sc;
        }
        res.best_weights.resize(st.n);
        const double hi_obs = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
        for (int i = 0; i < st.n; i++) {
            int c = ct.cell_of[i];
            double mult = (X_init[c] > 0.0) ? W_best[c] / X_init[c] : 1.0;
            res.best_weights[i] = std::clamp(st.weights[i] * mult, lo, hi_obs);
        }
    } else {
        res.best_weights.assign(st.n, 0.0);
    }

    // Obs expansion: w_i = d_i × X[c]/X_init[c], hard clamp.
    // sum(X[c]) = n by construction (Sinkhorn sweeps + bisection preserve mass).
    // Normalization post-clamp would violate bounds (same reasoning as raking.cpp fix).
    const double hi_obs = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
    for (int i = 0; i < st.n; i++) {
        int c = ct.cell_of[i];
        double mult = (X_init[c] > 0.0) ? X[c] / X_init[c] : 1.0;
        st.weights[i] = std::clamp(st.weights[i] * mult, lo, hi_obs);
    }
    return res;
}

} // namespace lbw
```

### Step 2.3: Add sinkhorn.cpp to Makevars AND Makevars.in

Read: `cat src/Makevars src/Makevars.in`

In BOTH files, add `sinkhorn.cpp` to `PKG_SOURCES`. Also verify `calib_validate.cpp` is present (it should be).

### Step 2.4: Build gate
```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```
Expected: `* DONE`. Fix any compile errors before wiring.

---

## Task 3 — Wire sinkhorn into c_api.cpp + r_bridge.cpp

### Step 3.1: c_api.cpp

Add includes at the top:
```cpp
#include "sinkhorn.hpp"
```

Read the SINKHORN stub: `grep -n "case RK_ALG_SINKHORN" src/c_api.cpp`

Replace only the `RK_ALG_SINKHORN` case (keep CHEBYSHEV/GREG/GRAKE as stubs):
```cpp
case RK_ALG_SINKHORN: {
    auto sres = lbw::sinkhorn_solve(st);
    if (result) {
        result->status                       = sres.status;
        result->iterations                   = sres.iterations;
        result->max_error                    = sres.max_error;
        result->convergence_metric           = sres.convergence_metric;
        result->convergence_rule             = sres.convergence_rule;
        result->convergence_tol              = sres.convergence_tol;
        result->convergence_iter             = sres.convergence_iter;
        result->convergence_objective        = sres.convergence_objective;
        result->convergence_minimized_metric = sres.convergence_minimized_metric;
        result->best_error                   = sres.best_error;
        result->best_iter                    = sres.best_iter;
        result->algorithm_used               = static_cast<int>(alg);
        std::strncpy(result->message, sres.message, sizeof(result->message) - 1);
    }
    return sres.status;
}
case RK_ALG_CHEBYSHEV:
case RK_ALG_GREG:
case RK_ALG_GRAKE: {
    if (result) {
        result->status = RK_ERR_BADARG;
        snprintf(result->message, sizeof(result->message),
                 "method not yet implemented in this build");
        result->algorithm_used = static_cast<int>(algorithm);  // named unimplemented alg, not AUTO
    }
    return RK_ERR_BADARG;
}
```

### Step 3.2: r_bridge.cpp

Add include: `#include "sinkhorn.hpp"`

Find the stub early-exit block (`grep -n "sinkhorn\|Stub methods" src/r_bridge.cpp | head -5`). Remove "sinkhorn" from it:
```cpp
/* Stub methods (Plans D/E): error immediately. */
if (strcmp(method_str, "chebyshev") == 0 ||
    strcmp(method_str, "greg")      == 0 ||
    strcmp(method_str, "grake")     == 0) {
    Rf_error("leafblower: invalid arguments \342\200\224 method not yet implemented in this build");
}
```

Then add the `"sinkhorn"` dispatch branch. Find the existing branches pattern (`grep -n "strcmp.*method_str.*sinkhorn\|} else if.*sinkhorn\|} else if.*raking" src/r_bridge.cpp | head -5`). Add after the raking branch or auto branch:
```cpp
} else if (strcmp(method_str, "sinkhorn") == 0) {
    auto res = lbw::sinkhorn_solve(st);
    pack_solver_result(res);
    res_best_error  = res.best_error;
    res_best_iter   = res.best_iter;
    if (!res.best_weights.empty())
        res_best_weights = std::move(res.best_weights);
    else
        res_best_weights.assign(st.n, 0.0);
    res_status = res.status;
```

Close the else-if chain properly (check how it ends in the existing code).

### Step 3.3: Build gate
```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

### Step 3.4: Smoke test
```bash
Rscript -e '
  set.seed(42)
  data <- data.frame(
    a = factor(sample(c("1","2","3"), 500, replace=TRUE)),
    b = factor(sample(c("1","2"), 500, replace=TRUE))
  )
  target <- list(a=c("1"=1/3,"2"=1/3,"3"=1/3), b=c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight=5, method="sinkhorn",
                           max_iterations=300, attach_weights=FALSE)
  r <- attr(w, "result")
  cat(sprintf("status=%d max_error=%.4e kl=%.4e iter=%d\n",
              r$status, r$max_error, r$kl, r$iterations))
  stopifnot(r$status == 0, r$max_error < 0.01)
' 2>&1 | grep -v "^Welcome\|Working"
```

### Step 3.5: Full regression
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 355 (A1 still SKIPs in worktree — no parquet).

### Step 3.6: Commit
```bash
git add src/sinkhorn.hpp src/sinkhorn.cpp src/c_api.cpp src/r_bridge.cpp \
        src/Makevars src/Makevars.in
git commit -m "$(cat <<'EOF'
feat(sinkhorn): KL Bregman Dykstra solver — true KL minimum

sinkhorn_solve: alternating Sinkhorn sweeps + bisection KL projection
onto capacity box with log-domain Dykstra corrections (spec §4).
Monotone KL decrease per iteration. All 6 metrics computed same as
raking/ieppa. Wired in c_api.cpp and r_bridge.cpp; CHEBYSHEV/GREG/GRAKE
stubs retained. A1 test GREEN when stepstone data available.
EOF
)"
bd close <ticket-id>
```

---

## Task 4 — Verify A1 GREEN

Run from main repo (has stepstone parquet):
```bash
cd /home/dd/Gemini/leafblower
Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | grep -E "A1|FAIL|PASS|SKIP" | head -10
```
Expected: A1 PASSES — sinkhorn KL ≤ iEPPA KL at best_iter (3.008e-3).

---

## Final Verification

- [ ] `Rscript -e 'devtools::test()' 2>&1 | tail -3` → FAIL 0, PASS ≥ 355
- [ ] T1a no longer includes "sinkhorn" in stub loop
- [ ] `grep "case RK_ALG_SINKHORN" src/c_api.cpp` → dispatches to sinkhorn_solve
- [ ] `grep "sinkhorn" src/r_bridge.cpp | grep "Rf_error"` → 0 (stub error removed)
- [ ] A1 GREEN from main repo
