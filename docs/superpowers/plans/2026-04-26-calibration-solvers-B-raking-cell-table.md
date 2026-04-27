# Calibration Solvers — Plan B: Raking Cell-Table Migration

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Migrate `method="raking"` from obs-level O(K×n) IPF to cell-table O(K×M_cell) IPF. Same algorithm, same convergence framework, ~263× speedup on stepstone. Hard bounds enforced per-obs via post-exit clamp. Raking default metric changed to `kl+improvement` (raking IS a KL minimizer for feasible problems).

**Architecture:** `raking_solve(CalibState& st)` gains a `CellTable ct` built at entry. The inner loop works on `X[c]` (cell masses, M_cell≈6k) instead of `w[i]` (obs weights, n≈1.58M). All metrics computed cell-level. Post-exit obs expansion: `w_i = d_i × X[c]/X_init[c]` then hard clamp. harvest.R adds ieppa-style kl default for raking. Dead code (`sum_weights_ilp`, `compute_errRp`) ticketed for cleanup but not removed.

**Spec:** `docs/superpowers/specs/2026-04-25-calibration-solvers-design.md` §8

**Baseline:** FAIL 0 | PASS 346 | SKIP 4

---

## File Structure

| File | Action |
|---|---|
| `src/raking.cpp` | Add `compute_errRp_ct`; replace obs-level loops with cell-level; add post-exit obs expansion |
| `R/harvest.R` | Add raking kl+improvement default (same pattern as ieppa) |
| `data-raw/gen_raking_obs_ref.R` | Script to generate obs-level raking reference fixture |
| `tests/testthat/fixtures/raking_obs_reference_stepstone.rds` | Obs-level max_err reference (generated before migration) |
| `tests/testthat/test-calibration-solvers.R` | Add A5 test (RED → GREEN) |

---

## Key data structures

**Current obs-level (replaced):**
- `w[i]` — obs weights, size n
- `q[i]` — Dykstra box correction per obs
- `q_hyp` — scalar Dykstra hyperplane correction

**New cell-level:**
- `X[c]` — cell masses, size M_cell (`X[c] = Σ_{i∈c} st.weights[i]` initially)
- `X_init[c]` — initial cell masses (for obs expansion at exit)
- `p[c]` — Dykstra box correction per cell (replaces `q[i]`)
- `q_hyp` — scalar, same concept but shift denominator is `ct.M_cell` not `st.n`
- `X_prev[c]` — previous cell masses for pct/l1 metrics
- `W_best[c]` — best-iterate cell snapshot

**Hyperplane shift formula for cells:**

Obs-level: `shift = (n - s) / n` — n obs × shift = (n-s) ✓  
Cell-level: `shift = (n - s) / M_cell` — M_cell cells × shift = (n-s) ✓  

Use `ct.M_cell` (not `st.n`) as denominator throughout.

**CellTable fields used:**
- `ct.M_cell`, `ct.cell_of[i]`, `ct.n_per_cell[c]`, `ct.g_per_cell[k][c]`

---

## Task 1 — Generate obs-level reference fixture

**Files:** `data-raw/gen_raking_obs_ref.R`, `tests/testthat/fixtures/raking_obs_reference_stepstone.rds`

This must run BEFORE the migration so we have a ground truth for the correctness check.

- [ ] **Step 1.1: Create `data-raw/gen_raking_obs_ref.R`**

```r
#!/usr/bin/env Rscript
# Generates tests/testthat/fixtures/raking_obs_reference_stepstone.rds
# Run BEFORE Plan B migration to capture obs-level raking reference.
suppressPackageStartupMessages({
  library(leafblower); library(arrow); library(jsonlite)
})
data <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
tgt_raw <- jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json")
target <- lapply(tgt_raw, function(x) { v <- unlist(x); v / sum(v) })

cat("Running obs-level raking to record reference...\n")
t0 <- proc.time()["elapsed"]
suppressWarnings(
  w <- leafblower::harvest(data, target, method="raking",
                           max_weight=5, max_iterations=500, attach_weights=FALSE)
)
elapsed_obs <- proc.time()["elapsed"] - t0

r <- attr(w, "result")
ref <- list(
  max_error     = r$max_error,
  best_error    = r$best_error,
  iterations    = r$iterations,
  status        = r$status,
  elapsed_obs   = elapsed_obs,
  version       = as.character(packageVersion("leafblower")),
  date          = Sys.Date()
)
cat(sprintf("  max_error=%.4e  elapsed=%.1fs  iterations=%d\n",
            ref$max_error, ref$elapsed_obs, ref$iterations))
saveRDS(ref, "tests/testthat/fixtures/raking_obs_reference_stepstone.rds")
cat("Saved.\n")
```

- [ ] **Step 1.2: Run to generate fixture** (must run from main repo, not worktree)
```bash
cd /home/dd/Gemini/leafblower
Rscript data-raw/gen_raking_obs_ref.R 2>&1
ls -la tests/testthat/fixtures/raking_obs_reference_stepstone.rds
```

- [ ] **Step 1.3: Write A5 test (RED state)**

Append to `tests/testthat/test-calibration-solvers.R`:
```r
test_that("A5: raking cell-table: correct + speedup vs obs-level reference", {
  skip_on_cran()
  skip_if_not_installed("arrow")
  skip_if_not_installed("jsonlite")
  skip_if(!file.exists("benchmarks/stepstone_fulldata_bench_data.parquet"),
          "stepstone benchmark data not available")
  ref_path <- test_path("fixtures/raking_obs_reference_stepstone.rds")
  skip_if(!file.exists(ref_path), "obs-level raking reference fixture not generated")
  ref <- readRDS(ref_path)

  library(arrow); library(jsonlite)
  data <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
  tgt_raw <- jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json")
  target <- lapply(tgt_raw, function(x) { v <- unlist(x); v / sum(v) })

  t0 <- proc.time()["elapsed"]
  suppressWarnings(
    w_cell <- leafblower::harvest(data, target, method="raking",
                                  max_weight=5, max_iterations=500, attach_weights=FALSE)
  )
  elapsed_cell <- proc.time()["elapsed"] - t0

  r_cell <- attr(w_cell, "result")

  # Convergence
  expect_equal(r_cell$status, 0L, label="raking cell-table must converge")

  # Correctness: max_err within 1e-8 of obs-level reference (homogeneous d_i)
  expect_lt(abs(r_cell$max_error - ref$max_error), 1e-8,
            label="cell-table max_err must match obs-level to 1e-8")

  # Speedup: < 5% of obs-level elapsed (spec A5: elapsed_cell < elapsed_obs * 0.05)
  expect_lt(elapsed_cell, ref$elapsed_obs * 0.05,
            label=sprintf("cell-table must be < 5%% of obs-level elapsed (obs=%.1fs)", ref$elapsed_obs))

  # Hard bounds: all weights in [min_weight, max_weight]
  expect_true(all(w_cell >= 1/5 - 1e-10 & w_cell <= 5 + 1e-10),
              label="all weights within [1/5, 5] bounds")
})
```

Before writing the implementation, confirm the test WOULD fail with current code (if parquet were present, the speedup assertion would eventually pass but correctness might differ). In the worktree without parquet, it will SKIP — which is acceptable since the correctness of the RED state is documented.

- [ ] **Step 1.4: Commit fixture + test**
```bash
git add data-raw/gen_raking_obs_ref.R \
        tests/testthat/fixtures/raking_obs_reference_stepstone.rds \
        tests/testthat/test-calibration-solvers.R
git commit -m "test(A5): raking cell-table correctness + speedup red test

Reference fixture captures obs-level max_err on stepstone before migration.
A5 test verifies: cell-table max_err within 1e-8 of obs-level, elapsed < 5%
of obs-level, all weights within bounds."
```

---

## Task 2 — File: `src/raking.cpp`; `R/harvest.R`

**Ticket:** Create before starting:
```bash
bd create --title="feat(raking): migrate to cell-table backend" --type=feature --priority=1
bd create --title="chore(raking): remove dead sum_weights_ilp + compute_errRp after cell-table migration" --type=task --priority=4
```

### Step 2.1: Add cell_table.hpp include + compute_errRp_ct

Add to raking.cpp includes (after existing `#include` block):
```cpp
#include "cell_table.hpp"
```

After the existing `compute_errRp` static function (~line 80), add:
```cpp
// Cell-table errRp: O(K * M_cell). bucket pre-allocated to max_cats.
static double compute_errRp_ct(const CalibState& st,
                                const CellTable& ct,
                                const std::vector<double>& X,
                                std::vector<double>& bucket) {
    double W = 0.0;
    for (int c = 0; c < ct.M_cell; c++) W += X[c];
    double err = 0.0;
    for (int k = 0; k < st.K; k++) {
        std::fill(bucket.begin(), bucket.begin() + st.cat_counts[k], 0.0);
        for (int c = 0; c < ct.M_cell; c++) {
            int g = ct.g_per_cell[k][c];
            if (g >= 0 && g < st.cat_counts[k]) bucket[g] += X[c];
        }
        for (int j = 0; j < st.cat_counts[k]; j++) {
            double e = std::fabs(bucket[j] / W - st.targets[k][j]);
            if (e > err) err = e;
        }
    }
    return err;
}
```

### Step 2.2: Replace raking_solve entry and initialization

Find `RakingResult raking_solve(CalibState& st) {` (~line 83). Replace the entire function body (everything from the opening brace to the closing `}` before `} // namespace lbw`) with the new implementation below.

**IMPORTANT:** Read the full existing function first: `sed -n '83,477p' src/raking.cpp` to understand scope before replacing.

New function body (entry through loop setup):

```cpp
RakingResult raking_solve(CalibState& st) {
    static constexpr double kEmptyBucketThreshold = 1e-15;
    static constexpr int    kErrCheckInterval     = 10;
    static constexpr int    kMaxNoImprove         = 5;

    RakingResult res;
    res.status     = RK_ERR_NOCONV;
    res.iterations = 0;
    res.max_error  = 1.0;

    // Build cell table: O(n log n) one-time cost.
    CellTable ct;
    if (build_cell_table(st.n, st.K, st.group_ids, st.cat_counts,
                         st.weights, ct) != 0) {
        // RakingResult has no message field; caller gets RK_ERR_BADARG status.
        res.status = RK_ERR_BADARG;
        return res;
    }

    // Initial cell masses: X[c] = Σ_{i∈c} st.weights[i]
    std::vector<double> X(ct.M_cell, 0.0);
    for (int i = 0; i < st.n; i++)
        X[ct.cell_of[i]] += st.weights[i];
    std::vector<double> X_init(X);

    // Cell bounds: L_c = lo * n_per_cell[c], U_c = hi * n_per_cell[c]
    const double lo = st.min_weight;
    const double hi = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
    std::vector<double> L_cell(ct.M_cell), U_cell(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) {
        L_cell[c] = lo * static_cast<double>(ct.n_per_cell[c]);
        U_cell[c] = hi * static_cast<double>(ct.n_per_cell[c]);
    }

    // Dykstra corrections: p[c] box, q_hyp scalar hyperplane.
    std::vector<double> p(ct.M_cell, 0.0);
    double q_hyp = 0.0;

    bool is_infeasible = false;
    int max_cats = *std::max_element(st.cat_counts, st.cat_counts + st.K);
    std::vector<double> bucket(max_cats), scale(max_cats);

    // Descent monitor
    double min_errRp_window = std::numeric_limits<double>::infinity();
    int n_no_improve = 0;

    // pct/l1 tracking at cell level
    std::vector<double> X_prev(X);

    // Convergence rule state
    double prev_metric_for_rule = std::numeric_limits<double>::infinity();

    // Best-iterate tracking (cell-level snapshot)
    double best_metric_seen = std::numeric_limits<double>::infinity();
    int    best_iter_val    = 0;
    std::vector<double> W_best(ct.M_cell, 0.0);
```

### Step 2.3: Main loop

```cpp
    for (int iter = 1; iter <= st.inner_max_iter; iter++) {
        res.iterations = iter;

        // Cell-level cyclic IPF
        double W_total = 0.0;
        for (int c = 0; c < ct.M_cell; c++) W_total += X[c];

        for (int k = 0; k < st.K; k++) {
            std::fill(bucket.begin(), bucket.begin() + st.cat_counts[k], 0.0);
            for (int c = 0; c < ct.M_cell; c++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k]) bucket[g] += X[c];
            }
            std::fill(scale.begin(), scale.begin() + st.cat_counts[k], 1.0);
            for (int j = 0; j < st.cat_counts[k]; j++) {
                double Tkj = st.targets[k][j] * W_total;
                if (bucket[j] < kEmptyBucketThreshold * W_total) {
                    if (Tkj > 0.0) is_infeasible = true;
                } else {
                    scale[j] = Tkj / bucket[j];
                }
            }
            for (int c = 0; c < ct.M_cell; c++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k]) X[c] *= scale[g];
            }
        }

        // Dykstra box: X[c] = clamp(X[c] + p[c], L_cell[c], U_cell[c])
        for (int c = 0; c < ct.M_cell; c++) {
            double yc = X[c] + p[c];
            double Xc = std::clamp(yc, L_cell[c], U_cell[c]);
            p[c] = yc - Xc;
            X[c] = Xc;
        }

        // Dykstra hyperplane: sum(X) = n.
        // Shift is (n - s) / M_cell — M_cell cells × shift = (n-s), corrects total to n.
        {
            double s = 0.0;
            for (int c = 0; c < ct.M_cell; c++) { X[c] += q_hyp; s += X[c]; }
            double shift = (static_cast<double>(st.n) - s) / static_cast<double>(ct.M_cell);
            for (int c = 0; c < ct.M_cell; c++) X[c] += shift;
            q_hyp = -shift;
        }

        // Convergence check
        if (iter == 1 || iter % kErrCheckInterval == 0 || iter == st.inner_max_iter) {
            double errRp = compute_errRp_ct(st, ct, X, bucket);
            res.max_error = errRp;

            // BLOCK 1: MAX_ERR best-iterate
            if (st.convergence_cfg.metric == lbw::CalibMetric::MAX_ERR) {
                if (errRp < best_metric_seen) {
                    best_metric_seen = errRp;
                    best_iter_val    = iter;
                    W_best           = X;
                }
            }

            // pct_change + l1_weight (cell-level approximation: Σ|ΔX[c]|/n)
            double pct_change = 0.0, l1_weight_sum = 0.0;
            for (int c = 0; c < ct.M_cell; c++) {
                double diff = std::fabs(X[c] - X_prev[c]);
                double rel  = diff / std::max(X_prev[c], 1e-12);
                if (rel > pct_change) pct_change = rel;
                l1_weight_sum += diff;
            }
            double l1_weight = l1_weight_sum / static_cast<double>(st.n);

            // Extra metrics (gated, same as obs-level)
            const lbw::CalibMetric metric = st.convergence_cfg.metric;
            const auto& cfg_m = st.convergence_cfg;
            const bool about_to_converge =
                (cfg_m.absolute_tol > 0.0 && errRp < cfg_m.absolute_tol) ||
                [&]() {
                    double prev_copy = prev_metric_for_rule;
                    double active = (metric == lbw::CalibMetric::MAX_ERR)   ? errRp :
                                    (metric == lbw::CalibMetric::L1_WEIGHT) ? l1_weight : -1.0;
                    return active >= 0.0 && lbw::apply_rule(cfg_m.rule, active, prev_copy, cfg_m.pct_tol);
                }();
            const bool need_extra_metrics =
                (metric == lbw::CalibMetric::MEAN_ERR   ||
                 metric == lbw::CalibMetric::KL         ||
                 metric == lbw::CalibMetric::CHI2       ||
                 metric == lbw::CalibMetric::GRAKE_NORM ||
                 iter == st.inner_max_iter               ||
                 about_to_converge);

            double W_tot2 = 0.0;
            for (int c = 0; c < ct.M_cell; c++) W_tot2 += X[c];
            constexpr double kMetricEps = 1e-10;
            constexpr double kChi2Floor = 1.0;
            double mean_err_sum = 0.0, kl_max = 0.0, chi2_total = 0.0, grake_norm = 0.0;
            if (need_extra_metrics && W_tot2 > 0.0) {
                for (int k = 0; k < st.K; k++) {
                    const int nj = st.cat_counts[k];
                    std::fill(bucket.begin(), bucket.begin() + nj, 0.0);
                    for (int c = 0; c < ct.M_cell; c++) {
                        int g = ct.g_per_cell[k][c];
                        if (g >= 0 && g < nj) bucket[g] += X[c];
                    }
                    double max_k = 0.0, kl_k = 0.0;
                    for (int j = 0; j < nj; j++) {
                        double S_p   = bucket[j] / W_tot2;
                        double T     = st.targets[k][j];
                        double err   = std::fabs(S_p - T);
                        if (err > max_k) max_k = err;
                        if (T > 0.0)
                            kl_k += T * std::log((T + kMetricEps) / (S_p + kMetricEps));
                        double obs    = bucket[j];
                        double pop_kj = T * W_tot2;
                        chi2_total += (obs - pop_kj) * (obs - pop_kj) / (pop_kj + kChi2Floor);
                        double nm = std::fabs(obs - pop_kj) / (1.0 + std::fabs(pop_kj));
                        if (nm > grake_norm) grake_norm = nm;
                    }
                    mean_err_sum += max_k;
                    if (kl_k > kl_max) kl_max = kl_k;
                }
                // BLOCK 2: best-iterate for non-MAX_ERR
                if (st.convergence_cfg.metric != lbw::CalibMetric::MAX_ERR) {
                    const double mean_err_blk2 = (st.K > 0)
                        ? mean_err_sum / static_cast<double>(st.K) : 0.0;
                    const double curr_best = lbw::select_metric(
                        st.convergence_cfg.metric,
                        errRp, mean_err_blk2, kl_max, chi2_total, grake_norm, l1_weight);
                    if (std::isfinite(curr_best) && curr_best < best_metric_seen) {
                        best_metric_seen = curr_best;
                        best_iter_val    = iter;
                        W_best           = X;
                    }
                }
            }
            double mean_err = (st.K > 0) ? mean_err_sum / static_cast<double>(st.K) : 0.0;

            res.l1_weight_change = l1_weight;
            res.mean_error       = mean_err;
            res.kl               = kl_max;
            res.chi2             = chi2_total;
            res.grake_norm       = grake_norm;
            for (int c = 0; c < ct.M_cell; c++) X_prev[c] = X[c];

            // Descent monitor
            if (!std::isfinite(min_errRp_window)) {
                min_errRp_window = errRp; n_no_improve = 0;
            } else {
                const double rel_eps = 0.01 * min_errRp_window;
                const double eps = std::max(rel_eps, st.tol_abs);
                if (errRp < min_errRp_window - eps) {
                    min_errRp_window = errRp; n_no_improve = 0;
                } else {
                    n_no_improve++;
                }
            }

            if (st.verbose >= 1) {
                char msg[256];
                std::snprintf(msg, 256, "raking iter %d: errRp=%.2e", iter, errRp);
                st.log(msg);
            }

            // Convergence dispatch (identical to obs-level)
            {
                const auto& cfg = st.convergence_cfg;
                const double curr_metric = lbw::select_metric(
                    cfg.metric, errRp, mean_err, kl_max, chi2_total, grake_norm, l1_weight);
                bool converged_abs = (cfg.absolute_tol > 0.0) && (curr_metric < cfg.absolute_tol);
                const bool converged_pct = lbw::apply_rule(
                    cfg.rule, curr_metric, prev_metric_for_rule, cfg.pct_tol);
                bool have_pct = (cfg.pct_tol > 0.0), have_abs = (cfg.absolute_tol > 0.0);
                bool converged = false;
                if (have_pct && have_abs) {
                    converged = (cfg.stop_when == lbw::CalibStopWhen::ALL)
                                ? (converged_pct && converged_abs)
                                : (converged_pct || converged_abs);
                } else if (have_pct)  converged = converged_pct;
                else if (have_abs)    converged = converged_abs;
                else                  converged = (errRp < st.tol_abs);

                if (converged) {
                    res.status             = is_infeasible ? RK_ERR_INFEAS : RK_OK;
                    res.convergence_metric = static_cast<int>(cfg.metric);
                    res.convergence_rule   = static_cast<int>(cfg.rule);
                    res.convergence_tol    = cfg.pct_tol;
                    res.convergence_iter   = iter;
                    break;
                }
            }

            if (n_no_improve >= kMaxNoImprove) {
                res.status = RK_ERR_NOCONV;
                if (st.verbose >= 1) {
                    char msg[256];
                    std::snprintf(msg, 256,
                        "raking: errRp stalled for %d checks (last=%.2e, window_min=%.2e); "
                        "aborting at iter %d.", n_no_improve, errRp, min_errRp_window, iter);
                    st.log(msg);
                }
                break;
            }
        }
    }

    if (is_infeasible && res.status == RK_ERR_NOCONV)
        res.status = RK_ERR_INFEAS;
```

### Step 2.4: Exit block (post-loop Dykstra + obs expansion)

```cpp
    // Post-loop Dykstra finalizer at cell level.
    for (int fixup = 0; fixup < 20; fixup++) {
        bool box_ok = true;
        for (int c = 0; c < ct.M_cell; c++) {
            double yc = X[c] + p[c];
            double Xc = std::clamp(yc, L_cell[c], U_cell[c]);
            p[c] = yc - Xc;
            if (yc != Xc) box_ok = false;
            X[c] = Xc;
        }
        {
            double s = 0.0;
            for (int c = 0; c < ct.M_cell; c++) { X[c] += q_hyp; s += X[c]; }
            double shift = (static_cast<double>(st.n) - s) / static_cast<double>(ct.M_cell);
            for (int c = 0; c < ct.M_cell; c++) X[c] += shift;
            q_hyp = -shift;
        }
        if (box_ok) break;
    }

    // Best-iterate: normalize cell snapshot, then expand to obs.
    res.convergence_objective        = best_metric_seen;
    res.convergence_minimized_metric = static_cast<int>(st.convergence_cfg.metric);
    res.best_error = best_metric_seen;
    res.best_iter  = best_iter_val;
    if (std::isfinite(best_metric_seen)) {
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

    // Post-exit obs expansion: w_i = d_i × X[c]/X_init[c], hard clamp.
    const double hi_obs = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
    for (int i = 0; i < st.n; i++) {
        int c = ct.cell_of[i];
        double mult = (X_init[c] > 0.0) ? X[c] / X_init[c] : 1.0;
        st.weights[i] = std::clamp(st.weights[i] * mult, lo, hi_obs);
    }

    // Solver-contract normalization: sum(w) = n.
    double total_w = 0.0;
    for (int i = 0; i < st.n; i++) total_w += st.weights[i];
    if (std::isfinite(total_w) && total_w > 0.0) {
        const double norm = static_cast<double>(st.n) / total_w;
        for (int i = 0; i < st.n; i++) st.weights[i] *= norm;
    }

    return res;
}

} // namespace lbw
```

### Step 2.5: Note dead code

`sum_weights_ilp` and `compute_errRp` (obs-level) are now dead code. Do NOT remove them — the cleanup ticket is already created (Step 2 ticket creation). They will be removed in a follow-up.

### Step 2.6: Change raking default metric in R/harvest.R

Read: `grep -n "method.*ieppa\|conv\$metric.*kl\|ieppa.*kl" R/harvest.R | head -10`

Find the ieppa kl override block (added in Plan A Task 3):
```r
if (method == "ieppa" && is.null(convergence[["metric"]]) && ...)
  conv$metric <- "kl"
```

Add an analogous block for raking immediately after:
```r
# Raking is a KL minimizer — same default as ieppa.
if (method == "raking" &&
    is.null(convergence[["metric"]]) &&
    is.null(convergence[["criterion"]]) &&
    is.null(convergence[["improvement"]]) &&
    is.null(convergence[["pct"]]) &&
    is.null(convergence[["absolute"]])) {
  conv$metric <- "kl"
}
```

### Step 2.7: Build gate
```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```
Expected: `* DONE (leafblower)`. Fix any compile error before proceeding.

### Step 2.8: Full regression
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -5
```
Expected: FAIL 0, PASS ≥ 346.

If raking tests FAIL: read the failure. Do NOT adjust tolerances blindly.

### Step 2.9: Commit
```bash
git add src/raking.cpp R/harvest.R
git commit -m "$(cat <<'EOF'
feat(raking): migrate to cell-table backend + kl default

Cell-table O(K×M_cell) IPF replaces obs-level O(K×n). CellTable built
at entry; X[c] replaces w[i]. Dykstra box on cells (L_c/U_c). Hyperplane
shift = (n - s) / M_cell (not n). Post-exit obs expansion: w_i = d_i ×
X[c]/X_init[c] then hard clamp. All convergence metrics preserved.
compute_errRp_ct added; obs-level compute_errRp/sum_weights_ilp ticketed
for removal. harvest.R: raking now defaults to kl+improvement (same as
ieppa — raking IS a KL minimizer for feasible problems).
EOF
)"
```

---

## Task 3 — Final verification

- [ ] **Step 3.1: Full regression**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 346.

- [ ] **Step 3.2: Raking default metric verified**
```bash
Rscript -e '
  set.seed(1)
  data <- data.frame(a=factor(sample(c("1","2"),200,TRUE)))
  target <- list(a=c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight=3, method="raking", attach_weights=FALSE)
  r <- attr(w, "result")
  cat("metric:", r$convergence_used$metric, "\n")
  stopifnot(r$convergence_used$metric == "kl")
' 2>&1
```

- [ ] **Step 3.3: Smoke test**
```bash
Rscript -e '
  set.seed(1)
  data <- data.frame(a=factor(sample(c("1","2","3"),500,TRUE)),
                     b=factor(sample(c("1","2"),500,TRUE)))
  target <- list(a=c("1"=0.4,"2"=0.4,"3"=0.2), b=c("1"=0.6,"2"=0.4))
  w <- leafblower::harvest(data, target, max_weight=5, method="raking", attach_weights=FALSE)
  r <- attr(w, "result")
  cat(sprintf("status=%d max_error=%.4e iterations=%d\n", r$status, r$max_error, r$iterations))
  stopifnot(r$status == 0, r$max_error < 0.01)
' 2>&1
```

---

## Final Verification

- [ ] `grep "compute_errRp_ct\|cell_table.hpp" src/raking.cpp | wc -l` → > 0
- [ ] `grep "convergence_objective.*best_metric\|best_metric.*convergence_objective" src/raking.cpp` → present
- [ ] `grep "conv\$metric.*kl.*raking\|raking.*kl" R/harvest.R` → raking override present
- [ ] FAIL 0, PASS ≥ 346
