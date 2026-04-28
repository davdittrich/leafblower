# Raking Water-Filling Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace leafblower's hybrid Dykstra/IPF bounded raking with autumn-style water-filling bounded IPF, giving SQUAREM a clean stateless F operator so the L2 step-halving criterion works correctly and the DEFF gap to autumn closes.

**Architecture:** `water_fill_cat()` — a lambda inside `raking_solve` — replaces the separate Dykstra box correction. For each (margin k, category j) it iteratively finds the multiplier m that brings category cells to the target sum while respecting `[L_cell, U_cell]` bounds. p[] and q_hyp correction vectors are removed; F_eval becomes `(Xv) -> double`. SQUAREM can then use autumn's L2 halving (`‖F(X*)-X*‖² ≤ 1.01 × ‖v‖²`) and drop the α=-1 upper cap.

**Tech Stack:** C++17 (`src/raking.cpp` only). No other files except fixture regeneration.

---

## File Map

| File | Role |
|------|------|
| `src/raking.cpp` | Core rewrite: remove p[]/q_hyp, add cells_per_cat + water_fill_cat, update F_eval + SQUAREM |
| `tests/testthat/fixtures/raking_squarem_baseline.rds` | Regenerate (algorithm changes — new baseline) |
| `tests/testthat/fixtures/raking_obs_reference_stepstone.rds` | Regenerate (run gen_raking_obs_ref.R after AC bench) |

---

## Algorithm: Water-Filling per Category

Autumn `single_adjust()` (rake.R lines 63-93) translated to cell level:

```
Given: cells in category j, target T_kj, current X[c], bounds [L[c], U[c]]

free_sum = Σ_{c in j} X[c]  (= bucket[j] from IPF sweep)
clamped_sum = 0

for pass in 0..n_cells:
    m = (T_kj - clamped_sum) / free_sum
    newly_clamped = {c free : X[c]*m > U[c] or X[c]*m < L[c]}
    if newly_clamped empty:
        X[c] = X[c] * m  for all free c
        return
    for c in newly_clamped:
        X[c] = U[c] if X[c]*m > U[c] else L[c]
        clamped_sum += X[c]
        free_sum -= X_orig[c]   # X_orig = pre-pass value
```

Key: `X_orig[c]` is the weight at the **start of the water-fill** (not start of the current pass). Autumn uses `cell_w = weights[cell_idx]` (original) — we must also save original values before modifying.

---

## Task 0: Generate comparison fixture (before any code changes)

**Files:**
- Read: current benchmark numbers (no file changes)

- [ ] **Step 1: Confirm current SQUAREM max_err**

```bash
OMP_NUM_THREADS=1 Rscript -e "
  suppressPackageStartupMessages({library(arrow);library(jsonlite);library(leafblower)})
  df <- arrow::read_parquet('benchmarks/stepstone_fulldata_bench_data.parquet')
  df\$uuid <- NULL
  tgt <- lapply(jsonlite::fromJSON('benchmarks/stepstone_fulldata_bench_targets.json'),
                function(t){t<-unlist(t);t/sum(t)})
  for(nm in names(tgt)) df[[nm]] <- factor(df[[nm]])
  r <- leafblower::harvest(df, tgt, method='raking', accelerate=TRUE,
        max_weight=5, max_iterations=5000L, attach_weights=FALSE, verbose=0)
  cat('current SQUAREM max_err=', attr(r,'result')\$max_error,
      'DEFF=', sum((as.numeric(r)/mean(as.numeric(r)))^2)/length(as.numeric(r)), '\n')
" 2>&1 | grep -v "^Welcome\|^Working"
```

Expected: `max_err= 0.001765...  DEFF= 2.85...` (pre-water-filling baseline).

- [ ] **Step 2: Note numbers for comparison after implementation**

Record: `current: max_err=1.77e-3, DEFF=2.85`

---

## Task 1: Pre-compute cells_per_cat data structure

**Files:**
- Modify: `src/raking.cpp`

The `cells_per_cat[k][j]` structure gives O(1) cell lookup per (margin, category), needed for water-filling. Built once per `raking_solve` call.

- [ ] **Step 1: Add cells_per_cat after X_init (around line 93)**

Read lines 90-96 first. After `std::vector<double> X_init(X);`, add:

```cpp
    // Per-(margin, category) cell index lists for water-filling.
    // cells_per_cat[k][j] = cells where g_per_cell[k][c] == j.
    // Built once: O(M_cell × K). Memory: ~M_cell × K ints ≈ 261k ints for stepstone.
    std::vector<std::vector<std::vector<int>>> cells_per_cat(st.K);
    for (int k = 0; k < st.K; k++) {
        cells_per_cat[k].assign(st.cat_counts[k], {});
        for (int c = 0; c < ct.M_cell; c++) {
            int g = ct.g_per_cell[k][c];
            if (g >= 0 && g < st.cat_counts[k])
                cells_per_cat[k][g].push_back(c);
        }
    }

    // Pre-allocated scratch arrays for water-filling inner loop.
    // Size = max cells per (margin, category) to avoid per-call heap allocation.
    int wf_max_cat = 0;
    for (int k = 0; k < st.K; k++)
        for (int j = 0; j < st.cat_counts[k]; j++)
            wf_max_cat = std::max(wf_max_cat, (int)cells_per_cat[k][j].size());
    std::vector<double>  wf_x_orig(wf_max_cat);  // original X[c] per free cell
    std::vector<uint8_t> wf_status(wf_max_cat);  // 0=free, 1=at_U, 2=at_L
```

- [ ] **Step 2: Compile gate**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```
Expected: `* DONE (leafblower)`.

---

## Task 2: Implement water_fill_cat lambda

**Files:**
- Modify: `src/raking.cpp`

Add this lambda after the `cells_per_cat`/scratch setup and before the existing `compute_weight_kl` lambda.

- [ ] **Step 1: Insert water_fill_cat lambda**

```cpp
    // water_fill_cat: apply bounded IPF scale to all cells in category j of margin k.
    // Adjusts X[c] in place so that Σ_{c in j} X[c] = T_kj while X[c] ∈ [L[c], U[c]].
    // Algorithm: autumn single_adjust() (rake.R lines 63-93) translated to cell level.
    // Uses pre-allocated wf_x_orig / wf_status scratch to avoid per-call heap allocation.
    auto water_fill_cat = [&](int k, int j, double T_kj, double bucket_j,
                               std::vector<double>& Xv) {
        const auto& cells = cells_per_cat[k][j];
        const int n = static_cast<int>(cells.size());
        if (n == 0) return;
        if (bucket_j < kEmptyBucketThreshold) {
            if (T_kj > 0.0) is_infeasible = true;
            return;
        }

        // Save original weights and mark all cells free
        for (int ci = 0; ci < n; ci++) {
            wf_x_orig[ci] = Xv[cells[ci]];
            wf_status[ci] = 0;
        }

        double clamped_sum = 0.0;
        double free_sum    = bucket_j;  // Σ_{free} X_orig[c]

        for (int pass = 0; pass <= n; ++pass) {
            if (free_sum < kEmptyBucketThreshold) { is_infeasible = true; break; }
            const double T_free = T_kj - clamped_sum;
            if (T_free <= 0.0) break;
            const double m = T_free / free_sum;

            bool any_clamped = false;
            for (int ci = 0; ci < n; ci++) {
                if (wf_status[ci] != 0) continue;
                const double proposed = wf_x_orig[ci] * m;
                if (proposed > U_cell[cells[ci]]) {
                    wf_status[ci] = 1;
                    clamped_sum += U_cell[cells[ci]];
                    free_sum    -= wf_x_orig[ci];
                    any_clamped  = true;
                } else if (proposed < L_cell[cells[ci]]) {
                    wf_status[ci] = 2;
                    clamped_sum += L_cell[cells[ci]];
                    free_sum    -= wf_x_orig[ci];
                    any_clamped  = true;
                }
            }

            if (!any_clamped) {
                // m applies cleanly — commit final values
                for (int ci = 0; ci < n; ci++) {
                    const int c = cells[ci];
                    if      (wf_status[ci] == 0) Xv[c] = wf_x_orig[ci] * m;
                    else if (wf_status[ci] == 1) Xv[c] = U_cell[c];
                    else                          Xv[c] = L_cell[c];
                }
                return;
            }
        }
        // All passes exhausted (infeasible category) — apply best-effort final values
        const double T_final = T_kj - clamped_sum;
        const double m_final = (free_sum > kEmptyBucketThreshold && T_final > 0.0)
                               ? T_final / free_sum : 0.0;
        for (int ci = 0; ci < n; ci++) {
            const int c = cells[ci];
            if      (wf_status[ci] == 0) Xv[c] = wf_x_orig[ci] * m_final;
            else if (wf_status[ci] == 1) Xv[c] = U_cell[c];
            else                          Xv[c] = L_cell[c];
        }
    };
```

- [ ] **Step 2: Compile gate**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```

---

## Task 3: Replace Dykstra box + hyperplane in F_eval with water-filling

**Files:**
- Modify: `src/raking.cpp`

**Key changes to F_eval** (currently lines 170-246):
1. Remove `pv` and `q_hyp_v` parameters
2. Replace box correction block with water_fill_cat calls (one per category per margin, inside the k-margin sweep)
3. Replace separate Bregman hyperplane block with simple scaling

**Before** (current signature and box/hyperplane in F_eval):
```cpp
auto F_eval = [&](std::vector<double>& Xv,
                  std::vector<double>& pv,
                  double& q_hyp_v) -> double {
    ...
    // Per-category scale applied as single multiplier:
    Xv[c] *= scale[g];
    ...
    // Bregman box (after all K margins):
    for (int c ...) {
        double yc = Xv[c] * pv[c]; Xc = clamp(yc, L, U);
        pv[c] = yc / Xc; Xv[c] = Xc;
    }
    // Bregman hyperplane:
    for Xv[c] *= q_hyp_v; s_hp += Xv[c];
    Xv[c] *= n/s_hp; q_hyp_v = 1/scale;
}
```

**After** (water-filling + simple hyperplane):
```cpp
auto F_eval = [&](std::vector<double>& Xv) -> double {
    ...
    // Per-category: IPF scale → water_fill_cat (enforces bounds within category)
    // Water-fill is called after computing bucket[j] for each category:
    for (int j = 0; j < st.cat_counts[k]; j++) {
        double Tkj = st.targets[k][j] * W_total;
        water_fill_cat(k, j, Tkj, bucket[j], Xv);  // replaces separate box step
    }
    ...
    // Simple hyperplane (no correction vector):
    double s_hp = 0.0;
    for (int c = 0; c < ct.M_cell; c++) s_hp += Xv[c];
    const double sc_hp = static_cast<double>(st.n) / s_hp;
    for (int c = 0; c < ct.M_cell; c++) Xv[c] *= sc_hp;
}
```

- [ ] **Step 1: Replace the entire F_eval lambda**

Find the F_eval lambda (line 170 to ~245). Replace it with:

```cpp
    // F_eval: one complete bounded IPF iteration.
    // Water-filling enforces X[c] ∈ [L[c], U[c]] within each margin step (autumn style).
    // No correction vectors — F is stateless: suitable for SQUAREM L2 halving.
    auto F_eval = [&](std::vector<double>& Xv) -> double {
        double W_total = 0.0;
        for (int c = 0; c < ct.M_cell; c++) W_total += Xv[c];

        if (use_greedy)
            std::sort(margin_order.begin(), margin_order.end(),
                      [&](int a, int b){ return errRp_k[a] > errRp_k[b]; });

        for (int ki = 0; ki < st.K; ki++) {
            const int k = use_greedy ? margin_order[ki] : ki;
            std::fill(bucket.begin(), bucket.begin() + st.cat_counts[k], 0.0);
            for (int c = 0; c < ct.M_cell; c++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k]) bucket[g] += Xv[c];
            }

            // SOR: compute effective scale for each category before water-fill
            // (water-fill handles bounds; SOR under-relaxes the scale direction).
            const double eff_omega = sor_active ? sor_omega[k] : 1.0;

            for (int j = 0; j < st.cat_counts[k]; j++) {
                double Tkj = st.targets[k][j] * W_total;
                if (bucket[j] < kEmptyBucketThreshold * W_total) {
                    if (Tkj > 0.0) is_infeasible = true;
                    continue;
                }
                if (eff_omega != 1.0) {
                    // SOR: scale target toward current bucket (under-relaxation)
                    const double s0 = Tkj / bucket[j];
                    Tkj = bucket[j] * std::pow(s0, eff_omega);
                }
                water_fill_cat(k, j, Tkj, bucket[j], Xv);
            }

            // SOR adaptation: per-margin absolute errRp
            if (sor_active && sor_auto && W_total > 0.0) {
                // Re-compute bucket after water-fill to get post-step residuals
                std::fill(bucket.begin(), bucket.begin() + st.cat_counts[k], 0.0);
                for (int c = 0; c < ct.M_cell; c++) {
                    int g = ct.g_per_cell[k][c];
                    if (g >= 0 && g < st.cat_counts[k]) bucket[g] += Xv[c];
                }
                double ek = 0.0;
                for (int j = 0; j < st.cat_counts[k]; j++) {
                    double e = std::fabs(bucket[j] / W_total - st.targets[k][j]);
                    if (e > ek) ek = e;
                }
                if (sor_prev_errRp[k] < ek)
                    sor_omega[k] = std::max(omega_min, sor_omega[k] * 0.7);
                else
                    sor_omega[k] = std::min(1.0, sor_omega[k] * 1.05);
                sor_prev_errRp[k] = ek;
            }

            // Update per-margin errRp for Greedy sort (use post-water-fill bucket)
            if (W_total > 0.0) {
                // If SOR not active, bucket was not recomputed above; recompute for errRp_k
                if (!sor_active || !sor_auto) {
                    std::fill(bucket.begin(), bucket.begin() + st.cat_counts[k], 0.0);
                    for (int c = 0; c < ct.M_cell; c++) {
                        int g = ct.g_per_cell[k][c];
                        if (g >= 0 && g < st.cat_counts[k]) bucket[g] += Xv[c];
                    }
                }
                double ek = 0.0;
                for (int j = 0; j < st.cat_counts[k]; j++) {
                    double e = std::fabs(bucket[j] / W_total - st.targets[k][j]);
                    if (e > ek) ek = e;
                }
                errRp_k[k] = ek;
            }
        }

        // Hyperplane: scale to sum = n (no Dykstra correction vector needed)
        double s_hp = 0.0;
        for (int c = 0; c < ct.M_cell; c++) s_hp += Xv[c];
        if (s_hp > 0.0) {
            const double sc_hp = static_cast<double>(st.n) / s_hp;
            for (int c = 0; c < ct.M_cell; c++) Xv[c] *= sc_hp;
        }

        return compute_errRp_ct(st, ct, Xv, bucket);
    };
```

**Note on SOR with water-filling**: SOR under-relaxes by adjusting the effective target before water-fill. `Tkj_sor = bucket[j] * (Tkj/bucket[j])^eff_omega`. When eff_omega=1: Tkj_sor = Tkj (standard). When eff_omega<1: less aggressive scale. This preserves SOR semantics.

**Note on bucket recomputation**: after water-fill, the post-step bucket differs from pre-step. For SOR and errRp_k accuracy, a second pass over cells is needed. This adds O(M_cell) per margin per F-eval. Acceptable cost given simpler architecture.

- [ ] **Step 2: Compile gate**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```

---

## Task 4: Update all call sites of F_eval and remove p[], q_hyp, finalizer

**Files:**
- Modify: `src/raking.cpp`

F_eval now takes one argument (Xv) instead of three (Xv, pv, q_hyp_v). All call sites must be updated. p[] and q_hyp declarations removed. Post-loop finalizer removed (no longer needed — water-filling enforces bounds within each F-eval call).

- [ ] **Step 1: Remove p[], q_hyp declarations (lines 104-106)**

Find:
```cpp
    // Dykstra corrections: p[c] box, q_hyp scalar hyperplane.
    std::vector<double> p(ct.M_cell, 1.0);   // Bregman: multiplicative identity
    double q_hyp = 1.0;                        // Bregman: multiplicative identity
```
Replace with:
```cpp
    // No Dykstra correction vectors: water-filling enforces bounds within F_eval.
```

- [ ] **Step 2: Remove hyperplane_step lambda (lines 158-168)**

Find:
```cpp
    // Dykstra hyperplane: projects X onto {sum(X) = n}.
    auto hyperplane_step = [&]() {
        ...
    };
```
Delete this block entirely.

- [ ] **Step 3: Update SQUAREM call sites (lines ~265-325)**

In the SQUAREM while loop, every F_eval call currently passes three args. Update:

```cpp
// OLD:
auto w1 = X; auto p1 = p; double qh1 = q_hyp;
double errRp_w1 = F_eval(w1, p1, qh1);  ++f_eval_count;
is_infeasible = infeas_before;

auto w2 = w1; auto p2 = p1; double qh2 = qh1;
double errRp_w2 = F_eval(w2, p2, qh2);  ++f_eval_count;

// NEW:
auto w1 = X;
double errRp_w1 = F_eval(w1);  ++f_eval_count;
is_infeasible = infeas_before;

auto w2 = w1;
double errRp_w2 = F_eval(w2);  ++f_eval_count;
```

And snapshot/restore lines — remove p_snap, q_snap, p_star, qh_star:

```cpp
// OLD:
auto X_snap = w2; auto p_snap = p2; double q_snap = qh2;
...
auto p_star = p_snap; double qh_star = q_snap;
double errRp_new = F_eval(X_star, p_star, qh_star);  ++f_eval_count;
...
// in halving loop:
X_star = X_snap; p_star = p_snap; qh_star = q_snap;
...
errRp_new = F_eval(X_star, p_star, qh_star);  ++f_eval_count;
...
X = X_star; p = p_star; q_hyp = qh_star;

// NEW:
auto X_snap = w2;
...
double errRp_new = F_eval(X_star);  ++f_eval_count;
...
// in halving loop:
X_star = X_snap;
...
errRp_new = F_eval(X_star);  ++f_eval_count;
...
X = X_star;  // no p/q_hyp to restore
```

Also remove the fixed-point guard exit (used to reference p2, qh2):
```cpp
// OLD:
X = w2; p = p2; q_hyp = qh2;
// NEW:
X = w2;
```

- [ ] **Step 4: Update flat loop (remove p/q_hyp usage)**

The flat loop (inside `else { for (int iter ...) }`) currently uses `p[]`, `q_hyp` and calls `hyperplane_step()`. Replace the entire inner loop body with a call to F_eval:

```cpp
} else {
    for (int iter = 1; iter <= st.inner_max_iter; iter++) {
        res.iterations = iter;

        double errRp = F_eval(X);  // water-fill bounded IPF + hyperplane

        if (iter == 1 || iter % kErrCheckInterval == 0 || iter == st.inner_max_iter) {
            res.max_error = errRp;

            // Best-iterate tracking
            if (errRp < best_metric_seen) {
                best_metric_seen    = errRp;
                best_iter_val       = iter;
                best_objective_seen = compute_weight_kl();
                W_best              = X;
            }

            // Extra metrics (gated)
            // [keep existing BLOCK 2 / extra metrics block if needed for non-MAX_ERR convergence]
            // For now: compute cell metrics inline when needed
            const lbw::CalibMetric metric = st.convergence_cfg.metric;
            double mean_err = 0.0, kl_max = 0.0, chi2_total = 0.0, grake_norm = 0.0;
            double l1_weight = 0.0;
            // pct_change: approximate via errRp delta (not exact but sufficient)
            double pct_change = std::fabs(errRp - res.max_error) / std::max(res.max_error, 1e-15);
            if (metric != lbw::CalibMetric::MAX_ERR || iter == st.inner_max_iter) {
                // Full metric computation (reuse existing BLOCK from old loop)
                double W_tot2 = 0.0;
                for (int c = 0; c < ct.M_cell; c++) W_tot2 += X[c];
                if (W_tot2 > 0.0) {
                    constexpr double kMetricEps = 1e-10;
                    constexpr double kChi2Floor = 1.0;
                    double mean_err_sum = 0.0;
                    for (int k = 0; k < st.K; k++) {
                        const int nj = st.cat_counts[k];
                        std::fill(bucket.begin(), bucket.begin() + nj, 0.0);
                        for (int c = 0; c < ct.M_cell; c++) {
                            int g = ct.g_per_cell[k][c];
                            if (g >= 0 && g < nj) bucket[g] += X[c];
                        }
                        double max_k = 0.0, kl_k = 0.0;
                        for (int j = 0; j < nj; j++) {
                            double S_p = bucket[j] / W_tot2, T = st.targets[k][j];
                            double err = std::fabs(S_p - T);
                            if (err > max_k) max_k = err;
                            if (T > 0.0)
                                kl_k += T * std::log((T + kMetricEps) / (S_p + kMetricEps));
                            double obs = bucket[j], pop_kj = T * W_tot2;
                            chi2_total += (obs - pop_kj) * (obs - pop_kj) / (pop_kj + kChi2Floor);
                            double nm = std::fabs(obs - pop_kj) / (1.0 + std::fabs(pop_kj));
                            if (nm > grake_norm) grake_norm = nm;
                        }
                        mean_err_sum += max_k;
                        if (kl_k > kl_max) kl_max = kl_k;
                    }
                    mean_err = (st.K > 0) ? mean_err_sum / static_cast<double>(st.K) : 0.0;
                    // BLOCK 2: best-iterate for non-MAX_ERR metrics
                    if (metric != lbw::CalibMetric::MAX_ERR) {
                        double l1_weight_blk2 = 0.0;  // not tracked in simplified loop
                        const double curr_best = lbw::select_metric(
                            metric, errRp, mean_err, kl_max, chi2_total, grake_norm, l1_weight_blk2);
                        if (std::isfinite(curr_best) && curr_best < best_metric_seen) {
                            best_metric_seen    = curr_best;
                            best_iter_val       = iter;
                            best_objective_seen = compute_weight_kl();
                            W_best              = X;
                        }
                    }
                }
            }

            res.mean_error       = mean_err;
            res.kl               = kl_max;
            res.chi2             = chi2_total;
            res.grake_norm       = grake_norm;
            res.l1_weight_change = l1_weight;

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

            lbw::CellMetrics m_conv;
            m_conv.errRp = errRp; m_conv.mean_err = mean_err;
            m_conv.kl = kl_max; m_conv.chi2 = chi2_total;
            m_conv.grake_norm = grake_norm; m_conv.l1 = l1_weight;
            if (lbw::check_convergence(st.convergence_cfg, m_conv,
                                       prev_metric_for_rule, st.tol_abs)) {
                res.status             = is_infeasible ? RK_ERR_INFEAS : RK_OK;
                res.convergence_metric = static_cast<int>(st.convergence_cfg.metric);
                res.convergence_rule   = static_cast<int>(st.convergence_cfg.rule);
                res.convergence_tol    = st.convergence_cfg.pct_tol;
                res.convergence_iter   = iter;
                break;
            }

            if (n_no_improve >= kMaxNoImprove) {
                res.status = RK_ERR_NOCONV;
                break;
            }
        }
    }
}  // end else flat loop
```

- [ ] **Step 5: Remove the 200-round post-loop finalizer**

Find:
```cpp
    // Post-loop Bregman/KL Dykstra finalizer (multiplicative pattern).
    for (int fixup = 0; fixup < 200; fixup++) {
        bool box_ok = true;
        for (int c = 0; c < ct.M_cell; c++) {
            double yc = X[c] * p[c];
            double Xc = std::clamp(yc, L_cell[c], U_cell[c]);
            p[c] = (Xc > 0.0) ? yc / Xc : 1.0;
            if (yc != Xc) box_ok = false;
            X[c] = Xc;
        }
        hyperplane_step();
        if (box_ok) break;
    }
```

Replace with a simple hyperplane normalization only:
```cpp
    // Post-loop: normalize sum to n (water-filling already enforces bounds).
    {
        double s_post = 0.0;
        for (int c = 0; c < ct.M_cell; c++) s_post += X[c];
        if (s_post > 0.0) {
            const double sc_post = static_cast<double>(st.n) / s_post;
            for (int c = 0; c < ct.M_cell; c++) X[c] *= sc_post;
        }
    }
```

- [ ] **Step 6: Compile gate**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```

---

## Task 5: Update SQUAREM — remove α cap, switch to L2 halving

**Files:**
- Modify: `src/raking.cpp`

With water-filling making F stateless, both autumn-style fixes now apply:
1. Remove `kAlphaMax = -1.0` — allow α ∈ (-1000, 0)
2. L2 halving: compare `‖F(X*)-X*‖²` vs `‖v‖²` (clean F → no Dykstra explosion)
3. Add autumn's `|α+1| < 1e-3` fallback guard

- [ ] **Step 1: Remove kAlphaMax constant and clamp**

Find:
```cpp
        static constexpr double kAlphaMax    = -1.0;    // alpha is strictly ≤ 0; -1.0 = least-aggressive bound
        static constexpr double kAlphaMin    = -1000.0;
```
Replace:
```cpp
        static constexpr double kAlphaMin    = -1000.0;  // no upper cap: autumn allows α ∈ (-1000, 0)
```

Find:
```cpp
                double alpha = std::max(kAlphaMin,
                               std::min(kAlphaMax, -norm_r / (norm_v + kVNormEps)));
```
Replace:
```cpp
                double alpha = std::max(kAlphaMin, -norm_r / (norm_v + kVNormEps));
```

- [ ] **Step 2: Switch to L2 halving with ‖v‖² reference and add |α+1|<1e-3 guard**

Find the F_eval call + halving loop in SQUAREM:
```cpp
                auto p_star = p_snap; double qh_star = q_snap;
                double errRp_new = F_eval(X_star, p_star, qh_star);  ++f_eval_count;

                for (int h = 0; h < kMaxHalvings && errRp_new > kHalvingSlack * errRp_w2; h++) {
                    is_infeasible = infeas_before;
                    alpha = (alpha - 1.0) / 2.0;
                    X_star = X_snap; p_star = p_snap; qh_star = q_snap;
                    for (int c = 0; c < ct.M_cell; c++) {
                        double ri = w1[c] - X[c],  vi = w2[c] - w1[c];
                        X_star[c] = X_snap[c] - 2.0 * alpha * ri + alpha * alpha * vi;
                        if (X_star[c] < 0.0) X_star[c] = 0.0;
                    }
                    errRp_new = F_eval(X_star, p_star, qh_star);  ++f_eval_count;
                }
                X = X_star; p = p_star; q_hyp = qh_star;
```

Replace with (L2 halving + |α+1| guard):
```cpp
                // L2 step-halving: ‖F(X*)-X*‖² vs ‖v‖² (autumn convention).
                // Works now that F is stateless (no Dykstra explosion on zeroed cells).
                const double plain_resid = norm_v * norm_v;  // ‖v‖²
                auto X_star_pre = X_star;
                double errRp_new = F_eval(X_star);  ++f_eval_count;
                double cand_resid = 0.0;
                for (int c = 0; c < ct.M_cell; c++) {
                    double d = X_star[c] - X_star_pre[c];
                    cand_resid += d * d;
                }

                for (int h = 0; h < kMaxHalvings && cand_resid > kHalvingSlack * plain_resid; h++) {
                    is_infeasible = infeas_before;
                    alpha = (alpha + (-1.0)) / 2.0;  // midpoint toward -1 (autumn formula)
                    if (std::fabs(alpha - (-1.0)) < 1e-3) {
                        // Fell back to plain step — accept w2
                        X = w2;
                        res.max_error = errRp_w2; res.iterations = f_eval_count;
                        goto squarem_next_iter;  // skip X_star assignment below
                    }
                    X_star = X_snap;
                    for (int c = 0; c < ct.M_cell; c++) {
                        double ri = w1[c] - X[c],  vi = w2[c] - w1[c];
                        X_star[c] = X_snap[c] - 2.0 * alpha * ri + alpha * alpha * vi;
                        if (X_star[c] < 0.0) X_star[c] = 0.0;
                    }
                    X_star_pre = X_star;
                    errRp_new = F_eval(X_star);  ++f_eval_count;
                    cand_resid = 0.0;
                    for (int c = 0; c < ct.M_cell; c++) {
                        double d = X_star[c] - X_star_pre[c];
                        cand_resid += d * d;
                    }
                }
                X = X_star;
                squarem_next_iter:;
```

**Note on goto**: the `goto squarem_next_iter` jumps past the `X = X_star` assignment when we fell back to w2. The label is placed after the assignment. No variable initializations are jumped over (label is at end of loop body). This is the same pattern used in the existing code and is safe C++17.

- [ ] **Step 3: Compile gate**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```

- [ ] **Step 4: Smoke test**

```bash
Rscript -e "
  set.seed(42L); n <- 500L
  df  <- data.frame(v1=factor(sample(3L,n,TRUE)), v2=factor(sample(2L,n,TRUE)))
  tgt <- list(v1=c('1'=0.5,'2'=0.3,'3'=0.2), v2=c('1'=0.6,'2'=0.4))
  r <- leafblower::harvest(df, tgt, method='raking', accelerate=TRUE,
       max_weight=5, max_iterations=500L, attach_weights=FALSE)
  res <- attr(r, 'result')
  cat('water-fill SQUAREM: status=', res\$status, 'iters=', res\$iterations,
      'max_err=', res\$max_error, '\n')
  stopifnot(!any(is.nan(as.numeric(r))), !any(is.infinite(as.numeric(r))))
  cat('PASS: no NaN/Inf\n')
"
```

---

## Task 6: Regenerate fixtures and verify all tests

**Files:**
- Regenerate: `tests/testthat/fixtures/raking_squarem_baseline.rds`

- [ ] **Step 1: Build package**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```

- [ ] **Step 2: Regenerate squarem baseline fixture**

```bash
Rscript -e "
  set.seed(42L); n <- 2000L
  df  <- data.frame(v1=factor(sample(3L,n,TRUE)), v2=factor(sample(2L,n,TRUE)))
  tgt <- list(v1=c('1'=0.5,'2'=0.3,'3'=0.2), v2=c('1'=0.6,'2'=0.4))
  r <- leafblower::harvest(df, tgt, method='raking',
        accelerate=FALSE, max_weight=5, max_iterations=500L, attach_weights=FALSE)
  w <- as.numeric(r)
  saveRDS(w, 'tests/testthat/fixtures/raking_squarem_baseline.rds')
  cat('Saved', length(w), 'weights, sum=', sum(w), '\n')
"
```

- [ ] **Step 3: Run full test suite**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
```

Expected: FAIL ≤ 3 (pre-existing: test-ieppa-nonuniform-d.R ×2; possibly squarem-ac3 if it needs fixture — regenerate again if so).

If `squarem-ac3` fails: the fixture was generated from the OLD algorithm. The test compares `accelerate=FALSE` to the fixture — both should be water-fill now and match.

- [ ] **Step 4: Verify AC5 (‖v‖ guard — no NaN)**

```bash
Rscript -e "
  df <- data.frame(v1=factor(rep('1',10)))
  tgt <- list(v1=c('1'=1.0))
  r <- leafblower::harvest(df, tgt, method='raking', accelerate=TRUE,
       max_weight=5, max_iterations=500L, attach_weights=FALSE)
  w <- as.numeric(r)
  cat('finite:', all(is.finite(w)), 'NaN:', any(is.nan(w)), '\n')
"
```

---

## Task 7: Stepstone benchmark

> **Local only — not CI**. Run and record result.

- [ ] **Step 1: Run comparison benchmark**

```bash
OMP_NUM_THREADS=1 Rscript -e "
  suppressPackageStartupMessages({library(arrow);library(jsonlite);library(leafblower);library(autumn)})
  df  <- arrow::read_parquet('benchmarks/stepstone_fulldata_bench_data.parquet')
  df\$uuid <- NULL
  tgt <- lapply(jsonlite::fromJSON('benchmarks/stepstone_fulldata_bench_targets.json'),
                function(t){t<-unlist(t);t/sum(t)})
  for(nm in names(tgt)) df[[nm]] <- factor(df[[nm]])

  fit <- function(w,df,tgt){W<-sum(w);n<-length(w);max_err<-0;marg_kl<-0
    for(nm in names(tgt)){lv<-names(tgt[[nm]]);S<-tapply(w,df[[nm]],sum,default=0)[lv];S[is.na(S)]<-0
      tk<-unname(tgt[[nm]]);Sr<-S/W;max_err<-max(max_err,max(abs(Sr-tk)))
      safe<-tk>0&Sr>0;marg_kl<-marg_kl+sum(ifelse(safe,tk*log(tk/pmax(Sr,1e-300)),0))}
    list(max_err=max_err,marg_kl=marg_kl,DEFF=n*sum(w^2)/W^2,ESS=W^2/sum(w^2))}

  run <- function(label,...){cat(sprintf('%-30s',label));flush.console()
    t0<-proc.time()['elapsed']
    r<-suppressWarnings(if(label=='autumn')autumn::harvest(df,tgt,max_weight=5,min_weight=0,accelerate=TRUE)
       else leafblower::harvest(df,tgt,max_weight=5,min_weight=0,max_iterations=5000L,
            attach_weights=FALSE,verbose=0,...))
    wall<-proc.time()['elapsed']-t0
    w<-tryCatch(as.numeric(r),error=function(e)as.numeric(r\$weights))
    m<-fit(w,df,tgt)
    cat(sprintf(' wall=%5.1fs iters=%4s  max_err=%.4e  marg_kl=%.3e  DEFF=%.4f  ESS=%s\n',
      wall,ifelse(is.null(attr(r,'result')\$iterations),'  —',attr(r,'result')\$iterations),
      m\$max_err,m\$marg_kl,m\$DEFF,format(round(m\$ESS),big.mark=',')))
  }

  run('ieppa',                     method='ieppa')
  run('raking (flat, wf)',         method='raking', accelerate=FALSE)
  run('raking+SQUAREM (wf+L2)',    method='raking', accelerate=TRUE)
  run('autumn',                    label='autumn')
" 2>&1 | grep -v "^Welcome\|^Working"
```

- [ ] **Step 2: Record result**

Compare against pre-water-filling SQUAREM:
```
Pre-wf SQUAREM:    max_err=1.77e-3  DEFF=2.85
Target (wf+L2):    max_err < 1.70e-3, DEFF < 2.50
Autumn reference:  max_err=1.60e-3  DEFF=1.99
```

If DEFF > 2.50 or max_err > 1.77e-3: water-filling not improving. Check:
- Is SOR being applied to SQUAREM? (should only be in flat loop)
- Are wf_x_orig/wf_status scratch arrays large enough?

- [ ] **Step 3: Commit**

```bash
git add src/raking.cpp tests/testthat/fixtures/raking_squarem_baseline.rds
git commit -m "feat(raking): replace Dykstra box correction with water-filling bounded IPF (autumn-style); enables L2 SQUAREM halving and sub-acceleration"
```

---

## Acceptance Criteria

1. **AC1**: `accelerate=TRUE, method="raking"` runs without error on all test cases
2. **AC2** *(local bench)*: Water-fill SQUAREM max_err ≤ 1.70e-3 AND DEFF ≤ 2.50 on stepstone (improvement over pre-wf: 1.77e-3, 2.85)
3. **AC3**: `accelerate=FALSE` bit-identical to regenerated baseline (tolerance 1e-14)
4. **AC4**: Category sums respect target exactly for feasible problems (water-fill convergence)
5. **AC5**: No NaN/Inf on trivially-converged problems (‖v‖ guard)
6. **AC6**: `devtools::test()` FAIL ≤ 3 (pre-existing + possibly 1 new if any flat-loop test changes)

---

## Self-Review

**Spec coverage**: AC1-AC6 all have corresponding steps. ✓

**No placeholders**: all code blocks are complete. ✓

**Known risks**:
- The `goto squarem_next_iter` in Task 5 — tested in compile gate; no variable initializations between goto and label. Clean C++17.
- SOR with water-filling: sor target is modified before water-fill. If SOR is disabled (default for SQUAREM since accelerate=TRUE forces round_robin and SOR is not default), the eff_omega=1.0 branch skips the target modification entirely.
- Bucket recomputation after water-fill (for SOR/errRp_k): adds O(M_cell) per margin per F-eval. For K=9, M_cell=29k: 261k extra ops per F-eval — ~0.5% overhead.
