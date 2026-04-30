// raking.cpp — bounded raking via cyclic IPF with water-filling box projections.
//
// Algorithm: F_eval performs one complete IPF sweep with per-category water-filling
// (autumn single_adjust) that enforces X[c] ∈ [L_c, U_c] within each margin step.
// No Dykstra correction vectors; bounds are enforced inline by water_fill_cat.
//
// Both the flat loop and SRAA-m-accelerated path call F_eval directly.
// Water-filling convergence: cyclic KL projections onto bounded margin constraints
// converge to the bounded KL minimum (Csiszar-Tusnady 1984).
//
// References (classical IPF family):
//   Deming W. E. & Stephan F. F. (1940), "On a Least Squares Adjustment of
//     a Sampled Frequency Table", Ann. Math. Stat. 11, 427-444.
//   Csiszar I. (1975), "I-Divergence Geometry of Probability Distributions
//     and Minimization Problems", Ann. Probab. 3, 146-158.

#include "lbw_config.h"
#include "lbw_math.hpp"
#include "raking.hpp"
#include "cell_table.hpp"
#include "calib_dispatch.hpp"
#include "sraa.hpp"
#include "leafblower.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <limits>
#include <numeric>
#include <vector>

namespace lbw {

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


// Constrained raking solver: cyclic IPF with per-category water-filling box projections.
// Each margin step applies water_fill_cat() — KL projection onto
// {Σ_c Xv[c]=T_kj, L_c≤Xv[c]≤U_c} (Csiszar-Tusnady 1984; autumn single_adjust).
// No Dykstra correction vectors; stateless F enables SRAA-m Anderson Acceleration.
// inner_max_iter is the single iteration budget; outer_max_iter is unused.
RakingResult raking_solve(CalibState& st) {
    static constexpr double kAbsoluteZeroThreshold = 1e-15;  // bucket_j / free_sum is genuinely zero
    static constexpr double kRelativeZeroFraction  = 1e-15;  // bucket[j] < fraction * W_total
    static constexpr int    kErrCheckInterval     = 10;
    static constexpr int    kMaxNoImprove         = 5;

    RakingResult res;
    res.base.status     = RK_ERR_BUDGET;  // initial; overwritten by criterion/stall; remains if budget exhausted
    res.base.iterations = 0;
    res.base.max_error  = 1.0;

    // Build cell table: O(n log n) one-time cost.
    CellTable ct;
    if (build_cell_table(st.n, st.K, st.group_ids, st.cat_counts,
                         st.weights, ct) != 0) {
        // RakingResult has no message field; caller gets RK_ERR_BADARG status.
        res.base.status = RK_ERR_BADARG;
        return res;
    }

    // Initial cell masses: X[c] = Σ_{i∈c} st.weights[i]
    std::vector<double> X(ct.M_cell, 0.0);
    for (int i = 0; i < st.n; i++)
        X[ct.cell_of[i]] += st.weights[i];
    std::vector<double> X_init(X);

    // Per-(margin, category) cell index lists for water-filling.
    // cells_per_cat[k][j] = cells where g_per_cell[k][c] == j.
    // Built once: O(M_cell × K). Memory: ~M_cell × K ints.
    std::vector<std::vector<std::vector<int>>> cells_per_cat(st.K);
    for (int k = 0; k < st.K; k++) {
        cells_per_cat[k].assign(st.cat_counts[k], {});
        for (int c = 0; c < ct.M_cell; c++) {
            int g = ct.g_per_cell[k][c];
            if (g >= 0 && g < st.cat_counts[k])
                cells_per_cat[k][g].push_back(c);
        }
    }

    // Pre-allocated scratch for water-filling inner loop.
    int wf_max_cat = 0;
    for (int k = 0; k < st.K; k++)
        for (int j = 0; j < st.cat_counts[k]; j++)
            wf_max_cat = std::max(wf_max_cat, (int)cells_per_cat[k][j].size());
    std::vector<double>  wf_x_orig(wf_max_cat);
    std::vector<uint8_t> wf_status(wf_max_cat);

    // Cell bounds: L_c = lo * n_per_cell[c], U_c = hi * n_per_cell[c]
    const double lo = st.min_weight;
    const double hi = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
    std::vector<double> L_cell(ct.M_cell), U_cell(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) {
        L_cell[c] = lo * static_cast<double>(ct.n_per_cell[c]);
        U_cell[c] = hi * static_cast<double>(ct.n_per_cell[c]);
    }

    // No Dykstra correction vectors: water-filling enforces bounds within F_eval.

    bool is_infeasible = false;

    // Pre-computed reciprocals: n_per_cell is constant; avoids division per cell per super-step.
    std::vector<double> inv_n_per_cell(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++)
        inv_n_per_cell[c] = 1.0 / static_cast<double>(ct.n_per_cell[c]);

    int max_cats = *std::max_element(st.cat_counts, st.cat_counts + st.K);
    std::vector<double> bucket(max_cats);

    // Descent monitor — tracks solver loss function (weight KL for flat loop, errRp for SRAA-m).
    double min_loss_window = std::numeric_limits<double>::infinity();
    int n_no_improve = 0;

    // pct/l1 tracking at cell level
    std::vector<double> X_prev(X);

    // Convergence rule state
    double prev_metric_for_rule = std::numeric_limits<double>::infinity();

    // Best-iterate tracking (cell-level snapshot)
    double best_metric_seen    = std::numeric_limits<double>::infinity();
    double best_objective_seen = 0.0;  // weight KL at best_iter (A1 fix)
    int    best_iter_val    = 0;
    std::vector<double> W_best(ct.M_cell, 0.0);

    // Scratch buffers for compute_weight_kl vectorized log path.
    std::vector<double> kl_ratio_scratch(ct.M_cell, 0.0);
    std::vector<double> kl_weight_scratch(ct.M_cell, 0.0);

    // SOR: wire st.sor_cfg into raking's IPF step (same API as ieppa).
    // Apply only when bounds are active (oscillation risk).
    const bool sor_active  = st.sor_cfg.enabled &&
                              (st.min_weight > 0.0 || hi < 1e300);
    const bool sor_auto    = st.sor_cfg.auto_adapt;
    const double omega_min = st.sor_cfg.omega_min;    // default 0.3
    const double omega_init = st.sor_cfg.omega_init;  // default 1.0
    std::vector<double> sor_omega(st.K, omega_init);
    std::vector<double> sor_prev_errRp(st.K, std::numeric_limits<double>::infinity());

    // Greedy scheduler: per-margin residuals + sort order.
    // errRp_k[k] updated each iter during sweep using bucket[] sums.
    // Survives across iterations — first iter uses uniform priority (1/K init).
    std::vector<double> errRp_k(st.K, 1.0 / st.K);
    std::vector<int> margin_order(st.K);
    std::iota(margin_order.begin(), margin_order.end(), 0);
    const bool use_greedy = (st.scheduler.mode == SchedulerMode::GREEDY);
    // R8: greedy reordering breaks SRAA-m's fixed-point geometry; demote silently.
    bool use_greedy_effective = use_greedy;
    if (st.accelerate && use_greedy_effective) {
        use_greedy_effective = false;
        st.log("[raking] greedy scheduler disabled under SRAA-m acceleration; using round_robin");
    }

    // Weight-space KL objective: Σ_c X[c]*log(X[c]/X_init[c])/n
    // Distinct from m.kl (marginal KL) — this is what raking actually minimizes.
    auto compute_weight_kl = [&]() -> double {
        const double inv_n = 1.0 / static_cast<double>(st.n);
        int valid_count = 0;
        for (int c = 0; c < ct.M_cell; c++) {
            if (X_init[c] > 0.0 && X[c] > 0.0) {
                kl_ratio_scratch[valid_count] = X[c] / X_init[c];
                kl_weight_scratch[valid_count] = X[c];
                valid_count++;
            }
        }
        lbw::bulk_log(kl_ratio_scratch.data(), kl_ratio_scratch.data(), valid_count);
        double wkl = 0.0;
        for (int i = 0; i < valid_count; i++)
            wkl += kl_weight_scratch[i] * kl_ratio_scratch[i] * inv_n;
        return std::isfinite(wkl) ? wkl : 0.0;
    };

    // water_fill_cat: KL projection of category j (margin k) onto
    // {Σ_{c∈j} Xv[c] = T_kj, L_cell[c] ≤ Xv[c] ≤ U_cell[c]}.
    // Algorithm: autumn single_adjust() (rake.R:63-93) at cell level.
    // wf_x_orig[] and wf_status[] are pre-allocated scratch (re-used every call).
    auto water_fill_cat = [&](int k, int j, double T_kj, double bucket_j,
                               std::vector<double>& Xv) {
        const auto& cells = cells_per_cat[k][j];
        const int n = static_cast<int>(cells.size());
        if (n == 0) return;
        if (bucket_j < kAbsoluteZeroThreshold) {
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
            if (free_sum < kAbsoluteZeroThreshold) { is_infeasible = true; break; }
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
                // m applies cleanly — commit final values and return
                for (int ci = 0; ci < n; ci++) {
                    const int c = cells[ci];
                    if      (wf_status[ci] == 0) Xv[c] = wf_x_orig[ci] * m;
                    else if (wf_status[ci] == 1) Xv[c] = U_cell[c];
                    else                          Xv[c] = L_cell[c];
                }
                return;
            }
        }
        // All passes exhausted (infeasible category) — commit best-effort values
        const double T_final = T_kj - clamped_sum;
        const double m_final = (free_sum > kAbsoluteZeroThreshold && T_final > 0.0)
                               ? T_final / free_sum : 0.0;
        for (int ci = 0; ci < n; ci++) {
            const int c = cells[ci];
            if      (wf_status[ci] == 0) Xv[c] = wf_x_orig[ci] * m_final;
            else if (wf_status[ci] == 1) Xv[c] = U_cell[c];
            else                          Xv[c] = L_cell[c];
        }
    };

    // last_F_metrics: populated by F_eval via compute_cell_metrics; reused in SRAA-m convergence check.
    // Declared here (outer scope) so the [&] lambda capture includes it.
    lbw::CellMetrics last_F_metrics;

    // F_eval: one complete bounded IPF iteration using water-filling.
    // Water-filling enforces X[c] ∈ [L[c], U[c]] within each margin step.
    // No correction vectors — F is stateless: enables SRAA-m L2 halving.
    // Mathematical basis: cyclic KL projections onto bounded margin constraints
    // converge to the bounded KL minimum (Csiszar-Tusnady 1984).
    auto F_eval = [&](std::vector<double>& Xv) -> double {
        double W_total = 0.0;
        for (int c = 0; c < ct.M_cell; c++) W_total += Xv[c];

        if (use_greedy_effective)
            std::sort(margin_order.begin(), margin_order.end(),
                      [&](int a, int b){ return errRp_k[a] > errRp_k[b]; });

        for (int ki = 0; ki < st.K; ki++) {
            const int k = use_greedy_effective ? margin_order[ki] : ki;

            // Compute pre-water-fill bucket (needed for errRp_k and SOR)
            std::fill(bucket.begin(), bucket.begin() + st.cat_counts[k], 0.0);
            for (int c = 0; c < ct.M_cell; c++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k]) bucket[g] += Xv[c];
            }

            // errRp_k for Greedy: pre-water-fill residual (how bad margin is NOW).
            // Must be computed BEFORE water-fill — post-fill residual is always near zero.
            if (W_total > 0.0) {
                double ek = 0.0;
                for (int j = 0; j < st.cat_counts[k]; j++) {
                    double e = std::fabs(bucket[j] / W_total - st.targets[k][j]);
                    if (e > ek) ek = e;
                }
                errRp_k[k] = ek;
            }

            // SOR effective omega (under-relaxation toward current bucket)
            const double eff_omega = sor_active ? sor_omega[k] : 1.0;

            // Apply water-filling per category
            for (int j = 0; j < st.cat_counts[k]; j++) {
                double Tkj = st.targets[k][j] * W_total;
                if (bucket[j] < kRelativeZeroFraction * W_total) {
                    if (Tkj > 0.0) is_infeasible = true;
                    continue;
                }
                if (eff_omega != 1.0) {
                    // SOR: under-relax target — equivalent to pow(scale, omega) at cell level
                    const double s0 = Tkj / bucket[j];
                    Tkj = bucket[j] * std::pow(s0, eff_omega);
                }
                water_fill_cat(k, j, Tkj, bucket[j], Xv);
            }

            // SOR adaptation: post-water-fill residual for omega adjustment
            if (sor_active && sor_auto && W_total > 0.0) {
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
        }

        // Hyperplane: normalize sum to n (no correction vector needed;
        // water-fill preserves category sums, so this is near-no-op for feasible problems)
        double s_hp = 0.0;
        for (int c = 0; c < ct.M_cell; c++) s_hp += Xv[c];
        if (s_hp > 0.0) {
            const double sc_hp = static_cast<double>(st.n) / s_hp;
            for (int c = 0; c < ct.M_cell; c++) Xv[c] *= sc_hp;
        }

        // compute_cell_metrics is a strict superset of compute_errRp_ct at same O(K×M_cell) cost.
        // Stores all metrics in last_F_metrics — SRAA-m convergence check reuses them directly.
        last_F_metrics = lbw::compute_cell_metrics(st, ct, Xv, static_cast<double>(st.n), bucket);
        return last_F_metrics.errRp;
    };

    // SRAA-m state (replaces SRAA-m step-halving while-loop)
    lbw::SRAAState rk_sraa;
    if (st.accelerate) rk_sraa.init(ct.M_cell, lbw::kSRAAm);

    if (st.accelerate) {
        // SRAA-m replaced by SRAA-m. X_prev_sq (stall detection) removed — SRAA-m-specific.
        // Loop mirrors greenkhorn: sraa_step called once per outer iteration until budget exhausted.
        int f_eval_budget = st.inner_max_iter;
        int f_evals_used = 0;
        int rk_outer_stall_count = 0;
        while (f_evals_used + 1 <= f_eval_budget) {
            rk_sraa.F_cur = X;  // seed F_cur with current X before each sraa_step call
            auto r = lbw::sraa_step(F_eval, X, L_cell, U_cell, rk_sraa);
            f_evals_used += r.f_evals;
            res.base.max_error  = r.err_result;
            res.base.iterations = f_evals_used;

            // Best-iterate tracking
            if (r.err_result < best_metric_seen) {
                best_metric_seen    = r.err_result;
                best_iter_val       = f_evals_used;
                best_objective_seen = compute_weight_kl();
                W_best              = X;
            }

            // Outer revert: if quality has degraded significantly, revert to best
            {
                double curr_quality_rk = r.err_result;
                if (curr_quality_rk > best_metric_seen * (1.0 + lbw::kSRAAOuterSlack)) {
                    if (++rk_outer_stall_count >= lbw::kSRAAOuterStallWindow) {
                        X = W_best;               // revert; no S_flat rebuild — F_eval handles it
                        rk_sraa.clear();
                        rk_outer_stall_count = 0;
                    }
                } else { rk_outer_stall_count = 0; }
            }

            // Convergence check
            {
                lbw::CellMetrics m_conv = last_F_metrics;
                m_conv.errRp = r.err_result;
                if (lbw::check_convergence(st.convergence_cfg, m_conv,
                                           prev_metric_for_rule, st.tol_abs)) {
                    res.base.status             = RK_OK;
                    res.base.convergence_metric = static_cast<int>(st.convergence_cfg.metric);
                    res.base.convergence_rule   = static_cast<int>(st.convergence_cfg.rule);
                    res.base.convergence_tol    = st.convergence_cfg.pct_tol;
                    res.base.convergence_iter   = f_evals_used;
                    break;
                }
            }

            if (st.verbose >= 1) {
                char msg[256];
                std::snprintf(msg, 256, "raking[sraa] f_evals=%d errRp=%.2e aa=%d",
                              f_evals_used, r.err_result, (int)r.aa_accepted);
                st.log(msg);
            }
        }
    } else {
        for (int iter = 1; iter <= st.inner_max_iter; iter++) {
            res.base.iterations = iter;

            double errRp = F_eval(X);

            if (iter == 1 || iter % kErrCheckInterval == 0 || iter == st.inner_max_iter) {
                res.base.max_error = errRp;

                // Best-iterate tracking (MAX_ERR metric)
                if (st.convergence_cfg.metric == lbw::CalibMetric::MAX_ERR) {
                    if (errRp < best_metric_seen) {
                        best_metric_seen    = errRp;
                        best_iter_val       = iter;
                        best_objective_seen = compute_weight_kl();
                        W_best              = X;
                    }
                }

                // L1 weight change from previous check
                double l1_weight = 0.0;
                for (int c = 0; c < ct.M_cell; c++)
                    l1_weight += std::fabs(X[c] - X_prev[c]);
                l1_weight /= static_cast<double>(st.n);
                for (int c = 0; c < ct.M_cell; c++) X_prev[c] = X[c];

                // Extra metrics for non-MAX_ERR convergence
                const lbw::CalibMetric metric = st.convergence_cfg.metric;
                double mean_err = 0.0, kl_max = 0.0, chi2_total = 0.0, grake_norm = 0.0;
                const bool need_extra = (metric == lbw::CalibMetric::MEAN_ERR ||
                                         metric == lbw::CalibMetric::KL ||
                                         metric == lbw::CalibMetric::CHI2 ||
                                         metric == lbw::CalibMetric::GRAKE_NORM ||
                                         iter == st.inner_max_iter);
                if (need_extra) {
                    // R2: was ~40-line duplicate of lbw::compute_cell_metrics.
                    double W_tot2 = 0.0;
                    for (int c = 0; c < ct.M_cell; c++) W_tot2 += X[c];
                    if (W_tot2 > 0.0) {
                        const lbw::CellMetrics m = lbw::compute_cell_metrics(st, ct, X, W_tot2, bucket);
                        mean_err   = m.mean_err;
                        kl_max     = m.kl;
                        chi2_total = m.chi2;
                        grake_norm = m.grake_norm;
                        if (metric != lbw::CalibMetric::MAX_ERR) {
                            const double curr_best = lbw::select_metric(metric, m);
                            if (std::isfinite(curr_best) && curr_best < best_metric_seen) {
                                best_metric_seen    = curr_best;
                                best_iter_val       = iter;
                                best_objective_seen = compute_weight_kl();
                                W_best              = X;
                            }
                        }
                    }
                }

                res.base.mean_error       = mean_err;
                res.base.kl               = kl_max;
                res.base.chi2             = chi2_total;
                res.base.grake_norm       = grake_norm;
                res.base.l1_weight_change = l1_weight;

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
                    // Converged — do NOT override with INFEAS here. water_fill_cat may
                    // transiently set is_infeasible when cells temporarily hit U_cell during
                    // convergence. INFEAS only overrides on stall (post-loop check below).
                    res.base.status             = RK_OK;
                    res.base.convergence_metric = static_cast<int>(st.convergence_cfg.metric);
                    res.base.convergence_rule   = static_cast<int>(st.convergence_cfg.rule);
                    res.base.convergence_tol    = st.convergence_cfg.pct_tol;
                    res.base.convergence_iter   = iter;
                    break;
                }

                // Weight KL stall: monotone for water-filling IPF (Csiszar-Tusnady).
                // KL plateau ↔ constrained KL minimum — correct stall signal.
                // Guard: wkl ≤ tol_abs means effectively at optimum → converged (not stalled).
                const double wkl_flat = compute_weight_kl();
                if (wkl_flat <= st.tol_abs) {
                    res.base.status = RK_OK; res.base.convergence_iter = iter; break;
                }
                if (!std::isfinite(min_loss_window)) {
                    min_loss_window = wkl_flat; n_no_improve = 0;
                } else if (wkl_flat < min_loss_window * (1.0 - st.convergence_cfg.pct_tol)) {
                    min_loss_window = wkl_flat; n_no_improve = 0;
                } else {
                    n_no_improve++;
                }

                if (n_no_improve >= kMaxNoImprove) {
                    res.base.status = RK_ERR_STALL;
                    if (st.verbose >= 1) {
                        char msg[256];
                        std::snprintf(msg, 256,
                            "raking: errRp stalled for %d consecutive checks "
                            "(wkl=%.2e); aborting at iter %d.",
                            n_no_improve, wkl_flat, iter);
                        st.log(msg);
                    }
                    break;
                }
            }
        }
    }  // end else flat loop

    // Water-filling detects partial infeasibility (some categories can't reach targets
    // within bounds). With Dykstra this was silently masked by bound violations.
    // For stalled iterations: return STALL (status=5) + best weights rather than
    // hard-erroring — caller can use the best achievable calibration.

    // Note: is_infeasible (water_fill_cat capacity exhaustion) is intentionally
    // NOT promoted to RK_ERR_INFEAS here. Upper-bound capacity capping is normal
    // bounded calibration — the solver returns the best achievable constrained
    // optimum with binding capacity constraints, not an infeasibility error.
    // True structural INFEAS (lower bounds exceed target) is caught at validation.

    // Post-loop: normalize sum to n (water-filling already enforces bounds)
    {
        double s_post = 0.0;
        for (int c = 0; c < ct.M_cell; c++) s_post += X[c];
        if (s_post > 0.0) {
            const double sc_post = static_cast<double>(st.n) / s_post;
            for (int c = 0; c < ct.M_cell; c++) X[c] *= sc_post;
        }
    }

    // Best-iterate: normalize cell snapshot, then expand to obs.
    res.base.convergence_solver_objective = best_objective_seen;
    res.base.convergence_minimized_metric = static_cast<int>(st.convergence_cfg.metric);
    res.base.best_error = best_metric_seen;
    res.base.best_iter  = best_iter_val;
    if (std::isfinite(best_metric_seen)) {
        double s = 0.0;
        for (int c = 0; c < ct.M_cell; c++) s += W_best[c];
        if (s > 0.0) {
            const double sc = static_cast<double>(st.n) / s;
            for (int c = 0; c < ct.M_cell; c++) W_best[c] *= sc;
        }
        res.base.best_weights.resize(st.n);
        const double hi_obs = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
        for (int i = 0; i < st.n; i++) {
            int c = ct.cell_of[i];
            double mult = (X_init[c] > 0.0) ? W_best[c] / X_init[c] : 1.0;
            res.base.best_weights[i] = std::clamp(st.weights[i] * mult, lo, hi_obs);
        }
    } else {
        res.base.best_weights.assign(st.n, 0.0);
    }

    // Post-exit obs expansion: w_i = d_i × X[c]/X_init[c], hard clamp.
    const double hi_obs = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
    apply_obs_expansion(ct, X, X_init, st.n, lo, hi_obs, st.weights);

    return res;
}
} // namespace lbw
