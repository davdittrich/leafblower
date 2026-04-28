// raking.cpp — bounded raking via cyclic IPF with water-filling box projections.
//
// Algorithm: F_eval performs one complete IPF sweep with per-category water-filling
// (autumn single_adjust) that enforces X[c] ∈ [L_c, U_c] within each margin step.
// No Dykstra correction vectors; bounds are enforced inline by water_fill_cat.
//
// Both the flat loop and SQUAREM-accelerated path call F_eval directly.
// Water-filling convergence: cyclic KL projections onto bounded margin constraints
// converge to the bounded KL minimum (Csiszar-Tusnady 1984).
//
// References (classical IPF family):
//   Deming W. E. & Stephan F. F. (1940), "On a Least Squares Adjustment of
//     a Sampled Frequency Table", Ann. Math. Stat. 11, 427-444.
//   Csiszar I. (1975), "I-Divergence Geometry of Probability Distributions
//     and Minimization Problems", Ann. Probab. 3, 146-158.

#include "lbw_config.h"
#include "raking.hpp"
#include "cell_table.hpp"
#include "calib_dispatch.hpp"
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


// Constrained raking solver: cyclic IPF for marginal projections + Dykstra box correction.
// Marginal step: pure IPF (Bregman/multiplicative projection — Euclidean Dykstra corrections
// diverge on multiplicative projections and are not used here).
// Box step: Dykstra additive correction q[i] prevents cycling at the [lo,hi]^n boundary.
// inner_max_iter is the single iteration budget; outer_max_iter is unused.
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
    int max_cats = *std::max_element(st.cat_counts, st.cat_counts + st.K);
    std::vector<double> bucket(max_cats);

    // Descent monitor
    double min_errRp_window = std::numeric_limits<double>::infinity();
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

    // Weight-space KL objective: Σ_c X[c]*log(X[c]/X_init[c])/n
    // Distinct from m.kl (marginal KL) — this is what raking actually minimizes.
    auto compute_weight_kl = [&]() -> double {
        double wkl = 0.0;
        const double inv_n = 1.0 / static_cast<double>(st.n);
        for (int c = 0; c < ct.M_cell; c++) {
            if (X_init[c] > 0.0 && X[c] > 0.0)
                wkl += X[c] * std::log(X[c] / X_init[c]) * inv_n;
        }
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
        const double m_final = (free_sum > kEmptyBucketThreshold && T_final > 0.0)
                               ? T_final / free_sum : 0.0;
        for (int ci = 0; ci < n; ci++) {
            const int c = cells[ci];
            if      (wf_status[ci] == 0) Xv[c] = wf_x_orig[ci] * m_final;
            else if (wf_status[ci] == 1) Xv[c] = U_cell[c];
            else                          Xv[c] = L_cell[c];
        }
    };

    // F_eval: one complete bounded IPF iteration using water-filling.
    // Water-filling enforces X[c] ∈ [L[c], U[c]] within each margin step.
    // No correction vectors — F is stateless: enables SQUAREM L2 halving.
    // Mathematical basis: cyclic KL projections onto bounded margin constraints
    // converge to the bounded KL minimum (Csiszar-Tusnady 1984).
    auto F_eval = [&](std::vector<double>& Xv) -> double {
        double W_total = 0.0;
        for (int c = 0; c < ct.M_cell; c++) W_total += Xv[c];

        if (use_greedy)
            std::sort(margin_order.begin(), margin_order.end(),
                      [&](int a, int b){ return errRp_k[a] > errRp_k[b]; });

        for (int ki = 0; ki < st.K; ki++) {
            const int k = use_greedy ? margin_order[ki] : ki;

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
                if (bucket[j] < kEmptyBucketThreshold * W_total) {
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

        return compute_errRp_ct(st, ct, Xv, bucket);
    };

    if (st.accelerate) {
        static constexpr double kAlphaMin    = -1000.0;  // no upper cap: autumn allows α ∈ (-1000, 0)
        static constexpr double kVNormEps    = 1e-300;
        static constexpr double kVNormRel    = 1e-10;
        static constexpr double kHalvingSlack = 1.01;
        static constexpr int    kMaxHalvings  = 16;

        if (st.inner_max_iter >= 3) {
            int f_eval_count = 0;
            while (f_eval_count + 3 <= st.inner_max_iter) {
                // Save infeasibility state. Intermediate F_eval probes (w1, w2, halving)
                // may false-flag infeasibility on extrapolated iterates.
                // Only the final accepted F_eval contributes to infeasibility status.
                bool infeas_before = is_infeasible;

                auto w1 = X;
                double errRp_w1 = F_eval(w1);  ++f_eval_count;
                (void)errRp_w1;  // advances IPF side effects (errRp_k); value unused
                is_infeasible = infeas_before;

                auto w2 = w1;
                double errRp_w2 = F_eval(w2);  ++f_eval_count;
                is_infeasible = infeas_before;

                res.iterations = f_eval_count;

                double norm_r = 0.0, norm_v = 0.0, norm_w2 = 0.0;
                for (int c = 0; c < ct.M_cell; c++) {
                    double ri = w1[c] - X[c],  vi = w2[c] - w1[c];
                    norm_r  += ri * ri;
                    norm_v  += vi * vi;
                    norm_w2 += w2[c] * w2[c];
                }
                norm_r = std::sqrt(norm_r);
                norm_v = std::sqrt(norm_v);
                norm_w2 = std::sqrt(norm_w2);

                // Fixed-point guard: ‖v‖/‖w2‖ < threshold → already converged
                if (norm_v / (norm_w2 + kVNormEps) < kVNormRel) {
                    X = w2;
                    res.max_error        = errRp_w2;
                    res.status           = is_infeasible ? RK_ERR_INFEAS : RK_OK;
                    res.convergence_iter = f_eval_count;
                    break;
                }

                // CBB step: α = -‖r‖/‖v‖, capped at kAlphaMin only (no upper cap).
                // Autumn allows α ∈ (-1000, 0): sub-acceleration (α > -1) is valid.
                double alpha = std::max(kAlphaMin, -norm_r / (norm_v + kVNormEps));

                // Snapshot: state at w2, before extrapolation
                auto X_snap = w2;

                // Extrapolate X* = X_snap - 2α·r + α²·v; clamp to ≥ 0
                auto X_star = X_snap;
                for (int c = 0; c < ct.M_cell; c++) {
                    double ri = w1[c] - X[c],  vi = w2[c] - w1[c];
                    X_star[c] = X_snap[c] - 2.0 * alpha * ri + alpha * alpha * vi;
                    if (X_star[c] < 0.0) X_star[c] = 0.0;
                }

                // L2 step-halving: ‖F(X*)-X*‖² vs ‖v‖².
                // Works with water-filling F: hyperplane step is near-no-op (sum ≈ n),
                // so ‖F(X*)-X*‖ cleanly reflects IPF movement, not Dykstra explosion.
                const double plain_resid = norm_v * norm_v;  // ‖v‖²
                auto X_star_pre = X_star;
                double errRp_new = F_eval(X_star);  ++f_eval_count;
                double cand_resid = 0.0;
                for (int c = 0; c < ct.M_cell; c++) {
                    double d = X_star[c] - X_star_pre[c];
                    cand_resid += d * d;
                }

                // Step-halving with boolean flag (not goto) to ensure convergence
                // check runs for both accepted-extrapolation and fell-back-to-w2 paths.
                bool fell_back = false;
                for (int h = 0; h < kMaxHalvings && cand_resid > kHalvingSlack * plain_resid; h++) {
                    is_infeasible = infeas_before;  // discard probe's infeasibility
                    alpha = (alpha + (-1.0)) / 2.0;  // midpoint toward -1 (autumn formula)
                    if (std::fabs(alpha - (-1.0)) < 1e-3) {
                        // Fell back to plain step
                        X = w2; errRp_new = errRp_w2; fell_back = true; break;
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
                if (!fell_back) X = X_star;
                res.max_error  = errRp_new;
                res.iterations = f_eval_count;

                // Best-iterate tracking
                if (errRp_new < best_metric_seen) {
                    best_metric_seen    = errRp_new;
                    best_iter_val       = f_eval_count;
                    best_objective_seen = compute_weight_kl();
                    W_best              = X;
                }

                // Convergence criterion
                lbw::CellMetrics m_conv; m_conv.errRp = errRp_new;
                if (lbw::check_convergence(st.convergence_cfg, m_conv,
                                           prev_metric_for_rule, st.tol_abs)) {
                    res.status             = is_infeasible ? RK_ERR_INFEAS : RK_OK;
                    res.convergence_metric = static_cast<int>(st.convergence_cfg.metric);
                    res.convergence_rule   = static_cast<int>(st.convergence_cfg.rule);
                    res.convergence_tol    = st.convergence_cfg.pct_tol;
                    res.convergence_iter   = f_eval_count;
                    break;
                }

                if (st.verbose >= 1) {
                    char msg[256];
                    std::snprintf(msg, 256, "raking[sq] f_eval=%d errRp=%.2e alpha=%.4g",
                                  f_eval_count, errRp_new, alpha);
                    st.log(msg);
                }

                if (!std::isfinite(min_errRp_window)) {
                    min_errRp_window = errRp_new; n_no_improve = 0;
                } else {
                    const double eps = std::max(0.01 * min_errRp_window, st.tol_abs);
                    if (errRp_new < min_errRp_window - eps) {
                        min_errRp_window = errRp_new; n_no_improve = 0;
                    } else {
                        n_no_improve++;
                    }
                }
                if (n_no_improve >= kMaxNoImprove) { res.status = RK_ERR_NOCONV; break; }
            }
        }
    } else {
        for (int iter = 1; iter <= st.inner_max_iter; iter++) {
            res.iterations = iter;

            double errRp = F_eval(X);

            if (iter == 1 || iter % kErrCheckInterval == 0 || iter == st.inner_max_iter) {
                res.max_error = errRp;

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
                    double W_tot2 = 0.0;
                    for (int c = 0; c < ct.M_cell; c++) W_tot2 += X[c];
                    constexpr double kMetricEps = 1e-10;
                    constexpr double kChi2Floor = 1.0;
                    double mean_sum = 0.0;
                    if (W_tot2 > 0.0) {
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
                            mean_sum += max_k;
                            if (kl_k > kl_max) kl_max = kl_k;
                        }
                        mean_err = (st.K > 0) ? mean_sum / static_cast<double>(st.K) : 0.0;
                        if (metric != lbw::CalibMetric::MAX_ERR) {
                            const double curr_best = lbw::select_metric(
                                metric, errRp, mean_err, kl_max, chi2_total, grake_norm, l1_weight);
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

    if (is_infeasible && res.status == RK_ERR_NOCONV)
        res.status = RK_ERR_INFEAS;

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
    res.convergence_solver_objective = best_objective_seen;
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
    apply_obs_expansion(ct, X, X_init, st.n, lo, hi_obs, st.weights);

    return res;
}
} // namespace lbw
