// raking.cpp — classical bounded raking via cyclic IPF + Dykstra projections.
//
// Algorithm composition (no published proof for the hybrid; convergence is
// empirical, verified by the test suite):
//   1. Cyclic IPF marginal update (Deming & Stephan 1940; Csiszar 1975
//      proves I-projection convergence onto the intersection of margin-
//      constraint hyperplanes for feasible problems).
//   2. Additive Dykstra box projection onto [min_weight, max_weight]^n
//      (Boyle & Dykstra 1986). Euclidean corrections accumulate in q[i].
//   3. Additive Dykstra hyperplane projection onto {w : sum(w) = n}.
//      Euclidean corrections accumulate in q_hyp[i].
//
// The composition of multiplicative (IPF) and Euclidean (Dykstra)
// projections is not covered by Boyle-Dykstra's single-metric theorem. The
// inline comment below is explicit: Euclidean Dykstra corrections are NOT
// applied to the IPF step — they diverge against multiplicative updates.
// Convergence is therefore empirical-only for this hybrid.
//
// References (classical IPF family):
//   Deming W. E. & Stephan F. F. (1940), "On a Least Squares Adjustment of
//     a Sampled Frequency Table", Ann. Math. Stat. 11, 427-444.
//   Csiszar I. (1975), "I-Divergence Geometry of Probability Distributions
//     and Minimization Problems", Ann. Probab. 3, 146-158.
//   Boyle J. P. & Dykstra R. L. (1986), "A Method for Finding Projections
//     onto the Intersection of Convex Sets in Hilbert Spaces", Advances in
//     Order Restricted Statistical Inference, Springer.

#include "lbw_config.h"
#include "raking.hpp"
#include "cell_table.hpp"
#include "calib_dispatch.hpp"
#include "leafblower.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <limits>
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

    // Dykstra hyperplane: projects X onto {sum(X) = n}.
    // M_cell cells each get +shift where shift = (n-s)/M_cell → total correction = (n-s).
    auto hyperplane_step = [&]() {
        double s = 0.0;
        for (int c = 0; c < ct.M_cell; c++) { X[c] += q_hyp; s += X[c]; }
        double shift = (static_cast<double>(st.n) - s) / static_cast<double>(ct.M_cell);
        for (int c = 0; c < ct.M_cell; c++) X[c] += shift;
        q_hyp = -shift;
    };

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
        hyperplane_step();

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

            lbw::CellMetrics m_conv;
            m_conv.errRp      = errRp;
            m_conv.mean_err   = mean_err;
            m_conv.kl         = kl_max;
            m_conv.chi2       = chi2_total;
            m_conv.grake_norm = grake_norm;
            m_conv.l1         = l1_weight;
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
                if (st.verbose >= 1) {
                    char msg[256];
                    std::snprintf(msg, 256,
                        "raking: errRp stalled for %d consecutive checks (last=%.2e, window_min=%.2e); likely near-infeasible bounds. "
                        "aborting at iter %d.", n_no_improve, errRp, min_errRp_window, iter);
                    st.log(msg);
                }
                break;
            }
        }
    }

    if (is_infeasible && res.status == RK_ERR_NOCONV)
        res.status = RK_ERR_INFEAS;

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
        hyperplane_step();
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
    apply_obs_expansion(ct, X, X_init, st.n, lo, hi_obs, st.weights);

    return res;
}
} // namespace lbw
