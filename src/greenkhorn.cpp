#include "lbw_config.h"
#include "greenkhorn.hpp"
#include "calib_dispatch.hpp"
#include "calib_validate.hpp"
#include "cell_table.hpp"
#include "sraa.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <limits>
#include <numeric>

namespace lbw {

GreenkornResult greenkhorn_solve(CalibState& st) {
    GreenkornResult res;
    res.base.status = RK_ERR_NOCONV;

    CellTable ct;
    std::vector<double> X_init;
    double hi_eff;
    std::vector<double> L_cell, U_cell;
    std::vector<int> cat_offset;
    int n_cats_total;
    if (lbw::solver_setup_ct(st, ct, X_init, hi_eff, L_cell, U_cell,
                              cat_offset, n_cats_total, res) != RK_OK)
        return res;

    const int M = ct.M_cell;
    const int K = st.K;
    int max_cats = lbw::max_cats_count(K, st.cat_counts);

    // Build cells_per_cat[k][j] = list of cell indices in bucket (k,j)
    auto cells_per_cat = lbw::build_cells_per_cat(ct, K, st.cat_counts);

    // Working copy; clamp to capacity bounds (X_init kept for obs-expansion reference)
    std::vector<double> X(X_init);
    for (int c = 0; c < M; c++)
        X[c] = std::clamp(X[c], L_cell[c], U_cell[c]);

    // Total mass W (maintained incrementally)
    double W = 0.0;
    for (int c = 0; c < M; c++) W += X[c];

    // Bucket sums S_flat[k * S_stride + j]
    const int S_stride = max_cats;
    std::vector<double> S_flat(K * S_stride, 0.0);
    for (int k = 0; k < K; k++)
        for (int j = 0; j < st.cat_counts[k]; j++)
            for (int c : cells_per_cat[k][j]) S_flat[k*S_stride+j] += X[c];

    // Per-margin residuals
    std::vector<double> errRp(K, 0.0);
    auto compute_errRp_k = [&](int k) -> double {
        double e = 0.0;
        for (int j = 0; j < st.cat_counts[k]; j++) {
            double ach = (W > 0.0) ? S_flat[k*S_stride+j] / W : 0.0;
            e = std::max(e, std::abs(ach - st.targets[k][j]));
        }
        return e;
    };
    for (int k = 0; k < K; k++) errRp[k] = compute_errRp_k(k);

    // G8c: best-iterate tracking via BestIterTracker (replaces ad-hoc vars).
    // greenkhorn does not minimise KL directly, so best_objective stays ∞.
    BestIterTracker best;
    {
        double init_errRp = *std::max_element(errRp.begin(), errRp.end());
        best.update(init_errRp, std::numeric_limits<double>::infinity(), 0, X);
    }

    // Convergence state
    std::vector<double> bucket_scratch(max_cats, 0.0);
    double prev_metric = std::numeric_limits<double>::infinity();
    double first_errRp = -1.0;  // B5: captured at first convergence check for stall/budget classify
    const CalibConvergenceCfg& cfg = st.convergence_cfg;
    constexpr int    kErrCheckInterval   = 10;
    constexpr double kEmptyBucketThreshold = 1e-15;

    // Single Greenkhorn step on margin k_step
    // Updates X, W, S_flat, errRp
    auto greenkhorn_step = [&](int k_step) {
        if (W <= 0.0) return;
        for (int j = 0; j < st.cat_counts[k_step]; j++) {
            double S_kj = S_flat[k_step * S_stride + j];
            if (S_kj < kEmptyBucketThreshold * W) continue;
            double f = st.targets[k_step][j] * W / S_kj;
            for (int c : cells_per_cat[k_step][j]) {
                double X_old = X[c];
                double X_new = std::clamp(X[c] * f, L_cell[c], U_cell[c]);
                double delta = X_new - X_old;
                X[c] = X_new;
                W += delta;
                if (std::abs(delta) < 1e-300) continue;
                for (int k2 = 0; k2 < K; k2++) {
                    if (k2 == k_step) continue;
                    int g2 = ct.g_per_cell[k2][c];
                    if (g2 >= 0 && g2 < st.cat_counts[k2])
                        S_flat[k2 * S_stride + g2] += delta;
                }
            }
            S_flat[k_step * S_stride + j] = 0.0;
            for (int c : cells_per_cat[k_step][j]) S_flat[k_step*S_stride+j] += X[c];
        }
        for (int k = 0; k < K; k++) errRp[k] = compute_errRp_k(k);
    };

    // SRAA-m state (replaces SQUAREM)
    lbw::SRAAState grk_sraa;
    std::vector<double> Sv_sraa;
    if (st.accelerate && K > 0) {
        grk_sraa.init(M, lbw::kSRAAm);
        Sv_sraa.assign((size_t)K * S_stride, 0.0);
    }
    auto f_eval_sraa = [&](std::vector<double>& Xv) -> double {
        // Xv must be grk_sraa.F_cur or grk_sraa.scratch — NEVER the outer X
        double Wv = 0.0;
        std::fill(Sv_sraa.begin(), Sv_sraa.end(), 0.0);
        for (int c = 0; c < M; c++) {
            Wv += Xv[c];
            for (int k = 0; k < K; k++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k])
                    Sv_sraa[k * S_stride + g] += Xv[c];
            }
        }
        std::swap(X, Xv); std::swap(S_flat, Sv_sraa); std::swap(W, Wv);
        // K greenkhorn steps: each picks the current argmax-violation margin.
        // F(X) is deterministic (same X → same errRp → same argmax sequence),
        // so the fixed-point map is stationary and SRAA is valid.
        for (int ki = 0; ki < K; ki++) {
            int k_star = (int)(std::max_element(errRp.begin(), errRp.end()) - errRp.begin());
            greenkhorn_step(k_star);
        }
        std::swap(X, Xv); std::swap(S_flat, Sv_sraa); std::swap(W, Wv);
        // Wv now holds W_after_K_steps. If zero (all cells at lower bound = 0),
        // the AA input was degenerate; return +inf so the NaN guard in sraa_step rejects it.
        if (Wv <= 0.0) return std::numeric_limits<double>::infinity();
        return *std::max_element(errRp.begin(), errRp.end());
    };

    // Main loop
    int outer_stall_count = 0;
    for (int iter = 0; iter < st.inner_max_iter; iter++) {
        // W<=0 guard
        if (W <= 0.0) {
            res.base.status = RK_ERR_INFEAS;
            std::snprintf(res.message, sizeof(res.message),
                "greenkhorn: total mass W<=0 (all cells at zero bound)");
            break;
        }

        if (st.accelerate && K > 0) {
            grk_sraa.F_cur = X;  // seed F_cur with current X before each sraa_step call
            auto r = lbw::sraa_step(f_eval_sraa, X, L_cell, U_cell, grk_sraa);
            // Rebuild W and S_flat from the accepted X (f_eval_sraa swaps back on exit)
            W = 0.0;
            std::fill(S_flat.begin(), S_flat.end(), 0.0);
            for (int c = 0; c < M; c++) {
                W += X[c];
                for (int k = 0; k < K; k++) {
                    int g = ct.g_per_cell[k][c];
                    if (g >= 0 && g < st.cat_counts[k])
                        S_flat[k * S_stride + g] += X[c];
                }
            }
            for (int k = 0; k < K; k++) errRp[k] = compute_errRp_k(k);
            res.base.iterations += K * r.f_evals;  // B6: f_evals=1 (plain) or 2 (AA attempted)
        } else {
            // Pure Greenkhorn: single margin per step
            int k_star = (int)(std::max_element(errRp.begin(), errRp.end()) - errRp.begin());
            greenkhorn_step(k_star);
            res.base.iterations = iter + 1;
        }

        // Best-iterate
        double curr_max = *std::max_element(errRp.begin(), errRp.end());
        if (curr_max < best.best_metric) {
            best.update(curr_max, std::numeric_limits<double>::infinity(),
                        res.base.iterations, X);
        } else if (st.accelerate && K > 0 &&
                   curr_max > best.best_metric * (1.0 + lbw::kSRAAOuterSlack)) {
            if (++outer_stall_count >= lbw::kSRAAOuterStallWindow) {
                X = best.best_weights;        // revert to outer-quality best
                grk_sraa.clear();             // restart AA history
                outer_stall_count = 0;
                // Rebuild W, S_flat, errRp from reverted X
                W = 0.0; std::fill(S_flat.begin(), S_flat.end(), 0.0);
                for (int c = 0; c < M; c++) {
                    W += X[c];
                    for (int k = 0; k < K; k++) {
                        int g = ct.g_per_cell[k][c];
                        if (g >= 0 && g < st.cat_counts[k])
                            S_flat[k * S_stride + g] += X[c];
                    }
                }
                for (int k = 0; k < K; k++) errRp[k] = compute_errRp_k(k);
                // After revert to best, errRp matches best.best_metric — no update needed.
            }
        } else { outer_stall_count = 0; }

        // Convergence check every kErrCheckInterval iters
        if ((iter + 1) % kErrCheckInterval == 0 || iter == st.inner_max_iter - 1) {
            lbw::CellMetrics m = lbw::compute_cell_metrics(st, ct, X, W, bucket_scratch);
            double current_errRp = *std::max_element(errRp.begin(), errRp.end());
            if (first_errRp < 0.0) first_errRp = current_errRp;  // B5: capture at first check
            bool converged = lbw::check_convergence(cfg, m, prev_metric, st.tol_abs);
            if (converged) {
                lbw::mark_converged(res, cfg, res.base.iterations);
                // B7: only overwrite best if convergence X is actually better
                if (current_errRp < best.best_metric) {
                    best.update(current_errRp, std::numeric_limits<double>::infinity(),
                                res.base.iterations, X);
                }
                break;
            }
        }
    }

    // G8c: write best-iterate fields from tracker.
    res.base.convergence_solver_objective = best.best_objective;  // ∞ for greenkhorn
    res.base.best_error = best.best_metric;
    res.base.best_iter  = best.best_iter;

    // Post-loop status
    if (res.base.status == RK_ERR_NOCONV) {
        double final_errRp = *std::max_element(errRp.begin(), errRp.end());
        // B5: compare against first_errRp (initial error), not prev_metric (last update = final error)
        res.base.status = (first_errRp > 0.0 && final_errRp < first_errRp * 0.999)
            ? RK_ERR_BUDGET : RK_ERR_STALL;
        if (res.base.status == RK_ERR_STALL) res.base.stall_kind = 2;  // KL plateau (Sinkhorn-class)
        std::snprintf(res.message, sizeof(res.message),
            "greenkhorn: %s after %d steps; best max_err=%.4e",
            res.base.status == RK_ERR_BUDGET ? "budget exhausted" : "stall",
            res.base.iterations, res.base.best_error);
    }

    // Weight reconstruction from best.best_weights (cell-level X snapshot at best iter)
    res.base.best_weights.resize(st.n);
    for (int i = 0; i < st.n; i++) {
        int c = ct.cell_of[i];
        res.base.best_weights[i] = (X_init[c] > 0.0)
            ? st.weights[i] * best.best_weights[c] / X_init[c]
            : st.weights[i];
    }
    res.base.max_error = best.best_metric;

    return res;
}

} // namespace lbw
