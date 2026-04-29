#include "lbw_config.h"
#include "greenkhorn.hpp"
#include "calib_dispatch.hpp"
#include "cell_table.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <limits>
#include <numeric>

namespace lbw {

GreenkornResult greenkhorn_solve(CalibState& st) {
    GreenkornResult res;
    res.status = RK_ERR_NOCONV;

    // Build cell table
    CellTable ct;
    {
        std::vector<const int32_t*> gids(st.K);
        for (int k = 0; k < st.K; k++) gids[k] = st.group_ids[k];
        if (build_cell_table(st.n, st.K, gids.data(), st.cat_counts, st.weights, ct) != RK_OK) {
            res.status = RK_ERR_INFEAS;
            std::snprintf(res.message, sizeof(res.message),
                "greenkhorn: cell table build failed");
            return res;
        }
    }

    const int M = ct.M_cell;
    const int K = st.K;
    int max_cats = 0;
    for (int k = 0; k < K; k++) max_cats = std::max(max_cats, st.cat_counts[k]);

    // Build cells_per_cat[k][j] = list of cell indices in bucket (k,j)
    std::vector<std::vector<std::vector<int>>> cells_per_cat(K);
    for (int k = 0; k < K; k++) {
        cells_per_cat[k].assign(st.cat_counts[k], {});
        for (int c = 0; c < M; c++) {
            int g = ct.g_per_cell[k][c];
            if (g >= 0 && g < st.cat_counts[k]) cells_per_cat[k][g].push_back(c);
        }
    }

    // Initialize cell masses from design weights
    std::vector<double> X(M, 0.0);
    for (int i = 0; i < st.n; i++) X[ct.cell_of[i]] += st.weights[i];
    const std::vector<double> X_init = X;

    // Capacity bounds per cell
    std::vector<double> L_cell(M), U_cell(M);
    for (int c = 0; c < M; c++) {
        L_cell[c] = st.min_weight * ct.n_per_cell[c];
        U_cell[c] = st.max_weight * ct.n_per_cell[c];
        X[c] = std::clamp(X[c], L_cell[c], U_cell[c]);
    }

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

    // Best-iterate tracking
    double best_errRp = *std::max_element(errRp.begin(), errRp.end());
    res.best_error = best_errRp;
    res.best_iter  = 0;
    std::vector<double> X_best = X;

    // Convergence state
    std::vector<double> bucket_scratch(max_cats, 0.0);
    double prev_metric = std::numeric_limits<double>::infinity();
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

    // F_eval: K greedy steps, sorted ONCE at entry (stationary for SQUAREM)
    auto F_eval = [&](std::vector<double>& X_in,
                      std::vector<double>& S_in,
                      double& W_in) {
        // Swap global state to X_in/S_in/W_in temporarily
        std::swap(X, X_in); std::swap(S_flat, S_in); std::swap(W, W_in);
        // Sort margins by errRp at round entry (fixed order for this F_eval call)
        std::vector<int> order(K);
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(),
            [&](int a, int b) { return errRp[a] > errRp[b]; });
        for (int ki = 0; ki < K; ki++) greenkhorn_step(order[ki]);
        // Swap back
        std::swap(X, X_in); std::swap(S_flat, S_in); std::swap(W, W_in);
    };

    // SQUAREM scratch buffers
    std::vector<double> sq_w1, sq_w2, sq_X_snap, sq_S1, sq_S2, sq_Ssnap;
    double sq_W1 = 0.0, sq_W2 = 0.0, sq_Wsnap = 0.0;
    if (st.accelerate) {
        sq_w1.assign(M, 0.0); sq_w2.assign(M, 0.0); sq_X_snap.assign(M, 0.0);
        sq_S1.assign(K*S_stride, 0.0); sq_S2.assign(K*S_stride, 0.0);
        sq_Ssnap.assign(K*S_stride, 0.0);
    }

    // Main loop
    for (int iter = 0; iter < st.inner_max_iter; iter++) {
        // W<=0 guard
        if (W <= 0.0) {
            res.status = RK_ERR_INFEAS;
            std::snprintf(res.message, sizeof(res.message),
                "greenkhorn: total mass W<=0 (all cells at zero bound)");
            break;
        }

        if (st.accelerate && K > 0) {
            // SQUAREM: one super-step = two F_eval calls + extrapolation
            sq_X_snap = X; sq_Ssnap = S_flat; sq_Wsnap = W;
            sq_w1 = X;  sq_S1 = S_flat; sq_W1 = W;
            F_eval(sq_w1, sq_S1, sq_W1);
            sq_w2 = sq_w1; sq_S2 = sq_S1; sq_W2 = sq_W1;
            F_eval(sq_w2, sq_S2, sq_W2);

            // CBB alpha (obs-level norms)
            double r_sq = 0.0, v_sq = 0.0;
            for (int c = 0; c < M; c++) {
                double inv_nc = (ct.n_per_cell[c] > 0) ? 1.0/ct.n_per_cell[c] : 0.0;
                double r_c = sq_w1[c] - sq_X_snap[c];
                double v_c = sq_w2[c] - 2.0*sq_w1[c] + sq_X_snap[c];
                r_sq += r_c * r_c * inv_nc;
                v_sq += v_c * v_c * inv_nc;
            }
            double alpha = (v_sq > 0.0) ? -std::sqrt(r_sq / v_sq) : -1.0;
            alpha = std::max(alpha, -4.0);

            // Extrapolate X_star and clamp
            std::vector<double> sq_X_star(M);
            for (int c = 0; c < M; c++) {
                double r_c = sq_w1[c] - sq_X_snap[c];
                double v_c = sq_w2[c] - 2.0*sq_w1[c] + sq_X_snap[c];
                sq_X_star[c] = std::clamp(sq_X_snap[c] - 2.0*alpha*r_c + alpha*alpha*v_c,
                                           L_cell[c], U_cell[c]);
            }
            // Rebuild S_star
            std::vector<double> sq_Sstar(K*S_stride, 0.0);
            double sq_Wstar = 0.0;
            for (int c = 0; c < M; c++) sq_Wstar += sq_X_star[c];
            for (int k = 0; k < K; k++)
                for (int j = 0; j < st.cat_counts[k]; j++)
                    for (int c : cells_per_cat[k][j]) sq_Sstar[k*S_stride+j] += sq_X_star[c];
            // Stabilization step
            F_eval(sq_X_star, sq_Sstar, sq_Wstar);

            // Compare: accept X_star or fall back to w2
            double err_star = 0.0;
            for (int k = 0; k < K; k++)
                for (int j = 0; j < st.cat_counts[k]; j++) {
                    double ach = (sq_Wstar > 0.0) ? sq_Sstar[k*S_stride+j]/sq_Wstar : 0.0;
                    err_star = std::max(err_star, std::abs(ach - st.targets[k][j]));
                }
            double err_w2 = *std::max_element(errRp.begin(), errRp.end());
            if (err_star <= err_w2 * 1.01) {
                X = sq_X_star; S_flat = sq_Sstar; W = sq_Wstar;
            } else {
                X = sq_w2; S_flat = sq_S2; W = sq_W2;
            }
            for (int k = 0; k < K; k++) errRp[k] = compute_errRp_k(k);
            res.iterations += K * 3;  // 3 F_eval calls per super-step
        } else {
            // Pure Greenkhorn: single margin per step
            int k_star = (int)(std::max_element(errRp.begin(), errRp.end()) - errRp.begin());
            greenkhorn_step(k_star);
            res.iterations = iter + 1;
        }

        // Best-iterate
        double curr_max = *std::max_element(errRp.begin(), errRp.end());
        if (curr_max < best_errRp) {
            best_errRp = curr_max;
            res.best_error = best_errRp;
            res.best_iter  = res.iterations;
            X_best = X;
        }

        // Convergence check every kErrCheckInterval iters
        if ((iter + 1) % kErrCheckInterval == 0 || iter == st.inner_max_iter - 1) {
            lbw::CellMetrics m = lbw::compute_cell_metrics(st, ct, X, W, bucket_scratch);
            bool converged = lbw::check_convergence(cfg, m, prev_metric, st.tol_abs);
            if (converged) {
                res.status = RK_OK;
                res.convergence_iter = res.iterations;
                X_best = X;
                break;
            }
        }
    }

    // Post-loop status
    if (res.status == RK_ERR_NOCONV) {
        double final_errRp = *std::max_element(errRp.begin(), errRp.end());
        res.status = (final_errRp < prev_metric * 0.999) ? RK_ERR_BUDGET : RK_ERR_STALL;
        std::snprintf(res.message, sizeof(res.message),
            "greenkhorn: %s after %d steps; best max_err=%.4e",
            res.status == RK_ERR_BUDGET ? "budget exhausted" : "stall",
            res.iterations, res.best_error);
    }

    // Weight reconstruction from X_best
    res.best_weights.resize(st.n);
    for (int i = 0; i < st.n; i++) {
        int c = ct.cell_of[i];
        res.best_weights[i] = (X_init[c] > 0.0)
            ? st.weights[i] * X_best[c] / X_init[c]
            : st.weights[i];
    }
    res.max_error = best_errRp;

    return res;
}

} // namespace lbw
