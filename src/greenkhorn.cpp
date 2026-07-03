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

    // CR-A7: greenkhorn drives convergence through compute_cell_metrics +
    // select_metric, which never populates CellMetrics.l1 (the obs-level weight
    // change, written only by raking/sinkhorn). metric="l1_weight" would feed a
    // phantom l1 == 0.0 into check_convergence and spuriously converge on the
    // first check interval. The greedy coordinate-wise loop keeps no per-check
    // prev-weight snapshot, so fall back to MAX_ERR to track the real objective
    // — mirrors chebyshev.cpp.
    if (st.convergence_cfg.metric == lbw::CalibMetric::L1_WEIGHT)
        st.convergence_cfg.metric = lbw::CalibMetric::MAX_ERR;

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

    // CXX.2 (leafblower-5fm8.5): DECOUPLE the per-iteration oscillation monitor
    // (errRp scale) from the REPORTED best-iterate (configured-metric scale).
    //  - `sraa_best`: cheap errRp-scale (max marginal-residual) monitor that drives
    //    the SRAA outer-stall revert. errRp IS the correct fast oscillation signal;
    //    computing full CellMetrics every iteration would be an O(K·M_cell) perf
    //    regression. Updated every iteration; always seeded with the initial iterate
    //    so the revert always has a valid fallback.
    //  - `best`: the REPORTED/SELECTED best-iterate feeding res.base.{best_error,
    //    max_error,best_iter,best_weights}; its scale is select_metric(cfg.metric).
    //
    // Cadence of the `best` update is metric-dependent, BY DESIGN:
    //  * cfg.metric == MAX_ERR (default): select_metric == max(errRp) == the
    //    per-iteration `curr_max` we already maintain cheaply. `best` is updated
    //    EVERY iteration from curr_max (NOT the interval-block fresh m.errRp), so the
    //    best-iterate is chosen on the per-iteration scale, incl. non-interval iters.
    //    NOTE (kxna.4): W is now re-anchored to exact ΣX at each interval check, which
    //    perturbs the subsequent incremental errRp at the ~1e-15 floor — so this path is
    //    NO LONGER byte-identical to the historical incremental-only-W code. The rebuilt
    //    W is the more-accurate value; R and Python rebuild identically, so parity holds.
    //  * cfg.metric != MAX_ERR (kl/chi2/…): the per-iteration errRp is the WRONG
    //    scale, so `best` is updated ONLY in the kErrCheckInterval block where the
    //    full CellMetrics `m` exists, via select_metric(cfg.metric, m) — mirroring
    //    sinkhorn.cpp. This is the bug fix: pre-CXX.2 `best.best_metric` was an errRp
    //    value mislabeled as the requested metric.
    // greenkhorn does not minimise KL directly, so best_objective stays ∞.
    const bool metric_is_max_err = (st.convergence_cfg.metric == CalibMetric::MAX_ERR);
    BestIterTracker best;
    BestIterTracker sraa_best;
    {
        double init_errRp = *std::max_element(errRp.begin(), errRp.end());
        sraa_best.update(init_errRp, std::numeric_limits<double>::infinity(), 0, X);
        // MAX_ERR default: seed `best` identically to the historical code (initial
        // iterate on the errRp==metric scale). Non-MAX_ERR: leave `best` at +inf so
        // the first interval check sets it on the correct metric scale.
        if (metric_is_max_err)
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
            const double f_base = st.targets[k_step][j] * W / S_kj;
            const double f = f_base;
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

        // CXX.2: per-iteration errRp oscillation monitor (sraa_best) → SRAA
        // outer-stall revert. The reported `best` mirrors this ONLY for MAX_ERR
        // (where select_metric == max(errRp) == curr_max); other metrics defer to
        // the interval block below.
        double curr_max = *std::max_element(errRp.begin(), errRp.end());
        if (curr_max < sraa_best.best_metric) {
            sraa_best.update(curr_max, std::numeric_limits<double>::infinity(),
                             res.base.iterations, X);
            // MAX_ERR default: track the per-iteration errRp-best EXACTLY as the
            // historical code did (byte-identical default path, incl. best_iter on
            // a non-interval iteration). select_metric(MAX_ERR)==curr_max, no recompute.
            if (metric_is_max_err)
                best.update(curr_max, std::numeric_limits<double>::infinity(),
                            res.base.iterations, X);
        } else if (st.accelerate && K > 0 &&
                   curr_max > sraa_best.best_metric * (1.0 + lbw::kSRAAOuterSlack)) {
            if (++outer_stall_count >= lbw::kSRAAOuterStallWindow) {
                X = sraa_best.best_weights;   // revert to errRp-best (oscillation monitor)
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
                // After revert to sraa_best, errRp matches sraa_best.best_metric — no update needed.
            }
        } else { outer_stall_count = 0; }

        // Convergence check every kErrCheckInterval iters
        if ((iter + 1) % kErrCheckInterval == 0 || iter == st.inner_max_iter - 1) {
            // CR-C4 (kxna.4): rebuild W from scratch instead of trusting the running
            // W += delta accumulated since iter 0. compute_cell_metrics divides by W;
            // FP drift in the incremental sum would corrupt every derived metric (and
            // the convergence decision). The block already pays O(M) here, the SRAA path
            // rebuilds identically, and re-anchoring W also bounds forward drift.
            W = 0.0;
            for (int c = 0; c < M; c++) W += X[c];
            lbw::CellMetrics m = lbw::compute_cell_metrics(st, ct, X, W, bucket_scratch);
            double current_errRp = *std::max_element(errRp.begin(), errRp.end());
            if (first_errRp < 0.0) first_errRp = current_errRp;  // B5: capture at first check
            // CXX.2: track best-iterate on the configured convergence metric
            // (matches every other solver via select_metric), not the errRp
            // fast proxy. For cfg.metric==MAX_ERR this equals current_errRp, so
            // default-metric behavior is unchanged; KL/CHI2/etc. now report and
            // select the iterate that is best under the requested metric.
            const double curr_metric = lbw::select_metric(cfg.metric, m);
            // CXX.2: for NON-MAX_ERR metrics, update the reported best-iterate on
            // the configured-metric scale at EVERY interval check (mirrors
            // sinkhorn.cpp). This is where the scale-correct `best` is built for
            // kl/chi2/… For MAX_ERR we SKIP it: the per-iteration block already
            // tracks max(errRp) on the metric scale, and the fresh m.errRp recompute
            // here can differ from the incremental errRp at ~1e-15 (kxna.4's W
            // re-anchoring perturbs both at that floor), so best_iter is chosen from
            // the incremental curr_max for consistency, not this interval-block value.
            if (!metric_is_max_err &&
                std::isfinite(curr_metric) && curr_metric < best.best_metric) {
                best.update(curr_metric, std::numeric_limits<double>::infinity(),
                            res.base.iterations, X);
            }
            if (lbw::check_convergence(cfg, m, prev_metric, st.tol_abs)) {
                lbw::mark_converged(res, cfg, res.base.iterations, st.tol_abs);
                break;
            }
        }
    }

    // CXX.2: for non-MAX_ERR metrics `best` is updated only at interval checks
    // (MAX_ERR is pre-seeded + per-iteration, so it is never empty). If a
    // non-default-metric solve broke out before any interval check fired (e.g. the
    // W<=0 guard on iter 0), `best` is still empty. Fall back to the always-seeded
    // errRp monitor so best_weights is never empty (avoids OOB in the reconstruction
    // loop below). best.best_iter < 0 ⇔ never updated.
    if (best.best_iter < 0)
        best = sraa_best;

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
            "greenkhorn: %s after %d steps; best err=%.4e",
            res.base.status == RK_ERR_BUDGET ? "budget exhausted" : "stall",
            res.base.iterations, res.base.best_error);
    }

    // CR-C2 (kxna.2): recompute honest cell metrics on the best-iterate masses so
    // the reported diagnostics reflect the RETURNED weights. Pre-fix
    // res.base.{kl,chi2,mean_error,grake_norm} were never assigned (default 0.0 =
    // phantom "perfect fit") via pack_solver_result. Field order mirrors
    // sinkhorn.cpp. NOTE: res.base.max_error is intentionally left on the
    // best-iterate configured-metric scale (best.best_metric) below — see CXX.2
    // (test-greenkhorn-best-metric.R); CR-C3 proposes changing it to errRp scale
    // but that conflicts with CXX.2 and is deferred pending resolution.
    {
        double W_best = 0.0;
        for (int c = 0; c < ct.M_cell; c++) W_best += best.best_weights[c];
        if (W_best > 0.0) {
            const lbw::CellMetrics mb =
                lbw::compute_cell_metrics(st, ct, best.best_weights, W_best, bucket_scratch);
            res.base.kl               = mb.kl;
            res.base.mean_error       = mb.mean_err;
            res.base.chi2             = mb.chi2;
            res.base.grake_norm       = mb.grake_norm;
            res.base.l1_weight_change = mb.l1;       // 0: greenkhorn keeps no prev-weight snapshot
        }
    }

    // Weight reconstruction from best.best_weights (cell-level X snapshot at best iter)
    res.base.best_weights.resize(st.n);
    for (int i = 0; i < st.n; i++) {
        int c = ct.cell_of[i];
        res.base.best_weights[i] = (X_init[c] > 0.0)
            ? st.weights[i] * best.best_weights[c] / X_init[c]
            : st.weights[i];
    }
    // CR-C1 (kxna.1): enforce the Σw=n + bounds_mode contract via the shared helper
    // (normalize→bounds order), matching oris_finalize. W is maintained purely
    // incrementally (W += delta) and never rescaled at exit, so without this
    // greenkhorn silently violates Σw=n whenever max_weight clamps bind.
    lbw::finalize_weights_buf(res.base.best_weights.data(), st.n, st, ct,
                              res.n_bounds_violated, res.n_bounds_clamped);

    // Unchanged (CXX.2 contract; CR-C3 deferred): best-iterate configured-metric scale.
    res.base.max_error = best.best_metric;

    return res;
}

} // namespace lbw
