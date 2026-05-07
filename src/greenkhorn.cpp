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

// ── T-H: Two-stage hierarchical Greenkhorn ───────────────────────────────────
//
// Glue layer: dispatches on p->hierarchical_enabled × hierarchical_mode.
// DOES NOT modify Greenkhorn math — all greedy IPF logic lives in
// greenkhorn_solve().
//
// Within-cell lambda signature: fn(weights, cell_id)
//   cell_id == -1 → Stage-1: full-data Greenkhorn on coarse margins only (spec §6).
//   cell_id >= 0  → Stage-2: cell-restricted Greenkhorn on fine margins only.
//
// Queue isolation: greenkhorn_solve() reconstructs X, W, S_flat, errRp, and
// cells_per_cat entirely from CalibState on every call. There is no global or
// static queue state. Within-cell calls receive a sub-CalibState with
// n_cell < N and K_fine < K, so the priority queue is seeded on cell-local
// residuals — not the full-data residuals from Stage-1 or prior Stage-2 calls.
//
// Lambda [&] guard: bool guards declared BEFORE auto fn = [&] per CLAUDE.md.
GreenkornResult greenkhorn_solve_hierarchical(CalibState& st, const rk_params_t* p) {
    // Early-out: no hierarchical params or disabled — single-stage passthrough.
    if (!p || p->hierarchical_enabled == 0 || !p->hierarchical_coarse_mask)
        return greenkhorn_solve(st);

    const int N = st.n;
    const int K = st.K;
    const int mode = p->hierarchical_mode;   // 0=refine, 1=exact
    const int  min_cell_n     = p->hierarchical_min_cell_n > 0 ? p->hierarchical_min_cell_n : 1;
    const double outer_tol    = p->hierarchical_outer_tol;
    const int  outer_iters    = p->hierarchical_outer_iterations > 0
                                    ? p->hierarchical_outer_iterations : 50;

    // Build coarse-cell partition (used for Stage-1 multiplier bookkeeping).
    lbw::CellPartition partition;
    {
        int rc = lbw::build_cell_partition(N, K, st.group_ids, st.cat_counts,
                                           p->hierarchical_coarse_mask,
                                           min_cell_n, partition);
        if (rc != RK_OK) {
            // Degenerate case (e.g. min_cell_n > N → cap = 0).
            // Validator already emitted Rf_warning; fall back to single-stage.
            return greenkhorn_solve(st);
        }
    }
    lbw::SparseMask sparse = lbw::build_sparse_mask(partition, min_cell_n);

    // Identify fine margins (coarse_mask[k]==0) and coarse margins (coarse_mask[k]==1).
    std::vector<int> fine_idx;          // indices into [0..K) of fine margins
    std::vector<int> coarse_idx;        // indices into [0..K) of coarse margins
    fine_idx.reserve(K);
    coarse_idx.reserve(K);
    for (int k = 0; k < K; k++) {
        if (p->hierarchical_coarse_mask[k]) coarse_idx.push_back(k);
        else                                fine_idx.push_back(k);
    }
    const int K_fine   = static_cast<int>(fine_idx.size());
    const int K_coarse = static_cast<int>(coarse_idx.size());

    // Per-cell fine-margin group_ids and targets pointers (cell-restricted copies).
    // Allocated once and reused per lambda call via shared buffers.
    // Buffer: fine group_ids per cell observation (size N × K_fine max).
    std::vector<std::vector<int32_t>> fine_gids_buf(K_fine, std::vector<int32_t>(N));
    std::vector<double>               fine_weights_buf(N);
    std::vector<const int32_t*>       fine_gids_ptrs(K_fine);
    std::vector<const double*>        fine_tgt_ptrs(K_fine);
    for (int fi = 0; fi < K_fine; fi++)
        fine_tgt_ptrs[fi] = st.targets[fine_idx[fi]];

    // Coarse-only sub-CalibState buffers for Stage-1 (T-M / spec §6).
    // Group_ids pointers alias the original arrays; no buffer copy required at Stage-1.
    // Stage-1 write-back: st.weights redirected to working_weights before the call;
    // greenkhorn_solve writes best_weights there directly — no separate staging buffer needed.
    std::vector<const int32_t*>       coarse_gids_ptrs(K_coarse);
    std::vector<const double*>        coarse_tgt_ptrs(K_coarse);
    std::vector<int>                  coarse_cats(K_coarse);
    for (int ci = 0; ci < K_coarse; ci++) {
        coarse_gids_ptrs[ci] = st.group_ids[coarse_idx[ci]];
        coarse_tgt_ptrs[ci]  = st.targets[coarse_idx[ci]];
        coarse_cats[ci]      = st.cat_counts[coarse_idx[ci]];
    }

    // stage1_multipliers: per-cell w_after/w_before ratio (length n_cells_total).
    // Computed after Stage-1 completes; needed by apply_sparse_inheritance.
    // Spec §6: cell-aggregate multiplier = Σ(weights[i in c]) / Σ(w_init[i in c]),
    // NOT a per-obs ratio (which would non-deterministically depend on obs order).
    std::vector<double> stage1_multipliers(partition.n_cells_total, 1.0);

    // Save original caller-owned weight buffer; redirect st.weights to our working copy.
    double* const original_weights_ptr = st.weights;

    // working_weights: shared working buffer used by both outer_iterate_strategy_a
    // and our fn. st.weights points at this buffer so greenkhorn_solve() writes here.
    std::vector<double> working_weights(st.weights, st.weights + N);
    st.weights = working_weights.data();  // redirect raw pointer

    // w_init: weight snapshot before Stage-1 (for sparse inheritance).
    std::vector<double> w_init(working_weights);

    // ── bool guards declared BEFORE [&] lambda (CLAUDE.md rule) ──
    bool stage1_done = false;

    // Within-cell callable for outer_iterate_strategy_a.
    // weights_ref is the same vector as working_weights (passed by outer loop).
    // cell_id == -1: full-data Stage-1 Greenkhorn on ALL margins.
    // cell_id >= 0:  cell-restricted Stage-2 Greenkhorn on fine margins only.
    //
    // Queue isolation verified: greenkhorn_solve() allocates cells_per_cat,
    // S_flat, errRp, and priority ordering locally inside each call. The sub-st
    // passed for cell_id>=0 has n=n_cell, K=K_fine, so the greedy argmax is
    // computed on cell-local residuals — independent of the full-data state.
    auto fn = [&](std::vector<double>& weights_ref, int cell_id) {
        // st.weights already points at weights_ref.data() (== working_weights.data()).
        (void)weights_ref;
        if (cell_id < 0) {
            // Stage-1: full-data Greenkhorn on coarse margins only (Stage-1 per spec §6).
            // NOTE: greenkhorn_solve() does NOT write back to sub.weights — it returns
            // calibrated obs-level weights in res.base.best_weights. Copy them into
            // working_weights (via st.weights) so subsequent Stage-2 calls see updated
            // starting weights and the outer loop convergence check is meaningful.
            if (K_coarse > 0) {
                CalibState sub  = st;
                sub.K           = K_coarse;
                sub.group_ids   = coarse_gids_ptrs.data();
                sub.targets     = coarse_tgt_ptrs.data();
                sub.cat_counts  = coarse_cats.data();
                auto res1 = greenkhorn_solve(sub);
                if (!res1.base.best_weights.empty() &&
                    static_cast<int>(res1.base.best_weights.size()) == N) {
                    std::copy(res1.base.best_weights.begin(),
                              res1.base.best_weights.end(), st.weights);
                }
            }
            // Record per-cell Stage-1 multipliers for sparse inheritance.
            // Spec §6: cell-aggregate ratio Σw / Σw_init (NOT per-obs overwrite).
            std::vector<double> num(partition.n_cells_total, 0.0);
            std::vector<double> den(partition.n_cells_total, 0.0);
            for (int i = 0; i < N; i++) {
                int c = partition.cell_id_per_obs[i];
                num[c] += st.weights[i];
                den[c] += w_init[i];
            }
            for (int c = 0; c < partition.n_cells_total; c++)
                stage1_multipliers[c] = (den[c] > 0.0) ? num[c] / den[c] : 1.0;
            stage1_done = true;
        } else {
            // Stage-2: cell-restricted Greenkhorn on fine margins only.
            if (K_fine == 0) return;
            const std::vector<int>& obs = partition.obs_indices_by_cell[cell_id];
            const int n_cell = static_cast<int>(obs.size());
            if (n_cell == 0) return;

            // Pack cell weights and fine group_ids into contiguous buffers.
            for (int ii = 0; ii < n_cell; ii++) {
                fine_weights_buf[ii] = st.weights[obs[ii]];
                for (int fi = 0; fi < K_fine; fi++)
                    fine_gids_buf[fi][ii] = st.group_ids[fine_idx[fi]][obs[ii]];
            }
            for (int fi = 0; fi < K_fine; fi++)
                fine_gids_ptrs[fi] = fine_gids_buf[fi].data();

            // Sub-CalibState: inherits all solver config, restricted to cell obs + fine margins.
            // Queue isolation: greenkhorn_solve() rebuilds all internal state from sub.
            // The sub has n=n_cell (not N) and K=K_fine (not K), so cells_per_cat,
            // S_flat, and errRp are constructed exclusively from cell-local data.
            CalibState sub   = st;
            sub.n            = n_cell;
            sub.K            = K_fine;
            sub.weights      = fine_weights_buf.data();
            sub.group_ids    = fine_gids_ptrs.data();
            sub.targets      = fine_tgt_ptrs.data();

            // cat_counts for fine margins.
            std::vector<int> fine_cats(K_fine);
            for (int fi = 0; fi < K_fine; fi++)
                fine_cats[fi] = st.cat_counts[fine_idx[fi]];
            sub.cat_counts = fine_cats.data();

            // NOTE: greenkhorn_solve() does NOT write back to sub.weights — calibrated
            // weights are in res.base.best_weights. Copy them into fine_weights_buf
            // so we can scatter them back to the shared working buffer.
            auto sub_res = greenkhorn_solve(sub);
            if (!sub_res.base.best_weights.empty() &&
                static_cast<int>(sub_res.base.best_weights.size()) == n_cell) {
                for (int ii = 0; ii < n_cell; ii++)
                    fine_weights_buf[ii] = sub_res.base.best_weights[ii];
            }
            // Write calibrated weights back to shared working buffer.
            for (int ii = 0; ii < n_cell; ii++)
                st.weights[obs[ii]] = fine_weights_buf[ii];
        }
    };

    GreenkornResult res;

    if (mode == 0) {
        // Strategy A (refine): outer_iterate_strategy_a alternates Stage-1 / Stage-2.
        // working_weights is both passed to fn and holds the live weights.
        lbw::OuterResult out = lbw::outer_iterate_strategy_a(
            fn, working_weights,
            partition, sparse, w_init, stage1_multipliers,
            st.convergence_cfg, outer_tol, outer_iters, N);

        res.base.status        = out.status;
        res.base.iterations    = out.iterations_used;
        res.base.max_error     = out.residual_final;
        res.hier_n_cells_total         = partition.n_cells_total;
        res.hier_n_cells_skipped       = partition.n_cells_skipped;
        res.hier_outer_iterations_used = out.iterations_used;
        res.hier_outer_residual_final  = out.residual_final;
        res.hier_levels_used           = 1;
        // Budget-exit: last_iterate_weights is in working_weights (fn writes there).
    } else {
        // Strategy B (exact): orthogonality already validated at entry.
        // Stage-1: full-data Greenkhorn on coarse margins only (spec §6).
        // greenkhorn_solve() returns calibrated obs-level weights in best_weights,
        // not in st.weights. Copy best_weights into working_weights for continuity.
        GreenkornResult s1{};
        if (K_coarse > 0) {
            CalibState sub  = st;
            sub.K           = K_coarse;
            sub.group_ids   = coarse_gids_ptrs.data();
            sub.targets     = coarse_tgt_ptrs.data();
            sub.cat_counts  = coarse_cats.data();
            s1 = greenkhorn_solve(sub);
            if (!s1.base.best_weights.empty() &&
                static_cast<int>(s1.base.best_weights.size()) == N) {
                std::copy(s1.base.best_weights.begin(),
                          s1.base.best_weights.end(), working_weights.begin());
            }
        }
        // Cell-aggregate Stage-1 multiplier (spec §6).
        {
            std::vector<double> num(partition.n_cells_total, 0.0);
            std::vector<double> den(partition.n_cells_total, 0.0);
            for (int i = 0; i < N; i++) {
                int c = partition.cell_id_per_obs[i];
                num[c] += working_weights[i];
                den[c] += w_init[i];
            }
            for (int c = 0; c < partition.n_cells_total; c++)
                stage1_multipliers[c] = (den[c] > 0.0) ? num[c] / den[c] : 1.0;
        }
        // Stage-2: one pass per non-sparse cell on fine margins only.
        for (int c = 0; c < partition.n_cells_total; c++) {
            if (!sparse.is_sparse_cell[c])
                fn(working_weights, c);
        }
        // Sparse cells inherit the cell-aggregate Stage-1 multiplier.
        lbw::apply_sparse_inheritance(working_weights, w_init, partition, sparse,
                                      stage1_multipliers);
        lbw::enforce_sigmaw_eq_n(working_weights, N);

        res.base.status        = s1.base.status;
        res.base.iterations    = s1.base.iterations;
        res.base.max_error     = s1.base.max_error;
        res.hier_n_cells_total         = partition.n_cells_total;
        res.hier_n_cells_skipped       = partition.n_cells_skipped;
        res.hier_outer_iterations_used = 1;
        res.hier_outer_residual_final  = 0.0;
        res.hier_levels_used           = 1;
    }

    // Write final weights back to caller's buffer and restore st.weights.
    std::copy(working_weights.begin(), working_weights.end(), original_weights_ptr);
    st.weights = original_weights_ptr;
    // Populate res.base.best_weights so callers (r_bridge, c_api) that read
    // res.base.best_weights get the final calibrated obs-level weights.
    // (greenkhorn_solve() stores calibrated weights in best_weights, not in
    //  st.weights; the hierarchical wrapper accumulates them in working_weights.)
    res.base.best_weights.assign(working_weights.begin(), working_weights.end());
    return res;
}

} // namespace lbw
