// Cold post-loop finalization for ORIS, split out of oris.cpp (uu8r.2).
// Called once per oris_solve() after the homotopy level loop exits: obs
// expansion, bounds enforcement, ALM projection, best-iterate fallback, and
// final diagnostics. Not on the hot path. Body is byte-identical to the former
// file-static definition in oris.cpp.

#include "oris_internal.hpp"
#include "leafblower.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <set>
#include <utility>
#include <vector>

namespace lbw {

void oris_finalize(
    CalibState&                               st,
    ORISResult&                               res,
    const CellTable&                          ct,
    std::vector<double>&                      X,
    const std::vector<double>&                X_init,
    const std::vector<double>&                L_cell,
    const std::vector<double>&                U_cell,
    bool                                      alm_active,
    double                                    capacity_mu_adaptive,
    const std::vector<double>&                lambda_cell,
    const BestIterTracker&                    best,
    bool                                      absolute_tol_fired,
    const std::set<std::pair<int,int>>&       structural_infeas_pairs,
    double                                    sor_min_omega,
    int                                       sor_n_damped,
    double                                    sor_omega_mean,
    int                                       sor_any_latched,
    int                                       sor_n_pinned_fb,
    int                                       sor_n_warmup_fb,
    int                                       sor_n_conv_fb,
    int                                       sor_n_resid_grew,
    int                                       sor_n_monotone_cd,
    const std::vector<std::pair<int,double>>& probe_samples
) {
    constexpr double kInfeasStallRatio = 10.0;

    // ════════════════════ Post-loop: ALM projection + obs expansion + bounds ════════════════════
    // ALM final projection: hard-clamp X[c] into [L_cell, U_cell], then a small
    // bounded rescale loop to recover sum(X)=n while staying within bounds.
    if (alm_active) {
        constexpr int    kMaxRescaleIters = 3;
        constexpr double kRescaleTol      = 1e-12;
        for (int c = 0; c < ct.M_cell; c++) X[c] = std::clamp(X[c], L_cell[c], U_cell[c]);
        double prev_total = 0.0;
        for (int iter = 0; iter < kMaxRescaleIters; iter++) {
            double total = 0.0;
            for (int c = 0; c < ct.M_cell; c++) total += X[c];
            if (std::abs(total - prev_total) < kRescaleTol * st.n || total <= 0.0) break;
            prev_total = total;
            const double rescale = static_cast<double>(st.n) / total;
            for (int c = 0; c < ct.M_cell; c++) X[c] = std::clamp(X[c] * rescale, L_cell[c], U_cell[c]);
        }
        double final_total = 0.0;
        for (int c = 0; c < ct.M_cell; c++) final_total += X[c];
        res.alm_sum_drift = std::abs(final_total - static_cast<double>(st.n));
        // ohi0: force cell-aggregate Σ_c X = n. The bounded rescale loop above
        // cannot reach it when cells are mass-saturated (infeasible), leaving a
        // residual. Apply one unclamped scale AFTER recording alm_sum_drift so the
        // infeasibility signal is preserved. No clamp: on infeasible inputs X may
        // exceed [L_cell,U_cell] — the honest cell-aggregate outcome. Does not
        // change obs weights: the obs-level Σw=n renorm below absorbs the uniform
        // X scaling.
        if (final_total > 0.0) {
            const double rescale_final = static_cast<double>(st.n) / final_total;
            for (int c = 0; c < ct.M_cell; c++) X[c] *= rescale_final;
        }
        if (res.alm_sum_drift > 1e-6 * st.n && st.verbose >= 1) {
            char msg[256];
            std::snprintf(msg, sizeof(msg), "[oris_soft] final projection sum drift = %.2e", res.alm_sum_drift);
            st.log(msg);
        }

        // Populate ALM diagnostics.
        res.alm_capacity_mu_final = capacity_mu_adaptive;
        double max_dual = 0.0;
        for (int c = 0; c < ct.M_cell; c++)
            max_dual = std::max(max_dual, std::abs(lambda_cell[c]));
        res.alm_max_dual_norm = max_dual;
        // alm_n_growth_events populated incrementally above.
    }

    // Classify RK_ERR_NOCONV → BUDGET or STALL.
    // RK_ERR_BUDGET (4): metric improved at some point → increase max_iterations.
    // RK_ERR_STALL  (5): metric never improved from initial → at constrained optimum.
    if (res.base.status == RK_ERR_NOCONV) {
        res.base.status = best.has_best() ? RK_ERR_BUDGET : RK_ERR_STALL;
        if (res.base.status == RK_ERR_STALL) {
            // leafblower-8eod: stall_kind driven by actual accelerate state at solver exit.
            // This path is NOT bijective with the user's accelerate flag (harvest.R heuristic
            // was wrong): oris emits status==5 for both SRAA (accelerate=true) and plain-BCD
            // (accelerate=false). Set stall_kind from st.accelerate here so harvest.R reads
            // the actual mechanism instead of the user input flag.
            res.base.stall_kind = st.accelerate ? 1 : 2;  // 1=wchange (SRAA), 2=kl (plain-BCD)
        }
    }

    // populate convergence diagnostics at solver exit.
    res.base.convergence_metric             = static_cast<int>(st.convergence_cfg.metric);
    res.base.convergence_rule               = static_cast<int>(st.convergence_cfg.rule);
    res.base.convergence_tol = absolute_tol_fired
        ? st.convergence_cfg.absolute_tol : st.convergence_cfg.pct_tol;
    // za9r: convergence_iter is pinned at the firing site (non-SRAA terminal
    // block / SRAA mark_converged). Do NOT clobber it with res.base.iterations
    // here; only reset to -1 when the solve did not converge.
    if (res.base.status != RK_OK) res.base.convergence_iter = -1;
    res.base.convergence_solver_objective   = best.best_objective;
    res.base.convergence_minimized_metric   = static_cast<int>(st.convergence_cfg.metric);

    // WU-E / G8b: expand best.best_weights (cell-level W ratio snapshot) to obs-level,
    // then finalize through the shared Σw=n + bounds_mode contract (same as the final
    // iterate below). If best.has_best() is false (solver exited before first check),
    // best_weights is all zeros.
    res.base.best_error = best.best_metric;
    res.base.best_iter  = best.best_iter;
    if (best.has_best()) {
        std::vector<double> best_weights_obs(st.n);
        for (int i = 0; i < st.n; i++) {
            best_weights_obs[i] = st.weights[i] * best.best_weights[ct.cell_of[i]];
        }
        // l6to: enforce Σw=n via the shared finalize_weights contract (normalize then
        // bounds_mode dispatch). The previous normalize-then-clamp broke Σw=n on
        // infeasible problems, and a post-clamp renormalize would re-break the bounds.
        // Cell mode: no clamp (cell-aggregate contract). Unit mode: per-cell water-fill
        // preserves Σw=n while enforcing [min_weight, max_weight].
        int best_nbv = 0, best_nbc = 0;
        finalize_weights_buf(best_weights_obs.data(), st.n, st, ct, best_nbv, best_nbc);
        res.base.best_weights = std::move(best_weights_obs);
    } else {
        res.base.best_weights.assign(st.n, 0.0);
    }

    // Expansion to observation weights.
    std::vector<double> mult(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) {
        mult[c] = (X_init[c] > 0.0) ? X[c] / X_init[c] : 0.0;
    }
    for (int i = 0; i < st.n; i++) {
        st.weights[i] = st.weights[i] * mult[ct.cell_of[i]];
    }

    // Solver-owned normalization (moved from wrapper 2026-04-24 per user directive)
    // + bounds_mode dispatch via the shared finalize_weights contract.
    // o2o6.1: the normalize (l7sg subnormal guard), cell-mode diagnostic-only
    // counting, and unit-mode per-cell water-fill core were duplicated here and in
    // calib_dispatch.hpp::finalize_weights_buf. Delegate to the single source.
    // Sanctioned order: normalize (Σw=n, applied AFTER expansion and BEFORE bounds
    // post-processing so unit-mode water-fill sees final-scale weights) → bounds
    // dispatch. Degenerate total_w leaves weights unchanged (status as set upstream).
    // The unit-mode gbib.1 Σ=n renorm + post-renorm clamp below stay inline: they
    // are final-iterate-specific (the best-iterate call at L149 deliberately omits
    // them), so folding them into the shared helper would alter best-iterate output.
    finalize_weights(st, ct, res.n_bounds_violated, res.n_bounds_clamped);

    if (st.bounds_mode != RK_BOUNDS_CELL) {
        // gbib.1: water-fill preserves Σ weights = n in the redistribute branch
        // (factor = 1 + excess/free_sum exactly absorbs `excess`), but the
        // pathological n_free==0 / free_sum<=0 break path (lines ~1930) leaves
        // the cell sum off by `excess` (all violators pinned at bounds, no free
        // obs to absorb). Likewise, kWaterFillMaxIter exhaustion can leave a
        // residual. Restore Σ=n (contract; see test-oris-bounds-mode L143).
        // Bounds may be re-broken by this rescale on STALL/BUDGET paths, but
        // that condition was already pathologically infeasible — Σ=n is the
        // stronger contract.
        double total_after = 0.0;
        for (int i = 0; i < st.n; i++) total_after += st.weights[i];
        if (std::isfinite(total_after) && total_after > 0.0) {
            const double renorm = static_cast<double>(st.n) / total_after;
            if (std::fabs(renorm - 1.0) > 1e-15) {
                for (int i = 0; i < st.n; i++) st.weights[i] *= renorm;
            }
        }
        // Clamp post-renorm: renorm can push water-fill-clamped obs above max_weight.
        // Drift invariant: after this clamp, Σweights may be slightly < n in the
        // pathological n_free==0 case (all weights pinned to bounds). The drift is
        // bounded by O(ε × n_pinned × max_weight) where ε = |renorm - 1|, typically
        // < 1e-10 × n. A second renorm after clamping is not applied — it would
        // re-violate bounds for the pinned cells, creating an oscillating fix-up loop.
        // Callers must tolerate Σ ≈ n (within double-precision rounding) in this path.
        for (int i = 0; i < st.n; i++) {
            st.weights[i] = std::max(st.min_weight, std::min(st.max_weight, st.weights[i]));
        }
    }

    // CR-A6 (mxcl.6): recompute the reported error on the FINAL, capacity-clamped
    // weights on the SRAA path. The SRAA sweep is unconstrained (capacity is
    // enforced only by the level-exit mass-preserving clamp), so res.base.max_error
    // was set to the per-step UNCONSTRAINED errRp (oris.cpp:877) — on a tight/
    // infeasible problem it reads ~2e-16 ("converged") while the bounds-clamped
    // returned weights miss margins badly. Mirror the non-SRAA path, which already
    // reports the constrained per-iter errRp (oris.cpp:1863). Aggregate the returned
    // obs weights back to cell masses so the recomputed errRp reflects exactly the
    // returned margins (bounds_mode-agnostic). Guarded to accelerate=TRUE so the
    // non-SRAA path stays byte-identical.
    if (st.accelerate) {
        std::vector<double> X_final(ct.M_cell, 0.0);
        for (int i = 0; i < st.n; i++) X_final[ct.cell_of[i]] += st.weights[i];
        double W_total = 0.0;
        for (int c = 0; c < ct.M_cell; c++) W_total += X_final[c];
        if (W_total > 0.0) {
            std::vector<double> S_scratch(lbw::max_cats_count(st.K, st.cat_counts));
            const lbw::CellMetrics cm =
                lbw::compute_cell_metrics(st, ct, X_final, W_total, S_scratch);
            res.base.max_error = cm.errRp;
            // CR-B13 (y2ks.13): best_error had the same honesty gap max_error had
            // before mxcl.6 — L122 copies best.best_metric, tracked on the PRE-clamp
            // accepted masses, so on a tight/infeasible BUDGET exit it reads ~2e-16
            // ("fixed point") while the returned bounds-clamped weights miss margins.
            // harvest.R:728,740 surfaces best_error verbatim in the BUDGET warning
            // LABELLED with convergence_used$metric, and harvest.R:754-757 divides it
            // by metric_prev_check (oris.cpp:1176, tracked as select_metric(metric,cm)).
            // Report the metric-FAMILY value of the returned solution — NOT cm.errRp,
            // which would mislabel units under the default metric (marginal_kl) and
            // corrupt the extrapolation ratio. Reuses the SAME cm (no second
            // aggregation); reduces to cm.errRp when metric=max_err.
            res.base.best_error = lbw::select_metric(st.convergence_cfg.metric, cm);
        } else {
            // Degenerate/NaN return (Σw ≤ 0): every achieved proportion is 0, so the
            // honest max margin error is the largest target. Do NOT keep the stale
            // unconstrained errRp (it would falsely read "converged").
            double max_t = 0.0;
            for (int k = 0; k < st.K; k++)
                for (int j = 0; j < st.cat_counts[k]; j++)
                    max_t = std::max(max_t, st.targets[k][j]);
            res.base.max_error = max_t;
            res.base.best_error = max_t;  // CR-B13: keep best_error honest too
        }
    }

    // PCT stall detection: pct_change < pct_tol (PCT converged) but max_error >> pct_tol
    // signals infeasible problem. Threshold 10x: well-posed problems have
    // errRp/pct_change ratio 1-5x; infeasible stalls show 100x+; 10x separates them.
    // Warning only — status unchanged for backward compatibility. Placed AFTER the
    // SRAA max_error recompute above so it sees the honest constrained error (the
    // stale unconstrained ~2e-16 would suppress the warning on infeasible SRAA runs).
    {
        const auto& cfg = st.convergence_cfg;
        if (res.base.status != RK_OK &&
            (cfg.metric == CalibMetric::MAX_ERR || cfg.metric == CalibMetric::MEAN_ERR) &&
            cfg.pct_tol > 0.0 &&
            res.base.max_error > kInfeasStallRatio * cfg.pct_tol &&
            st.log_fn != nullptr) {
            char msg[256];
            std::snprintf(msg, sizeof(msg),
                "PCT convergence stall: pct_change < %.3g but max_error=%.3g "
                "(%.0fx pct_tol). Possible contradictory or infeasible targets.",
                cfg.pct_tol, res.base.max_error,
                res.base.max_error / cfg.pct_tol);
            st.log_fn(msg, st.log_ctx);
        }
    }

    // B9: structural infeasibility must override RK_OK as well. SRAA convergence
    // (line ~1684) and the post-loop classifier (NOCONV→BUDGET/STALL above)
    // can each leave a "success" status set even though structural_infeas_pairs
    // is non-empty — those pairs are buckets that can never be satisfied.
    if (!structural_infeas_pairs.empty() &&
        (res.base.status == RK_OK ||
         res.base.status == RK_ERR_NOCONV ||
         res.base.status == RK_ERR_BUDGET ||
         res.base.status == RK_ERR_STALL)) {
        res.base.status = RK_ERR_INFEAS;
    }

    // kxna.20 (CR-C7c): re-gate RK_OK on the post-finalize RETURNED margin error in
    // bounds_mode="unit". The per-obs water-fill (finalize_weights at L161) can leave a
    // fully-pinned cell off-target after status was set RK_OK on the pre-finalize cell
    // iterate. Structural infeasibility already promoted to INFEAS above (wins); here a
    // remaining RK_OK with a drifted returned margin is demoted to STALL. Recompute the
    // returned errRp LOCALLY (both SRAA and non-SRAA unit paths) so the reported
    // res.base.max_error/best_error — set to metric-family values on the SRAA path
    // (L212/224) — are left untouched. compute_cell_metrics is a cold once-per-solve pass.
    if (st.bounds_mode == RK_BOUNDS_UNIT) {
        std::vector<double> X_ret(ct.M_cell, 0.0);
        for (int i = 0; i < st.n; i++) X_ret[ct.cell_of[i]] += st.weights[i];
        double W_ret = 0.0;
        for (int c = 0; c < ct.M_cell; c++) W_ret += X_ret[c];
        if (W_ret > 0.0) {
            std::vector<double> S_scratch(lbw::max_cats_count(st.K, st.cat_counts));
            const lbw::CellMetrics cmr =
                lbw::compute_cell_metrics(st, ct, X_ret, W_ret, S_scratch);
            lbw::regate_unit_status(res.base, st, cmr.errRp);
        }
    }

    if (st.verbose >= 1) {
        const char* status_label =
            (res.base.status == RK_OK) ? "converged" :
            (res.base.status == RK_ERR_NOCONV) ? "max_iter exhausted (NOCONV)" :
            (res.base.status == RK_ERR_INFEAS) ? "infeasible" : "error";
        char msg[256];
        std::snprintf(msg, sizeof(msg),
                      "ORIS %s in %d iters, errRp=%.3e",
                      status_label, res.base.iterations, res.base.max_error);
        st.log(msg);
    }

    if (st.verbose >= 1 && !structural_infeas_pairs.empty()) {
        char msg[256];
        size_t off = 0;
        off += std::snprintf(msg + off, sizeof(msg) - off,
                             "ORIS persistent infeasible cells: ");
        size_t idx = 0;
        const size_t total = structural_infeas_pairs.size();
        for (auto it = structural_infeas_pairs.begin();
             it != structural_infeas_pairs.end() && off < sizeof(msg) - 32;
             ++it, ++idx) {
            off += std::snprintf(msg + off, sizeof(msg) - off,
                                 "margin=%d cat=%d%s",
                                 it->first + 1,
                                 it->second + 1,
                                 (idx + 1 < total) ? ", " : "");
        }
        st.log(msg);
    }

    // store SOR diagnostics in result.
    res.sor_min_omega    = sor_min_omega;
    res.sor_n_damped     = sor_n_damped;
    res.sor_omega_mean   = sor_omega_mean;
    res.sor_any_latched  = sor_any_latched;
    res.sor_n_pinned_fb  = sor_n_pinned_fb;
    res.sor_n_warmup_fb  = sor_n_warmup_fb;
    res.sor_n_conv_fb    = sor_n_conv_fb;
    res.sor_n_resid_grew = sor_n_resid_grew;
    res.sor_n_monotone_cd = sor_n_monotone_cd;

    write_trajectory_csv(probe_samples);
}

}  // namespace lbw
