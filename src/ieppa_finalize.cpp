// Cold post-loop finalization for iEPPA, split out of ieppa.cpp (uu8r.2).
// Called once per ieppa_solve() after the homotopy level loop exits: obs
// expansion, bounds enforcement, ALM projection, best-iterate fallback, and
// final diagnostics. Not on the hot path. Body is byte-identical to the former
// file-static definition in ieppa.cpp.

#include "ieppa_internal.hpp"
#include "leafblower.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <set>
#include <utility>
#include <vector>

namespace lbw {

void ieppa_finalize(
    CalibState&                               st,
    IEPPAResult&                              res,
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
            std::snprintf(msg, sizeof(msg), "[ieppa_soft] final projection sum drift = %.2e", res.alm_sum_drift);
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
            // was wrong): ieppa emits status==5 for both SRAA (accelerate=true) and plain-BCD
            // (accelerate=false). Set stall_kind from st.accelerate here so harvest.R reads
            // the actual mechanism instead of the user input flag.
            res.base.stall_kind = st.accelerate ? 1 : 2;  // 1=wchange (SRAA), 2=kl (plain-BCD)
        }
    }

    // PCT stall detection: pct_change < pct_tol (PCT converged) but max_error >> pct_tol
    // signals infeasible problem. Threshold 10x: well-posed problems have
    // errRp/pct_change ratio 1-5x; infeasible stalls show 100x+; 10x separates them.
    // Warning only — status unchanged for backward compatibility.
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

    // Solver-owned normalization (moved from wrapper 2026-04-24 per user directive).
    // Applied AFTER expansion and BEFORE bounds_mode post-processing so unit-mode
    // water-fill sees final-scale weights and can strictly enforce
    // [min_weight, max_weight]. Contract: if total_w == 0 (degenerate: all zero
    // X_init[c] or all-zero mult[c]), weights remain unchanged (all zero) and
    // solver status is left as set upstream. X[c] is NOT rescaled — it is dead
    // after expansion (water-fill uses (void)target_sum below and redistributes
    // within each cell, so cell-aggregate scale is irrelevant).
    double total_w = 0.0;
    for (int i = 0; i < st.n; i++) total_w += st.weights[i];
    // l7sg: guard against subnormal total_w (~1e-310), which would overflow
    // norm = n/total_w to +inf. Below kMinSafeTotalWeight, treat as degenerate-zero
    // (weights left unchanged, status as set upstream).
    if (std::isfinite(total_w) && total_w > kMinSafeTotalWeight) {
        const double norm = static_cast<double>(st.n) / total_w;
        for (int i = 0; i < st.n; i++) st.weights[i] *= norm;
    }

    if (st.bounds_mode == RK_BOUNDS_CELL) {
        // Cell mode: count per-obs bound violations for diagnostic only.
        // Do NOT clamp — cell-mode contract is X[c] <= U_cell[c] (aggregate),
        // not w[i] <= max_weight. Clamping individual weights here distorts
        // marginals when d[i] varies within a cell (non-uniform design weights).
        int violations = 0;
        for (int i = 0; i < st.n; i++) {
            if (st.weights[i] > st.max_weight || st.weights[i] < st.min_weight)
                violations++;
        }
        res.n_bounds_violated = violations;
        res.n_bounds_clamped  = 0;
    } else {
        // Unit mode: per-cell water-filling.
        // Water-fill redistributes excess within each cell, preserving the
        // post-normalize cell sum X[c] exactly via three-way classification
        // in the scan (violator / pinned / free). Only strictly-free obs
        // enter free_sum, so factor = 1 + excess/free_sum_really_free
        // distributes excess in full.
        // Build cells_of_obs (list of obs indices per cell) in one pass.
        std::vector<std::vector<int>> cells_of_obs(ct.M_cell);
        for (int i = 0; i < st.n; i++) cells_of_obs[ct.cell_of[i]].push_back(i);

        const int kWaterFillMaxIter = std::max(50, st.K * 10);
        int total_clamped = 0;
        for (int c = 0; c < ct.M_cell; c++) {
            const auto& idxs = cells_of_obs[c];
            if (idxs.empty()) continue;

            for (int it = 0; it < kWaterFillMaxIter; it++) {
                double excess = 0.0;
                double free_sum = 0.0;
                int    n_free = 0;
                bool   any_violation = false;
                // Running counter: increment total_clamped at each normal-path
                // clamp (strict > max / < min). Pinned weights stay at bound
                // exactly, so the next-iter scan (strict-inequality here) and
                // the redistribute guard (strict < max && > min at "free"
                // branch below) both exclude them — no double count.
                // Pathological re-clamps (n_free==0, budget exhausted) are
                // redundant with these increments and do not re-count.
                for (int i : idxs) {
                    if (st.weights[i] > st.max_weight) {
                        excess += st.weights[i] - st.max_weight;
                        st.weights[i] = st.max_weight;
                        any_violation = true;
                        total_clamped++;
                    } else if (st.weights[i] < st.min_weight) {
                        excess -= st.min_weight - st.weights[i];
                        st.weights[i] = st.min_weight;
                        any_violation = true;
                        total_clamped++;
                    } else if (st.weights[i] == st.max_weight || st.weights[i] == st.min_weight) {
                        // Pinned from prior iter (set exactly via direct
                        // assignment above). Excluded from free_sum so
                        // factor = 1 + excess/free_sum_really_free
                        // distributes excess in full; cell-sum conservation
                        // holds exactly. FP equality is safe here because
                        // pinned obs was assigned, not computed. See
                        // leafblower-6s1o for the pre-fix under-distribution.
                    } else {
                        free_sum += st.weights[i];
                        n_free++;
                    }
                }
                if (!any_violation) break;
                if (n_free == 0 || free_sum <= 0.0) {
                    // Pathological: no room to redistribute. All violators
                    // already pinned by the scan above (line 585/589); cell
                    // sum may deviate from target by accumulated excess.
                    break;
                }
                // Redistribute excess proportionally over free observations.
                double factor = 1.0 + excess / free_sum;
                for (int i : idxs) {
                    if (st.weights[i] > st.min_weight && st.weights[i] < st.max_weight) {
                        st.weights[i] *= factor;
                    }
                }
                // Budget-exhaustion case (it == kWaterFillMaxIter - 1): if
                // factor-redistribution newly pushes a free obs above max,
                // the scan in the NEXT iter would clamp it — but there is no
                // next iter. In practice iter count is always enough because
                // water-fill converges geometrically for feasible problems.
            }
        }
        res.n_bounds_clamped = total_clamped;

        // gbib.1: water-fill preserves Σ weights = n in the redistribute branch
        // (factor = 1 + excess/free_sum exactly absorbs `excess`), but the
        // pathological n_free==0 / free_sum<=0 break path (lines ~1930) leaves
        // the cell sum off by `excess` (all violators pinned at bounds, no free
        // obs to absorb). Likewise, kWaterFillMaxIter exhaustion can leave a
        // residual. Restore Σ=n (contract; see test-ieppa-bounds-mode L143).
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

    if (st.verbose >= 1) {
        const char* status_label =
            (res.base.status == RK_OK) ? "converged" :
            (res.base.status == RK_ERR_NOCONV) ? "max_iter exhausted (NOCONV)" :
            (res.base.status == RK_ERR_INFEAS) ? "infeasible" : "error";
        char msg[256];
        std::snprintf(msg, sizeof(msg),
                      "iEPPA %s in %d iters, errRp=%.3e",
                      status_label, res.base.iterations, res.base.max_error);
        st.log(msg);
    }

    if (st.verbose >= 1 && !structural_infeas_pairs.empty()) {
        char msg[256];
        size_t off = 0;
        off += std::snprintf(msg + off, sizeof(msg) - off,
                             "iEPPA persistent infeasible cells: ");
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
    res.sor_min_omega = sor_min_omega;
    res.sor_n_damped  = sor_n_damped;

    write_trajectory_csv(probe_samples);
}

}  // namespace lbw
