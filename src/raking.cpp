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
#include "calib_dispatch.hpp"
#include "leafblower.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <limits>
#include <vector>

namespace lbw {

// Sum weights with 4-way ILP unroll.
// Separate accumulators break the loop-carried dependency chain, letting
// the compiler pipeline four additions in parallel. The tail loop handles n % 4.
static double sum_weights_ilp(const std::vector<double>& w, int n) {
    double W = 0.0, W1 = 0.0, W2 = 0.0, W3 = 0.0;
    const int n4 = n & ~3;
    for (int i = 0; i < n4; i += 4) {
        W  += w[i];   W1 += w[i+1];
        W2 += w[i+2]; W3 += w[i+3];
    }
    for (int i = n4; i < n; ++i) W += w[i];
    return W + W1 + W2 + W3;
}

// Compute errRp = max_k max_j |S_kj/W - tau_kj|
// O(n*K): single O(n) bucket accumulation pass per margin.
// bucket must be pre-allocated to at least max_cats elements by the caller;
// it is filled and reused across margins to avoid per-call heap allocation.
static double compute_errRp(const CalibState& st,
                              const std::vector<double>& w,
                              std::vector<double>& bucket) {
    double W = sum_weights_ilp(w, st.n);

    double err = 0.0;
    for (int k = 0; k < st.K; k++) {
        std::fill(bucket.begin(), bucket.begin() + st.cat_counts[k], 0.0);
        for (int i = 0; i < st.n; i++) {
            int g = st.group_ids[k][i];
            if (g >= 0) bucket[g] += w[i];
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
    static constexpr double kEmptyBucketThreshold   = 1e-15;   // relative threshold: bucket[j] < 1e-15*W → treat as empty, skip IPF scale
    static constexpr int    kErrCheckInterval        = 10;      // Check convergence every N inner iterations.
                                                                 // compute_errRp costs K O(n) passes — nearly as expensive as a full sweep.
                                                                 // Every-10 reduces that overhead by 90% at the cost of ≤9 extra IPF iters.
                                                                 // Exception: always check on iter 1 to catch problems that converge immediately.

    RakingResult res;
    res.status = RK_ERR_NOCONV;
    res.iterations = 0;
    res.max_error = 1.0;

    std::vector<double> w(st.weights, st.weights + st.n);

    std::vector<double> q(st.n, 0.0);
    // Dykstra hyperplane correction: at any fixed point of the hyperplane
    // projection, `q_hyp[i] = w[i] - w_proj = -shift` is identical for all i
    // (shift depends only on the total sum, not on i). Store scalar, not vector.
    double q_hyp = 0.0;

    double lo = st.min_weight;
    // 1e300 not numeric_limits::max(): prevents overflow in w[i] *= scale[g]
    // when bucket[g] is tiny (scale[g] can be large).
    double hi = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
    bool is_infeasible = false;

    int max_cats = *std::max_element(st.cat_counts, st.cat_counts + st.K);
    std::vector<double> bucket(max_cats), scale(max_cats);

    // Descent monitor state. Detects stalled convergence: no net progress over
    // a sliding window of kMaxNoImprove consecutive error-checks. Uses a
    // window-minimum comparison to avoid false positives on normal oscillation
    // (errRp often wobbles within noise even while the running minimum still
    // decreases). Firing condition: the window minimum has not improved by
    // more than a relative tolerance across kMaxNoImprove checks.
    constexpr int kMaxNoImprove = 5;  // 5 * kErrCheckInterval = 50 stalled iters
    double min_errRp_window = std::numeric_limits<double>::infinity();
    int n_no_improve = 0;

    // WU-B: W_prev tracks obs-level weights at the last convergence check for pct_change.
    std::vector<double> w_prev(w);

    // WU-D: prev_metric_for_rule tracks the active metric value from the previous
    // convergence check, used by IMPROVEMENT and PLATEAU CalibRule dispatch.
    double prev_metric_for_rule = std::numeric_limits<double>::infinity();

    // WU-E: best-iterate tracking (obs-level snapshot at min observed errRp).
    // Tracks errRp regardless of active convergence criterion.
    double best_errRp_seen = std::numeric_limits<double>::infinity();
    int    best_iter_val   = 0;
    std::vector<double> w_best(st.n, 0.0);

    for (int iter = 1; iter <= st.inner_max_iter; iter++) {
        res.iterations = iter;

        // Marginal projections: pure cyclic IPF (no clamp, no Euclidean correction).
        // Euclidean Dykstra corrections are incompatible with multiplicative IPF steps.
        for (int k = 0; k < st.K; k++) {
            // Bucket accumulation for IPF scale computation
            // W sum separated from scatter-add so the compiler can vectorise it.
            double W = sum_weights_ilp(w, st.n);

            // Bucket scatter-add: write aliases prevent vectorisation.
            std::fill(bucket.begin(), bucket.begin() + st.cat_counts[k], 0.0);
            for (int i = 0; i < st.n; i++) {
                int g = st.group_ids[k][i];
                if (g >= 0) bucket[g] += w[i];
            }

            // IPF scale factors
            std::fill(scale.begin(), scale.begin() + st.cat_counts[k], 1.0);
            for (int j = 0; j < st.cat_counts[k]; j++) {
                double Tkj = st.targets[k][j] * W;
                if (bucket[j] < kEmptyBucketThreshold * W) {
                    if (Tkj > 0.0) is_infeasible = true;
                } else {
                    scale[j] = Tkj / bucket[j];
                }
            }

            // IPF step — NO CLAMP. g==-1 (NA) entries pass through unchanged.
            for (int i = 0; i < st.n; i++) {
                int g = st.group_ids[k][i];
                if (g >= 0) w[i] *= scale[g];
            }
        }

        // Box projection [lo, hi]^n with Dykstra correction (mean=1 scale).
        // q[i] accumulates overshoot from previous box clamps.
#if defined(_OPENMP) || LBW_HAS_OMP_SIMD
#pragma omp simd
#endif
        for (int i = 0; i < st.n; i++) {
            double yi = w[i] + q[i];
            double wc = std::clamp(yi, lo, hi);
            q[i] = yi - wc;
            w[i] = wc;
        }

        // Dykstra hyperplane projection: {w : sum(w) = n}
        // q_hyp is a scalar: the hyperplane correction is uniform across i.
        {
            double s = 0.0;
#if defined(_OPENMP) || LBW_HAS_OMP_SIMD
#pragma omp simd reduction(+:s)
#endif
            for (int i = 0; i < st.n; i++) {
                w[i] += q_hyp;  // apply prior correction uniformly
                s += w[i];
            }
            double shift = (static_cast<double>(st.n) - s) / static_cast<double>(st.n);
#if defined(_OPENMP) || LBW_HAS_OMP_SIMD
#pragma omp simd
#endif
            for (int i = 0; i < st.n; i++) w[i] += shift;
            q_hyp = -shift;  // w_pre_proj - w_post_proj = (w + q_hyp_old) - (w + q_hyp_old + shift) = -shift
        }

        // Convergence check: run every kErrCheckInterval iters and on the final iter.
        if (iter == 1 || iter % kErrCheckInterval == 0 || iter == st.inner_max_iter) {
            double errRp = compute_errRp(st, w, bucket);
            res.max_error = errRp;

            // WU-E: update best-iterate snapshot (tracks errRp regardless of active criterion).
            if (errRp < best_errRp_seen) {
                best_errRp_seen = errRp;
                best_iter_val   = iter;
                w_best          = w;  // obs-level snapshot
            }

            // WU-B: pct_change (max-relative, L∞) and l1_weight (Σ|Δw|/n, L1/n).
            double pct_change = 0.0;
            double l1_weight_sum = 0.0;
            for (int i = 0; i < st.n; i++) {
                double diff = std::fabs(w[i] - w_prev[i]);
                double rel = diff / std::max(w_prev[i], 1e-12);
                if (rel > pct_change) pct_change = rel;
                l1_weight_sum += diff;
            }
            double l1_weight = l1_weight_sum / static_cast<double>(st.n);

            // WU-B: alternative metrics.
            //
            // Gate: skip the O(K*n) scatter-adds when the active stopping
            // criterion does not use these metrics. Always compute on the final
            // iteration (iter == st.inner_max_iter) so that calib_result is
            // fully populated on exit.
            const lbw::CalibMetric metric = st.convergence_cfg.metric;
            const auto& cfg_m = st.convergence_cfg;
            // Force full metrics when convergence is about to fire (same fix as ieppa).
            const bool about_to_converge =
                (cfg_m.absolute_tol > 0.0 && errRp < cfg_m.absolute_tol) ||
                (metric == lbw::CalibMetric::MAX_ERR &&
                 lbw::apply_rule(cfg_m.rule, errRp, prev_metric_for_rule, cfg_m.pct_tol)) ||
                (metric == lbw::CalibMetric::L1_WEIGHT &&
                 lbw::apply_rule(cfg_m.rule, l1_weight, prev_metric_for_rule, cfg_m.pct_tol));
            const bool need_extra_metrics =
                (metric == lbw::CalibMetric::MEAN_ERR   ||
                 metric == lbw::CalibMetric::KL         ||
                 metric == lbw::CalibMetric::CHI2       ||
                 metric == lbw::CalibMetric::GRAKE_NORM ||
                 iter == st.inner_max_iter               ||
                 about_to_converge);

            double W_total = sum_weights_ilp(w, st.n);
            constexpr double kMetricEps = 1e-10;
            constexpr double kChi2Floor = 1.0;
            double mean_err_sum = 0.0;
            double kl_max       = 0.0;
            double chi2_total   = 0.0;
            double grake_norm   = 0.0;  // WU-D: max over margins of |S_kj - T_kj*W| / (1 + |T_kj*W|)
            if (need_extra_metrics && W_total > 0.0) {
                for (int k = 0; k < st.K; k++) {
                    const int nj = st.cat_counts[k];
                    std::fill(bucket.begin(), bucket.begin() + nj, 0.0);
                    for (int i = 0; i < st.n; i++) {
                        int g = st.group_ids[k][i];
                        if (g >= 0 && g < nj) bucket[g] += w[i];
                    }
                    double max_k = 0.0;
                    double kl_k  = 0.0;
                    for (int j = 0; j < nj; j++) {
                        double S_p = bucket[j] / W_total;
                        double T   = st.targets[k][j];
                        double err = std::fabs(S_p - T);
                        if (err > max_k) max_k = err;
                        if (T > 0.0) {
                            kl_k += T * std::log((T + kMetricEps) / (S_p + kMetricEps));
                        }
                        double obs     = bucket[j];
                        double pop_kj  = T * W_total;
                        chi2_total += (obs - pop_kj) * (obs - pop_kj) / (pop_kj + kChi2Floor);
                        // WU-D: grake_norm = max_kj |S_kj - pop_kj| / (1 + |pop_kj|)
                        double nm = std::fabs(obs - pop_kj) / (1.0 + std::fabs(pop_kj));
                        if (nm > grake_norm) grake_norm = nm;
                    }
                    mean_err_sum += max_k;
                    if (kl_k > kl_max) kl_max = kl_k;
                }
            }
            double mean_err = (st.K > 0) ? (mean_err_sum / static_cast<double>(st.K)) : 0.0;

            // Store metrics and update w_prev (unconditional — intermediate checks
            // store 0 for gated metrics; final iter always populates all fields).
            res.l1_weight_change = l1_weight;    // Σ|Δw|/n (L1/n, scale-free)
            res.mean_error       = mean_err;
            res.kl               = kl_max;
            res.chi2             = chi2_total;
            res.grake_norm       = grake_norm;  // WU-D
            for (int i = 0; i < st.n; i++) w_prev[i] = w[i];

            // First check: no baseline yet, just record and reset counter.
            // On subsequent checks, require relative improvement of 1% of the
            // current window minimum (floored at tol_abs so a valid convergence
            // does not trip the monitor right before the convergence check).
            if (!std::isfinite(min_errRp_window)) {
                min_errRp_window = errRp;
                n_no_improve = 0;
            } else {
                const double rel_eps = 0.01 * min_errRp_window;
                const double eps = std::max(rel_eps, st.tol_abs);
                if (errRp < min_errRp_window - eps) {
                    min_errRp_window = errRp;
                    n_no_improve = 0;
                } else {
                    n_no_improve++;
                }
            }

            if (st.verbose >= 1) {
                char msg[256];
                std::snprintf(msg, 256, "raking iter %d: errRp=%.2e", iter, errRp);
                st.log(msg);
            }

            // WU-D: active-criterion + CalibRule dispatch + stop_when ANY/ALL logic.
            // The active metric value is selected by cfg.metric; the stopping rule
            // is determined by cfg.rule (THRESHOLD/IMPROVEMENT/PLATEAU).
            // Absolute-tol path uses the metric value directly vs absolute_tol.
            // PCT-tol path applies the CalibRule to the metric value.
            {
                const auto& cfg = st.convergence_cfg;

                // Select the active metric value (shared helper).
                const double curr_metric = lbw::select_metric(
                    cfg.metric, errRp, mean_err, kl_max, chi2_total, grake_norm, l1_weight);

                // Absolute-tol convergence: metric < absolute_tol.
                bool converged_abs = (cfg.absolute_tol > 0.0) && (curr_metric < cfg.absolute_tol);

                // PCT-tol convergence: apply CalibRule to curr_metric vs prev_metric_for_rule.
                // apply_rule takes prev by VALUE; update prev_metric_for_rule after the call.
                const bool converged_pct = lbw::apply_rule(
                    cfg.rule, curr_metric, prev_metric_for_rule, cfg.pct_tol);
                prev_metric_for_rule = curr_metric;

                bool have_pct = (cfg.pct_tol > 0.0);
                bool have_abs = (cfg.absolute_tol > 0.0);
                bool converged = false;
                if (have_pct && have_abs) {
                    converged = (cfg.stop_when == lbw::CalibStopWhen::ALL)
                                ? (converged_pct && converged_abs)
                                : (converged_pct || converged_abs);
                } else if (have_pct) {
                    converged = converged_pct;
                } else if (have_abs) {
                    converged = converged_abs;
                } else {
                    // Fallback: legacy tol_abs check.
                    converged = (errRp < st.tol_abs);
                }

                if (converged) {
                    res.status             = is_infeasible ? RK_ERR_INFEAS : RK_OK;
                    res.convergence_metric = static_cast<int>(cfg.metric);
                    res.convergence_rule   = static_cast<int>(cfg.rule);
                    res.convergence_tol    = cfg.pct_tol;
                    res.convergence_iter   = iter;
                    break;
                }
            }

            if (n_no_improve >= kMaxNoImprove) {
                res.status = RK_ERR_NOCONV;
                if (st.verbose >= 1) {
                    char msg[256];
                    std::snprintf(msg, 256,
                                  "raking: errRp stalled for %d consecutive checks "
                                  "(last=%.2e, window_min=%.2e); likely near-infeasible bounds. "
                                  "Aborting at iter %d.",
                                  n_no_improve, errRp, min_errRp_window, iter);
                    st.log(msg);
                }
                break;
            }
        }
    }

    // Infeasibility detected during iteration: override NOCONV with INFEAS.
    // Truly infeasible problems (empty bucket with positive target) can never
    // converge — the empty bucket always contributes τ_kj > 0 to errRp,
    // so errRp never drops below tol_abs and the convergence break never fires.
    // Check the flag here and return the correct status code.
    if (is_infeasible && res.status == RK_ERR_NOCONV)
        res.status = RK_ERR_INFEAS;

    // Post-loop Dykstra finalizer: alternate box+hyperplane until box-feasible.
    // The main loop exits after the hyperplane step; the resulting shift can push
    // weights fractionally above hi or below lo. Continue the Dykstra cycle
    // (box then hyperplane) until all weights are within [lo, hi].
    // At true convergence this terminates in 1 iteration (shift ~ floating-point
    // rounding). At the end, sum(w) = n, so harvest.R's /mean(weights) is a no-op.
    for (int fixup = 0; fixup < 20; fixup++) {
        bool box_ok = true;
        for (int i = 0; i < st.n; i++) {
            double yi = w[i] + q[i];
            double wc = std::clamp(yi, lo, hi);
            q[i] = yi - wc;
            if (yi != wc) box_ok = false;
            w[i] = wc;
        }
        // Hyperplane step restores sum(w) = n regardless of box changes.
        // q_hyp is scalar (uniform correction); applied uniformly to all w[i].
        double s = 0.0;
#if defined(_OPENMP) || LBW_HAS_OMP_SIMD
#pragma omp simd reduction(+:s)
#endif
        for (int i = 0; i < st.n; i++) {
            w[i] += q_hyp;
            s += w[i];
        }
        double shift = (static_cast<double>(st.n) - s) / static_cast<double>(st.n);
#if defined(_OPENMP) || LBW_HAS_OMP_SIMD
#pragma omp simd
#endif
        for (int i = 0; i < st.n; i++) w[i] += shift;
        q_hyp = -shift;
        if (box_ok) break;
    }

    // WU-E: finalize best-iterate result (obs-level, sum-normalize to n, no clamping).
    res.best_error = best_errRp_seen;
    res.best_iter  = best_iter_val;
    if (std::isfinite(best_errRp_seen)) {
        double s = 0.0;
        for (int i = 0; i < st.n; i++) s += w_best[i];
        if (s > 0.0) {
            const double scale = static_cast<double>(st.n) / s;
            for (int i = 0; i < st.n; i++) w_best[i] *= scale;
        }
        res.best_weights = std::move(w_best);
    } else {
        res.best_weights.assign(st.n, 0.0);
    }

    for (int i = 0; i < st.n; i++) st.weights[i] = w[i];

    // Solver-owned normalization (moved from wrapper 2026-04-24 per user directive).
    // Defensive: the preceding hyperplane finalizer (lines above) enforces
    // sum(w) = n unconditionally on every iteration via shift = (n - s)/n, so
    // this block is a no-op at both CONV and NOCONV. Kept for solver-contract
    // self-containment against future refactors that might weaken the invariant.
    // total_w == 0 is pathological (all-zero input weights); leave unchanged.
    double total_w = 0.0;
    for (int i = 0; i < st.n; i++) total_w += st.weights[i];
    if (std::isfinite(total_w) && total_w > 0.0) {
        const double norm = static_cast<double>(st.n) / total_w;
        for (int i = 0; i < st.n; i++) st.weights[i] *= norm;
    }

    return res;
}

} // namespace lbw
