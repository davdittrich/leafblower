#include "lbw_config.h"
#include "raking.hpp"
#include "leafblower.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
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
    std::vector<double> q_hyp(st.n, 0.0);  // Dykstra correction for hyperplane {w: sum(w)=n}

    double lo = st.min_weight;
    // 1e300 not numeric_limits::max(): prevents overflow in w[i] *= scale[g]
    // when bucket[g] is tiny (scale[g] can be large).
    double hi = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
    bool is_infeasible = false;

    int max_cats = *std::max_element(st.cat_counts, st.cat_counts + st.K);
    std::vector<double> bucket(max_cats), scale(max_cats);

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
        // q_hyp[i] accumulates overshoot from previous hyperplane projections.
        {
            for (int i = 0; i < st.n; i++) w[i] += q_hyp[i];
            double s = 0.0;
            for (int i = 0; i < st.n; i++) s += w[i];
            double shift = (static_cast<double>(st.n) - s) / static_cast<double>(st.n);
            for (int i = 0; i < st.n; i++) {
                double w_proj = w[i] + shift;
                q_hyp[i] = w[i] - w_proj;  // Dykstra correction: pre - post
                w[i] = w_proj;
            }
        }

        // Convergence check: run every kErrCheckInterval iters and on the final iter.
        if (iter == 1 || iter % kErrCheckInterval == 0 || iter == st.inner_max_iter) {
            double errRp = compute_errRp(st, w, bucket);
            res.max_error = errRp;

            if (st.verbose >= 1) {
                char msg[256];
                std::snprintf(msg, 256, "iEPPA iter %d: errRp=%.2e", iter, errRp);
                st.log(msg);
            }

            if (errRp < st.tol_abs) {
                res.status = is_infeasible ? RK_ERR_INFEAS : RK_OK;
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
        for (int i = 0; i < st.n; i++) w[i] += q_hyp[i];
        double s = 0.0;
        for (int i = 0; i < st.n; i++) s += w[i];
        double shift = (static_cast<double>(st.n) - s) / static_cast<double>(st.n);
        for (int i = 0; i < st.n; i++) {
            double w_proj = w[i] + shift;
            q_hyp[i] = w[i] - w_proj;
            w[i] = w_proj;
        }
        if (box_ok) break;
    }

    for (int i = 0; i < st.n; i++) st.weights[i] = w[i];

    return res;
}

} // namespace lbw
