#include "lbw_config.h"
#include "ieppa.hpp"
#include "leafblower.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <vector>

namespace lbw {

// Compute errRp = max_k max_j |S_kj/W - tau_kj|
// O(n*K): single O(n) bucket accumulation pass per margin.
// bucket must be pre-allocated to at least max_cats elements by the caller;
// it is filled and reused across margins to avoid per-call heap allocation.
static double compute_errRp(const CalibState& st,
                              const std::vector<double>& w,
                              std::vector<double>& bucket) {
    double W = 0.0, W1 = 0.0, W2 = 0.0, W3 = 0.0;
    int i4e = st.n & ~3;
    for (int i = 0; i < i4e; i += 4) {
        W  += w[i];   W1 += w[i+1];
        W2 += w[i+2]; W3 += w[i+3];
    }
    for (int i = i4e; i < st.n; ++i) W += w[i];
    W += W1 + W2 + W3;

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
IEPPAResult ieppa_solve(CalibState& st) {
    static constexpr double kEmptyBucketThreshold   = 1e-15;   // relative threshold: bucket[j] < 1e-15*W → treat as empty, skip IPF scale
    static constexpr double kWeightCollapseThreshold = 1e-300;  // weights collapsed: skip norm
    static constexpr int    kMaxFixupIterations      = 20;      // post-convergence fixup cap
    static constexpr int    kErrCheckInterval        = 10;      // Check convergence every N iterations instead of every 1.
                                                                 // compute_errRp costs K O(n) passes — nearly as expensive as a full sweep.
                                                                 // Every-10 reduces that overhead by 90% at the cost of ≤9 extra IPF iters.

    IEPPAResult res;
    res.status = RK_ERR_NOCONV;
    res.iterations = 0;
    res.max_error = 1.0;

    std::vector<double> w(st.weights, st.weights + st.n);

    std::vector<double> q(st.n, 0.0);

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
            double W = 0.0, W1 = 0.0, W2 = 0.0, W3 = 0.0;
            int ni = st.n, i4 = ni & ~3;
            for (int i = 0; i < i4; i += 4) {
                W  += w[i];   W1 += w[i+1];
                W2 += w[i+2]; W3 += w[i+3];
            }
            for (int i = i4; i < ni; ++i) W += w[i];
            W += W1 + W2 + W3;

            // Bucket scatter-add: write aliases prevent vectorisation.
            std::fill(bucket.begin(), bucket.begin() + st.cat_counts[k], 0.0);
            for (int i = 0; i < ni; i++) {
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

        // Normalize to mean=1 so box bounds match the scale R returns.
        // harvest.R divides by mean(weights) after C returns, so the constraint is
        // w_i/mean(w) <= max_weight. Working at mean=1 makes max(w) == max/mean(w).
        {
            double Wsum = 0.0;
            for (int i = 0; i < st.n; i++) Wsum += w[i];
            double wm = Wsum / st.n;
            if (wm > kWeightCollapseThreshold) {
                // Rescale q[] proportionally to w[]: q[i] represents Dykstra overshoot
                // in the same unit as w[i]. After renormalizing w[i] /= wm, the corrected
                // iterate y[i] = w[i] + q[i] must shift by the same factor to keep
                // the Dykstra fixed-point invariant intact.
#if defined(_OPENMP) || LBW_HAS_OMP_SIMD
#pragma omp simd
#endif
                for (int i = 0; i < st.n; i++) { w[i] /= wm; q[i] /= wm; }
            }
        }

        // Box projection [lo, hi]^n with Dykstra correction (mean=1 scale).
        // q[i] accumulates overshoot from previous box clamps.
#if defined(_OPENMP) || LBW_HAS_OMP_SIMD
#pragma omp simd
#endif
        for (int i = 0; i < st.n; i++) {
            double yi = w[i] + q[i];
            double wc = std::max(lo, std::min(hi, yi));
            q[i] = yi - wc;
            w[i] = wc;
        }

        // Convergence check: run every kErrCheckInterval iters and on the final iter.
        if (iter % kErrCheckInterval == 0 || iter == st.inner_max_iter) {
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

    // Final normalization-and-clamp fixup: harvest.R divides by mean(weights) after
    // C returns, so the effective constraint is max(w)/mean(w) <= max_weight.
    // The box projection inside the loop works in mean~=1 space but clamping reduces
    // mean slightly, causing post-R-normalization to push max fractionally above hi.
    // Fix: iterate renormalize→reclamp until the fixed point max(w)<=hi*mean(w) holds.
    bool fixup_converged = false;
    for (int fixup = 0; fixup < kMaxFixupIterations; fixup++) {
        double Wsum = 0.0;
        for (int i = 0; i < st.n; i++) Wsum += w[i];
        double wm = (Wsum > kWeightCollapseThreshold) ? Wsum / st.n : 1.0;
        bool changed = false;
        for (int i = 0; i < st.n; i++) {
            w[i] /= wm;  // normalize to mean=1
            double wc = std::max(lo, std::min(hi, w[i]));
            if (wc != w[i]) { w[i] = wc; changed = true; }
        }
        if (!changed) { fixup_converged = true; break; }
    }
    if (!fixup_converged)
        st.log("iEPPA: fixup loop did not reach fixed point in 20 iterations; "
               "weights may exceed max_weight by floating-point rounding");

    for (int i = 0; i < st.n; i++) st.weights[i] = w[i];

    return res;
}

} // namespace lbw
