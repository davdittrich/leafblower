#include "ieppa.hpp"
#include "leafblower.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <vector>

namespace lbw {

// Compute errRp = max_k max_j |S_kj/W - tau_kj|
// O(n*K): single O(n) bucket accumulation pass per margin.
static double compute_errRp(const CalibState& st,
                              const std::vector<double>& w) {
    double W = 0.0;
    for (int i = 0; i < st.n; i++) W += w[i];

    double err = 0.0;
    for (int k = 0; k < st.K; k++) {
        std::vector<double> bucket(st.cat_counts[k], 0.0);
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

// Bregman distance D(w, prox) = sum_i [w_i*log(w_i/p_i) - w_i + p_i]
static double bregman_dist(const std::vector<double>& w,
                            const std::vector<double>& prox) {
    double d = 0.0;
    for (int i = 0; i < (int)w.size(); i++) {
        double wi = w[i], pi = prox[i];
        if (wi > 1e-300 && pi > 1e-300)
            d += wi * std::log(wi / pi) - wi + pi;
        else
            d += pi;
    }
    return d;
}

// One full BCD sweep over K margins. Mutates w in-place.
// Clamp is applied on the mean-normalized scale so bounds are consistent
// with harvest.R's post-return mean=1 normalization.
static void bcd_sweep(CalibState& st, std::vector<double>& w,
                      bool& infeas_flag) {
    for (int k = 0; k < st.K; k++) {
        double W = 0.0;
        std::vector<double> bucket(st.cat_counts[k], 0.0);
        for (int i = 0; i < st.n; i++) {
            W += w[i];
            int g = st.group_ids[k][i];
            if (g >= 0) bucket[g] += w[i];
        }

        std::vector<double> scale(st.cat_counts[k], 1.0);
        for (int j = 0; j < st.cat_counts[k]; j++) {
            double Tkj = st.targets[k][j] * W;
            if (bucket[j] < 1e-15 * W) {
                if (Tkj > 0.0) infeas_flag = true;
            } else {
                scale[j] = Tkj / bucket[j];
            }
        }

        for (int i = 0; i < st.n; i++) {
            int g = st.group_ids[k][i];
            if (g >= 0) {
                double wi = w[i] * scale[g];
                wi = std::max(st.min_weight, std::min(st.max_weight, wi));
                w[i] = wi;
            }
        }
    }
}

IEPPAResult ieppa_solve(CalibState& st) {
    IEPPAResult res;
    res.status = RK_ERR_NOCONV;
    res.iterations = 0;
    res.max_error = 1.0;

    std::vector<double> w(st.n);
    for (int i = 0; i < st.n; i++) w[i] = st.weights[i];

    std::vector<double> w_prox(w);
    double normU = st.max_weight;
    bool infeas_flag = false;

    double last_outer_errRp = 1.0;
    double tolRp = 1.0;

    for (int outer = 1; outer <= st.outer_max_iter; outer++) {
        res.iterations = outer;
        double tolRb = 1.0 / std::pow((double)outer, 1.1);

        double errRp_inner = 1.0;
        for (int inner = 0; inner < st.inner_max_iter; inner++) {
            bcd_sweep(st, w, infeas_flag);

            errRp_inner = compute_errRp(st, w);
            double D = bregman_dist(w, w_prox);
            double breg_crit = D / (1.0 + normU);

            if (errRp_inner < tolRp && breg_crit < tolRb) break;
        }

        double errRp = compute_errRp(st, w);
        res.max_error = errRp;

        if (st.verbose >= 1) {
            char msg[256];
            std::snprintf(msg, 256, "iEPPA outer iter %d: errRp=%.2e, tolRp=%.2e",
                     outer, errRp, tolRp);
            st.log(msg);
        }

        if (errRp < st.tol_abs) {
            res.status = infeas_flag ? RK_ERR_INFEAS : RK_OK;
            break;
        }

        w_prox = w;
        last_outer_errRp = errRp;
        tolRp = std::max(1e-6, last_outer_errRp / 1.5);
    }

    if (infeas_flag && res.status == RK_OK) res.status = RK_ERR_INFEAS;

    // Project to feasible set: iteratively normalize to mean=1 then clamp,
    // until the fixed point max(w)/mean(w) <= max_weight is stable.
    // Needed because harvest.R divides by mean(weights) after C returns,
    // so we must ensure max(w)/mean(w) <= max_weight on the scale C returns.
    {
        double lo = (st.min_weight > 0.0) ? st.min_weight : 0.0;
        double hi = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
        for (int proj = 0; proj < 200; proj++) {
            double wsum2 = 0.0;
            for (int i = 0; i < st.n; i++) wsum2 += w[i];
            double wm = (wsum2 > 1e-300) ? wsum2 / st.n : 1.0;
            for (int i = 0; i < st.n; i++) w[i] /= wm;
            bool changed = false;
            for (int i = 0; i < st.n; i++) {
                double wc = std::max(lo, std::min(hi, w[i]));
                if (wc != w[i]) { w[i] = wc; changed = true; }
            }
            if (!changed) break;
        }
    }

    for (int i = 0; i < st.n; i++) st.weights[i] = w[i];

    res.max_error = compute_errRp(st, std::vector<double>(st.weights, st.weights + st.n));

    return res;
}

} // namespace lbw
