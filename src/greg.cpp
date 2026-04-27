#include "greg.hpp"
#include "calib_linalg.hpp"
#include "calib_dispatch.hpp"
#include "leafblower.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <numeric>
#include <limits>

namespace lbw {

GregResult greg_solve(CalibState& st) {
    GregResult res;
    res.status = RK_ERR_NOCONV;

    CellTable ct;
    if (build_cell_table(st.n, st.K, st.group_ids, st.cat_counts, st.weights, ct) != 0) {
        res.status = RK_ERR_BADARG;
        std::strncpy(res.message, "build_cell_table failed", sizeof(res.message) - 1);
        return res;
    }
    res.M_cell = ct.M_cell;

    std::vector<double> X_init(ct.M_cell, 0.0);
    for (int i = 0; i < st.n; i++) X_init[ct.cell_of[i]] += st.weights[i];

    const double lo = st.min_weight;
    const double hi = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
    std::vector<double> L_cell(ct.M_cell), U_cell(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) {
        L_cell[c] = lo * static_cast<double>(ct.n_per_cell[c]);
        U_cell[c] = hi * static_cast<double>(ct.n_per_cell[c]);
    }

    // cat_offset[k] = starting index for margin k in the n_cats_total vector
    std::vector<int> cat_offset(st.K);
    int n_cats_total = 0;
    for (int k = 0; k < st.K; k++) { cat_offset[k] = n_cats_total; n_cats_total += st.cat_counts[k]; }

    if (n_cats_total > kNCatsTotalMax) {
        res.status = RK_ERR_BADARG;
        std::snprintf(res.message, sizeof(res.message),
                      "n_cats_total=%d exceeds limit %d; use method='ieppa'",
                      n_cats_total, kNCatsTotalMax);
        return res;
    }

    std::vector<bool> fixed_lo(ct.M_cell, false), fixed_hi(ct.M_cell, false);
    std::vector<double> X(X_init);

    const int max_cats = *std::max_element(st.cat_counts, st.cat_counts + st.K);
    // Hoist per-iteration work vectors to avoid heap churn (D_eff/N/b reused each Newton step)
    std::vector<double> bucket_b(max_cats);
    std::vector<double> D_eff(ct.M_cell);
    std::vector<double> N(static_cast<size_t>(n_cats_total) * static_cast<size_t>(n_cats_total));
    std::vector<double> b(static_cast<size_t>(n_cats_total));
    const double n_total = static_cast<double>(st.n);

    static constexpr int kMaxNewtonIters = 10;
    static constexpr double kEps = 1e-10;

    for (int newton_iter = 0; newton_iter < kMaxNewtonIters; newton_iter++) {
        res.iterations = newton_iter + 1;

        std::fill(D_eff.begin(), D_eff.end(), 0.0);
        for (int c = 0; c < ct.M_cell; c++)
            if (!fixed_lo[c] && !fixed_hi[c] && X_init[c] > kEps)
                D_eff[c] = X_init[c];

        if (compute_normal_equations(ct, D_eff.data(), N.data(),
                                     cat_offset.data(), st.K,
                                     static_cast<size_t>(n_cats_total)) != RK_OK) {
            res.status = RK_ERR_BADARG; return res;
        }

        // b[k][j] = T_kj * n - sum_{c in (k,j)} X[c]  (marginal defect, fixed n)
        for (int k = 0; k < st.K; k++) {
            std::fill(bucket_b.begin(), bucket_b.begin() + st.cat_counts[k], 0.0);
            for (int c = 0; c < ct.M_cell; c++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k]) bucket_b[g] += X[c];
            }
            for (int j = 0; j < st.cat_counts[k]; j++)
                b[static_cast<size_t>(cat_offset[k]) + static_cast<size_t>(j)] =
                    st.targets[k][j] * n_total - bucket_b[j];
        }

        // LDLT solve: N * lambda = b
        if (ldlt_factor_inplace(N.data(), static_cast<size_t>(n_cats_total), 1e-10) != RK_OK) {
            res.status = RK_ERR_BADARG; return res;
        }
        ldlt_solve(N.data(), static_cast<size_t>(n_cats_total), b.data());
        const std::vector<double>& lambda = b;

        // Newton update: X_new[c] = X_init[c] * (1 + sum_k lambda[k, g_k[c]])
        bool any_clamped = false;
        for (int c = 0; c < ct.M_cell; c++) {
            if (fixed_lo[c]) { X[c] = L_cell[c]; continue; }
            if (fixed_hi[c]) { X[c] = U_cell[c]; continue; }
            if (X_init[c] <= kEps) continue;
            double sum_lam = 0.0;
            for (int k = 0; k < st.K; k++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k])
                    sum_lam += lambda[static_cast<size_t>(cat_offset[k]) + static_cast<size_t>(g)];
            }
            double X_new = X_init[c] * (1.0 + sum_lam);
            if (X_new < L_cell[c] - 1e-10) {
                X[c] = L_cell[c]; fixed_lo[c] = true; any_clamped = true;
            } else if (X_new > U_cell[c] + 1e-10) {
                X[c] = U_cell[c]; fixed_hi[c] = true; any_clamped = true;
            } else {
                X[c] = std::clamp(X_new, L_cell[c], U_cell[c]);
            }
        }

        if (!any_clamped) {
            res.status = RK_OK;
            res.convergence_iter = newton_iter + 1;
            break;
        }
    }

    // Compute 6 metrics at exit
    double W = 0.0;
    for (int c = 0; c < ct.M_cell; c++) W += X[c];
    {
        auto m = lbw::compute_cell_metrics(st, ct, X, W, bucket_b);
        res.max_error  = m.errRp;
        res.kl         = m.kl;
        res.chi2       = m.chi2;
        res.mean_error = m.mean_err;
        res.grake_norm = m.grake_norm;
        res.convergence_solver_objective = lbw::select_solver_objective(RK_ALG_GREG, m);
        res.best_error = m.chi2;
    }

    // Obs expansion + clamp
    const double hi_obs = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
    apply_obs_expansion(ct, X, X_init, st.n, lo, hi_obs, st.weights);

    res.best_weights.resize(st.n);
    std::copy(st.weights, st.weights + st.n, res.best_weights.begin());
    return res;
}

} // namespace lbw
