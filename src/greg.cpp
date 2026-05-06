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
    // GREG defaults differ from CalibResult — override here to preserve existing behavior.
    res.base.status                       = RK_ERR_NOCONV;
    res.base.convergence_metric           = static_cast<int>(CalibMetric::CHI2);
    res.base.convergence_rule             = 0;
    res.base.convergence_tol             = 0.0;
    res.base.convergence_iter             = 1;
    res.base.convergence_minimized_metric = static_cast<int>(CalibMetric::CHI2);
    res.base.convergence_solver_objective = std::numeric_limits<double>::infinity();
    res.base.best_iter                    = 1;

    CellTable ct;
    std::vector<double> X_init;
    double hi_eff;
    std::vector<double> L_cell, U_cell;
    std::vector<int> cat_offset;
    int n_cats_total;
    if (lbw::solver_setup_ct(st, ct, X_init, hi_eff, L_cell, U_cell,
                              cat_offset, n_cats_total, res) != RK_OK)
        return res;
    res.M_cell = ct.M_cell;
    const double lo = st.min_weight;
    const double hi = hi_eff;

    std::vector<bool> fixed_lo(ct.M_cell, false), fixed_hi(ct.M_cell, false);
    std::vector<double> X(X_init);

    const int max_cats = lbw::max_cats_count(st.K, st.cat_counts);
    // Hoist per-iteration work vectors to avoid heap churn (D_eff/N/b reused each Newton step)
    std::vector<double> bucket_b(max_cats);
    std::vector<double> D_eff(ct.M_cell);
    // G3: compute normal equations directly into N_factored; eliminate the N copy.
    std::vector<double> N_factored(static_cast<size_t>(n_cats_total) * static_cast<size_t>(n_cats_total));
    std::vector<double> b(static_cast<size_t>(n_cats_total));
    const double n_total = static_cast<double>(st.n);

    static constexpr int kMaxNewtonIters = 50;
    static constexpr double kEps = 1e-10;

    std::vector<bool> prev_fixed_lo(ct.M_cell, false), prev_fixed_hi(ct.M_cell, false);
    bool need_refactor = true; // R1: refactor only when active set changes

    for (int newton_iter = 0; newton_iter < kMaxNewtonIters; newton_iter++) {
        res.base.iterations = newton_iter + 1;

        if (need_refactor) {
            std::fill(D_eff.begin(), D_eff.end(), 0.0);
            for (int c = 0; c < ct.M_cell; c++)
                if (!fixed_lo[c] && !fixed_hi[c] && X_init[c] > kEps)
                    D_eff[c] = X_init[c];

            if (compute_normal_equations(ct, D_eff.data(), N_factored.data(),
                                         cat_offset.data(), st.K,
                                         static_cast<size_t>(n_cats_total)) != RK_OK) {
                res.base.status = RK_ERR_BADARG; return res;
            }
            // Tikhonov ridge: add τI to normal-equations matrix before Cholesky.
            if (st.ridge_lambda > 0.0) {
                for (size_t j = 0; j < static_cast<size_t>(n_cats_total); j++)
                    N_factored[j * static_cast<size_t>(n_cats_total) + j] += st.ridge_lambda;
            }
            if (cholesky_factor_inplace(N_factored.data(), static_cast<size_t>(n_cats_total), 1e-10) != RK_OK) {
                res.base.status = RK_ERR_BADARG; return res;
            }
        }

        // b[k][j] = T_kj * n - sum_{c in (k,j)} X[c]  (marginal defect; b changes with X every iter)
        for (int k = 0; k < st.K; k++) {
            lbw::aggregate_to_margin(ct, X, k, st.cat_counts[k], bucket_b.data());
            for (int j = 0; j < st.cat_counts[k]; j++)
                b[static_cast<size_t>(cat_offset[k]) + static_cast<size_t>(j)] =
                    st.targets[k][j] * n_total - bucket_b[j];
        }

        // LDLT solve using cached factored matrix (recomputed only on active-set change)
        cholesky_solve(N_factored.data(), static_cast<size_t>(n_cats_total), b.data());
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
            const double cell_tol = 1e-10 * std::max(1.0, X_init[c]);
            if (X_new < L_cell[c] - cell_tol) {
                X[c] = L_cell[c]; fixed_lo[c] = true; any_clamped = true;
            } else if (X_new > U_cell[c] + cell_tol) {
                X[c] = U_cell[c]; fixed_hi[c] = true; any_clamped = true;
            } else {
                X[c] = std::clamp(X_new, L_cell[c], U_cell[c]);
            }
        }

        // R1: invalidate factorization cache only when active set changed
        need_refactor = (fixed_lo != prev_fixed_lo || fixed_hi != prev_fixed_hi);
        prev_fixed_lo = fixed_lo;
        prev_fixed_hi = fixed_hi;

        // C6: KKT release pass — drop constraints whose multiplier has wrong sign
        for (int c = 0; c < ct.M_cell; ++c) {
            if (X_init[c] <= 1e-10) continue;
            double grad_c = (X[c] - X_init[c]) / X_init[c];
            if (fixed_lo[c] && grad_c < -1e-8) {
                fixed_lo[c] = false; need_refactor = true; any_clamped = true;
            }
            if (fixed_hi[c] && grad_c > 1e-8) {
                fixed_hi[c] = false; need_refactor = true; any_clamped = true;
            }
        }

        if (!any_clamped) {
            res.base.status = RK_OK;
            res.base.convergence_iter = newton_iter + 1;
            break;
        }
    }

    if (res.base.status == RK_ERR_NOCONV) {
        std::snprintf(res.message, sizeof(res.message),
                      "greg: no convergence after %d Newton steps; active set budget exceeded",
                      kMaxNewtonIters);
    }

    // Compute 6 metrics at exit
    double W = 0.0;
    for (int c = 0; c < ct.M_cell; c++) W += X[c];
    {
        auto m = lbw::compute_cell_metrics(st, ct, X, W, bucket_b);
        res.base.max_error  = m.errRp;
        res.base.kl         = m.kl;
        res.base.chi2       = m.chi2;
        res.base.mean_error = m.mean_err;
        res.base.grake_norm = m.grake_norm;
        res.base.convergence_solver_objective = lbw::select_solver_objective(RK_ALG_GREG, m);
        res.base.best_error = m.chi2;
    }

    // Obs expansion + clamp
    const double hi_obs = lbw::resolve_hi(st);
    apply_obs_expansion(ct, X, X_init, st.n, lo, hi_obs, st.weights);

    res.base.best_weights.resize(st.n);
    std::copy(st.weights, st.weights + st.n, res.base.best_weights.begin());
    return res;
}

} // namespace lbw
