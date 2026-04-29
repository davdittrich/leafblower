#include "lbw_config.h"
#include "logit_calib.hpp"
#include "calib_dispatch.hpp"
#include "calib_linalg.hpp"
#include "cell_table.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <numeric>

namespace lbw {

LogitCalibResult logit_calibrate(CalibState& st) {
    LogitCalibResult res;
    res.status = RK_ERR_NOCONV;

    // Input validation: logit with max_weight=Inf blows up
    if (!std::isfinite(st.max_weight) || st.max_weight <= 0.0) {
        res.status = RK_ERR_BADARG;
        std::snprintf(res.message, sizeof(res.message),
            "logit: max_weight must be finite and positive");
        return res;
    }
    if (st.min_weight >= st.max_weight) {
        res.status = RK_ERR_BADARG;
        std::snprintf(res.message, sizeof(res.message),
            "logit: min_weight (%.4g) >= max_weight (%.4g)",
            st.min_weight, st.max_weight);
        return res;
    }

    // Build cell table
    CellTable ct;
    {
        std::vector<const int32_t*> gids(st.K);
        for (int k = 0; k < st.K; k++) gids[k] = st.group_ids[k];
        if (build_cell_table(st.n, st.K, gids.data(), st.cat_counts, st.weights, ct) != RK_OK) {
            res.status = RK_ERR_INFEAS;
            std::snprintf(res.message, sizeof(res.message),
                "logit: cell table build failed");
            return res;
        }
    }

    const int M = ct.M_cell;
    const int K = st.K;

    // cells_per_cat[k][j] = list of cell indices in bucket (k,j)
    int max_cats = 0;
    for (int k = 0; k < K; k++) max_cats = std::max(max_cats, st.cat_counts[k]);
    std::vector<std::vector<std::vector<int>>> cells_per_cat(K);
    for (int k = 0; k < K; k++) {
        cells_per_cat[k].assign(st.cat_counts[k], {});
        for (int c = 0; c < M; c++) {
            int g = ct.g_per_cell[k][c];
            if (g >= 0 && g < st.cat_counts[k]) cells_per_cat[k][g].push_back(c);
        }
    }

    // cat_offset and nct
    int nct = 0;
    std::vector<int> cat_offset(K);
    for (int k = 0; k < K; k++) { cat_offset[k] = nct; nct += st.cat_counts[k]; }

    // X_init[c] from design weights
    std::vector<double> X_init(M, 0.0);
    for (int i = 0; i < st.n; i++) X_init[ct.cell_of[i]] += st.weights[i];

    // Capacity bounds per cell
    std::vector<double> L_cell(M), U_cell(M);
    for (int c = 0; c < M; c++) {
        L_cell[c] = st.min_weight * ct.n_per_cell[c];
        U_cell[c] = st.max_weight * ct.n_per_cell[c];
    }

    // Newton state
    std::vector<double> lambda(nct, 0.0);   // dual variables; init=0 (midpoint warm-start)
    std::vector<double> w(M, 0.0);          // calibrated cell masses
    std::vector<double> D_eff(M, 0.0);      // Newton weights = (U-L)*sigma*(1-sigma)
    std::vector<double> N(nct * nct, 0.0);  // normal equations matrix
    std::vector<double> b(nct, 0.0);        // residuals / Newton RHS
    std::vector<double> bucket_scratch(max_cats, 0.0);
    std::vector<double> w_best;

    double initial_resid = std::numeric_limits<double>::infinity();
    double best_resid    = std::numeric_limits<double>::infinity();
    double prev_metric   = std::numeric_limits<double>::infinity();
    const CalibConvergenceCfg& cfg = st.convergence_cfg;

    const int kMaxNewtonIters = std::min(50, st.inner_max_iter);

    for (int iter = 0; iter < kMaxNewtonIters; iter++) {
        res.iterations = iter + 1;

        // (1) Compute w[c] and D_eff[c] from lambda
        for (int c = 0; c < M; c++) {
            double z = 0.0;
            for (int k = 0; k < K; k++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k]) z += lambda[cat_offset[k] + g];
            }
            z = std::clamp(z, -700.0, 700.0);  // prevent exp overflow
            double sig   = 1.0 / (1.0 + std::exp(-z));
            double range = U_cell[c] - L_cell[c];
            w[c]     = L_cell[c] + range * sig;        // bounds by construction, no clamp
            D_eff[c] = range * sig * (1.0 - sig);      // D-S 1992 Newton weight
        }

        // (2) Residuals b[cat_offset[k]+j] = tau*n - sum_{c in bucket} w[c]
        std::fill(b.begin(), b.end(), 0.0);
        for (int k = 0; k < K; k++) {
            for (int j = 0; j < st.cat_counts[k]; j++) {
                double target = st.targets[k][j] * static_cast<double>(st.n);
                double S_kj = 0.0;
                for (int c : cells_per_cat[k][j]) S_kj += w[c];
                b[cat_offset[k] + j] = target - S_kj;
            }
        }

        // Track initial residual for post-loop status
        double max_b = 0.0;
        for (double bj : b) max_b = std::max(max_b, std::abs(bj));
        if (iter == 0) initial_resid = max_b;
        if (max_b < best_resid) {
            best_resid = max_b;
            w_best = w;
            res.best_iter = iter + 1;
        }

        // (3) Build N = A*diag(D_eff)*A^T and solve N*delta_lambda = b
        std::fill(N.begin(), N.end(), 0.0);
        if (compute_normal_equations(ct, D_eff.data(), N.data(),
                                      cat_offset.data(), K,
                                      static_cast<size_t>(nct)) != RK_OK) {
            res.status = RK_ERR_INFEAS;
            std::snprintf(res.message, sizeof(res.message),
                "logit: singular normal equations (degenerate bounds - L=U cells)");
            break;
        }
        if (ldlt_factor_inplace(N.data(), static_cast<size_t>(nct), 1e-10) != RK_OK) {
            res.status = RK_ERR_INFEAS;
            std::snprintf(res.message, sizeof(res.message),
                "logit: LDLT factorization failed (degenerate bounds)");
            break;
        }
        ldlt_solve(N.data(), static_cast<size_t>(nct), b.data());  // b = delta_lambda

        // (4) Update lambda
        for (int j = 0; j < nct; j++) lambda[j] += b[j];

        // (5) Convergence check via shared infrastructure
        double W_total = 0.0;
        for (int c = 0; c < M; c++) W_total += w[c];
        lbw::CellMetrics m = lbw::compute_cell_metrics(st, ct, w, W_total, bucket_scratch);
        bool converged = lbw::check_convergence(cfg, m, prev_metric, st.tol_abs);
        if (converged) {
            res.status = RK_OK;
            res.convergence_iter = iter + 1;
            w_best = w;
            break;
        }
    }

    // Post-loop status
    if (res.status == RK_ERR_NOCONV) {
        res.status = (best_resid < initial_resid * 0.999) ? RK_ERR_BUDGET : RK_ERR_STALL;
        std::snprintf(res.message, sizeof(res.message),
            "logit: %s after %d Newton steps",
            res.status == RK_ERR_BUDGET ? "budget exhausted" : "stall",
            res.iterations);
    }

    // Final metrics
    if (w_best.empty()) w_best = w;
    double W_best = 0.0;
    for (int c = 0; c < M; c++) W_best += w_best[c];
    lbw::CellMetrics m_best = lbw::compute_cell_metrics(st, ct, w_best, W_best, bucket_scratch);
    res.max_error  = m_best.errRp;
    res.best_error = m_best.errRp;
    res.mean_error = m_best.mean_err;
    res.kl         = m_best.kl;
    res.chi2       = m_best.chi2;
    res.grake_norm = m_best.grake_norm;

    // Weight reconstruction
    res.best_weights.resize(st.n);
    for (int i = 0; i < st.n; i++) {
        int c = ct.cell_of[i];
        res.best_weights[i] = (X_init[c] > 0.0)
            ? st.weights[i] * w_best[c] / X_init[c]
            : st.weights[i];
    }

    return res;
}

} // namespace lbw
