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
    // Logit defaults differ from CalibResult — override here to preserve existing behavior.
    res.base.status                       = RK_ERR_NOCONV;
    res.base.convergence_metric           = static_cast<int>(CalibMetric::CHI2);
    res.base.convergence_rule             = 0;
    res.base.convergence_tol             = 0.0;
    res.base.convergence_iter             = 1;
    res.base.convergence_minimized_metric = static_cast<int>(CalibMetric::CHI2);
    res.base.convergence_solver_objective = std::numeric_limits<double>::infinity();
    res.base.best_iter                    = 1;

    // Input validation: logit with max_weight=Inf blows up
    if (!std::isfinite(st.max_weight) || st.max_weight <= 0.0) {
        res.base.status = RK_ERR_BADARG;
        std::snprintf(res.message, sizeof(res.message),
            "logit: max_weight must be finite and positive");
        return res;
    }
    if (st.min_weight >= st.max_weight) {
        res.base.status = RK_ERR_BADARG;
        std::snprintf(res.message, sizeof(res.message),
            "logit: min_weight (%.4g) >= max_weight (%.4g)",
            st.min_weight, st.max_weight);
        return res;
    }

    CellTable ct;
    std::vector<double> X_init;
    double hi_eff;
    std::vector<double> L_cell, U_cell;
    std::vector<int> cat_offset;
    int nct;
    if (lbw::solver_setup_ct(st, ct, X_init, hi_eff, L_cell, U_cell,
                              cat_offset, nct, res) != RK_OK)
        return res;

    const int M = ct.M_cell;
    const int K = st.K;

    // cells_per_cat[k][j] = list of cell indices in bucket (k,j)
    int max_cats = lbw::max_cats_count(K, st.cat_counts);
    auto cells_per_cat = lbw::build_cells_per_cat(ct, K, st.cat_counts);

    // Newton state
    std::vector<double> lambda(nct, 0.0);   // dual variables; initialized below from design-weight logit-inverse
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
    constexpr double kDeffFloor = 1e-6;  // floor prevents D_eff→0 when sig saturates
    constexpr double kInitSigmaEps = 1e-4;  // clips sigma_target to [eps, 1-eps] bounding z_target
    constexpr int    kMaxHalvings = 10;   // max Armijo halvings; alpha_min = 2^-10 ≈ 1e-3
    constexpr double kArmijoC     = 0.01; // Armijo sufficient-decrease constant
    constexpr double kMaxDeltaZ   = 2.0;  // max z-shift per Newton step (norm guard)
    constexpr double kArmijoHalving = 0.5;
    constexpr double kLambdaInitRejectAbs = 10.0;  // reject lambda_0 if any component exceeds this

    // Layer 2: design-weight initialization — place lambda_0 in convergence basin
    // Solve (AA^T)lambda_0 = Az_target where z_target[c] = logit(sigma_target[c])
    // and sigma_target[c] = clip((X_init[c]-L[c])/(U[c]-L[c]), eps, 1-eps)
    {
        std::vector<double> z_target(M, 0.0);
        for (int c = 0; c < M; c++) {
            double range = U_cell[c] - L_cell[c];
            if (range < 1e-12) { z_target[c] = 0.0; continue; }  // degenerate: L==U
            double sig0 = (X_init[c] - L_cell[c]) / range;
            sig0 = std::clamp(sig0, kInitSigmaEps, 1.0 - kInitSigmaEps);
            z_target[c] = std::log(sig0 / (1.0 - sig0));
        }
        // b_init[cat_offset[k]+j] = sum of z_target over bucket (k,j)
        std::vector<double> b_init(nct, 0.0);
        for (int k = 0; k < K; k++)
            for (int j = 0; j < st.cat_counts[k]; j++)
                for (int c : cells_per_cat[k][j])
                    b_init[cat_offset[k] + j] += z_target[c];
        // Normal equations with D_eff=1 (uniform: purely geometric initialization)
        std::vector<double> D_ones(M, 1.0);
        std::vector<double> N_init(nct * nct, 0.0);
        if (compute_normal_equations(ct, D_ones.data(), N_init.data(),
                                     cat_offset.data(), K, (size_t)nct) == RK_OK &&
            ldlt_factor_inplace(N_init.data(), (size_t)nct, 1e-10) == RK_OK) {
            ldlt_solve(N_init.data(), (size_t)nct, b_init.data());
            // Sanity: ill-conditioned AA^T (redundant margins) can produce huge lambda_0
            // which saturates all cells and makes the main Newton step degenerate.
            // Reject if any component exceeds 10 (z_c shifts >10*K would saturate sigma).
            double max_lambda_init = 0.0;
            for (double lj : b_init) max_lambda_init = std::max(max_lambda_init, std::abs(lj));
            if (max_lambda_init <= kLambdaInitRejectAbs) lambda = std::move(b_init);
        }
        // if factor fails or lambda_0 is ill-conditioned: lambda stays zero (Armijo handles it)
    }

    // Armijo scratch buffers (pre-allocated to avoid per-halving allocation)
    std::vector<double> w_trial(M, 0.0);        // trial weights for Armijo
    std::vector<double> b_trial(nct, 0.0);      // trial residuals for Armijo
    std::vector<double> lambda_trial(nct, 0.0); // trial lambda for Armijo

    for (int iter = 0; iter < kMaxNewtonIters; iter++) {
        res.base.iterations = iter + 1;

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
            D_eff[c] = std::max(kDeffFloor * range, range * sig * (1.0 - sig));
        }

        // Early-exit: if >50% of cells saturated (|z|>650), lambda is too large to recover
        if (iter == 0) {
            int n_saturated = 0;
            for (int c = 0; c < M; c++) {
                double z = 0.0;
                for (int k = 0; k < K; k++) {
                    int g = ct.g_per_cell[k][c];
                    if (g >= 0 && g < st.cat_counts[k]) z += lambda[cat_offset[k] + g];
                }
                if (std::fabs(z) > 650.0) ++n_saturated;
            }
            if (n_saturated > M / 2) {
                res.base.status = RK_ERR_NOCONV;
                std::snprintf(res.message, sizeof(res.message),
                    "logit: >50%% of cells saturated (|z|>650) — lambda too large; reduce bounds");
                return res;
            }
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
            res.base.best_iter = iter + 1;
        }

        // Capture ||b_current||² before ldlt_solve overwrites b with delta_lambda
        double resid_sq_0 = 0.0;
        double max_b_mag  = 0.0;
        for (double bj : b) { resid_sq_0 += bj * bj; max_b_mag = std::max(max_b_mag, std::abs(bj)); }

        // (3) Build N = A*diag(D_eff)*A^T and solve N*delta_lambda = b
        std::fill(N.begin(), N.end(), 0.0);
        if (compute_normal_equations(ct, D_eff.data(), N.data(),
                                      cat_offset.data(), K,
                                      static_cast<size_t>(nct)) != RK_OK) {
            res.base.status = RK_ERR_INFEAS;
            std::snprintf(res.message, sizeof(res.message),
                "logit: singular normal equations (degenerate bounds - L=U cells)");
            break;
        }
        // Adaptive eps: bound null-space amplification from rank-deficient AA^T*D
        // (overlapping margins can make N nearly singular even with D_eff > 0).
        // Choosing eps = max_b / kMaxDeltaZ ensures null-space components of delta_lambda
        // remain O(kMaxDeltaZ), keeping alpha_max bounded away from zero.
        double eps_ldlt = std::max(1e-10, max_b_mag / kMaxDeltaZ);
        if (ldlt_factor_inplace(N.data(), static_cast<size_t>(nct), eps_ldlt) != RK_OK) {
            res.base.status = RK_ERR_INFEAS;
            std::snprintf(res.message, sizeof(res.message),
                "logit: LDLT factorization failed (degenerate bounds)");
            break;
        }
        ldlt_solve(N.data(), static_cast<size_t>(nct), b.data());  // b = delta_lambda

        // (4) Armijo line search with step-norm guard
        // Norm guard: cap alpha so no cell z-coord shifts more than kMaxDeltaZ
        // (b = delta_lambda after ldlt_solve)
        double max_delta_z = 0.0;
        for (int c = 0; c < M; c++) {
            double dz = 0.0;
            for (int k = 0; k < K; k++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k]) dz += std::abs(b[cat_offset[k] + g]);
            }
            max_delta_z = std::max(max_delta_z, dz);
        }
        double alpha = (max_delta_z > 0.0) ? std::min(1.0, kMaxDeltaZ / max_delta_z) : 1.0;

        bool armijo_improved = false;
        for (int halv = 0; halv < kMaxHalvings; halv++) {
            for (int j = 0; j < nct; j++) lambda_trial[j] = lambda[j] + alpha * b[j];
            // Recompute w_trial from lambda_trial (logit link)
            for (int c = 0; c < M; c++) {
                double z = 0.0;
                for (int k = 0; k < K; k++) {
                    int g = ct.g_per_cell[k][c];
                    if (g >= 0 && g < st.cat_counts[k]) z += lambda_trial[cat_offset[k] + g];
                }
                z = std::clamp(z, -700.0, 700.0);
                double sig = 1.0 / (1.0 + std::exp(-z));
                double range = U_cell[c] - L_cell[c];
                w_trial[c] = L_cell[c] + range * sig;
            }
            // Compute b_trial residuals from w_trial
            std::fill(b_trial.begin(), b_trial.end(), 0.0);
            for (int k = 0; k < K; k++)
                for (int j = 0; j < st.cat_counts[k]; j++) {
                    double target = st.targets[k][j] * static_cast<double>(st.n);
                    double S_kj = 0.0;
                    for (int c : cells_per_cat[k][j]) S_kj += w_trial[c];
                    b_trial[cat_offset[k] + j] = target - S_kj;
                }
            double resid_sq_trial = 0.0;
            for (double bj : b_trial) resid_sq_trial += bj * bj;
            if (resid_sq_trial < resid_sq_0 * (1.0 - kArmijoC * alpha)) {
                armijo_improved = true;
                break;
            }
            alpha *= kArmijoHalving;
        }
        if (!armijo_improved && st.verbose >= 1) {
            Rprintf("[logit] Newton step: Armijo exhausted (alpha=%.2e), accepting best available\n",
                alpha);
        }
        for (int j = 0; j < nct; j++) lambda[j] += alpha * b[j];

        // Recompute w from updated lambda so convergence check sees post-step weights
        for (int c = 0; c < M; c++) {
            double z = 0.0;
            for (int k = 0; k < K; k++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k]) z += lambda[cat_offset[k] + g];
            }
            z = std::clamp(z, -700.0, 700.0);
            double sig = 1.0 / (1.0 + std::exp(-z));
            w[c] = L_cell[c] + (U_cell[c] - L_cell[c]) * sig;
        }

        // (5) Convergence check via shared infrastructure
        double W_total = 0.0;
        for (int c = 0; c < M; c++) W_total += w[c];
        lbw::CellMetrics m = lbw::compute_cell_metrics(st, ct, w, W_total, bucket_scratch);
        bool converged = lbw::check_convergence(cfg, m, prev_metric, st.tol_abs);
        if (converged) {
            lbw::mark_converged(res, cfg, iter + 1);
            w_best = w;
            break;
        }
    }

    // Post-loop status
    if (res.base.status == RK_ERR_NOCONV) {
        res.base.status = (best_resid < initial_resid * 0.999) ? RK_ERR_BUDGET : RK_ERR_STALL;
        std::snprintf(res.message, sizeof(res.message),
            "logit: %s after %d Newton steps",
            res.base.status == RK_ERR_BUDGET ? "budget exhausted" : "stall",
            res.base.iterations);
    }

    // Final metrics
    if (w_best.empty()) w_best = w;
    double W_best = 0.0;
    for (int c = 0; c < M; c++) W_best += w_best[c];
    lbw::CellMetrics m_best = lbw::compute_cell_metrics(st, ct, w_best, W_best, bucket_scratch);
    res.base.max_error  = m_best.errRp;
    res.base.best_error = m_best.errRp;
    res.base.mean_error = m_best.mean_err;
    res.base.kl         = m_best.kl;
    res.base.chi2       = m_best.chi2;
    res.base.grake_norm = m_best.grake_norm;

    // Weight reconstruction
    res.base.best_weights.resize(st.n);
    for (int i = 0; i < st.n; i++) {
        int c = ct.cell_of[i];
        res.base.best_weights[i] = (X_init[c] > 0.0)
            ? st.weights[i] * w_best[c] / X_init[c]
            : st.weights[i];
    }

    return res;
}

} // namespace lbw
