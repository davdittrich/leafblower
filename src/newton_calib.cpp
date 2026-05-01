#include "newton_calib.hpp"
#include "leafblower.h"
#include "calib_dispatch.hpp"
#include "calib_linalg.hpp"
#include "cell_table.hpp"
#include "calib_validate.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cstring>
#include <cstdio>
#include <vector>

namespace lbw {

NewtonCalibResult newton_calibrate(CalibState& st) {
    NewtonCalibResult res;
    res.base.status = RK_ERR_NOCONV;

    // ── 1. Preamble: cell-table build + validation ──────────────────────────
    CellTable ct;
    std::vector<double> X_init, L_cell, U_cell;
    std::vector<int> cat_offset;
    double hi_eff = 0.0;
    int n_cats_total = 0;
    if (solver_setup_ct(st, ct, X_init, hi_eff, L_cell, U_cell,
                        cat_offset, n_cats_total, res) != RK_OK) {
        return res;
    }

    const int n = st.n;
    const int K = st.K;

    // ── 2. Reference-category elimination ───────────────────────────────────
    // lam_off[k] = Σ_{k'<k} (cat_counts[k'] - 1)
    // Free dual variables: j=1..cat_counts[k]-1 per margin k.
    // Category j=0 is the reference (λ fixed to 0).
    std::vector<int> lam_off(K + 1, 0);
    for (int k = 0; k < K; ++k)
        lam_off[k + 1] = lam_off[k] + std::max(0, st.cat_counts[k] - 1);
    const int n_lam = lam_off[K];
    res.n_lambda = n_lam;

    if (n_lam <= 0) {
        // All margins have only 1 category — trivially calibrated.
        res.base.status     = RK_OK;
        res.base.max_error  = 0.0;
        res.base.iterations = 0;
        return res;
    }
    if (n_lam > 4096) {
        res.base.status = RK_ERR_BADARG;
        std::snprintf(res.message, sizeof(res.message),
            "newton_kl: n_lambda=%d exceeds 4096 (sum of free categories)", n_lam);
        return res;
    }

    // ── 3. State allocation ──────────────────────────────────────────────────
    std::vector<double> lam(n_lam, 0.0);
    std::vector<double> H(static_cast<size_t>(n_lam) * n_lam, 0.0);
    std::vector<double> G(n_lam, 0.0);
    std::vector<double> delta(n_lam, 0.0);
    // Targets vector T[a] = st.targets[k][j] for a = lam_off[k] + j - 1, j>=1
    std::vector<double> T(n_lam, 0.0);
    for (int k = 0; k < K; ++k)
        for (int j = 1; j < st.cat_counts[k]; ++j)
            T[lam_off[k] + j - 1] = st.targets[k][j];

    const int max_iter = (st.outer_max_iter > 0) ? st.outer_max_iter : 50;
    const double tol_abs = (st.tol_abs > 0.0) ? st.tol_abs : 1e-6;

    // ── 4. Dual objective: g(λ) = log Z(λ) - Σ_a T[a]*λ[a] ─────────────────
    // g is convex; Newton minimizes it. At optimum ∇g = 0 ↔ margins calibrated.
    auto eval_dual = [&](const std::vector<double>& lam_v, double& Z_out) -> double {
        double Z = 0.0;
        for (int i = 0; i < n; ++i) {
            double u = 0.0;
            for (int k = 0; k < K; ++k) {
                const int j = st.group_ids[k][i];
                if (j > 0) u += lam_v[lam_off[k] + j - 1];
            }
            Z += st.weights[i] * std::exp(u);
        }
        Z_out = Z;
        double Tlam = 0.0;
        for (int a = 0; a < n_lam; ++a) Tlam += T[a] * lam_v[a];
        return std::log(Z) - Tlam;
    };

    double Z_curr = 0.0;
    double g_curr = eval_dual(lam, Z_curr);

    // ── 5. Newton iterations ─────────────────────────────────────────────────
    int iter = 0;
    for (iter = 0; iter < max_iter; ++iter) {

        // 5a. Zero accumulators for this step.
        std::fill(G.begin(), G.end(), 0.0);
        std::fill(H.begin(), H.end(), 0.0);
        double Z = 0.0;

        // 5b. Single obs-level pass: accumulate Z, G, H (upper triangle).
        for (int i = 0; i < n; ++i) {
            double u = 0.0;
            for (int k = 0; k < K; ++k) {
                const int j = st.group_ids[k][i];
                if (j > 0) u += lam[lam_off[k] + j - 1];
            }
            const double f_i = st.weights[i] * std::exp(u);
            Z += f_i;

            // Gradient accumulation
            for (int k = 0; k < K; ++k) {
                const int j = st.group_ids[k][i];
                if (j > 0) G[lam_off[k] + j - 1] += f_i;
            }

            // Hessian accumulation (upper triangle only)
            for (int k1 = 0; k1 < K; ++k1) {
                const int j1 = st.group_ids[k1][i];
                if (j1 <= 0) continue;
                const int a = lam_off[k1] + j1 - 1;
                for (int k2 = k1; k2 < K; ++k2) {
                    const int j2 = st.group_ids[k2][i];
                    if (j2 <= 0) continue;
                    const int b = lam_off[k2] + j2 - 1;
                    H[a * n_lam + b] += f_i;
                }
            }
        }
        Z_curr = Z;

        // 5c. Normalize G and H by Z.
        const double inv_Z = 1.0 / Z;
        for (int a = 0; a < n_lam; ++a) G[a] *= inv_Z;
        for (size_t idx = 0; idx < H.size(); ++idx) H[idx] *= inv_Z;

        // 5d. Subtract outer product G⊗G from H (Schur complement, upper triangle).
        for (int a = 0; a < n_lam; ++a)
            for (int b = a; b < n_lam; ++b)
                H[a * n_lam + b] -= G[a] * G[b];

        // 5e. Mirror upper triangle to lower.
        for (int a = 0; a < n_lam; ++a)
            for (int b = a + 1; b < n_lam; ++b)
                H[b * n_lam + a] = H[a * n_lam + b];

        // 5f. Subtract targets to get gradient ∇g = G - T.
        for (int a = 0; a < n_lam; ++a) G[a] -= T[a];

        // 5g. Convergence check: ||∇g||_∞ < tol_abs.
        double dual_gap = 0.0;
        for (int a = 0; a < n_lam; ++a)
            dual_gap = std::max(dual_gap, std::fabs(G[a]));
        res.dual_gap = dual_gap;

        if (dual_gap < tol_abs) {
            mark_converged(res, st.convergence_cfg, iter);
            break;
        }

        // 5h. Newton direction: solve H·δ = G (descent direction: λ -= α·δ).
        std::copy(G.begin(), G.end(), delta.begin());
        if (ldlt_factor_inplace(H.data(), static_cast<size_t>(n_lam), 1e-12) != RK_OK) {
            res.base.status = RK_ERR_BADARG;
            std::snprintf(res.message, sizeof(res.message),
                "newton_kl: LDLT factorization failed at iter %d", iter);
            return res;
        }
        ldlt_solve(H.data(), static_cast<size_t>(n_lam), delta.data());

        // 5i. Armijo backtracking line search.
        // g is convex; we minimize, so descent direction is -δ (λ -= α·δ).
        // Sufficient decrease: g(λ - α·δ) ≤ g(λ) - c1·α·(G·δ).
        // G·δ = G^T H^{-1} G ≥ 0 (H PSD), so this is a reduction condition.
        double g_dot_d = 0.0;
        for (int a = 0; a < n_lam; ++a) g_dot_d += G[a] * delta[a];

        double alpha = 1.0;
        std::vector<double> lam_trial(n_lam);
        double Z_trial = 0.0, g_trial = 0.0;
        constexpr double c1 = 1e-4;
        bool accepted = false;
        for (int ls = 0; ls < 30; ++ls) {
            for (int a = 0; a < n_lam; ++a) lam_trial[a] = lam[a] - alpha * delta[a];
            g_trial = eval_dual(lam_trial, Z_trial);
            if (g_trial <= g_curr - c1 * alpha * g_dot_d) {
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }
        if (!accepted) {
            // Line search exhausted — take a tiny step anyway to avoid stall
            alpha = 1.0 / (1 << 30);
            for (int a = 0; a < n_lam; ++a) lam_trial[a] = lam[a] - alpha * delta[a];
            g_trial = eval_dual(lam_trial, Z_trial);
        }
        res.line_alpha = alpha;

        // 5j. Update λ and dual objective.
        double step2 = 0.0;
        for (int a = 0; a < n_lam; ++a) {
            const double d = alpha * delta[a];
            lam[a] -= d;
            step2 += d * d;
        }
        res.step_norm = std::sqrt(step2);
        g_curr  = g_trial;
        Z_curr  = Z_trial;

        // Secondary stop: step too small to make progress.
        if (res.step_norm < 1e-14) {
            mark_converged(res, st.convergence_cfg, iter);
            break;
        }
    }
    res.base.iterations = iter + 1;

    // ── 6. Recover obs weights: w_i = weights[i] * exp(u_i) / Z * n ────────
    std::vector<double> w(n);
    const double scale = static_cast<double>(n) / Z_curr;
    int n_violated = 0;
    for (int i = 0; i < n; ++i) {
        double u = 0.0;
        for (int k = 0; k < K; ++k) {
            const int j = st.group_ids[k][i];
            if (j > 0) u += lam[lam_off[k] + j - 1];
        }
        const double w_i = st.weights[i] * std::exp(u) * scale;
        w[i] = w_i;
        if (w_i > st.max_weight || w_i < st.min_weight) ++n_violated;
    }

    // 5% bounds-violation fallback.
    const double frac_violated = static_cast<double>(n_violated) / n;
    if (frac_violated > 0.05) {
        res.base.status = RK_ERR_NOCONV;
        std::snprintf(res.message, sizeof(res.message),
            "newton_kl: %.1f%% of obs violate [min,max] weight bounds",
            frac_violated * 100.0);
    }

    // ── 7. Write calibrated weights and compute max_error ───────────────────
    for (int i = 0; i < n; ++i) st.weights[i] = w[i];

    double max_err = 0.0;
    for (int k = 0; k < K; ++k) {
        double total = 0.0;
        std::vector<double> achieved(st.cat_counts[k], 0.0);
        for (int i = 0; i < n; ++i) {
            const int j = st.group_ids[k][i];
            if (j >= 0 && j < st.cat_counts[k]) {
                achieved[j] += w[i];
                total        += w[i];
            }
        }
        if (total <= 0.0) continue;
        const double inv_tot = 1.0 / total;
        for (int j = 0; j < st.cat_counts[k]; ++j) {
            const double err = std::fabs(achieved[j] * inv_tot - st.targets[k][j]);
            if (err > max_err) max_err = err;
        }
    }
    res.base.max_error    = max_err;
    res.base.best_weights = w;

    return res;
}

} // namespace lbw
