#include "newton_calib.hpp"
#include "leafblower.h"
#include "calib_dispatch.hpp"
#include "calib_linalg.hpp"
#include "cell_table.hpp"
#include "calib_validate.hpp"
#include <cmath>
#include <algorithm>
#include <limits>
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

    // LM damping state. Init 1e-3 (mildly damped); adaptive μ via Marquardt
    // gain ratio: ρ>0.75 → μ/=3 (model good), ρ<0.25 → μ×10. Init 1.0 was too
    // damped on well-conditioned (stepstone K=9): max_err stuck at ~1.0 after
    // 33 iters — μ couldn't relax fast enough. 1e-3 lets Newton step at full
    // length on iter 0 when H is well-conditioned.
    const double lm_mu_min = 1e-12;
    const double lm_mu_max = 1e12;
    double lm_mu = 1e-3;

    // ── 4. Dual objective with log-sum-exp stabilization ───────────────────
    // g(λ) = log Z(λ) - T·λ where Z(λ) = Σ d_i exp(u_i). For K=20 with skewed
    // targets |u_i| reaches 60+ → raw exp(u_i) overflows/underflows, Z→0/∞,
    // false-convergence (NaN compares short-circuit gap=0). Stable form:
    // u_max = max_i u_i; Z = exp(u_max)·Σ d_i exp(u_i − u_max);
    // log Z = u_max + log Σ d_i exp(u_i − u_max).
    // All ratios (G/Z, H/Z, w_i = exp(u_i)/Z·n) are u_max-shift invariant.
    auto compute_u = [&](const std::vector<double>& lam_v, int i) -> double {
        double u = 0.0;
        for (int k = 0; k < K; ++k) {
            const int j = st.group_ids[k][i];
            if (j > 0) u += lam_v[lam_off[k] + j - 1];
        }
        return u;
    };
    auto compute_u_max = [&](const std::vector<double>& lam_v) -> double {
        double u_max = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < n; ++i) {
            const double u = compute_u(lam_v, i);
            if (u > u_max) u_max = u;
        }
        return u_max;
    };
    // eval_dual: returns g(λ); Z_out = Σ d_i exp(u_i − u_max); u_max_out = max u.
    auto eval_dual = [&](const std::vector<double>& lam_v,
                         double& Z_stable_out, double& u_max_out) -> double {
        const double u_max = compute_u_max(lam_v);
        double Z_stable = 0.0;
        for (int i = 0; i < n; ++i) {
            Z_stable += st.weights[i] * std::exp(compute_u(lam_v, i) - u_max);
        }
        Z_stable_out = Z_stable;
        u_max_out    = u_max;
        double Tlam = 0.0;
        for (int a = 0; a < n_lam; ++a) Tlam += T[a] * lam_v[a];
        return u_max + std::log(Z_stable) - Tlam;
    };

    double Z_curr = 0.0;     // = Z_stable at current λ (u_max factored out)
    double u_max_curr = 0.0; // current λ's u_max for stable recovery in step 6
    double g_curr = eval_dual(lam, Z_curr, u_max_curr);

    // ── 5. Newton iterations ─────────────────────────────────────────────────
    // run_newton_inner: encapsulates the per-ε inner Newton loop.
    //   T_eps       — target vector to converge against (pass current T for now)
    //   max_iter_inner — iteration budget for this inner call
    // Captures outer-scope state by reference: lam, lm_mu, Z_curr, u_max_curr,
    //   res, H, G, delta, n_lam, n, K, lam_off, st, g_curr, tol_abs, lm_mu_min,
    //   lm_mu_max, compute_u, compute_u_max, eval_dual.
    // Returns true  on converged-or-step-stalled exit (res.base.status == RK_OK).
    // Returns false on BADARG/NOCONV (res.base.status >= 1; caller halts outer).
    // Note: the BADARG (LDLT 3× fail) exit path now restores best-iterate λ
    // before returning, matching the NOCONV path. Pre-refactor returned with
    // unrestored λ on this path; no test exercises it, but the change brings
    // the two failure paths into structural parity.
    // Writes res.base.iterations to LOCAL iter count; outer caller reads it
    // BEFORE the next inner call (which would overwrite it).
    auto run_newton_inner = [&](const std::vector<double>& T_eps,
                                int max_iter_inner) -> bool {
        int consecutive_failed = 0;  // consecutive line-search failures; cap at 3
        // Best-iterate tracking: rank-deficient near-optimum can drift λ to ∞ even
        // while Armijo accepts (model under-predicts decrease). Save λ at lowest
        // gap seen; restore before recovery if the loop exits in a worse state.
        std::vector<double> lam_best(lam);
        double best_gap     = std::numeric_limits<double>::infinity();
        double best_Z       = Z_curr;
        double best_u_max   = u_max_curr;
        int    best_iter_id = 0;
        int iter = 0;
        for (iter = 0; iter < max_iter_inner; ++iter) {

            // 5a. Zero accumulators for this step.
            std::fill(G.begin(), G.end(), 0.0);
            std::fill(H.begin(), H.end(), 0.0);
            double Z = 0.0;

            // 5b. LSE prep: scan u_max once, then accumulate using exp(u_i − u_max).
            // u_max factor cancels in all ratios (G/Z, H/Z) so this is a pure
            // numerical-stability shift. Cost: one extra O(n·K) pass per Newton iter.
            const double u_max_step = compute_u_max(lam);
            for (int i = 0; i < n; ++i) {
                const double u = compute_u(lam, i);
                const double f_i = st.weights[i] * std::exp(u - u_max_step);
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
            Z_curr     = Z;            // = Σ d_i exp(u_i − u_max_step)
            u_max_curr = u_max_step;   // shift used during this step's accumulation

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

            // 5e2. Mean diagonal — additive floor so damping is effective even when
            //      H has rank-collapsed coordinates (pure μ·diag would be zero there).
            double d_floor = 0.0;
            for (int a = 0; a < n_lam; ++a) d_floor += H[a * n_lam + a];
            d_floor /= std::max(n_lam, 1);

            // 5f. Subtract targets to get gradient ∇g = G - T_eps.
            for (int a = 0; a < n_lam; ++a) G[a] -= T_eps[a];

            // 5g. Convergence check: ||∇g||_∞ < tol_abs.
            double dual_gap = 0.0;
            for (int a = 0; a < n_lam; ++a)
                dual_gap = std::max(dual_gap, std::fabs(G[a]));
            res.dual_gap = dual_gap;

            // Track best λ (lowest gap seen) for end-of-loop fallback recovery.
            if (dual_gap < best_gap) {
                best_gap     = dual_gap;
                best_Z       = Z;
                best_u_max   = u_max_step;
                best_iter_id = iter;
                lam_best     = lam;
            }

            if (dual_gap < tol_abs) {
                res.lm_mu_final = lm_mu;
                mark_converged(res, st.convergence_cfg, iter);
                break;
            }

            // 5h. LM damped Newton direction: H_damp·δ = G, λ -= α·δ.
            // Save pre-damp H to compute δᵀHδ for the Marquardt gain ratio.
            std::vector<double> H_pre(H);  // n_lam×n_lam copy — ≤80×80 doubles = trivial.

            bool solve_ok = false;
            for (int retry = 0; retry < 3 && !solve_ok; ++retry) {
                // Restore H (LDLT factored in-place, so retry needs the original).
                std::copy(H_pre.begin(), H_pre.end(), H.begin());
                // Scale-invariant LM damping with additive floor — handles rank collapse.
                for (int a = 0; a < n_lam; ++a) {
                    H[a * n_lam + a] = std::max(H[a * n_lam + a] * (1.0 + lm_mu),
                                                 lm_mu * d_floor);
                }
                std::copy(G.begin(), G.end(), delta.begin());
                if (ldlt_factor_inplace(H.data(), static_cast<size_t>(n_lam), 1e-12) != RK_OK) {
                    lm_mu = std::min(lm_mu_max, lm_mu * 10.0);
                    continue;
                }
                ldlt_solve(H.data(), static_cast<size_t>(n_lam), delta.data());
                solve_ok = true;
            }
            if (!solve_ok) {
                res.base.status = RK_ERR_BADARG;
                std::snprintf(res.message, sizeof(res.message),
                    "newton_kl: LDLT failed 3x even at lm_mu=%.2e", lm_mu);
                res.lm_mu_final = lm_mu;
                // Best-iterate restoration before early exit.
                if (best_gap < res.dual_gap * 0.99) {
                    lam        = lam_best;
                    Z_curr     = best_Z;
                    u_max_curr = best_u_max;
                    res.dual_gap        = best_gap;
                    res.base.iterations = best_iter_id + 1;
                } else {
                    res.base.iterations = iter + 1;
                }
                return false;
            }

            // δᵀ H_pre δ using pre-damp (PSD) Hessian for the gain-ratio denominator.
            // Clamp ≥ 0: under severe rank deficiency (severe-skew K=20, many empty
            // cells), Schur-complement subtraction can produce tiny-negative diagonal
            // due to FP rounding. PSD is the mathematical truth; clamping avoids
            // contaminating ρ with spurious numerical violations.
            double delta_H_delta = 0.0;
            for (int a = 0; a < n_lam; ++a)
                for (int b = 0; b < n_lam; ++b)
                    delta_H_delta += delta[a] * H_pre[a * n_lam + b] * delta[b];
            if (delta_H_delta < 0.0) delta_H_delta = 0.0;

            // G·δ = G^T H_damp^{-1} G ≥ 0 (descent curvature).
            double g_dot_d = 0.0;
            for (int a = 0; a < n_lam; ++a) g_dot_d += G[a] * delta[a];

            // 5i. Armijo backtracking line search.
            // Sufficient decrease: g(λ - α·δ) ≤ g(λ) - c1·α·(G·δ).
            double alpha = 1.0;
            std::vector<double> lam_trial(n_lam);
            double Z_trial = 0.0, g_trial = 0.0, u_max_trial = 0.0;
            constexpr double c1 = 1e-4;
            bool accepted = false;
            for (int ls = 0; ls < 30; ++ls) {
                for (int a = 0; a < n_lam; ++a) lam_trial[a] = lam[a] - alpha * delta[a];
                g_trial = eval_dual(lam_trial, Z_trial, u_max_trial);
                if (g_trial <= g_curr - c1 * alpha * g_dot_d) { accepted = true; break; }
                alpha *= 0.5;
            }
            res.line_alpha = alpha;

            // 5j. Three-way Armijo outcome.
            if (!accepted) {
                // FAILED — line search exhausted; back off μ, do NOT advance λ.
                lm_mu = std::min(lm_mu_max, lm_mu * 10.0);
                if (++consecutive_failed >= 3) {
                    res.base.status = RK_ERR_NOCONV;
                    res.lm_mu_final = lm_mu;
                    break;
                }
                continue;  // re-enter iter with new lm_mu, same λ
            }
            consecutive_failed = 0;

            // Marquardt gain ratio: actual / predicted reduction.
            // Guard denominator — predicted can be tiny/negative near singularity.
            const double predicted = alpha * g_dot_d
                                     - 0.5 * alpha * alpha * delta_H_delta;
            const double rho = (g_curr - g_trial) / std::max(predicted, 1e-300);

            if (alpha >= 0.999 && rho > 0.75) {
                // FULL_ACCEPT with good model quality — tighten damping.
                lm_mu = std::max(lm_mu_min, lm_mu / 3.0);
            } else if (rho < 0.25) {
                // Poor quadratic model (accepted Armijo but low gain) — back off.
                lm_mu = std::min(lm_mu_max, lm_mu * 10.0);
            }
            // else: BACKTRACKED or moderate ρ — keep lm_mu unchanged.

            // 5k. Update λ and dual objective.
            double step2 = 0.0;
            for (int a = 0; a < n_lam; ++a) {
                const double d = alpha * delta[a];
                lam[a] -= d;
                step2 += d * d;
            }
            res.step_norm = std::sqrt(step2);
            g_curr     = g_trial;
            Z_curr     = Z_trial;
            u_max_curr = u_max_trial;

            // Secondary stop: step too small to make progress.
            if (res.step_norm < 1e-14) {
                res.lm_mu_final = lm_mu;
                mark_converged(res, st.convergence_cfg, iter);
                break;
            }
        }
        res.base.iterations = iter + 1;

        // ── 5l. Best-iterate restoration ──────────────────────────────────────
        // If the loop's final λ has a worse gap than the best one seen, restore.
        // Defensive against drift: rank-deficient Hessian + accepted-by-Armijo can
        // walk λ past the optimum into unbounded-dual regions; gap stops decreasing
        // but g keeps decreasing (descent direction is wrong but still descent).
        // Compare conservatively (factor 1.01) so a stale best at iter 0 doesn't
        // overrule a slightly-noisy iter-N+ result.
        if (best_gap < res.dual_gap * 0.99) {
            lam        = lam_best;
            Z_curr     = best_Z;
            u_max_curr = best_u_max;
            res.dual_gap         = best_gap;
            res.base.iterations  = best_iter_id + 1;
        }

        res.lm_mu_final = lm_mu;
        return res.base.status == RK_OK;
    };

    run_newton_inner(T, max_iter);

    // ── 6. Recover obs weights via stable form ─────────────────────────────
    // w_i = d_i · exp(u_i) / Z_real · n
    //     = d_i · exp(u_i − u_max) / Z_stable · n   (u_max factor cancels)
    // Z_curr and u_max_curr are the LSE-stable values from the last accepted
    // iter (set in step 5b for the convergence-break path, or in step 5j
    // after Armijo accept). They correspond to the current `lam`.
    std::vector<double> w(n);
    const double scale = static_cast<double>(n) / Z_curr;
    int n_violated = 0;
    for (int i = 0; i < n; ++i) {
        const double u = compute_u(lam, i);
        const double w_i = st.weights[i] * std::exp(u - u_max_curr) * scale;
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
    res.lm_mu_final       = lm_mu;

    return res;
}

} // namespace lbw
