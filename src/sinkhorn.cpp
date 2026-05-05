#include "sinkhorn.hpp"
#include "lbw_math.hpp"
#include "calib_dispatch.hpp"
#include "calib_validate.hpp"
#include "leafblower.h"
#include <cmath>
#include <cstring>
#include <algorithm>
#include <limits>

namespace lbw {

// Bisect on log-scale multiplier μ to project X onto capacity box preserving KL geometry.
// Finds μ s.t. Σ_c clamp(X[c]*exp(a[c]+μ), L[c], U[c]) = target_mass.
// f(μ) = Σ clamp(X*exp(a+μ), L, U) - target is strictly increasing → bisection valid.
static bool bisect_capacity(const std::vector<double>& X,
                             const std::vector<double>& a,
                             const std::vector<double>& L,
                             const std::vector<double>& U,
                             int M_cell,
                             double target_mass,
                             double& mu_out,
                             std::vector<double>& X_proj)
{
    double sum_L = 0.0, sum_U = 0.0;
    for (int c = 0; c < M_cell; c++) { sum_L += L[c]; sum_U += U[c]; }
    if (sum_L > target_mass + 1e-9 || sum_U < target_mass - 1e-9) return false;

    auto f = [&](double mu) -> double {
        double s = 0.0;
        for (int c = 0; c < M_cell; c++)
            s += std::clamp(X[c] * std::exp(a[c] + mu), L[c], U[c]);
        return s - target_mass;
    };

    // Data-derived bracket: guarantees f(lo) <= 0 and f(hi) >= 0.
    // f(mu) = sum_c clamp(X[c]*exp(a[c]+mu), L[c], U[c]) - target_mass.
    // Each c is fully clamped to L[c] when mu <= log(L[c])-log(X[c])-a[c],
    // and fully clamped to U[c] when mu >= log(U[c])-log(X[c])-a[c].
    // Taking lo = min over c of the lower thresholds (minus margin) clamps
    // every term to L[c] => f(lo) = sum_L - target <= 0 (caller ensured
    // sum_L <= target).  Symmetric for hi.  Replaces the broken
    // lo*=2.0 expansion (doubling a negative makes it MORE negative,
    // and the loop capped at -800 truncates problems needing mu < -800).
    constexpr double kEpsLU = 1e-300;
    double lo = std::numeric_limits<double>::infinity();
    double hi = -std::numeric_limits<double>::infinity();
    for (int c = 0; c < M_cell; c++) {
        const double Xc = std::max(X[c], kEpsLU);
        const double lo_c = std::log(std::max(L[c], kEpsLU)) - std::log(Xc) - a[c];
        const double hi_c = std::log(std::max(U[c], kEpsLU)) - std::log(Xc) - a[c];
        if (lo_c < lo) lo = lo_c;
        if (hi_c > hi) hi = hi_c;
    }
    lo -= 1.0;  // safety margin
    hi += 1.0;
    if (!std::isfinite(lo) || !std::isfinite(hi) || lo >= hi) return false;
    if (f(lo) > 0.0 || f(hi) < 0.0) return false;  // bracket failed — infeasible
    for (int i = 0; i < 80; i++) {
        double mid = 0.5 * (lo + hi);
        if (f(mid) < 0.0) lo = mid; else hi = mid;
        if (hi - lo < 1e-12) break;
    }
    mu_out = 0.5 * (lo + hi);
    for (int c = 0; c < M_cell; c++)
        X_proj[c] = std::clamp(X[c] * std::exp(a[c] + mu_out), L[c], U[c]);
    return true;
}

// Fast variant: accepts pre-computed exp_a[c]=exp(a[c]) so bisection eval is
// exp_a[c]*exp(mu) — one scalar exp per step instead of M_cell vector exps.
static bool bisect_capacity_fast(
    const std::vector<double>& X,
    const double* exp_a_data,
    const std::vector<double>& L,
    const std::vector<double>& U,
    int M_cell, double target_mass,
    double& mu_out, std::vector<double>& X_proj)
{
    double sum_L = 0.0, sum_U = 0.0;
    for (int c = 0; c < M_cell; c++) { sum_L += L[c]; sum_U += U[c]; }
    if (sum_L > target_mass + 1e-9 || sum_U < target_mass - 1e-9) return false;

    auto f = [&](double mu_val) -> double {
        const double exp_mu = std::exp(mu_val);
        double s = 0.0;
        for (int c = 0; c < M_cell; c++)
            s += std::clamp(X[c] * exp_a_data[c] * exp_mu, L[c], U[c]);
        return s - target_mass;
    };

    // Data-derived bracket: in this fast path the per-cell coefficient is
    // X[c]*exp_a_data[c]; clamp threshold is mu = log(L_or_U[c]) - log(coef).
    // See bisect_capacity for derivation; same fix replaces the broken
    // lo*=2.0 negative-doubling expansion.
    constexpr double kEpsLU = 1e-300;
    double lo = std::numeric_limits<double>::infinity();
    double hi = -std::numeric_limits<double>::infinity();
    for (int c = 0; c < M_cell; c++) {
        const double coef = std::max(X[c] * exp_a_data[c], kEpsLU);
        const double lo_c = std::log(std::max(L[c], kEpsLU)) - std::log(coef);
        const double hi_c = std::log(std::max(U[c], kEpsLU)) - std::log(coef);
        if (lo_c < lo) lo = lo_c;
        if (hi_c > hi) hi = hi_c;
    }
    lo -= 1.0;
    hi += 1.0;
    if (!std::isfinite(lo) || !std::isfinite(hi) || lo >= hi) return false;
    if (f(lo) > 0.0 || f(hi) < 0.0) return false;
    for (int i = 0; i < 80; i++) {
        double mid = 0.5 * (lo + hi);
        if (f(mid) < 0.0) lo = mid; else hi = mid;
        if (hi - lo < 1e-12) break;
    }
    mu_out = 0.5 * (lo + hi);
    const double exp_mu_out = std::exp(mu_out);
    for (int c = 0; c < M_cell; c++)
        X_proj[c] = std::clamp(X[c] * exp_a_data[c] * exp_mu_out, L[c], U[c]);
    return true;
}

SinkhornResult sinkhorn_solve(CalibState& st) {
    static constexpr int    kErrCheckInterval = 10;
    static constexpr double kAmax             = 30.0;  // exp(30)≈1e13 >> max practical weight ratio

    SinkhornResult res;
    res.base.status               = RK_ERR_NOCONV;
    res.base.convergence_rule     = 0;    // sinkhorn: no improvement criterion
    res.base.convergence_tol      = 0.001;
    res.base.convergence_solver_objective = std::numeric_limits<double>::infinity();

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
    std::vector<double> X(X_init);  // working copy; X_init kept for obs-expansion
    const double lo = st.min_weight;
    const double hi = hi_eff;

    // log-domain Dykstra correction for capacity box
    std::vector<double> a(ct.M_cell, 0.0);
    // exp_a[c] = exp(a[c]); a init=0 → exp=1. Updated via bulk_scaled_exp after each a[] write.
    std::vector<double> exp_a(ct.M_cell, 1.0);

    int max_cats = lbw::max_cats_count(st.K, st.cat_counts);

    // G8c: best-iterate tracking via BestIterTracker (replaces ad-hoc vars).
    // Fixes 0.0 sentinel bug: best_objective_seen was init'd to 0.0, making
    // "no best recorded" indistinguishable from "objective == 0".
    BestIterTracker best;

    // Scratch buffers for compute_weight_kl.
    std::vector<double> kl_ratio_buf(ct.M_cell);
    std::vector<double> kl_weight_buf(ct.M_cell);

    std::vector<double> bucket(max_cats);
    std::vector<double> scale(max_cats);
    std::vector<double> X_proj(ct.M_cell);
    std::vector<double> X_prev(X);
    double prev_metric_for_rule = std::numeric_limits<double>::infinity();
    // j98p.1: track whether Dykstra accumulator a[] is still all-zero.
    // Short-circuit is only safe when a[] is zero; non-zero a[] means bisection
    // must run even when X is in bounds to maintain the fixed-point invariant.
    bool a_is_zero = true;

    for (int iter = 1; iter <= st.inner_max_iter; iter++) {
        res.base.iterations = iter;

        double W_total = 0.0;
        for (int c = 0; c < ct.M_cell; c++) W_total += X[c];

        for (int k = 0; k < st.K; k++) {
            lbw::aggregate_to_margin(ct, X, k, st.cat_counts[k], bucket.data());
            std::fill(scale.begin(), scale.begin() + st.cat_counts[k], 1.0);
            for (int j = 0; j < st.cat_counts[k]; j++) {
                if (bucket[j] < 1e-300) continue;
                double ratio = st.targets[k][j] * W_total / bucket[j];
                if (ratio <= 0.0) continue;
                scale[j] = ratio;
            }
            for (int c = 0; c < ct.M_cell; c++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k]) X[c] *= scale[g];
            }
        }

        // j98p.2: target_mass must be st.n, not Σ X[c]. Post-sweep X[c] values
        // drift below n when capacity bounds clamp cells; using Σ X[c] as the
        // bisection target propagates that drift into mu, corrupting subsequent
        // Bregman projections. When no clamping occurs Σ X[c] == n exactly, so
        // this fix is a no-op on the unclamped path.
        const double target_mass = static_cast<double>(st.n);

        // Short-circuit: if X already within capacity bounds AND a[] is all-zero,
        // bisection is a no-op (projection is identity; Dykstra correction adds 0).
        // j98p.1: when a[] carries accumulated corrections the bisection still
        // needs to run even when X is in bounds, to maintain the fixed-point.
        bool needs_projection = !a_is_zero;
        if (!needs_projection) {
            for (int c = 0; c < ct.M_cell; c++) {
                if (X[c] < L_cell[c] - 1e-12 || X[c] > U_cell[c] + 1e-12) {
                    needs_projection = true;
                    break;
                }
            }
        }
        double mu = 0.0;
        if (needs_projection) {
            if (!bisect_capacity_fast(X, exp_a.data(), L_cell, U_cell, ct.M_cell, target_mass, mu, X_proj)) {
                res.base.status = RK_ERR_INFEAS;
                break;
            }
            // Dykstra correction: accumulate log-space adjustment.
            // Sign convention: a[] is added INSIDE the bisection eval (X*exp(a+mu)),
            // so a[] must store log(X_pre/X_post) so that next iter's X*exp(a) reconstructs
            // the pre-clamp value the previous bisection saw. exp(+a) un-projects;
            // hence log(X) - log(X_proj), the OPPOSITE of the textbook
            // a += log(X_proj) - log(X) form (which assumes exp(-a) un-projects).
            // Verified by fixed-point trace: at convergence the update yields zero. (B11)
            for (int c = 0; c < ct.M_cell; c++) {
                if (X[c] > 1e-300 && X_proj[c] > 1e-300)
                    a[c] += std::log(X[c]) - std::log(X_proj[c]);
                a[c] = std::clamp(a[c], -kAmax, kAmax);
            }
            a_is_zero = false;  // j98p.1: a[] now carries non-zero corrections
            lbw::bulk_scaled_exp(1.0, a.data(), exp_a.data(), ct.M_cell);
            // 773f.5: apply projection result to X only on projection path. No-op copy eliminated.
            for (int c = 0; c < ct.M_cell; c++) X[c] = X_proj[c];
        }
        // B10: On non-projection path, X unchanged — do NOT zero a[]. Accumulated correction from prior
        // projected iterations is valid history; zeroing would corrupt subsequent bisect_capacity calls.

        if (iter == 1 || iter % kErrCheckInterval == 0 || iter == st.inner_max_iter) {
            // B10: recompute W after sweeps+bisection — X mass may have shifted
            double W_current = 0.0;
            for (int c = 0; c < ct.M_cell; ++c) W_current += X[c];
            const double W = W_current;
            auto m = lbw::compute_cell_metrics(st, ct, X, W, bucket);

            double l1_sum = 0.0;
            for (int c = 0; c < ct.M_cell; c++)
                l1_sum += std::fabs(X[c] - X_prev[c]);
            m.l1 = l1_sum / static_cast<double>(st.n);

            res.base.max_error        = m.errRp;
            res.base.kl               = m.kl;
            res.base.mean_error       = m.mean_err;
            res.base.chi2             = m.chi2;
            res.base.grake_norm       = m.grake_norm;
            res.base.l1_weight_change = m.l1;
            for (int c = 0; c < ct.M_cell; c++) X_prev[c] = X[c];

            const double curr_best = lbw::select_metric(
                st.convergence_cfg.metric, m);
            if (std::isfinite(curr_best) && curr_best < best.best_metric) {
                best.update(curr_best,
                            lbw::compute_weight_kl(X, X_init, ct.M_cell, st.n,
                                                   kl_ratio_buf.data(), kl_weight_buf.data()),
                            iter, X);
            }

            if (lbw::check_convergence(st.convergence_cfg, m, prev_metric_for_rule, st.tol_abs)) {
                lbw::mark_converged(res, st.convergence_cfg, iter);
                break;
            }
        }
    }

    // G8c: write best-iterate fields from tracker.
    res.base.convergence_solver_objective = best.best_objective;
    res.base.convergence_minimized_metric = static_cast<int>(st.convergence_cfg.metric);
    res.base.best_error = best.best_metric;
    res.base.best_iter  = best.best_iter;

    if (best.has_best()) {
        std::vector<double> w_snap = best.best_weights;  // cell-level X snapshot
        double s = 0.0;
        for (int c = 0; c < ct.M_cell; c++) s += w_snap[c];
        if (s > 0.0) {
            const double sc = static_cast<double>(st.n) / s;
            for (int c = 0; c < ct.M_cell; c++) w_snap[c] *= sc;
        }
        res.base.best_weights.resize(st.n);
        const double hi_obs = lbw::resolve_hi(st);
        for (int i = 0; i < st.n; i++) {
            int c = ct.cell_of[i];
            double mult = (X_init[c] > 0.0) ? w_snap[c] / X_init[c] : 1.0;
            res.base.best_weights[i] = std::clamp(st.weights[i] * mult, lo, hi_obs);
        }
    } else {
        res.base.best_weights.assign(st.n, 0.0);
    }

    // Obs expansion: w_i = d_i × X[c]/X_init[c], hard clamp.
    // sum(X[c])=n preserved by Sinkhorn+bisection → no normalization needed.
    const double hi_obs = lbw::resolve_hi(st);
    apply_obs_expansion(ct, X, X_init, st.n, lo, hi_obs, st.weights);
    return res;
}

} // namespace lbw
