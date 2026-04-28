#include "sinkhorn.hpp"
#include "calib_dispatch.hpp"
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

    double lo = -50.0, hi = 50.0;
    while (lo > -500.0 && f(lo) > 0.0) lo *= 2.0;
    while (hi < 500.0  && f(hi) < 0.0) hi *= 2.0;
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

SinkhornResult sinkhorn_solve(CalibState& st) {
    static constexpr int    kErrCheckInterval = 10;
    static constexpr double kAmax             = 30.0;  // exp(30)≈1e13 >> max practical weight ratio

    SinkhornResult res;
    res.status = RK_ERR_NOCONV;

    CellTable ct;
    if (build_cell_table(st.n, st.K, st.group_ids, st.cat_counts, st.weights, ct) != 0) {
        res.status = RK_ERR_BADARG;
        std::strncpy(res.message, "build_cell_table failed", sizeof(res.message) - 1);
        return res;
    }
    res.M_cell = ct.M_cell;

    std::vector<double> X(ct.M_cell, 0.0);
    for (int i = 0; i < st.n; i++) X[ct.cell_of[i]] += st.weights[i];
    const std::vector<double> X_init(X);

    const double lo = st.min_weight;
    const double hi = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
    std::vector<double> L_cell(ct.M_cell), U_cell(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) {
        L_cell[c] = lo * static_cast<double>(ct.n_per_cell[c]);
        U_cell[c] = hi * static_cast<double>(ct.n_per_cell[c]);
    }

    // log-domain Dykstra correction for capacity box
    std::vector<double> a(ct.M_cell, 0.0);

    int max_cats = *std::max_element(st.cat_counts, st.cat_counts + st.K);

    double best_metric_seen    = std::numeric_limits<double>::infinity();
    double best_objective_seen = 0.0;  // weight KL at best_iter (A1 fix)
    int    best_iter_val    = 0;
    std::vector<double> W_best(ct.M_cell, 0.0);

    // Weight-space KL objective: Σ_c X[c]*log(X[c]/X_init[c])/n
    // Distinct from m.kl (marginal KL) — this is what sinkhorn actually minimizes.
    auto compute_weight_kl = [&]() -> double {
        double wkl = 0.0;
        const double inv_n = 1.0 / static_cast<double>(st.n);
        for (int c = 0; c < ct.M_cell; c++) {
            if (X_init[c] > 0.0 && X[c] > 0.0)
                wkl += X[c] * std::log(X[c] / X_init[c]) * inv_n;
        }
        return std::isfinite(wkl) ? wkl : 0.0;
    };

    std::vector<double> bucket(max_cats);
    std::vector<double> scale(max_cats);
    std::vector<double> X_proj(ct.M_cell);
    std::vector<double> X_prev(X);
    double prev_metric_for_rule = std::numeric_limits<double>::infinity();

    for (int iter = 1; iter <= st.inner_max_iter; iter++) {
        res.iterations = iter;

        double W_total = 0.0;
        for (int c = 0; c < ct.M_cell; c++) W_total += X[c];

        for (int k = 0; k < st.K; k++) {
            std::fill(bucket.begin(), bucket.begin() + st.cat_counts[k], 0.0);
            for (int c = 0; c < ct.M_cell; c++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k]) bucket[g] += X[c];
            }
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

        double target_mass = 0.0;
        for (int c = 0; c < ct.M_cell; c++) target_mass += X[c];

        // Short-circuit: if X already within capacity bounds, bisection is a no-op
        // (projection is identity; Dykstra correction adds 0).
        bool needs_projection = false;
        for (int c = 0; c < ct.M_cell; c++) {
            if (X[c] < L_cell[c] - 1e-12 || X[c] > U_cell[c] + 1e-12) {
                needs_projection = true;
                break;
            }
        }
        double mu = 0.0;
        if (needs_projection) {
            if (!bisect_capacity(X, a, L_cell, U_cell, ct.M_cell, target_mass, mu, X_proj)) {
                res.status = RK_ERR_INFEAS;
                break;
            }
        } else {
            X_proj = X;
        }
        if (needs_projection) {
            for (int c = 0; c < ct.M_cell; c++) {
                if (X[c] > 1e-300 && X_proj[c] > 1e-300)
                    a[c] += std::log(X[c]) - std::log(X_proj[c]);
                a[c] = std::clamp(a[c], -kAmax, kAmax);
            }
        }
        for (int c = 0; c < ct.M_cell; c++) X[c] = X_proj[c];

        if (iter == 1 || iter % kErrCheckInterval == 0 || iter == st.inner_max_iter) {
            const double W = W_total;  // preserved by Sinkhorn sweeps + bisection
            auto m = lbw::compute_cell_metrics(st, ct, X, W, bucket);

            double l1_sum = 0.0;
            for (int c = 0; c < ct.M_cell; c++)
                l1_sum += std::fabs(X[c] - X_prev[c]);
            m.l1 = l1_sum / static_cast<double>(st.n);

            res.max_error        = m.errRp;
            res.kl               = m.kl;
            res.mean_error       = m.mean_err;
            res.chi2             = m.chi2;
            res.grake_norm       = m.grake_norm;
            res.l1_weight_change = m.l1;
            for (int c = 0; c < ct.M_cell; c++) X_prev[c] = X[c];

            const double curr_best = lbw::select_metric(
                st.convergence_cfg.metric,
                m.errRp, m.mean_err, m.kl, m.chi2, m.grake_norm, m.l1);
            if (std::isfinite(curr_best) && curr_best < best_metric_seen) {
                // === BEST-ITER UPDATE ===
                best_metric_seen    = curr_best;
                best_iter_val       = iter;
                best_objective_seen = compute_weight_kl();
                W_best              = X;
                // === END BEST-ITER UPDATE ===
            }

            if (lbw::check_convergence(st.convergence_cfg, m, prev_metric_for_rule, st.tol_abs)) {
                res.status             = RK_OK;
                res.convergence_metric = static_cast<int>(st.convergence_cfg.metric);
                res.convergence_rule   = static_cast<int>(st.convergence_cfg.rule);
                res.convergence_tol    = st.convergence_cfg.pct_tol;
                res.convergence_iter   = iter;
                break;
            }
        }
    }

    res.convergence_solver_objective = best_objective_seen;
    res.convergence_minimized_metric = static_cast<int>(st.convergence_cfg.metric);
    res.best_error = best_metric_seen;
    res.best_iter  = best_iter_val;

    if (std::isfinite(best_metric_seen)) {
        double s = 0.0;
        for (int c = 0; c < ct.M_cell; c++) s += W_best[c];
        if (s > 0.0) {
            const double sc = static_cast<double>(st.n) / s;
            for (int c = 0; c < ct.M_cell; c++) W_best[c] *= sc;
        }
        res.best_weights.resize(st.n);
        const double hi_obs = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
        for (int i = 0; i < st.n; i++) {
            int c = ct.cell_of[i];
            double mult = (X_init[c] > 0.0) ? W_best[c] / X_init[c] : 1.0;
            res.best_weights[i] = std::clamp(st.weights[i] * mult, lo, hi_obs);
        }
    } else {
        res.best_weights.assign(st.n, 0.0);
    }

    // Obs expansion: w_i = d_i × X[c]/X_init[c], hard clamp.
    // sum(X[c])=n preserved by Sinkhorn+bisection → no normalization needed.
    const double hi_obs = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
    apply_obs_expansion(ct, X, X_init, st.n, lo, hi_obs, st.weights);
    return res;
}

} // namespace lbw
