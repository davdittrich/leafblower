#include "chebyshev.hpp"
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

ChebyshevResult chebyshev_ipm(CalibState& st, LpVariant variant)
{
    static constexpr int    kMaxIpm    = 500;   // hard cap; user controls via max_iterations
    static constexpr double kSigma     = 0.1;    // centering parameter
    static constexpr double kTolMu     = 1e-6;   // complementarity gap convergence threshold
    static constexpr double kEps       = 1e-14;  // strict interior buffer
    static constexpr double kEpsLdlt   = 1e-10;  // LDLT perturbation
    static constexpr double kStepScale = 0.99;   // line search safety factor
    ChebyshevResult res;
    res.status = RK_ERR_NOCONV;

    // L1_WEIGHT is not computable by IPM (no prev-weight reference).
    // Fall back to MAX_ERR so the improvement rule tracks the actual objective.
    if (st.convergence_cfg.metric == lbw::CalibMetric::L1_WEIGHT)
        st.convergence_cfg.metric = lbw::CalibMetric::MAX_ERR;

    CellTable ct;
    if (build_cell_table(st.n, st.K, st.group_ids, st.cat_counts, st.weights, ct) != 0) {
        res.status = RK_ERR_BADARG;
        std::strncpy(res.message, "build_cell_table failed", sizeof(res.message) - 1);
        return res;
    }
    res.M_cell = ct.M_cell;

    const double lo = st.min_weight;
    const double hi = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
    std::vector<double> L_cell(ct.M_cell), U_cell(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) {
        L_cell[c] = lo * static_cast<double>(ct.n_per_cell[c]);
        U_cell[c] = hi * static_cast<double>(ct.n_per_cell[c]);
    }

    std::vector<int> cat_offset(st.K);
    int nct = 0;
    for (int k = 0; k < st.K; k++) { cat_offset[k] = nct; nct += st.cat_counts[k]; }

    if (nct > kNCatsTotalMax) {
        res.status = RK_ERR_BADARG;
        std::snprintf(res.message, sizeof(res.message),
                      "n_cats_total=%d exceeds limit %d", nct, kNCatsTotalMax);
        return res;
    }

    const double n_d = static_cast<double>(st.n);

    // w_kj scaling
    std::vector<double> w_kj(nct);
    for (int k = 0; k < st.K; k++)
        for (int j = 0; j < st.cat_counts[k]; j++) {
            int m = cat_offset[k] + j;
            w_kj[m] = (variant == LpVariant::CHEBYSHEV)
                      ? n_d : 1.0 + st.targets[k][j] * n_d;
        }

    std::vector<double> Tgt(nct);
    for (int k = 0; k < st.K; k++)
        for (int j = 0; j < st.cat_counts[k]; j++)
            Tgt[cat_offset[k]+j] = st.targets[k][j] * n_d;

    // T_flat[m] = pure target proportion (without n_d factor).
    // Slack updates use T_flat[m]*W_current to track the actual calibration
    // objective max_m |S[m]/W - T[m]| rather than |S[m]/n_d - T[m]|.
    std::vector<double> T_flat(nct);
    for (int m = 0; m < nct; m++) T_flat[m] = Tgt[m] / n_d;

    // Initial cell masses
    std::vector<double> X_init(ct.M_cell, 0.0);
    for (int i = 0; i < st.n; i++) X_init[ct.cell_of[i]] += st.weights[i];

    // Marginal sum helper
    auto compute_S = [&](const std::vector<double>& Xv, std::vector<double>& Sv) {
        std::fill(Sv.begin(), Sv.end(), 0.0);
        for (int k = 0; k < st.K; k++)
            for (int c = 0; c < ct.M_cell; c++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k]) Sv[cat_offset[k]+g] += Xv[c];
            }
    };

    // Warm start: strictly interior
    std::vector<double> X(ct.M_cell);
    {
        double max_gap = 0.0;
        for (int c = 0; c < ct.M_cell; c++) max_gap = std::max(max_gap, U_cell[c]-L_cell[c]);
        double eps_shift = std::max(1e-8 * max_gap, 1e-10);
        for (int c = 0; c < ct.M_cell; c++) {
            double gap = U_cell[c] - L_cell[c];
            if (gap < 2.0*eps_shift) X[c] = 0.5*(L_cell[c]+U_cell[c]);
            else X[c] = std::clamp(X_init[c], L_cell[c]+eps_shift, U_cell[c]-eps_shift);
        }
    }

    std::vector<double> S(nct);
    compute_S(X, S);

    // Initial delta: strictly above current max violation (use W_init for consistent scaling)
    double W_init_pre = 0.0;
    for (int c = 0; c < ct.M_cell; c++) W_init_pre += X[c];
    double delta = 0.0;
    for (int m = 0; m < nct; m++)
        delta = std::max(delta, std::fabs(S[m]-T_flat[m]*W_init_pre) / w_kj[m]);
    delta = 1.1*delta + 1e-8;

    // Primal slacks
    std::vector<double> s_lo(ct.M_cell), s_hi(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) {
        s_lo[c] = std::max(X[c]-L_cell[c], kEps);
        s_hi[c] = std::max(U_cell[c]-X[c], kEps);
    }
    std::vector<double> s_up(nct), s_dn(nct);
    for (int m = 0; m < nct; m++) {
        s_up[m] = std::max(T_flat[m]*W_init_pre+w_kj[m]*delta-S[m], kEps);
        s_dn[m] = std::max(S[m]-T_flat[m]*W_init_pre+w_kj[m]*delta, kEps);
    }
    double s_delta = delta;

    // Dual variables: initialize at μ/s (central path)
    double mu = 1.0;
    std::vector<double> y_lo(ct.M_cell), y_hi(ct.M_cell), y_up(nct), y_dn(nct);
    double y_delta = mu / s_delta;
    for (int c = 0; c < ct.M_cell; c++) { y_lo[c]=mu/s_lo[c]; y_hi[c]=mu/s_hi[c]; }
    for (int m = 0; m < nct; m++) { y_up[m]=mu/s_up[m]; y_dn[m]=mu/s_dn[m]; }

    // Count complement pairs for μ computation
    const int n_comp = 2*ct.M_cell + 2*nct + 1;  // s_lo,s_hi,s_up,s_dn,s_delta

    // Hoisted work vectors
    const int max_cats = *std::max_element(st.cat_counts, st.cat_counts+st.K);
    std::vector<double> D_eff(ct.M_cell), D_marg(nct);
    std::vector<double> N0((size_t)nct*(size_t)nct);
    std::vector<double> u_vec(nct), v_vec(nct), rhs_v(nct), w_sol(nct);
    std::vector<double> dX(ct.M_cell);
    std::vector<double> dS_up(nct), dS_dn(nct);
    std::vector<double> dY_lo(ct.M_cell), dY_hi(ct.M_cell), dY_up(nct), dY_dn(nct);
    std::vector<double> dlambda(nct);   // hoisted: Sherman-Morrison result
    std::vector<double> delta_S(nct);   // hoisted: ΔS[m] = Σ_c A_mc*ΔX[c]
    std::vector<double> bucket_tmp(max_cats);
    double best_delta = delta;
    // Track best X by actual calibration error (errRp) rather than LP delta.
    // LP delta can hit the floor (kEps) due to numerical degeneration while the
    // primal is far from the LP optimum; errRp directly measures calibration quality.
    double best_errRp = std::numeric_limits<double>::infinity();
    std::vector<double> X_best(X);
    // Convergence rule state — uses CalibState cfg (not hardcoded tol)
    double prev_metric_for_rule = std::numeric_limits<double>::infinity();

    const int max_ipm = std::min(kMaxIpm, st.inner_max_iter);
    for (int iter = 0; iter < max_ipm; iter++) {
        res.iterations = iter+1;
        compute_S(X, S);

        // Recompute μ from actual complementarity
        {
            double comp = y_delta * s_delta;
            for (int c = 0; c < ct.M_cell; c++) comp += y_lo[c]*s_lo[c] + y_hi[c]*s_hi[c];
            for (int m = 0; m < nct; m++) comp += y_up[m]*s_up[m] + y_dn[m]*s_dn[m];
            mu = comp / n_comp;
        }

        // Per-iteration metrics needed for convergence dispatch.
        // O(K*nct) pass — cheap vs LDLT O(nct³). Computes all 6 metrics from current S[].
        double W = 0.0;
        for (int c = 0; c < ct.M_cell; c++) W += X[c];
        CellMetrics cm;
        if (W > 1e-300) cm = lbw::compute_cell_metrics(st, ct, X, W, bucket_tmp);
        double errRp = cm.errRp, mean_err = cm.mean_err;
        double kl_max = cm.kl, chi2_total = cm.chi2, grake_norm = cm.grake_norm;
        double l1_weight = 0.0;  // not tracked per IPM step (no prev weights)

        // Track X with the best actual calibration error seen so far.
        if (errRp < best_errRp) { best_errRp = errRp; res.best_iter = iter+1; X_best = X; }

        // Convergence dispatch — uses CalibState cfg (metric+rule+tol), same as all other solvers.
        // apply_rule updates prev_metric_for_rule in-place (tracks improvement across iterations).
        {
            const auto& cfg = st.convergence_cfg;
            const double curr_metric = lbw::select_metric(
                cfg.metric, errRp, mean_err, kl_max, chi2_total, grake_norm, l1_weight);
            bool converged_abs = (cfg.absolute_tol > 0.0) && (curr_metric < cfg.absolute_tol);
            const bool converged_pct = lbw::apply_rule(
                cfg.rule, curr_metric, prev_metric_for_rule, cfg.pct_tol);
            // apply_rule updates prev_metric_for_rule in-place — no separate assignment needed
            bool have_pct = (cfg.pct_tol > 0.0), have_abs = (cfg.absolute_tol > 0.0);
            // IPM convergence: μ → 0 is the correct criterion (not improvement on errRp).
            // CalibRule (improvement/plateau) is designed for iterative projection methods;
            // for IPM it fires prematurely while the algorithm is still making progress.
            // Primary: μ < kTolMu. Secondary: user absolute_tol if set.
            bool converged = (mu < kTolMu);
            if (have_abs) converged = converged || converged_abs;

            if (converged) {
                res.status             = RK_OK;
                res.convergence_metric = static_cast<int>(cfg.metric);
                res.convergence_rule   = static_cast<int>(cfg.rule);
                res.convergence_tol    = cfg.pct_tol;
                res.convergence_iter   = iter+1;
                break;
            }
        }

        const double sigma_mu = kSigma * mu;

        // D_eff: 1/D_eff[c] = sum of all barrier weights for cell c
        for (int c = 0; c < ct.M_cell; c++) {
            double inv_D = y_lo[c]/s_lo[c] + y_hi[c]/s_hi[c];
            for (int k = 0; k < st.K; k++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k]) {
                    int m = cat_offset[k]+g;
                    inv_D += y_up[m]/s_up[m] + y_dn[m]/s_dn[m];
                }
            }
            D_eff[c] = (inv_D > 1e-300) ? 1.0/inv_D : 1e300;
        }

        // D_marg[m] = effective margin weight after Schur complement
        for (int m = 0; m < nct; m++)
            D_marg[m] = 1.0 / (y_up[m]/s_up[m] + y_dn[m]/s_dn[m] + 1e-300);

        // N_0 = A * D_eff * A^T
        if (compute_normal_equations(ct, D_eff.data(), N0.data(),
                                     cat_offset.data(), st.K,
                                     static_cast<size_t>(nct)) != RK_OK) {
            res.status = RK_ERR_BADARG; return res;
        }

        if (W < 1e-300) { res.status = RK_ERR_INFEAS; return res; }

        // RHS WITH complementarity centering (correct formula):
        // rhs[m] = -(S[m] - T_flat[m]*W)  +  D_marg[m]*( rmu_up/s_up - rmu_dn/s_dn )
        // Using T_flat[m]*W (= T[m]*W) as the dynamic target so the Newton step
        // drives S[m]/W → T[m] directly rather than S[m] → T[m]*n_d.
        // where rmu_up = sigma*mu - s_up*y_up  (centering residual)
        for (int m = 0; m < nct; m++) {
            double rmu_up = sigma_mu - s_up[m]*y_up[m];
            double rmu_dn = sigma_mu - s_dn[m]*y_dn[m];
            rhs_v[m] = -(S[m] - T_flat[m]*W)
                       + D_marg[m] * (rmu_up/s_up[m] - rmu_dn/s_dn[m]);
        }

        // δ stationarity: r_δ = 1 - Σ_m w_m*(y_up+y_dn) - y_delta
        // Margin centering contribution to Δδ: Σ_m w_m*(rmu_up/s_up + rmu_dn/s_dn)
        double rmu_delta = sigma_mu - s_delta*y_delta;
        double r_delta_stat = 1.0, margin_delta_center = 0.0;
        r_delta_stat -= y_delta;
        for (int m = 0; m < nct; m++) {
            r_delta_stat -= w_kj[m]*(y_up[m]+y_dn[m]);
            double rmu_up_m = sigma_mu - s_up[m]*y_up[m];
            double rmu_dn_m = sigma_mu - s_dn[m]*y_dn[m];
            margin_delta_center += w_kj[m]*(rmu_up_m/s_up[m] + rmu_dn_m/s_dn[m]);
        }

        // Sherman-Morrison: u[m]=w_kj[m], Theta = second derivative of barrier w.r.t. δ
        double Theta = y_delta / s_delta;
        for (int m = 0; m < nct; m++)
            Theta += w_kj[m]*w_kj[m] * (y_up[m]/s_up[m] + y_dn[m]/s_dn[m]);
        double alpha_sm = (Theta > 1e-300) ? 1.0/Theta : 0.0;
        for (int m = 0; m < nct; m++) u_vec[m] = w_kj[m];

        // LDLT factor N_0
        if (ldlt_factor_inplace(N0.data(), static_cast<size_t>(nct), kEpsLdlt) != RK_OK) {
            res.status = RK_ERR_BADARG; return res;
        }

        // Two back-solves (Sherman-Morrison)
        std::copy(rhs_v.begin(), rhs_v.end(), w_sol.begin());
        ldlt_solve(N0.data(), static_cast<size_t>(nct), w_sol.data());
        std::copy(u_vec.begin(), u_vec.end(), v_vec.begin());
        ldlt_solve(N0.data(), static_cast<size_t>(nct), v_vec.data());

        double utv = 0.0, utw = 0.0;
        for (int m = 0; m < nct; m++) { utv += u_vec[m]*v_vec[m]; utw += u_vec[m]*w_sol[m]; }
        double sm_denom = 1.0 + alpha_sm*utv;
        double sm_coeff = (std::fabs(sm_denom) > 1e-300) ? (alpha_sm*utw/sm_denom) : 0.0;
        for (int m = 0; m < nct; m++) dlambda[m] = w_sol[m] - sm_coeff*v_vec[m];

        // ΔX[c] = D_eff[c] * Σ_k Δλ[m_k]
        std::fill(dX.begin(), dX.end(), 0.0);
        for (int k = 0; k < st.K; k++)
            for (int c = 0; c < ct.M_cell; c++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k])
                    dX[c] += dlambda[cat_offset[k]+g];
            }
        for (int c = 0; c < ct.M_cell; c++) dX[c] *= D_eff[c];

        // Δδ: full formula including δ-stationarity residual and margin centering
        double w_dot_dlambda = 0.0;
        for (int m = 0; m < nct; m++) w_dot_dlambda += w_kj[m]*dlambda[m];
        double d_delta = alpha_sm * (
            rmu_delta/s_delta           // δ complementarity centering
          + margin_delta_center         // margin complementarity contribution to Δδ
          - r_delta_stat                // δ stationarity residual
          - (y_delta/s_delta)*w_dot_dlambda  // Schur coupling
        );

        // Compute ΔS from ΔX in O(K*M_cell)
        std::fill(delta_S.begin(), delta_S.end(), 0.0);
        for (int k = 0; k < st.K; k++)
            for (int c = 0; c < ct.M_cell; c++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k])
                    delta_S[cat_offset[k]+g] += dX[c];
            }
        for (int m = 0; m < nct; m++) {
            dS_up[m] = w_kj[m]*d_delta - delta_S[m];   // Δs_up = w·Δδ − ΔS
            dS_dn[m] = delta_S[m] + w_kj[m]*d_delta;   // Δs_dn = ΔS + w·Δδ
        }

        // Dual Newton steps (correct, not reset!)
        for (int m = 0; m < nct; m++) {
            double rmu_up = sigma_mu - s_up[m]*y_up[m];
            double rmu_dn = sigma_mu - s_dn[m]*y_dn[m];
            dY_up[m] = (rmu_up - y_up[m]*dS_up[m]) / s_up[m];
            dY_dn[m] = (rmu_dn - y_dn[m]*dS_dn[m]) / s_dn[m];
        }
        for (int c = 0; c < ct.M_cell; c++) {
            double rmu_lo = sigma_mu - s_lo[c]*y_lo[c];
            double rmu_hi = sigma_mu - s_hi[c]*y_hi[c];
            dY_lo[c] = (rmu_lo - y_lo[c]*dX[c]) / s_lo[c];
            dY_hi[c] = (rmu_hi + y_hi[c]*dX[c]) / s_hi[c];   // Δs_hi = -ΔX
        }
        double dY_delta = (rmu_delta - y_delta*d_delta) / s_delta;

        // Separate primal and dual line searches
        double alpha_p = 1.0, alpha_d = 1.0;
        for (int c = 0; c < ct.M_cell; c++) {
            if (dX[c] > 0.0)  alpha_p = std::min(alpha_p, kStepScale*s_hi[c]/dX[c]);
            if (dX[c] < 0.0)  alpha_p = std::min(alpha_p, -kStepScale*s_lo[c]/dX[c]);
            if (dY_lo[c] < 0.0) alpha_d = std::min(alpha_d, -kStepScale*y_lo[c]/dY_lo[c]);
            if (dY_hi[c] < 0.0) alpha_d = std::min(alpha_d, -kStepScale*y_hi[c]/dY_hi[c]);
        }
        for (int m = 0; m < nct; m++) {
            if (dS_up[m] < 0.0) alpha_p = std::min(alpha_p, -kStepScale*s_up[m]/dS_up[m]);
            if (dS_dn[m] < 0.0) alpha_p = std::min(alpha_p, -kStepScale*s_dn[m]/dS_dn[m]);
            if (dY_up[m] < 0.0) alpha_d = std::min(alpha_d, -kStepScale*y_up[m]/dY_up[m]);
            if (dY_dn[m] < 0.0) alpha_d = std::min(alpha_d, -kStepScale*y_dn[m]/dY_dn[m]);
        }
        if (d_delta < 0.0) alpha_p = std::min(alpha_p, -kStepScale*s_delta/d_delta);
        if (dY_delta < 0.0) alpha_d = std::min(alpha_d, -kStepScale*y_delta/dY_delta);
        alpha_p = std::max(alpha_p, 1e-10);
        alpha_d = std::max(alpha_d, 1e-10);

        // Update primal
        for (int c = 0; c < ct.M_cell; c++) X[c] += alpha_p*dX[c];
        delta += alpha_p*d_delta;
        if (delta < kEps) delta = kEps;
        s_delta = delta;


        for (int c = 0; c < ct.M_cell; c++) {
            s_lo[c] = std::max(X[c]-L_cell[c], kEps);
            s_hi[c] = std::max(U_cell[c]-X[c], kEps);
        }
        compute_S(X, S);
        // Use T_flat[m]*W_current as the dynamic target so that the LP tracks
        // the actual calibration objective max_m |S[m]/W - T[m]| rather than
        // max_m |S[m]/n_d - T[m]|. W drifts when bounds are asymmetric; fixing
        // the target to T[m]*n_d causes slacks to diverge when W != n_d.
        {
            double W_upd = 0.0;
            for (int c = 0; c < ct.M_cell; c++) W_upd += X[c];
            for (int m = 0; m < nct; m++) {
                s_up[m] = std::max(T_flat[m]*W_upd + w_kj[m]*delta - S[m], kEps);
                s_dn[m] = std::max(S[m] - T_flat[m]*W_upd + w_kj[m]*delta, kEps);
            }
        }

        // Update dual (Newton step, not reset!)
        y_delta += alpha_d*dY_delta;
        if (y_delta < kEps) y_delta = kEps;
        for (int c = 0; c < ct.M_cell; c++) {
            y_lo[c] += alpha_d*dY_lo[c]; if (y_lo[c] < kEps) y_lo[c] = kEps;
            y_hi[c] += alpha_d*dY_hi[c]; if (y_hi[c] < kEps) y_hi[c] = kEps;
        }
        for (int m = 0; m < nct; m++) {
            y_up[m] += alpha_d*dY_up[m]; if (y_up[m] < kEps) y_up[m] = kEps;
            y_dn[m] += alpha_d*dY_dn[m]; if (y_dn[m] < kEps) y_dn[m] = kEps;
        }

        // Guard: cap complementarity products s*y to prevent dual explosion.
        // When duals blow up (mu increases by >100x), reset them to maintain
        // the central path mu_target = current_complementarity / n_comp.
        // This prevents runaway dual variables while preserving primal progress.
        {
            double comp_new = y_delta*s_delta;
            for (int c = 0; c < ct.M_cell; c++) comp_new += y_lo[c]*s_lo[c] + y_hi[c]*s_hi[c];
            for (int m = 0; m < nct; m++) comp_new += y_up[m]*s_up[m] + y_dn[m]*s_dn[m];
            double mu_new = comp_new / n_comp;
            if (mu_new > 100.0 * mu) {
                // Dual explosion: re-center at current mu (prevents runaway duals).
                // Using current mu (not sigma*mu) is intentional: sigma*mu resets too aggressively
                // and causes subsequent steps to overshoot in the dual direction.
                y_delta = mu / s_delta;
                for (int c = 0; c < ct.M_cell; c++) { y_lo[c] = mu/s_lo[c]; y_hi[c] = mu/s_hi[c]; }
                for (int m = 0; m < nct; m++) { y_up[m] = mu/s_up[m]; y_dn[m] = mu/s_dn[m]; }
            }
        }

        if (delta < best_delta) { best_delta = delta; }
    }

    // Populate metrics
    res.convergence_objective = best_delta;
    res.best_error = best_errRp;  // actual calibration error at best_iter
    res.convergence_minimized_metric = static_cast<int>(
        variant == LpVariant::CHEBYSHEV ? CalibMetric::MAX_ERR : CalibMetric::GRAKE_NORM);
    res.convergence_metric = res.convergence_minimized_metric;

    // Use X_best (X at lowest errRp iteration) for final metrics and weights.
    // The final IPM iterate may have drifted due to numerical degeneration; the
    // best-errRp iterate gives the most accurate calibration.
    const std::vector<double>& X_out = X_best;

    // All 6 metrics
    double W_final = 0.0;
    for (int c = 0; c < ct.M_cell; c++) W_final += X_out[c];
    {
        auto m2 = lbw::compute_cell_metrics(st, ct, X_out, W_final, bucket_tmp);
        res.max_error  = m2.errRp;
        res.kl         = m2.kl;
        res.chi2       = m2.chi2;
        res.mean_error = m2.mean_err;
        res.grake_norm = m2.grake_norm;
    }

    // Obs expansion using X_out (best-errRp iterate)
    const double hi_obs = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
    for (int i = 0; i < st.n; i++) {
        int c = ct.cell_of[i];
        double mult = (X_init[c] > 1e-10) ? X_out[c]/X_init[c] : 1.0;
        st.weights[i] = std::clamp(st.weights[i]*mult, lo, hi_obs);
    }
    res.best_weights.resize(st.n);
    std::copy(st.weights, st.weights+st.n, res.best_weights.begin());
    return res;
}

} // namespace lbw
