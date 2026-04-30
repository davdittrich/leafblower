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

ChebyshevResult chebyshev_ipm(
    CalibState& st,
    LpVariant   variant,
    const std::vector<double>& w_warm_obs,
    double      delta_warm)
{
    static constexpr int    kMaxIpm    = 500;   // hard cap; K=9 overlapping margins don't converge
                                                // regardless of budget — Mehrotra needs LP stability
                                                // improvements for ill-conditioned systems (follow-on)
    static constexpr double kSigma     = 0.1;    // centering parameter
    static constexpr double kTolMu     = 1e-6;   // complementarity gap convergence threshold
    static constexpr double kEps       = 1e-14;  // strict interior buffer
    static constexpr double kEpsLdlt   = 1e-10;  // LDLT perturbation
    static constexpr double kStepScale = 0.99;   // line search safety factor
    static constexpr int    kInfeasPersistence = 5;  // consecutive negative-slack iters before INFEAS
    static constexpr double kWarmStartRelEps       = 1e-8;   // fractional shift off bound for strict-interior warm start
    static constexpr double kWarmStartAbsEps        = 1e-10; // absolute floor when gap is tiny
    static constexpr double kPrimalMachinePrecConv  = 1e-8;  // Mehrotra: accept when best errRp at machine precision
    (void)variant;  // GRAKE removed; parameter retained for ABI stability
    ChebyshevResult res;
    // Chebyshev defaults differ from CalibResult — override here to preserve existing behavior.
    res.base.status                       = RK_ERR_NOCONV;
    res.base.convergence_rule             = 0;
    res.base.convergence_tol             = 0.0;
    res.base.best_iter                    = 1;
    res.base.convergence_solver_objective = std::numeric_limits<double>::infinity();

    // L1_WEIGHT is not computable by IPM (no prev-weight reference).
    // Fall back to MAX_ERR so the improvement rule tracks the actual objective.
    if (st.convergence_cfg.metric == lbw::CalibMetric::L1_WEIGHT)
        st.convergence_cfg.metric = lbw::CalibMetric::MAX_ERR;

    CellTable ct;
    if (build_cell_table(st.n, st.K, st.group_ids, st.cat_counts, st.weights, ct) != 0) {
        res.base.status = RK_ERR_BADARG;
        std::strncpy(res.message, "build_cell_table failed", sizeof(res.message) - 1);
        return res;
    }
    res.M_cell = ct.M_cell;

    const double lo = st.min_weight;
    const double hi = lbw::resolve_hi(st);
    std::vector<double> L_cell, U_cell;
    lbw::compute_cell_bounds(ct, lo, hi, L_cell, U_cell);

    std::vector<int> cat_offset(st.K);
    int nct = 0;
    for (int k = 0; k < st.K; k++) { cat_offset[k] = nct; nct += st.cat_counts[k]; }

    if (nct > kNCatsTotalMax) {
        res.base.status = RK_ERR_BADARG;
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
            w_kj[m] = n_d;
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

    // Reference category elimination: drop last category per multi-cat margin
    // to break the normalization degeneracy (schur_nu=0 for sum-to-1 targets).
    // Used only for the verbose schur_nu diagnostic at verbose >= 2, iter == 0.
    int nct_red_count = 0;
    for (int k = 0; k < st.K; k++)
        if (st.cat_counts[k] >= 2) nct_red_count++;
    const int nct_red = nct - nct_red_count;

    std::vector<int> full_to_red(nct, -1);
    {
        int nr = 0;
        for (int k = 0; k < st.K; k++)
            for (int j = 0; j < st.cat_counts[k]; j++) {
                int m = cat_offset[k] + j;
                if (!(st.cat_counts[k] >= 2 && j == st.cat_counts[k]-1))
                    full_to_red[m] = nr++;
            }
    }

    // Initial cell masses — cold start from current obs weights
    std::vector<double> X_init(ct.M_cell, 0.0);
    for (int i = 0; i < st.n; i++) X_init[ct.cell_of[i]] += st.weights[i];

    {
        rk_result_t tmp_res = {};
        if (calib_validate_preentry(ct, st, &tmp_res, X_init.data(), nct) != RK_OK) {
            res.base.status = tmp_res.status;
            std::strncpy(res.message, tmp_res.message, sizeof(res.message) - 1);
            return res;
        }
    }

    // Warm-start override: aggregate ieppa obs-level weights → cell masses,
    // then apply mass-preserving clamp (clamp → rescale → reclamp).
    if (!w_warm_obs.empty() && static_cast<int>(w_warm_obs.size()) == st.n) {
        std::vector<double> X_warm(ct.M_cell, 0.0);
        for (int i = 0; i < st.n; i++)
            X_warm[ct.cell_of[i]] += w_warm_obs[i];

        double total_pre = 0.0, total_post = 0.0;
        for (int c = 0; c < ct.M_cell; c++) total_pre  += X_warm[c];
        for (int c = 0; c < ct.M_cell; c++)
            X_warm[c] = std::clamp(X_warm[c], L_cell[c], U_cell[c]);
        for (int c = 0; c < ct.M_cell; c++) total_post += X_warm[c];
        if (total_post > 0.0 && total_pre > 0.0) {
            double scale = total_pre / total_post;
            for (int c = 0; c < ct.M_cell; c++)
                X_warm[c] = std::clamp(X_warm[c] * scale, L_cell[c], U_cell[c]);
        }
        X_warm.swap(X_init);  // O(1) — X_init now holds warm-start masses
    }
    // delta_warm is reserved for future use; do NOT apply directly as LP delta
    // (units differ from max_m |S[m]/W - T[m]| / w_kj[m]).
    (void)delta_warm;

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
        double eps_shift = std::max(kWarmStartRelEps * max_gap, kWarmStartAbsEps);
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
    const int max_cats = lbw::max_cats_count(st.K, st.cat_counts);
    std::vector<double> D_eff(ct.M_cell), D_marg(nct);
    // Full nct system (for δ Sherman-Morrison — must use full space to preserve E1/E2)
    std::vector<double> N0((size_t)nct * (size_t)nct);
    std::vector<double> u_vec(nct), v_vec(nct);
    std::vector<double> dlambda(nct);
    std::vector<double> dX(ct.M_cell);
    std::vector<double> dS_up(nct), dS_dn(nct);
    std::vector<double> dY_lo(ct.M_cell), dY_hi(ct.M_cell), dY_up(nct), dY_dn(nct);
    std::vector<double> delta_S(nct);   // full nct — reference margins included
    std::vector<double> bucket_tmp(max_cats);
    // Mehrotra predictor-corrector workspace
    std::vector<double> D_jac(nct);                       // Jacobi scaling diagonal
    std::vector<double> rhs_A(nct), rhs_B(nct);          // Phase A / Phase B RHS
    std::vector<double> dlambda_A(nct), dlambda_B(nct);  // Phase A / Phase B solutions
    std::vector<double> dX_A(ct.M_cell), dX_B(ct.M_cell);
    std::vector<double> dS_up_A(nct), dS_dn_A(nct);      // Phase A slack steps
    double best_delta = delta;
    // Track best X by actual calibration error (errRp) rather than LP delta.
    // LP delta can hit the floor (kEps) due to numerical degeneration while the
    // primal is far from the LP optimum; errRp directly measures calibration quality.
    double best_errRp = std::numeric_limits<double>::infinity();
    std::vector<double> X_best(X);
    // Convergence rule state — uses CalibState cfg (not hardcoded tol)
    double prev_metric_for_rule = std::numeric_limits<double>::infinity();
    int slack_violations = 0;  // consecutive iterations with negative s_up or s_dn

    const int max_ipm = std::min(kMaxIpm, st.inner_max_iter);
    for (int iter = 0; iter < max_ipm; iter++) {
        res.base.iterations = iter+1;
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
        if (W < 1e-300) { res.base.status = RK_ERR_INFEAS; return res; }
        CellMetrics cm = lbw::compute_cell_metrics(st, ct, X, W, bucket_tmp);
        // Only save X_best when X is fully finite — Mehrotra on ill-conditioned K>=9 systems
        // can produce NaN X values; saving NaN as best would corrupt the final result.
        if (cm.errRp < best_errRp && std::isfinite(cm.errRp)) {
            bool x_ok = true;
            for (int c = 0; c < ct.M_cell && x_ok; c++) x_ok = std::isfinite(X[c]);
            if (x_ok) { best_errRp = cm.errRp; res.base.best_iter = iter+1; X_best = X; }
        }
        // If X has gone NaN (numerical drift in ill-conditioned system), stop iterating.
        // Use whatever X_best was saved before NaN propagated.
        if (!std::isfinite(cm.errRp)) break;
        double errRp = cm.errRp, mean_err = cm.mean_err;
        double kl_max = cm.kl, chi2_total = cm.chi2, grake_norm = cm.grake_norm;
        double l1_weight = 0.0;  // not tracked per IPM step (no prev weights)

        // Convergence dispatch — uses CalibState cfg (metric+rule+tol), same as all other solvers.
        // apply_rule updates prev_metric_for_rule in-place (tracks improvement across iterations).
        {
            const auto& cfg = st.convergence_cfg;
            const double curr_metric = lbw::select_metric(cfg.metric, cm);
            bool converged_abs = (cfg.absolute_tol > 0.0) && (curr_metric < cfg.absolute_tol);
            const bool converged_pct = lbw::apply_rule(
                cfg.rule, curr_metric, prev_metric_for_rule, cfg.pct_tol);
            // apply_rule updates prev_metric_for_rule in-place — no separate assignment needed
            bool have_pct = (cfg.pct_tol > 0.0), have_abs = (cfg.absolute_tol > 0.0);
            // IPM convergence: μ → 0 is the correct criterion (not improvement on errRp).
            // CalibRule (improvement/plateau) is designed for iterative projection methods;
            // for IPM it fires prematurely while the algorithm is still making progress.
            // Primary: μ < kTolMu. Secondary: user absolute_tol if set.
            // Tertiary: best_errRp < 1e-8 — Mehrotra drives primal to machine precision while
            // μ stays large (degenerate complementarity); accept when best calibration is perfect.
            // Guard iter>0: warm-start may already have perfect errRp on iter 0; require a step.
            // Primary: complementarity gap. Tertiary: best_errRp < 1e-8 fires when warm-start
            // brings errRp to machine precision while μ stays large (degenerate complementarity).
            // The X_best NaN guard above ensures best_errRp < 1e-8 only fires on valid solutions.
            // Guard iter>0: warm-start may already have perfect errRp on iter 0; require a step.
            bool converged = (mu < kTolMu) || (iter > 0 && best_errRp < kPrimalMachinePrecConv);
            if (have_abs) converged = converged || converged_abs;

            if (converged) {
                res.base.status             = RK_OK;
                res.base.convergence_metric = static_cast<int>(cfg.metric);
                res.base.convergence_rule   = static_cast<int>(cfg.rule);
                res.base.convergence_tol    = cfg.pct_tol;
                res.base.convergence_iter   = iter+1;
                break;
            }
        }

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

        // δ stationarity residual (centering-independent part)
        double r_delta_stat = 1.0;
        r_delta_stat -= y_delta;
        for (int m = 0; m < nct; m++)
            r_delta_stat -= w_kj[m]*(y_up[m]+y_dn[m]);

        // Sherman-Morrison: u[m]=w_kj[m], Theta = second derivative of barrier w.r.t. δ
        double Theta = y_delta / s_delta;
        for (int m = 0; m < nct; m++)
            Theta += w_kj[m]*w_kj[m] * (y_up[m]/s_up[m] + y_dn[m]/s_dn[m]);
        double alpha_sm = (Theta > 1e-300) ? 1.0/Theta : 0.0;
        for (int m = 0; m < nct; m++) u_vec[m] = w_kj[m];

        // N_0 = A * D_eff * A^T (full nct×nct) — rebuilt fresh each iteration
        if (lbw::compute_normal_equations(ct, D_eff.data(), N0.data(),
                                          cat_offset.data(), st.K,
                                          static_cast<size_t>(nct)) != RK_OK) {
            res.base.status = RK_ERR_BADARG; return res;
        }

        // Jacobi diagonal preconditioning: D_jac[j] = 1/sqrt(max(N[j*nct+j], 1e-12))
        // Scale N in-place: N_scaled[i][j] = D_jac[i]*N[i][j]*D_jac[j]
        for (int j = 0; j < nct; j++)
            D_jac[j] = 1.0 / std::sqrt(std::max(N0[(size_t)j*nct+j], 1e-12));
        for (int i = 0; i < nct; i++)
            for (int j = 0; j < nct; j++)
                N0[(size_t)i*nct+j] *= D_jac[i] * D_jac[j];

        // LDLT factor scaled N_0 once per iteration
        if (ldlt_factor_inplace(N0.data(), static_cast<size_t>(nct), kEpsLdlt) != RK_OK) {
            res.base.status = RK_ERR_BADARG; return res;
        }
        res.n_factorizations++;

        // ν diagnostic at verbose >= 2, first iteration only.
        // Uses reduced N (reference categories eliminated) to avoid degenerate near-zero schur_nu.
        // schur_nu = D_nu - e_red^T * N_red^{-1} * e_red.
        if (st.verbose >= 2 && iter == 0 && nct_red > 0) {
            std::vector<double> N_red((size_t)nct_red * (size_t)nct_red);
            std::vector<double> e_red(nct_red, 0.0), w_e_red(nct_red);
            if (lbw::compute_normal_equations_reduced(ct, D_eff.data(), N_red.data(),
                                                      cat_offset.data(), st.K,
                                                      static_cast<size_t>(nct_red),
                                                      full_to_red.data()) == RK_OK &&
                ldlt_factor_inplace(N_red.data(), static_cast<size_t>(nct_red), kEpsLdlt) == RK_OK) {
                for (int c = 0; c < ct.M_cell; c++)
                    for (int k = 0; k < st.K; k++) {
                        int g = ct.g_per_cell[k][c];
                        if (g < 0 || g >= st.cat_counts[k]) continue;
                        int m = cat_offset[k] + g;
                        int nr = full_to_red[m];
                        if (nr >= 0) e_red[nr] += D_eff[c];
                    }
                std::copy(e_red.begin(), e_red.end(), w_e_red.begin());
                ldlt_solve(N_red.data(), static_cast<size_t>(nct_red), w_e_red.data());
                double D_nu = 0.0;
                for (int c = 0; c < ct.M_cell; c++) D_nu += D_eff[c];
                double eTw_e = 0.0;
                for (int nr = 0; nr < nct_red; nr++) eTw_e += e_red[nr] * w_e_red[nr];
                const double schur_nu = D_nu - eTw_e;
                char msg[128];
                std::snprintf(msg, sizeof(msg), "chebyshev: schur_nu=%.4e (iter 0)", schur_nu);
                st.log(msg);
            }
        }

        // Mehrotra predictor-corrector; Jacobi preconditioning applied above.
        // n_comp = 2*ct.M_cell + 2*nct + 1 > 0 always.
        double sigma_mu;
        if (n_comp == 0) {
            // Degenerate: Jacobi-preconditioned single-step with fixed σ
            sigma_mu = kSigma * mu;
            const double rmu_delta_f = sigma_mu - s_delta*y_delta;
            double margin_delta_center_f = 0.0;
            for (int m = 0; m < nct; m++) {
                double rmu_up = sigma_mu - s_up[m]*y_up[m];
                double rmu_dn = sigma_mu - s_dn[m]*y_dn[m];
                rhs_A[m] = -(S[m] - T_flat[m]*W)
                           + D_marg[m] * (rmu_up/s_up[m] - rmu_dn/s_dn[m]);
                margin_delta_center_f += w_kj[m]*(rmu_up/s_up[m] + rmu_dn/s_dn[m]);
            }
            // Scale, solve, unscale (using Jacobi-preconditioned N0)
            for (int m = 0; m < nct; m++) rhs_A[m] *= D_jac[m];
            ldlt_solve(N0.data(), static_cast<size_t>(nct), rhs_A.data());
            for (int m = 0; m < nct; m++) dlambda_A[m] = D_jac[m] * rhs_A[m];
            // SM correction
            for (int m = 0; m < nct; m++) v_vec[m] = D_jac[m] * u_vec[m];
            ldlt_solve(N0.data(), static_cast<size_t>(nct), v_vec.data());
            for (int m = 0; m < nct; m++) v_vec[m] *= D_jac[m];
            double utv_f = 0.0, utw_f = 0.0;
            for (int m = 0; m < nct; m++) { utv_f += u_vec[m]*v_vec[m]; utw_f += u_vec[m]*dlambda_A[m]; }
            double sm_denom_f = 1.0 + alpha_sm*utv_f;
            double sm_coeff_f = (std::fabs(sm_denom_f) > 1e-300) ? (alpha_sm*utw_f/sm_denom_f) : 0.0;
            for (int m = 0; m < nct; m++) dlambda[m] = dlambda_A[m] - sm_coeff_f * v_vec[m];
            // dX from dlambda
            std::fill(dX.begin(), dX.end(), 0.0);
            for (int k = 0; k < st.K; k++)
                for (int c = 0; c < ct.M_cell; c++) {
                    int g = ct.g_per_cell[k][c];
                    if (g >= 0 && g < st.cat_counts[k])
                        dX[c] += dlambda[cat_offset[k]+g];
                }
            for (int c = 0; c < ct.M_cell; c++) dX[c] = D_eff[c] * dX[c];
            double w_dot_dlambda_f = 0.0;
            for (int m = 0; m < nct; m++) w_dot_dlambda_f += w_kj[m]*dlambda[m];
            double d_delta_f = alpha_sm * (
                rmu_delta_f/s_delta
              + margin_delta_center_f
              - r_delta_stat
              - (y_delta/s_delta)*w_dot_dlambda_f
            );
            std::fill(delta_S.begin(), delta_S.end(), 0.0);
            for (int k = 0; k < st.K; k++)
                for (int c = 0; c < ct.M_cell; c++) {
                    int g = ct.g_per_cell[k][c];
                    if (g >= 0 && g < st.cat_counts[k])
                        delta_S[cat_offset[k]+g] += dX[c];
                }
            for (int m = 0; m < nct; m++) {
                dS_up[m] = w_kj[m]*d_delta_f - delta_S[m];
                dS_dn[m] = delta_S[m] + w_kj[m]*d_delta_f;
            }
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
                dY_hi[c] = (rmu_hi + y_hi[c]*dX[c]) / s_hi[c];
            }
            double dY_delta_f = (rmu_delta_f - y_delta*d_delta_f) / s_delta;
            // Line search and update
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
            if (d_delta_f < 0.0) alpha_p = std::min(alpha_p, -kStepScale*s_delta/d_delta_f);
            if (dY_delta_f < 0.0) alpha_d = std::min(alpha_d, -kStepScale*y_delta/dY_delta_f);
            for (int c = 0; c < ct.M_cell; c++) X[c] += alpha_p*dX[c];
            delta += alpha_p*d_delta_f;
            if (delta < kEps) delta = kEps;
            s_delta = delta;
            for (int c = 0; c < ct.M_cell; c++) {
                s_lo[c] = std::max(X[c]-L_cell[c], kEps);
                s_hi[c] = std::max(U_cell[c]-X[c], kEps);
            }
            compute_S(X, S);
            {
                double W_upd = 0.0;
                for (int c = 0; c < ct.M_cell; c++) W_upd += X[c];
                int viol_this_iter = 0;
                for (int m = 0; m < nct; m++) {
                    double raw_sup = T_flat[m]*W_upd + w_kj[m]*delta - S[m];
                    double raw_sdn = S[m] - T_flat[m]*W_upd + w_kj[m]*delta;
                    if (raw_sup < 0.0 || raw_sdn < 0.0) viol_this_iter++;
                    s_up[m] = std::max(raw_sup, kEps);
                    s_dn[m] = std::max(raw_sdn, kEps);
                }
                if (viol_this_iter > 0) {
                    slack_violations++;
                    if (slack_violations > kInfeasPersistence) {
                        std::snprintf(res.message, sizeof(res.message),
                                      "chebyshev: %d consecutive iters with negative slacks — INFEAS",
                                      slack_violations);
                        res.base.status = RK_ERR_INFEAS;
                        break;
                    }
                } else { slack_violations = 0; }
            }
            y_delta += alpha_d*dY_delta_f; if (y_delta < kEps) y_delta = kEps;
            for (int c = 0; c < ct.M_cell; c++) {
                y_lo[c] += alpha_d*dY_lo[c]; if (y_lo[c] < kEps) y_lo[c] = kEps;
                y_hi[c] += alpha_d*dY_hi[c]; if (y_hi[c] < kEps) y_hi[c] = kEps;
            }
            for (int m = 0; m < nct; m++) {
                y_up[m] += alpha_d*dY_up[m]; if (y_up[m] < kEps) y_up[m] = kEps;
                y_dn[m] += alpha_d*dY_dn[m]; if (y_dn[m] < kEps) y_dn[m] = kEps;
            }
            {
                double comp_new = y_delta*s_delta;
                for (int c = 0; c < ct.M_cell; c++) comp_new += y_lo[c]*s_lo[c] + y_hi[c]*s_hi[c];
                for (int m = 0; m < nct; m++) comp_new += y_up[m]*s_up[m] + y_dn[m]*s_dn[m];
                double mu_new = comp_new / n_comp;
                if (mu_new > 100.0 * mu) {
                    y_delta = mu / s_delta;
                    for (int c = 0; c < ct.M_cell; c++) { y_lo[c] = mu/s_lo[c]; y_hi[c] = mu/s_hi[c]; }
                    for (int m = 0; m < nct; m++) { y_up[m] = mu/s_up[m]; y_dn[m] = mu/s_dn[m]; }
                }
            }
            if (delta < best_delta) { best_delta = delta; }
            continue;
        } else {
            // ── Phase A: affine predictor (σ = 0, no centering) ──────────────────
            // RHS_A[m] = -(S[m] - T_flat[m]*W) + D_marg[m]*(-s_up*y_up/s_up + s_dn*y_dn/s_dn)
            //          = -(S[m] - T_flat[m]*W) - D_marg[m]*(y_up - y_dn)
            // (rmu_up = 0 - s_up*y_up with σ=0; contribution = -y_up; similarly for dn)
            for (int m = 0; m < nct; m++)
                rhs_A[m] = -(S[m] - T_flat[m]*W)
                           + D_marg[m] * (-y_up[m] + y_dn[m]);

            // δ centering contribution to δ step (σ=0 → rmu_delta = -s_delta*y_delta)
            // margin_delta_center = Σ w[m]*(rmu_up/s_up + rmu_dn/s_dn); with σ=0: -y_up-y_dn
            double margin_delta_center_A = 0.0;
            for (int m = 0; m < nct; m++)
                margin_delta_center_A += w_kj[m]*(-y_up[m] - y_dn[m]);
            double rmu_delta_A = -s_delta*y_delta;

            // Scale RHS_A: rhs_A_scaled[j] = D_jac[j] * rhs_A[j]
            for (int m = 0; m < nct; m++) rhs_A[m] *= D_jac[m];
            // Solve Phase A (scaled)
            ldlt_solve(N0.data(), static_cast<size_t>(nct), rhs_A.data());
            // Unscale: dlambda_A[j] = D_jac[j] * rhs_A_scaled_solved[j]
            for (int m = 0; m < nct; m++) dlambda_A[m] = D_jac[m] * rhs_A[m];

            // SM correction for Phase A (δ direction)
            // v_vec = N_0^{-1} · (D_jac * u_vec)
            for (int m = 0; m < nct; m++) v_vec[m] = D_jac[m] * u_vec[m];
            ldlt_solve(N0.data(), static_cast<size_t>(nct), v_vec.data());
            for (int m = 0; m < nct; m++) v_vec[m] *= D_jac[m];

            double utv = 0.0, utw_A = 0.0;
            for (int m = 0; m < nct; m++) { utv += u_vec[m]*v_vec[m]; utw_A += u_vec[m]*dlambda_A[m]; }
            double sm_denom = 1.0 + alpha_sm*utv;
            double sm_coeff_A = (std::fabs(sm_denom) > 1e-300) ? (alpha_sm*utw_A/sm_denom) : 0.0;
            for (int m = 0; m < nct; m++) dlambda_A[m] -= sm_coeff_A * v_vec[m];

            // dx_A, d_delta_A from dlambda_A
            std::fill(dX_A.begin(), dX_A.end(), 0.0);
            for (int k = 0; k < st.K; k++)
                for (int c = 0; c < ct.M_cell; c++) {
                    int g = ct.g_per_cell[k][c];
                    if (g >= 0 && g < st.cat_counts[k])
                        dX_A[c] += dlambda_A[cat_offset[k]+g];
                }
            for (int c = 0; c < ct.M_cell; c++) dX_A[c] = D_eff[c] * dX_A[c];

            double w_dot_dlambda_A = 0.0;
            for (int m = 0; m < nct; m++) w_dot_dlambda_A += w_kj[m]*dlambda_A[m];
            double d_delta_A = alpha_sm * (
                rmu_delta_A/s_delta
              + margin_delta_center_A
              - r_delta_stat
              - (y_delta/s_delta)*w_dot_dlambda_A
            );

            // delta_S_A[m] = Σ_c A[m,c]*dX_A[c]
            std::fill(delta_S.begin(), delta_S.end(), 0.0);
            for (int k = 0; k < st.K; k++)
                for (int c = 0; c < ct.M_cell; c++) {
                    int g = ct.g_per_cell[k][c];
                    if (g >= 0 && g < st.cat_counts[k])
                        delta_S[cat_offset[k]+g] += dX_A[c];
                }
            for (int m = 0; m < nct; m++) {
                dS_up_A[m] = w_kj[m]*d_delta_A - delta_S[m];
                dS_dn_A[m] = delta_S[m] + w_kj[m]*d_delta_A;
            }

            // α_aff: max α s.t. (x+α·dx_A ≥ 0) AND (s+α·ds_A ≥ 0) with 0.99 damping
            double alpha_aff = 1.0;
            for (int c = 0; c < ct.M_cell; c++) {
                if (dX_A[c] > 0.0)  alpha_aff = std::min(alpha_aff, kStepScale*s_hi[c]/dX_A[c]);
                if (dX_A[c] < 0.0)  alpha_aff = std::min(alpha_aff, -kStepScale*s_lo[c]/dX_A[c]);
            }
            for (int m = 0; m < nct; m++) {
                if (dS_up_A[m] < 0.0) alpha_aff = std::min(alpha_aff, -kStepScale*s_up[m]/dS_up_A[m]);
                if (dS_dn_A[m] < 0.0) alpha_aff = std::min(alpha_aff, -kStepScale*s_dn[m]/dS_dn_A[m]);
            }
            if (d_delta_A < 0.0) alpha_aff = std::min(alpha_aff, -kStepScale*s_delta/d_delta_A);

            // μ_aff = clamp(dot(x+α_aff·dx_A, s+α_aff·ds_A) / n_comp, 0, μ*100)
            double comp_aff = y_delta*(s_delta + alpha_aff*d_delta_A);
            for (int c = 0; c < ct.M_cell; c++) {
                comp_aff += y_lo[c]*(s_lo[c] + alpha_aff*dX_A[c])
                          + y_hi[c]*(s_hi[c] - alpha_aff*dX_A[c]);
            }
            for (int m = 0; m < nct; m++) {
                comp_aff += y_up[m]*(s_up[m] + alpha_aff*dS_up_A[m])
                          + y_dn[m]*(s_dn[m] + alpha_aff*dS_dn_A[m]);
            }
            double mu_aff = std::clamp(comp_aff / (double)n_comp, 0.0, mu * 100.0);

            // σ = clamp((μ_aff/μ)³, 1e-8, 1.0); cap at kSigma when predictor stalls
            // (alpha_aff → 0 when primal is near-optimal → sigma → 1 → stall)
            double ratio = (mu > 1e-300) ? (mu_aff / mu) : 1.0;
            double sigma = std::clamp(ratio*ratio*ratio, 1e-8, 1.0);
            if (alpha_aff < 1e-4) sigma = std::min(sigma, kSigma);  // predictor stall guard
            sigma_mu = sigma * mu;

            // ── Phase B: corrector (centering + second-order) ─────────────────────
            // rhs_B[m] = rhs_A_original[m] + σμ/s_up - σμ/s_dn
            //          + D_marg[m]*( σμ/s_up - σμ/s_dn - dx_A·ds_up_A/s_up + dx_A·ds_dn_A/s_dn )
            // Combining: second-order correction = -dS_up_A[m]*y_up_affine_step/s_up
            // Full Mehrotra corrector for rhs_B:
            // rhs_B[m] = -(S-TW) + D_marg*(σμ - s_up*y_up - dS_up_A*dY_up_A)/s_up
            //                     - D_marg*(σμ - s_dn*y_dn - dS_dn_A*dY_dn_A)/s_dn
            // where dY_up_A[m] = (-s_up[m]*y_up[m] - y_up[m]*dS_up_A[m]) / s_up[m]
            //                   = (-y_up[m]*(s_up[m]+dS_up_A[m])) / s_up[m]  (σ=0 affine)
            // So dS_up_A * dY_up_A = dS_up_A * (-y_up*(s_up+dS_up_A))/s_up
            // Simplified: Mehrotra second-order term = -dS_up_A[m]*dY_up_A_unscaled
            // where dY_up_A_unscaled = (0 - y_up*dS_up_A)/s_up  (σ=0)
            // => second_order_up[m] = -dS_up_A[m] * (-y_up[m]*dS_up_A[m]/s_up[m])
            //                       = dS_up_A[m]^2 * y_up[m] / s_up[m]
            // Full corrector rhs:
            // rhs_B[m] = -(S-TW) + D_marg*(σμ - s_up*y_up + dS_up_A^2*y_up/s_up - (-σμ + s_dn*y_dn - dS_dn_A^2*y_dn/s_dn)) / ...
            // Compact form matching the spec's "rhs_B = rhs_A + σ*μ*e - dx_A·ds_A":
            // The "dx_A·ds_A" term in the Schur-complement reduced system is:
            // Δ_m = D_marg[m] * (dS_up_A[m]*(y_up[m]*dS_up_A[m]/s_up[m])
            //                  - dS_dn_A[m]*(y_dn[m]*dS_dn_A[m]/s_dn[m]))
            // Mehrotra second-order correction: -Δs_A*Δy_A for each complementarity pair.
            // With σ=0 affine: Δy_up_A = (-s_up*y_up - y_up*Δs_up_A)/s_up
            //   → -Δs_up_A*Δy_up_A = y_up*Δs_up_A + y_up*Δs_up_A²/s_up
            for (int m = 0; m < nct; m++) {
                double corr_up = y_up[m]*dS_up_A[m] + y_up[m]*dS_up_A[m]*dS_up_A[m]/s_up[m];
                double corr_dn = y_dn[m]*dS_dn_A[m] + y_dn[m]*dS_dn_A[m]*dS_dn_A[m]/s_dn[m];
                double rmu_up_B = sigma_mu - s_up[m]*y_up[m] + corr_up;
                double rmu_dn_B = sigma_mu - s_dn[m]*y_dn[m] + corr_dn;
                rhs_B[m] = -(S[m] - T_flat[m]*W)
                           + D_marg[m] * (rmu_up_B/s_up[m] - rmu_dn_B/s_dn[m]);
            }

            // δ Phase B centering with second-order correction
            // Δy_delta_A = (-s_delta*y_delta - y_delta*d_delta_A)/s_delta
            // -d_delta_A*Δy_delta_A = y_delta*d_delta_A + y_delta*d_delta_A²/s_delta
            double corr_delta = y_delta*d_delta_A + y_delta*d_delta_A*d_delta_A/s_delta;
            double rmu_delta_B = sigma_mu - s_delta*y_delta + corr_delta;
            double margin_delta_center_B = 0.0;
            for (int m = 0; m < nct; m++) {
                double corr_up = y_up[m]*dS_up_A[m] + y_up[m]*dS_up_A[m]*dS_up_A[m]/s_up[m];
                double corr_dn = y_dn[m]*dS_dn_A[m] + y_dn[m]*dS_dn_A[m]*dS_dn_A[m]/s_dn[m];
                double rmu_up_B = sigma_mu - s_up[m]*y_up[m] + corr_up;
                double rmu_dn_B = sigma_mu - s_dn[m]*y_dn[m] + corr_dn;
                margin_delta_center_B += w_kj[m]*(rmu_up_B/s_up[m] + rmu_dn_B/s_dn[m]);
            }

            // Scale and solve Phase B — REUSE same factored N0 (no refactor)
            for (int m = 0; m < nct; m++) rhs_B[m] *= D_jac[m];
            ldlt_solve(N0.data(), static_cast<size_t>(nct), rhs_B.data());
            for (int m = 0; m < nct; m++) dlambda_B[m] = D_jac[m] * rhs_B[m];

            // SM correction for Phase B (same v_vec from Phase A)
            double utw_B = 0.0;
            for (int m = 0; m < nct; m++) utw_B += u_vec[m]*dlambda_B[m];
            double sm_coeff_B = (std::fabs(sm_denom) > 1e-300) ? (alpha_sm*utw_B/sm_denom) : 0.0;
            for (int m = 0; m < nct; m++) dlambda_B[m] -= sm_coeff_B * v_vec[m];

            // dx_B from dlambda_B
            std::fill(dX_B.begin(), dX_B.end(), 0.0);
            for (int k = 0; k < st.K; k++)
                for (int c = 0; c < ct.M_cell; c++) {
                    int g = ct.g_per_cell[k][c];
                    if (g >= 0 && g < st.cat_counts[k])
                        dX_B[c] += dlambda_B[cat_offset[k]+g];
                }
            for (int c = 0; c < ct.M_cell; c++) dX_B[c] = D_eff[c] * dX_B[c];

            double w_dot_dlambda_B = 0.0;
            for (int m = 0; m < nct; m++) w_dot_dlambda_B += w_kj[m]*dlambda_B[m];
            double d_delta_B = alpha_sm * (
                rmu_delta_B/s_delta
              + margin_delta_center_B
              - r_delta_stat
              - (y_delta/s_delta)*w_dot_dlambda_B
            );

            // Alias Phase B result to the standard dlambda/dX/d_delta variables
            for (int m = 0; m < nct; m++) dlambda[m] = dlambda_B[m];
            for (int c = 0; c < ct.M_cell; c++) dX[c] = dX_B[c];

            // Compute delta_S from dX_B
            std::fill(delta_S.begin(), delta_S.end(), 0.0);
            for (int k = 0; k < st.K; k++)
                for (int c = 0; c < ct.M_cell; c++) {
                    int g = ct.g_per_cell[k][c];
                    if (g >= 0 && g < st.cat_counts[k])
                        delta_S[cat_offset[k]+g] += dX_B[c];
                }
            for (int m = 0; m < nct; m++) {
                dS_up[m] = w_kj[m]*d_delta_B - delta_S[m];
                dS_dn[m] = delta_S[m] + w_kj[m]*d_delta_B;
            }

            // Dual Newton steps using Phase B centering residuals (same second-order correction)
            for (int m = 0; m < nct; m++) {
                double corr_up = y_up[m]*dS_up_A[m] + y_up[m]*dS_up_A[m]*dS_up_A[m]/s_up[m];
                double corr_dn = y_dn[m]*dS_dn_A[m] + y_dn[m]*dS_dn_A[m]*dS_dn_A[m]/s_dn[m];
                double rmu_up_B = sigma_mu - s_up[m]*y_up[m] + corr_up;
                double rmu_dn_B = sigma_mu - s_dn[m]*y_dn[m] + corr_dn;
                dY_up[m] = (rmu_up_B - y_up[m]*dS_up[m]) / s_up[m];
                dY_dn[m] = (rmu_dn_B - y_dn[m]*dS_dn[m]) / s_dn[m];
            }
            for (int c = 0; c < ct.M_cell; c++) {
                double corr_lo = y_lo[c]*dX_A[c] + y_lo[c]*dX_A[c]*dX_A[c]/s_lo[c];
                double corr_hi = y_hi[c]*dX_A[c] + y_hi[c]*dX_A[c]*dX_A[c]/s_hi[c];
                // Δs_lo = +ΔX (primal slack increases with X), Δs_hi = -ΔX
                double rmu_lo_B = sigma_mu - s_lo[c]*y_lo[c] + corr_lo;
                double rmu_hi_B = sigma_mu - s_hi[c]*y_hi[c] + corr_hi;
                dY_lo[c] = (rmu_lo_B - y_lo[c]*dX[c]) / s_lo[c];
                dY_hi[c] = (rmu_hi_B + y_hi[c]*dX[c]) / s_hi[c];
            }
            double dY_delta = (rmu_delta_B - y_delta*d_delta_B) / s_delta;

            // Step lengths (0.99 damping)
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
            if (d_delta_B < 0.0) alpha_p = std::min(alpha_p, -kStepScale*s_delta/d_delta_B);
            if (dY_delta  < 0.0) alpha_d = std::min(alpha_d, -kStepScale*y_delta/dY_delta);

            // Update primal
            for (int c = 0; c < ct.M_cell; c++) X[c] += alpha_p*dX[c];
            delta += alpha_p*d_delta_B;
            if (delta < kEps) delta = kEps;
            s_delta = delta;

            for (int c = 0; c < ct.M_cell; c++) {
                s_lo[c] = std::max(X[c]-L_cell[c], kEps);
                s_hi[c] = std::max(U_cell[c]-X[c], kEps);
            }
            compute_S(X, S);
            {
                double W_upd = 0.0;
                for (int c = 0; c < ct.M_cell; c++) W_upd += X[c];
                int viol_this_iter = 0;
                for (int m = 0; m < nct; m++) {
                    double raw_sup = T_flat[m]*W_upd + w_kj[m]*delta - S[m];
                    double raw_sdn = S[m] - T_flat[m]*W_upd + w_kj[m]*delta;
                    if (raw_sup < 0.0 || raw_sdn < 0.0) viol_this_iter++;
                    s_up[m] = std::max(raw_sup, kEps);
                    s_dn[m] = std::max(raw_sdn, kEps);
                }
                if (viol_this_iter > 0) {
                    slack_violations++;
                    if (slack_violations > kInfeasPersistence) {
                        std::snprintf(res.message, sizeof(res.message),
                                      "chebyshev: %d consecutive iters with negative slacks — INFEAS",
                                      slack_violations);
                        res.base.status = RK_ERR_INFEAS;
                        break;
                    }
                } else {
                    slack_violations = 0;
                }
            }

            // Update dual
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

            // Dual explosion guard
            {
                double comp_new = y_delta*s_delta;
                for (int c = 0; c < ct.M_cell; c++) comp_new += y_lo[c]*s_lo[c] + y_hi[c]*s_hi[c];
                for (int m = 0; m < nct; m++) comp_new += y_up[m]*s_up[m] + y_dn[m]*s_dn[m];
                double mu_new = comp_new / n_comp;
                if (mu_new > 100.0 * mu) {
                    y_delta = mu / s_delta;
                    for (int c = 0; c < ct.M_cell; c++) { y_lo[c] = mu/s_lo[c]; y_hi[c] = mu/s_hi[c]; }
                    for (int m = 0; m < nct; m++) { y_up[m] = mu/s_up[m]; y_dn[m] = mu/s_dn[m]; }
                }
            }

            if (delta < best_delta) { best_delta = delta; }
            continue;  // skip the old single-step code below (never reached)
        }

        // n_comp > 0 always; the continue above is always taken. Unreachable.
        if (delta < best_delta) { best_delta = delta; }
    }

    // Populate metrics
    res.base.convergence_solver_objective = best_delta;
    res.base.best_error = best_errRp;  // actual calibration error at best_iter
    res.base.convergence_minimized_metric = static_cast<int>(CalibMetric::MAX_ERR);
    res.base.convergence_metric = res.base.convergence_minimized_metric;

    // Use X_best (X at lowest errRp iteration) for final metrics and weights.
    // The final IPM iterate may have drifted due to numerical degeneration; the
    // best-errRp iterate gives the most accurate calibration.
    const std::vector<double>& X_out = X_best;

    // All 6 metrics
    double W_final = 0.0;
    for (int c = 0; c < ct.M_cell; c++) W_final += X_out[c];
    if (W_final > 1e-300) {
        auto m2 = lbw::compute_cell_metrics(st, ct, X_out, W_final, bucket_tmp);
        res.base.max_error  = m2.errRp;
        res.base.kl         = m2.kl;
        res.base.chi2       = m2.chi2;
        res.base.mean_error = m2.mean_err;
        res.base.grake_norm = m2.grake_norm;
    }
    // else: X_best collapsed — leave metrics at default 0.0 (not NaN)

    // Obs expansion using X_out (best-errRp iterate)
    const double hi_obs = lbw::resolve_hi(st);
    for (int i = 0; i < st.n; i++) {
        int c = ct.cell_of[i];
        double mult = (X_init[c] > 1e-10) ? X_out[c]/X_init[c] : 1.0;
        st.weights[i] = std::clamp(st.weights[i]*mult, lo, hi_obs);
    }
    res.base.best_weights.resize(st.n);
    std::copy(st.weights, st.weights+st.n, res.base.best_weights.begin());
    return res;
}

} // namespace lbw
