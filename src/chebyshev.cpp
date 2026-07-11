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
#include "trajectory.hpp"  // STUDY-BRANCH-ONLY-DO-NOT-MERGE

namespace lbw {

ChebyshevResult chebyshev_ipm(
    CalibState& st,
    const std::vector<double>& w_warm_obs)
{
    static constexpr int    kMaxIpm    = 500;   // hard cap; K=9 overlapping margins don't converge
                                                // regardless of budget — Mehrotra needs LP stability
                                                // improvements for ill-conditioned systems (follow-on)
    static constexpr double kSigma     = 0.1;    // centering parameter
    static constexpr double kTolMu     = 1e-6;   // complementarity gap convergence threshold
    static constexpr double kEps       = 1e-14;  // strict interior buffer
    static constexpr double kEpsCholesky   = 1e-10;  // LDLT perturbation
    static constexpr double kStepScale = 0.99;   // line search safety factor
    static constexpr int    kInfeasPersistence = 5;  // consecutive negative-slack iters before INFEAS
    static constexpr double kWarmStartRelEps       = 1e-8;   // fractional shift off bound for strict-interior warm start
    static constexpr double kWarmStartAbsEps        = 1e-10; // absolute floor when gap is tiny
    static constexpr double kEpsChebyshev           = 1e-10; // obs-expansion guard: skip cells with near-zero X_init
    static constexpr double kPrimalMachinePrecConv  = 1e-8;  // Mehrotra: accept when best errRp at machine precision
    ChebyshevResult res;
    // Chebyshev defaults differ from CalibResult — override here to preserve existing behavior.
    res.base.status                       = RK_ERR_NOCONV;
    res.base.convergence_rule             = 0;
    res.base.convergence_tol             = 0.0;
    res.base.best_iter                    = 1;
    res.base.convergence_solver_objective = std::numeric_limits<double>::infinity();

    // L1_WEIGHT is not computable by IPM (no prev-weight reference).
    // Fall back to MAX_ERR so the improvement rule tracks the actual objective.
    // CR-C16 (kxna.16): reject a non-positive inner iteration budget up front, BEFORE
    // any mutation of the caller-owned st (side-effect-free reject). Otherwise
    // max_ipm = min(kMaxIpm, inner_max_iter) <= 0 runs the IPM loop zero times and
    // returns the interior-shifted init as X_best with best_error=inf, NOCONV and an
    // empty message. Matches the public-wrapper guard (max_iterations>=1, kxna.15) and
    // logit's low-level BADARG for the same misconfiguration (kxna.21).
    if (st.inner_max_iter < 1) {
        res.base.status = RK_ERR_BADARG;
        std::snprintf(res.message, sizeof(res.message),
            "chebyshev: inner_max_iter (%d) must be >= 1", st.inner_max_iter);
        return res;
    }

    if (st.convergence_cfg.metric == lbw::CalibMetric::L1_WEIGHT)
        st.convergence_cfg.metric = lbw::CalibMetric::MAX_ERR;

    CellTable ct;
    std::vector<double> X_init;
    double hi_eff;
    std::vector<double> L_cell, U_cell;
    std::vector<int> cat_offset;
    int nct;
    if (lbw::solver_setup_ct(st, ct, X_init, hi_eff, L_cell, U_cell,
                              cat_offset, nct, res) != RK_OK)
        return res;
    res.M_cell = ct.M_cell;

    const double n_d = static_cast<double>(st.n);

    // Margin scaling w_kj[m] is the run-constant scalar n_d for every (k,j)
    // (folded from a per-margin vector — xc1s.14). Used directly as n_d below.

    std::vector<double> Tgt(nct);
    for (int k = 0; k < st.K; k++)
        for (int j = 0; j < st.cat_counts[k]; j++)
            Tgt[cat_offset[k]+j] = st.targets[k][j] * n_d;

    // T_flat[m] = pure target proportion (without n_d factor).
    // Slack updates use T_flat[m]*W_current to track the actual calibration
    // objective max_m |S[m]/W - T[m]| rather than |S[m]/n_d - T[m]|.
    std::vector<double> T_flat(nct);
    for (int m = 0; m < nct; m++) T_flat[m] = Tgt[m] / n_d;

    // Reference category elimination: drop first category (j=0) per multi-cat
    // margin to break the normalization degeneracy (schur_nu=0 for sum-to-1
    // targets). Used only for the verbose schur_nu diagnostic at verbose >= 2,
    // iter == 0.
    // Reference category = j=0 (first) per ANOVA convention; matches newton_calib.
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
                if (!(st.cat_counts[k] >= 2 && j == 0))
                    full_to_red[m] = nr++;
            }
    }

    // Design cell masses (Σ_{i∈c} st.weights[i]) — the exit-expansion denominator.
    // Populated only on the warm-start path, where the swap below overwrites X_init
    // with warm masses; empty ⇒ no warm start ⇒ X_init still holds the design masses.
    std::vector<double> X_design_agg;

    // Warm-start override: aggregate oris obs-level weights → cell masses,
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
        X_warm.swap(X_init);  // O(1) — X_init ← warm masses, X_warm ← design masses
        X_design_agg = std::move(X_warm);  // keep design masses for the exit denom
    }

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
        delta = std::max(delta, std::fabs(S[m]-T_flat[m]*W_init_pre) / n_d);
    delta = 1.1*delta + 1e-8;

    // Primal slacks
    std::vector<double> s_lo(ct.M_cell), s_hi(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) {
        s_lo[c] = std::max(X[c]-L_cell[c], kEps);
        s_hi[c] = std::max(U_cell[c]-X[c], kEps);
    }
    std::vector<double> s_up(nct), s_dn(nct);
    for (int m = 0; m < nct; m++) {
        s_up[m] = std::max(T_flat[m]*W_init_pre+n_d*delta-S[m], kEps);
        s_dn[m] = std::max(S[m]-T_flat[m]*W_init_pre+n_d*delta, kEps);
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
    std::vector<bool> locked(ct.M_cell, false);  // cells with inv_D≈0: step forced to 0

    // Reduced normal equations + reduced Newton vectors (reference-cat elim → schur_nu > 0)
    std::vector<double> N_red((size_t)nct_red * (size_t)nct_red);
    std::vector<double> v_red(nct_red);
    std::vector<double> e_red(nct_red), w_e_red(nct_red), tmp_red(nct_red);
    std::vector<double> rhs_A_red(nct_red), rhs_B_red(nct_red);
    std::vector<double> dlambda_A_red(nct_red), dlambda_B_red(nct_red);
    std::vector<double> D_jac_red(nct_red);

    // red_to_full inverse map (built once; mirrors full_to_red built at L88–97)
    std::vector<int> red_to_full(nct_red);
    {
        int nr = 0;
        for (int m = 0; m < nct; m++) if (full_to_red[m] >= 0) red_to_full[nr++] = m;
    }

    // The reduced Sherman-Morrison update vector u_red[nr] == w_kj[red_to_full[nr]] == n_d
    // and the squared margin scale w_kj_sq[m] == n_d*n_d are both run-constant scalars
    // (folded from vectors — xc1s.14); used directly as n_d / n_d*n_d below.

    static constexpr double kSchurNuMin = 1e-8;  // tight; D_nu = O(n_d·M_cell)

    // Full-nct outputs (unchanged: duals, slacks, primal step on cells/δ)
    std::vector<double> dX(ct.M_cell);
    std::vector<double> dS_up(nct), dS_dn(nct);
    std::vector<double> dY_lo(ct.M_cell), dY_hi(ct.M_cell), dY_up(nct), dY_dn(nct);
    std::vector<double> delta_S(nct);   // full nct — reference margins included
    std::vector<double> bucket_tmp(max_cats);
    std::vector<double> dX_A(ct.M_cell), dX_B(ct.M_cell);
    std::vector<double> dS_up_A(nct), dS_dn_A(nct);      // Phase A slack steps
    double best_delta = delta;
    // Track best X by actual calibration error (errRp) rather than LP delta.
    // LP delta can hit the floor (kEps) due to numerical degeneration while the
    // primal is far from the LP optimum; errRp directly measures calibration quality.
    double best_errRp = std::numeric_limits<double>::infinity();
    std::vector<double> X_best(X);
    int slack_violations = 0;  // consecutive iterations with negative s_up or s_dn

    // STUDY-BRANCH-ONLY-DO-NOT-MERGE: generic trajectory probe state
    const std::vector<int> traj_probe_targets = lbw::traj_parse_iters();
    std::deque<int> traj_probe_queue(traj_probe_targets.begin(), traj_probe_targets.end());
    std::vector<std::pair<int,double>> traj_probe_samples;
    const int max_ipm = std::min(kMaxIpm, st.inner_max_iter);
    for (int iter = 0; iter < max_ipm; iter++) {
        res.base.iterations = iter+1;
        // CR-H9 (xc1s.9): the loop-top compute_S(X, S) was removed as redundant. The
        // primal X IS updated each iteration (Mehrotra step below), but the compute_S at
        // the END of the loop body (post-primal-update) already refreshes S from the new X,
        // and nothing writes X between there and the next loop top — so S is already current
        // on entry. Pre-loop compute_S (~line 151) seeds S for the first iteration.

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
        if (W < 1e-300) { res.base.status = RK_ERR_INFEAS; goto finalize; }
        CellMetrics cm = lbw::compute_cell_metrics(st, ct, X, W, bucket_tmp);
        // Only save X_best when X is fully finite — Mehrotra on ill-conditioned K>=9 systems
        // can produce NaN X values; saving NaN as best would corrupt the final result.
        if (cm.errRp < best_errRp && std::isfinite(cm.errRp)) {
            bool x_ok = true;
            for (int c = 0; c < ct.M_cell && x_ok; c++) x_ok = x_ok && std::isfinite(X[c]);
            if (x_ok) { best_errRp = cm.errRp; res.base.best_iter = iter+1; X_best = X; }
        }
        // If X has gone NaN (numerical drift in ill-conditioned system), stop iterating.
        // Use whatever X_best was saved before NaN propagated.
        if (!std::isfinite(cm.errRp)) break;
        // CR-H9 (xc1s.9): six per-iteration locals (errRp/mean_err/kl_max/chi2_total/
        // grake_norm/l1_weight) were copied out of cm and never read — removed; cm is
        // consulted directly (via select_metric) below.

        // Convergence dispatch — uses CalibState cfg (metric + absolute_tol).
        {
            const auto& cfg = st.convergence_cfg;
            const double curr_metric = lbw::select_metric(cfg.metric, cm);
            lbw::traj_record(traj_probe_queue, iter+1, curr_metric, traj_probe_samples);  // STUDY-BRANCH-ONLY-DO-NOT-MERGE
            bool converged_abs = (cfg.absolute_tol > 0.0) && (curr_metric < cfg.absolute_tol);
            // CR-H9 (xc1s.9): removed the orphaned apply_rule()/converged_pct result +
            // prev_metric_for_rule — the rule result was never read (the pct branch below
            // keys on best_errRp, not the rule) and the state had no other consumer. The
            // IPM-vs-CalibRule rationale is documented just below.
            bool have_pct = (cfg.pct_tol > 0.0), have_abs = (cfg.absolute_tol > 0.0);
            // IPM convergence: μ → 0 is the correct criterion (not improvement on errRp).
            // CalibRule (improvement/plateau) is designed for iterative projection methods;
            // for IPM it fires prematurely while the algorithm is still making progress.
            // Primary:      μ < kTolMu — complementarity gap closed.
            // Secondary:    user absolute_tol if set.
            // Tertiary:     best_errRp < 1e-8 — Mehrotra drives primal to machine precision while
            //               μ stays large (degenerate complementarity).
            // Quarternary:  best_errRp < pct_tol — IPM stalled but primal quality meets user
            //               tolerance. For tight K≥3 problems μ may not reach 1e-6 within kMaxIpm
            //               while the calibration objective (max marginal error) is already met.
            bool converged = (mu < kTolMu) || (iter > 0 && best_errRp < kPrimalMachinePrecConv);
            if (have_abs) converged = converged || converged_abs;
            if (!converged && have_pct && iter > 0 && std::isfinite(best_errRp))
                converged = (best_errRp < cfg.pct_tol);

            if (converged) {
                lbw::mark_converged(res, cfg, iter+1, st.tol_abs);
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
            if (inv_D > 1e-300) {
                D_eff[c]  = 1.0 / inv_D;
                locked[c] = false;
            } else {
                D_eff[c]  = 0.0;   // locked: cell cannot move this iteration
                locked[c] = true;
            }
        }

        // D_marg[m] = effective margin weight after Schur complement
        for (int m = 0; m < nct; m++)
            D_marg[m] = 1.0 / (y_up[m]/s_up[m] + y_dn[m]/s_dn[m] + 1e-300);

        // δ stationarity residual (centering-independent part)
        // r_delta_stat: invariant across Phase A/B — no dual update between phases.
        double r_delta_stat_acc = 1.0 - y_delta;
        for (int m = 0; m < nct; m++)
            r_delta_stat_acc -= n_d*(y_up[m]+y_dn[m]);
        const double r_delta_stat = r_delta_stat_acc;

        // Sherman-Morrison: the reduced update vector is the constant scalar n_d
        // (== w_kj[red_to_full[nr]]); Theta = 1/alpha_sm. The squared margin scale
        // n_d*n_d (== w_kj_sq[m]) is a run-constant per multi-cat margin.
        double Theta = y_delta / s_delta;
        for (int m = 0; m < nct; m++)
            Theta += n_d*n_d * (y_up[m]/s_up[m] + y_dn[m]/s_dn[m]);
        double alpha_sm = (Theta > 1e-300) ? 1.0/Theta : 0.0;

        // N_red = A_red * D_eff * A_red^T (reduced nct_red×nct_red) — rebuilt fresh each iteration
        if (lbw::compute_normal_equations_reduced(ct, D_eff.data(), N_red.data(),
                                                  cat_offset.data(), st.K,
                                                  static_cast<size_t>(nct_red),
                                                  full_to_red.data()) != RK_OK) {
            res.base.status = RK_ERR_BADARG; goto finalize;
        }

        // Jacobi diagonal preconditioning on N_red
        for (int j = 0; j < nct_red; j++)
            D_jac_red[j] = 1.0 / std::sqrt(std::max(N_red[(size_t)j*nct_red+j], 1e-12));
        for (int i = 0; i < nct_red; i++)
            for (int j = 0; j < nct_red; j++)
                N_red[(size_t)i*nct_red+j] *= D_jac_red[i] * D_jac_red[j];

        // LDLT factor scaled N_red once per iteration
        if (cholesky_factor_inplace(N_red.data(), static_cast<size_t>(nct_red), kEpsCholesky) != RK_OK) {
            res.base.status = RK_ERR_BADARG; goto finalize;
        }
        res.n_factorizations++;

        // Once-per-iteration ν workspace: e_red[nr], D_nu, r_nu
        std::fill(e_red.begin(), e_red.end(), 0.0);
        double D_nu = 0.0;
        for (int c = 0; c < ct.M_cell; c++) {
            D_nu += D_eff[c];
            for (int k = 0; k < st.K; k++) {
                int g = ct.g_per_cell[k][c];
                if (g < 0 || g >= st.cat_counts[k]) continue;
                int nr = full_to_red[cat_offset[k] + g];
                if (nr >= 0) e_red[nr] += D_eff[c];
            }
        }
        const double r_nu = W - n_d;   // W from L222–223

        // Mehrotra predictor-corrector; Jacobi preconditioning applied above.
        // n_comp = 2*ct.M_cell + 2*nct + 1 > 0 always.
        double sigma_mu;
        {
            // ── Phase A: affine predictor (σ = 0, no centering) ──────────────────
            // RHS in reduced space: reference categories dropped
            for (int nr = 0; nr < nct_red; nr++) {
                int m = red_to_full[nr];
                rhs_A_red[nr] = -(S[m] - T_flat[m]*W)
                                + D_marg[m] * (-y_up[m] + y_dn[m]);
            }

            // δ centering contribution (σ=0 → rmu_delta = -s_delta*y_delta)
            // margin_delta_center sums over full nct (δ stationarity uses all margins)
            double margin_delta_center_A = 0.0;
            for (int m = 0; m < nct; m++)
                margin_delta_center_A += n_d*(-y_up[m] - y_dn[m]);
            double rmu_delta_A = -s_delta*y_delta;

            // Solve N_red · dlambda_A_red = rhs_A_red (Jacobi-scaled)
            for (int nr = 0; nr < nct_red; nr++) rhs_A_red[nr] *= D_jac_red[nr];
            cholesky_solve(N_red.data(), static_cast<size_t>(nct_red), rhs_A_red.data());
            for (int nr = 0; nr < nct_red; nr++) dlambda_A_red[nr] = D_jac_red[nr] * rhs_A_red[nr];

            // δ-SM correction on dlambda_A_red (first SM)
            for (int nr = 0; nr < nct_red; nr++) tmp_red[nr] = D_jac_red[nr] * n_d;
            cholesky_solve(N_red.data(), static_cast<size_t>(nct_red), tmp_red.data());
            for (int nr = 0; nr < nct_red; nr++) v_red[nr] = D_jac_red[nr] * tmp_red[nr];
            double utv = 0.0, utw_A = 0.0;
            for (int nr = 0; nr < nct_red; nr++) {
                utv   += n_d*v_red[nr];
                utw_A += n_d*dlambda_A_red[nr];
            }
            double sm_denom = 1.0 + alpha_sm*utv;
            double sm_coeff_A = (std::fabs(sm_denom) > 1e-300) ? (alpha_sm*utw_A/sm_denom) : 0.0;
            for (int nr = 0; nr < nct_red; nr++) dlambda_A_red[nr] -= sm_coeff_A * v_red[nr];

            // ── ν correction (Phase A): third back-solve + Schur Δν ─────────────
            for (int nr = 0; nr < nct_red; nr++) tmp_red[nr] = D_jac_red[nr] * e_red[nr];
            cholesky_solve(N_red.data(), static_cast<size_t>(nct_red), tmp_red.data());
            for (int nr = 0; nr < nct_red; nr++) w_e_red[nr] = D_jac_red[nr] * tmp_red[nr];
            double ute_A = 0.0;
            for (int nr = 0; nr < nct_red; nr++) ute_A += n_d * w_e_red[nr];
            double sm_coeff_e_A = (std::fabs(sm_denom) > 1e-300) ? (alpha_sm * ute_A / sm_denom) : 0.0;
            for (int nr = 0; nr < nct_red; nr++) w_e_red[nr] -= sm_coeff_e_A * v_red[nr];

            double eTw_e_A = 0.0, eTdl_A = 0.0;
            for (int nr = 0; nr < nct_red; nr++) {
                eTw_e_A += e_red[nr] * w_e_red[nr];
                eTdl_A  += e_red[nr] * dlambda_A_red[nr];
            }
            const double schur_nu_A = D_nu - eTw_e_A;
            const double dnu_A = (schur_nu_A > kSchurNuMin) ? (-r_nu - eTdl_A) / schur_nu_A : 0.0;
            for (int nr = 0; nr < nct_red; nr++) dlambda_A_red[nr] -= dnu_A * w_e_red[nr];

            // Verbose iter-0 diagnostic (replaces removed diagnostic block)
            if (st.verbose >= 2 && iter == 0) {
                char msg[160];
                std::snprintf(msg, sizeof(msg),
                    "chebyshev: schur_nu=%.4e dnu=%.4e r_nu=%.4e (iter 0, Phase A)",
                    schur_nu_A, dnu_A, r_nu);
                st.log(msg);
            }

            // dX_A from dlambda_A_red + Δν · D_eff
            std::fill(dX_A.begin(), dX_A.end(), 0.0);
            for (int c = 0; c < ct.M_cell; c++) {
                if (locked[c]) continue;  // D_eff[c]==0: step stays 0 to avoid explosion
                double sum_dlam = 0.0;
                for (int k = 0; k < st.K; k++) {
                    int g = ct.g_per_cell[k][c];
                    if (g < 0 || g >= st.cat_counts[k]) continue;
                    int nr = full_to_red[cat_offset[k] + g];
                    if (nr >= 0) sum_dlam += dlambda_A_red[nr];
                }
                dX_A[c] = D_eff[c] * (sum_dlam + dnu_A);
            }

            // w_dot_dlambda_A AFTER ν correction
            double w_dot_dlambda_A = 0.0;
            for (int nr = 0; nr < nct_red; nr++)
                w_dot_dlambda_A += n_d * dlambda_A_red[nr];
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
                dS_up_A[m] = n_d*d_delta_A - delta_S[m];
                dS_dn_A[m] = delta_S[m] + n_d*d_delta_A;
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
            if (alpha_aff < 1e-4) sigma = std::max(sigma, kSigma);  // predictor stall guard
            sigma_mu = sigma * mu;

            // ── Phase B: corrector (centering + second-order) ─────────────────────
            // r_delta_stat cached at iter top; duals not mutated since Phase A.
            //
            // corr = -Δs_aff·Δy_aff (the Mehrotra 2nd-order term). The affine
            // complementarity row gives Δy_aff = -y - y·Δs_aff/s (FULL form,
            // incl. the -y from the -sy residual), so
            //   -Δs_aff·Δy_aff = y·Δs_aff + y·Δs_aff²/s.
            // BOTH terms are required. The linear y·Δs_aff piece is NOT a stray
            // predictor residual — dropping it (a previously proposed "fix",
            // leafblower-xiox) breaks the correction. Verified against a CLARABEL
            // LP reference: with this term the solver hits the Chebyshev optimum
            // to ~1e-11 on converged problems; removing it does not affect the
            // K≥9 stall (a separate conditioning issue, see header). Sign flips
            // for corr_hi below because Δs_hi = -ΔX.
            // RHS_B (reduced rows only — reference categories dropped)
            for (int nr = 0; nr < nct_red; nr++) {
                int m = red_to_full[nr];
                double corr_up = y_up[m]*dS_up_A[m] + y_up[m]*dS_up_A[m]*dS_up_A[m]/s_up[m];
                double corr_dn = y_dn[m]*dS_dn_A[m] + y_dn[m]*dS_dn_A[m]*dS_dn_A[m]/s_dn[m];
                double rmu_up_B = sigma_mu - s_up[m]*y_up[m] + corr_up;
                double rmu_dn_B = sigma_mu - s_dn[m]*y_dn[m] + corr_dn;
                rhs_B_red[nr] = -(S[m] - T_flat[m]*W)
                                + D_marg[m] * (rmu_up_B/s_up[m] - rmu_dn_B/s_dn[m]);
            }

            // δ centering (full nct; δ stationarity row spans ALL margins including reference)
            double corr_delta = y_delta*d_delta_A + y_delta*d_delta_A*d_delta_A/s_delta;
            double rmu_delta_B = sigma_mu - s_delta*y_delta + corr_delta;
            double margin_delta_center_B = 0.0;
            for (int m = 0; m < nct; m++) {
                double corr_up = y_up[m]*dS_up_A[m] + y_up[m]*dS_up_A[m]*dS_up_A[m]/s_up[m];
                double corr_dn = y_dn[m]*dS_dn_A[m] + y_dn[m]*dS_dn_A[m]*dS_dn_A[m]/s_dn[m];
                double rmu_up_B = sigma_mu - s_up[m]*y_up[m] + corr_up;
                double rmu_dn_B = sigma_mu - s_dn[m]*y_dn[m] + corr_dn;
                margin_delta_center_B += n_d*(rmu_up_B/s_up[m] + rmu_dn_B/s_dn[m]);
            }

            // Solve in reduced space — REUSE factored N_red (no refactor)
            for (int nr = 0; nr < nct_red; nr++) rhs_B_red[nr] *= D_jac_red[nr];
            cholesky_solve(N_red.data(), static_cast<size_t>(nct_red), rhs_B_red.data());
            for (int nr = 0; nr < nct_red; nr++) dlambda_B_red[nr] = D_jac_red[nr] * rhs_B_red[nr];

            // δ-SM correction on dlambda_B_red (same v_red and sm_denom from Phase A)
            double utw_B = 0.0;
            for (int nr = 0; nr < nct_red; nr++) utw_B += n_d*dlambda_B_red[nr];
            double sm_coeff_B = (std::fabs(sm_denom) > 1e-300) ? (alpha_sm*utw_B/sm_denom) : 0.0;
            for (int nr = 0; nr < nct_red; nr++) dlambda_B_red[nr] -= sm_coeff_B * v_red[nr];

            // ── ν correction (Phase B): reuse w_e_red (SM-corrected) and schur_nu from Phase A ──
            // e_red is fixed within iteration → w_e_red identical between phases.
            // T46: reuse eTw_e_A computed at line 411 (line invariant within iter).
            const double schur_nu_B = D_nu - eTw_e_A;
            double eTdl_B = 0.0;
            for (int nr = 0; nr < nct_red; nr++) eTdl_B += e_red[nr] * dlambda_B_red[nr];
            const double dnu_B = (schur_nu_B > kSchurNuMin) ? (-r_nu - eTdl_B) / schur_nu_B : 0.0;
            for (int nr = 0; nr < nct_red; nr++) dlambda_B_red[nr] -= dnu_B * w_e_red[nr];

            // dX_B from dlambda_B_red + Δν · D_eff
            std::fill(dX_B.begin(), dX_B.end(), 0.0);
            for (int c = 0; c < ct.M_cell; c++) {
                if (locked[c]) continue;  // D_eff[c]==0: step stays 0 to avoid explosion
                double sum_dlam = 0.0;
                for (int k = 0; k < st.K; k++) {
                    int g = ct.g_per_cell[k][c];
                    if (g < 0 || g >= st.cat_counts[k]) continue;
                    int nr = full_to_red[cat_offset[k] + g];
                    if (nr >= 0) sum_dlam += dlambda_B_red[nr];
                }
                dX_B[c] = D_eff[c] * (sum_dlam + dnu_B);
            }

            double w_dot_dlambda_B = 0.0;
            for (int nr = 0; nr < nct_red; nr++)
                w_dot_dlambda_B += n_d * dlambda_B_red[nr];
            double d_delta_B = alpha_sm * (
                rmu_delta_B/s_delta
              + margin_delta_center_B
              - r_delta_stat
              - (y_delta/s_delta)*w_dot_dlambda_B
            );

            // Alias Phase B result to dX/d_delta for downstream line search
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
                dS_up[m] = n_d*d_delta_B - delta_S[m];
                dS_dn[m] = delta_S[m] + n_d*d_delta_B;
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
                double corr_hi = -y_hi[c]*dX_A[c] + y_hi[c]*dX_A[c]*dX_A[c]/s_hi[c];
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
                    double raw_sup = T_flat[m]*W_upd + n_d*delta - S[m];
                    double raw_sdn = S[m] - T_flat[m]*W_upd + n_d*delta;
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
        }  // end Mehrotra predictor-corrector block
    }

finalize:
    // Single exit: early failure returns (RK_ERR_INFEAS/RK_ERR_BADARG) jump here so
    // best_weights is populated from X_best and obs-expansion runs on every path
    // (izql). Failure status set before the goto is preserved.
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

    // Obs expansion using X_out (best-errRp iterate). No per-obs clamp here:
    // Σw=n enforcement + bounds_mode post-processing are delegated to the shared
    // finalize_weights helper (canonical oris_finalize contract). This is the
    // only exit path (early failures goto finalize above), so Σw=n now holds on
    // every path (xl44).
    // Denominator MUST be the design cell mass Σ_{i∈c} st.weights[i]: st.weights
    // still hold design weights, so mult = X_out/X_design gives Σ_{i∈c} w_i =
    // X_out[c]. On the warm-start path X_init was overwritten with warm masses,
    // so use the preserved design aggregate there (X_init[c]/X_warm[c] would
    // return ≈ design weights, since a good warm start makes mult ≈ 1).
    const std::vector<double>& X_den = X_design_agg.empty() ? X_init : X_design_agg;
    for (int i = 0; i < st.n; i++) {
        int c = ct.cell_of[i];
        double mult = (X_den[c] > kEpsChebyshev) ? X_out[c]/X_den[c] : 1.0;
        st.weights[i] = st.weights[i] * mult;
    }
    lbw::finalize_weights(st, ct, res.n_bounds_violated, res.n_bounds_clamped);

    res.base.best_weights.resize(st.n);
    std::copy(st.weights, st.weights+st.n, res.base.best_weights.begin());
    lbw::traj_write_csv(traj_probe_samples, "errRp");  // STUDY-BRANCH-ONLY-DO-NOT-MERGE
    return res;
}

} // namespace lbw
