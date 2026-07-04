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
#ifndef LBW_NO_R
#  include <R_ext/Print.h>
#  include <R_ext/Lapack.h>
#  include <R_ext/RS.h>     // F77_CALL
#  ifndef FCONE
#    define FCONE
#  endif
#else
#  define Rprintf(...) fprintf(stderr, __VA_ARGS__)
// System LAPACK declaration for the Python build (eb79.18 min-norm warm-start).
extern "C" {
    void dsyevd_(const char* jobz, const char* uplo, int* n, double* A, int* lda,
                 double* w, double* work, int* lwork, int* iwork,
                 int* liwork, int* info);
}
#  define F77_CALL(x) x ## _
#  define FCONE
#endif

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
    // CR-C15b (kxna.21): defense-in-depth for the low-level/direct-c_api path — the
    // R/Python wrappers guard max_iterations>=1, but a direct caller with
    // inner_max_iter=0 makes kMaxNewtonIters=min(50,0)=0, the loop never runs, and
    // w_best stays value-initialized (all zeros) => degenerate all-zero weights + STALL.
    // Reject up front (consistent with the wrapper semantics and chebyshev kxna.16).
    if (st.inner_max_iter < 1) {
        res.base.status = RK_ERR_BADARG;
        std::snprintf(res.message, sizeof(res.message),
            "logit: inner_max_iter (%d) must be >= 1", st.inner_max_iter);
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

    // eb79.23: structural feasibility pre-check. Each margin bucket (k,j) can hold at most
    // Σ U_cell and at least Σ L_cell (cells are independently box-clampable), so its achievable
    // mass is exactly [Σ L_cell, Σ U_cell]. If the target count st.targets[k][j]*n lies STRICTLY
    // outside that interval, NO weighting can meet that margin — report INFEAS (relax bounds)
    // instead of iterating to the cap and mis-reporting BUDGET ("increase iterations", futile).
    // NECESSARY but not SUFFICIENT: a jointly-infeasible problem where every margin is
    // individually reachable is NOT caught here and correctly falls through to the Newton loop's
    // BUDGET/STALL (full joint feasibility is an LP, out of scope). Matches the sinkhorn
    // bisect_capacity_fast / oris structural_infeas convention. One-shot O(M·K).
    {
        constexpr double kFeasEps = 1e-9;  // relative FP margin: exact-corner is NOT flagged
        for (int k = 0; k < K; k++) {
            for (int j = 0; j < st.cat_counts[k]; j++) {
                double min_mass = 0.0, max_mass = 0.0;
                for (int c : cells_per_cat[k][j]) { min_mass += L_cell[c]; max_mass += U_cell[c]; }
                const double target_count = st.targets[k][j] * static_cast<double>(st.n);
                const bool infeasible =
                    (target_count > max_mass * (1.0 + kFeasEps)) ||
                    (target_count < min_mass * (1.0 - kFeasEps)) ||
                    (cells_per_cat[k][j].empty() && target_count > 0.0);
                if (infeasible) {
                    res.base.status = RK_ERR_INFEAS;
                    std::snprintf(res.message, sizeof(res.message),
                        "logit: infeasible — margin %d category %d: target %.6g outside "
                        "achievable [%.6g, %.6g] (relax min_weight/max_weight)",
                        k, j, target_count, min_mass, max_mass);
                    return res;
                }
            }
        }
    }

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
    const CalibConvergenceCfg& cfg = st.convergence_cfg;

    const int kMaxNewtonIters = std::min(50, st.inner_max_iter);
    constexpr double kDeffFloor = 1e-6;  // floor prevents D_eff→0 when sig saturates (sig→0 or 1)
    // Rationale: D_eff[c] = max(kDeffFloor*range, range*sig*(1-sig)).
    // At saturation (|z|→∞, sig→0 or 1), sig*(1-sig)→0 and D_eff→0 without this floor,
    // causing the Newton step b/D_eff to blow up. kDeffFloor=1e-6 scales with range
    // (cell capacity), so the absolute floor is proportional to the feasible weight
    // interval — correctly small for tight-bound cells and correctly large for wide ones.
    // Profiling on convergent configs (K=2..4, ranges 0.3–5.0, seeds 1/10/42):
    //   floor=1e-8*range, 1e-6*range (current), 1e-4*range, 1e-4/range (inverse)
    //   all produce identical iters and max_err — the floor is inactive on all
    //   tested convergent configs. On infeasible/near-infeasible tight-bound
    //   configs (where logit hits iter ceiling regardless), none of the variants
    //   changes the budget-exhausted outcome either.
    // Proposed 1e-4/range (inverse scaling) was REJECTED: it scales INVERSELY with
    // range, giving a LARGER absolute floor for tight-bound cells (small range),
    // which would ADD more damping — the opposite of the stated "over-damp" concern.
    // The current 1e-6*range is correct: tight-bound cells get a proportionally
    // small floor, wide-bound cells get a proportionally large one.
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
        // Normal equations with D_eff=1 (uniform: purely geometric initialization).
        // N_init = A·A^T is rank-deficient when margins are collinear/redundant
        // (rank < nct). eb79.18: solve (AA^T)λ0 = b_init via a MIN-NORM pseudo-inverse
        // (symmetric eigendecomposition, dsyevd), NOT ridge-Cholesky. The old ridge
        // solve put ≈z_target on EACH redundant dual, so z_c = Σ_k λ0[k] summed the
        // SAME z_target over the R redundant margins ⇒ z_c = R·z_target ⇒ σ→0 ⇒ every
        // cell collapsed to L at iter 0 ⇒ D_eff floor ⇒ non-descent Newton ⇒ STALL
        // (verified: eb79.18 WU1 trace — max|λ0|=2.303=|z_target|, w=[L,L], resid≈90).
        // The min-norm solution distributes z_target across the redundant duals (each
        // = z_target/R), so z_c = z_target and the start sits at the feasible design
        // point. Exact (residual-0) start holds iff z_target ∈ range(A^T); otherwise
        // this is the best L2 projection — still a bounded, well-conditioned start.
        std::vector<double> D_ones(M, 1.0);
        std::vector<double> N_init(static_cast<size_t>(nct) * nct, 0.0);
        if (compute_normal_equations(ct, D_ones.data(), N_init.data(),
                                     cat_offset.data(), K, (size_t)nct) == RK_OK) {
            // Symmetric eigendecomposition N_init = V diag(σ) V^T. dsyevd overwrites
            // N_init with V (column-major: V[:,i] = N_init[i*nct .. i*nct+nct-1]);
            // eigvals ascending in `eigvals`.
            std::vector<double> eigvals(nct, 0.0);
            int dsy_n = nct, dsy_lda = nct, dsy_info = 0;
            double wkopt = 0.0; int iwkopt = 0, lwork = -1, liwork = -1;
            F77_CALL(dsyevd)("V", "U", &dsy_n, N_init.data(), &dsy_lda, eigvals.data(),
                             &wkopt, &lwork, &iwkopt, &liwork, &dsy_info FCONE FCONE);
            if (dsy_info == 0) {
                lwork = static_cast<int>(wkopt); liwork = iwkopt;
                std::vector<double> work(std::max(1, lwork));
                std::vector<int>    iwork(std::max(1, liwork));
                F77_CALL(dsyevd)("V", "U", &dsy_n, N_init.data(), &dsy_lda, eigvals.data(),
                                 work.data(), &lwork, iwork.data(), &liwork, &dsy_info FCONE FCONE);
            }
            if (dsy_info == 0 && eigvals[nct - 1] > 0.0) {
                // Min-norm pseudo-inverse: λ0 = Σ_{σ_i > tol} (v_i^T b_init / σ_i) v_i.
                // tol = kEigRankRtol·σ_max separates the true spectrum from the
                // rank-deficiency null space (machine-noise eigenvalues ~1e-13·σ_max).
                // 1e-10 is a conservative double-precision rank cutoff (well above FP
                // noise, well below any genuine eigenvalue of a scaled A·A^T).
                constexpr double kEigRankRtol = 1e-10;
                const double tol = kEigRankRtol * eigvals[nct - 1];
                std::vector<double> lam0(nct, 0.0);
                for (int i = 0; i < nct; i++) {
                    if (eigvals[i] <= tol) continue;              // null-space direction
                    const double* v = &N_init[static_cast<size_t>(i) * nct];
                    double vtb = 0.0;
                    for (int r = 0; r < nct; r++) vtb += v[r] * b_init[r];
                    const double coef = vtb / eigvals[i];
                    for (int r = 0; r < nct; r++) lam0[r] += coef * v[r];
                }
                // Backstop: reject a still-pathological λ0 (min-norm should not trip
                // this on consistent problems) — lambda stays zero, Armijo handles it.
                double max_lambda_init = 0.0;
                for (double lj : lam0) max_lambda_init = std::max(max_lambda_init, std::abs(lj));
                if (max_lambda_init <= kLambdaInitRejectAbs) lambda = std::move(lam0);
            }
        }
        // if eigendecomposition fails or λ0 is ill-conditioned: lambda stays zero.
    }

    // Armijo scratch buffers (pre-allocated to avoid per-halving allocation)
    std::vector<double> w_trial(M, 0.0);        // trial weights for Armijo
    std::vector<double> sig_trial(M, 0.0);      // xc1s.2: trial sigmoids, reused on the accepted step
    std::vector<double> b_trial(nct, 0.0);      // trial residuals for Armijo
    std::vector<double> lambda_trial(nct, 0.0); // trial lambda for Armijo

    // xc1s.2: the logit link's per-cell z-accumulation
    // z_c = clamp(Σ_k λ[cat_offset[k] + g_per_cell[k][c]], ±700)  (±700 guards exp
    // overflow) was triplicated (iter-0 init, Armijo trial, post-step recompute).
    // Fold it into one closure parameterised on the active lambda vector.
    auto cell_z = [&](const std::vector<double>& lam, int c) -> double {
        double z = 0.0;
        for (int k = 0; k < K; k++) {
            int g = ct.g_per_cell[k][c];
            if (g >= 0 && g < st.cat_counts[k]) z += lam[cat_offset[k] + g];
        }
        return std::clamp(z, -700.0, 700.0);
    };

    // Cache for w[]+D_eff[]: L4 writes these after accepting a Newton step; L1 reads them
    // on iter>=1, skipping the redundant O(M*K) recompute. Invalidated only at iter==0.
    std::vector<double> w_cached(M, 0.0);
    std::vector<double> D_eff_cached(M, 0.0);
    bool wDeff_cache_valid = false;

    for (int iter = 0; iter < kMaxNewtonIters; iter++) {
        res.base.iterations = iter + 1;

        // (1) Compute w[c] and D_eff[c] from lambda
        if (wDeff_cache_valid) {
            // Cache valid: L4 of prior iter already computed w[]+D_eff[] from the accepted
            // lambda (same vector). Copy instead of recomputing O(M*K) inner products.
            for (int c = 0; c < M; c++) {
                w[c]     = w_cached[c];
                D_eff[c] = D_eff_cached[c];
            }
        } else {
            // iter==0 path: compute w[]+D_eff[] from scratch and fill the cache.
            // CR-C9 (kxna.9): the old ">50% cells saturated (|z|>650) → NOCONV" early-exit
            // was REMOVED as provably dead + misplaced. It ran ONLY here (cache invalid at
            // iter 0), so it inspected the INITIAL z only. z = A^T·λ0 is the orthogonal
            // projection of z_target onto range(A^T) (λ0 = truncated-pinv(AA^T)·A·z_target),
            // so ‖z‖₂ ≤ ‖z_target‖₂, and each z_target component is clipped to
            // ≤ log((1−1e-4)/1e-4) ≈ 9.21. Firing needs |z|>650 on >M/2 cells, i.e.
            // ‖z‖₂² > (M/2)·650² ≈ 2.1e5·M, contradicting ‖z‖₂² ≤ M·9.21² ≈ 85·M — impossible
            // for ANY K, M (a rejected warm-start leaves λ=0, so z=0). The check never saw the
            // LATER iterates where λ grows and cells actually saturate; re-wiring it to run
            // every iteration would spuriously fail TRANSIENT Newton overshoot the next step
            // corrects (and the old exit returned BEFORE best_weights reconstruction, so it
            // would have returned empty weights had it fired). Genuine saturation is already
            // handled: kDeffFloor keeps D_eff>0 (no divide-by-zero as sig→0/1), and a solution
            // pinned at the bounds lands in the STALL/BUDGET constrained-optimum classification.
            for (int c = 0; c < M; c++) {
                double z     = cell_z(lambda, c);
                double sig   = 1.0 / (1.0 + std::exp(-z));
                double range = U_cell[c] - L_cell[c];
                w[c]          = L_cell[c] + range * sig;
                D_eff[c]      = std::max(kDeffFloor * range, range * sig * (1.0 - sig));
                w_cached[c]     = w[c];
                D_eff_cached[c] = D_eff[c];
            }
            wDeff_cache_valid = true;
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

        // Capture ||b_current||² before cholesky_solve overwrites b with delta_lambda
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
        if (cholesky_factor_inplace(N.data(), static_cast<size_t>(nct), eps_ldlt) != RK_OK) {
            res.base.status = RK_ERR_INFEAS;
            std::snprintf(res.message, sizeof(res.message),
                "logit: LDLT factorization failed (degenerate bounds)");
            break;
        }
        cholesky_solve(N.data(), static_cast<size_t>(nct), b.data());  // b = delta_lambda

        // (4) Armijo line search with step-norm guard
        // Norm guard: cap alpha so no cell z-coord shifts more than kMaxDeltaZ
        // (b = delta_lambda after cholesky_solve)
        double max_delta_z = 0.0;
        for (int c = 0; c < M; c++) {
            // Cell logit coordinate is the SIGNED margin sum z_c = Σ_k λ[k,g]
            // (matches the z-accumulation at the Armijo recompute below), so the
            // per-step shift is Δz_c = Σ_k Δλ[k,g]. Accumulate signed, |·| once:
            // an abs-per-margin sum over-estimates |Δz_c| (triangle inequality)
            // and throttles alpha smaller than necessary on every K>=2 cell.
            double dz = 0.0;
            for (int k = 0; k < K; k++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k]) dz += b[cat_offset[k] + g];
            }
            max_delta_z = std::max(max_delta_z, std::fabs(dz));
        }
        double alpha = (max_delta_z > 0.0) ? std::min(1.0, kMaxDeltaZ / max_delta_z) : 1.0;

        // CR-C8 (kxna.8): remember the best TESTED trial so that if Armijo exhausts
        // all kMaxHalvings we accept the argmin-residual tested step, not the 11th
        // never-evaluated (further-halved) alpha the loop would otherwise leave.
        bool armijo_improved = false;
        double best_trial_alpha = alpha;
        double best_trial_resid_sq = std::numeric_limits<double>::infinity();
        for (int halv = 0; halv < kMaxHalvings; halv++) {
            for (int j = 0; j < nct; j++) lambda_trial[j] = lambda[j] + alpha * b[j];
            // Recompute w_trial (+ retain sig_trial) from lambda_trial (logit link).
            for (int c = 0; c < M; c++) {
                double z = cell_z(lambda_trial, c);
                double sig = 1.0 / (1.0 + std::exp(-z));
                double range = U_cell[c] - L_cell[c];
                sig_trial[c] = sig;
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
            if (resid_sq_trial < best_trial_resid_sq) {
                best_trial_resid_sq = resid_sq_trial;
                best_trial_alpha = alpha;
            }
            if (resid_sq_trial < resid_sq_0 * (1.0 - kArmijoC * alpha)) {
                armijo_improved = true;
                break;
            }
            alpha *= kArmijoHalving;
        }
        if (!armijo_improved) {
            // Armijo exhausted: apply the best TESTED trial, not the further-halved
            // (never-evaluated) alpha the loop's final halving would leave.
            alpha = best_trial_alpha;
            if (st.verbose >= 1)
                Rprintf("[logit] Newton step: Armijo exhausted, accepting best tested alpha=%.2e\n",
                        alpha);
        }
        // Apply the accepted step and refresh w+D_eff (+ cache) for the convergence
        // check and next iter's L1.
        if (armijo_improved) {
            // xc1s.2: adopt the accepted Armijo trial wholesale. lambda_trial is
            // lambda + alpha*b at the accepted alpha, and w_trial/sig_trial are its
            // logit weights — all computed together, so this is the self-consistent
            // post-step state. Reusing it skips a redundant O(M*K) z-accumulation +
            // O(M) exp per accepted step. Adopting lambda_trial (rather than a
            // separate `lambda += alpha*b`) keeps lambda, w and D_eff mutually
            // consistent; it equals the += form up to FP rounding of the SAME
            // quantity lambda+alpha*b — mathematically equivalent, not a regression.
            std::copy(lambda_trial.begin(), lambda_trial.end(), lambda.begin());
            for (int c = 0; c < M; c++) {
                double range = U_cell[c] - L_cell[c];
                double sig   = sig_trial[c];
                w[c]          = w_trial[c];
                D_eff[c]      = std::max(kDeffFloor * range, range * sig * (1.0 - sig));
                w_cached[c]     = w[c];
                D_eff_cached[c] = D_eff[c];
            }
        } else {
            // Armijo exhausted: alpha (best_trial_alpha) differs from the last
            // trial's alpha, so lambda_trial/w_trial/sig_trial are stale — apply
            // the step and recompute w+D_eff from the updated lambda.
            for (int j = 0; j < nct; j++) lambda[j] += alpha * b[j];
            for (int c = 0; c < M; c++) {
                double z     = cell_z(lambda, c);
                double sig   = 1.0 / (1.0 + std::exp(-z));
                double range = U_cell[c] - L_cell[c];
                w[c]          = L_cell[c] + range * sig;
                D_eff[c]      = std::max(kDeffFloor * range, range * sig * (1.0 - sig));
                w_cached[c]     = w[c];
                D_eff_cached[c] = D_eff[c];
            }
        }
        wDeff_cache_valid = true;

        // (5) Convergence check: gate on the ABSOLUTE-count residual (the quantity Newton
        // minimizes), recomputed FRESH from the post-step w[] (b was overwritten by
        // cholesky_solve into delta_lambda; the max_b tracker at iter-top is one iter lagged).
        double max_abs_resid = 0.0;
        for (int k = 0; k < K; k++) {
            for (int j = 0; j < st.cat_counts[k]; j++) {
                double target = st.targets[k][j] * static_cast<double>(st.n);
                double S_kj = 0.0;
                for (int c : cells_per_cat[k][j]) S_kj += w[c];
                max_abs_resid = std::max(max_abs_resid, std::fabs(target - S_kj));
            }
        }
        // CR-C7 (kxna.7): incorporate the POST-step residual into best-iterate
        // tracking. The iter-top update (pre-step max_b) lags by one iteration and
        // never sees the FINAL accepted step's improvement, so on BUDGET/STALL exits
        // that step was silently discarded (logit returns w_best) and the stale
        // best_resid biased the BUDGET-vs-STALL classification. max_abs_resid is the
        // fresh post-step absolute-count residual.
        if (max_abs_resid < best_resid) {
            best_resid = max_abs_resid;
            w_best = w;
            res.base.best_iter = iter + 1;
        }
        // eb79.22: absolute feasibility is SUFFICIENT for OK-convergence — margins met in
        // absolute-count space ⇒ solved. The prior gate ALSO AND-required the improvement/
        // plateau rule, which never fires while the residual is still monotonically shrinking;
        // cleanly-converging problems (esp. rank-deficient / collinear ones, which shrink to
        // machine-zero without plateauing) therefore reached the optimum yet ran to the
        // iteration cap and mis-reported BUDGET. eb79.16's scale-blind guard is PRESERVED: OK
        // still REQUIRES the absolute-count residual to be met, so a proportion-only match
        // under scale drift cannot certify success (the proportion metric — compute_cell_metrics
        // normalized by the current total W — is exactly what would falsely certify a stall as
        // max_error=0, so it must NOT gate OK). A metric plateau WITHOUT feasibility stays a
        // STALL (post-loop status below).
        bool converged = (max_abs_resid <= st.tol_abs * static_cast<double>(st.n));
        if (converged) {
            lbw::mark_converged(res, cfg, iter + 1, st.tol_abs * static_cast<double>(st.n));  // CR-C10b: governing threshold is tol_abs*n
            w_best = w;
            break;
        }
    }

    // Post-loop status
    if (res.base.status == RK_ERR_NOCONV) {
        res.base.status = (best_resid < initial_resid * 0.999) ? RK_ERR_BUDGET : RK_ERR_STALL;
        if (res.base.status == RK_ERR_STALL) res.base.stall_kind = 2;  // KL/metric plateau (BCD-class)
        std::snprintf(res.message, sizeof(res.message),
            "logit: %s after %d Newton steps",
            res.base.status == RK_ERR_BUDGET ? "budget exhausted" : "stall",
            res.base.iterations);
    }

    if (w_best.empty()) w_best = w;

    // CR-C7b (dtk8): reconstruct the obs-level weights from the best cell iterate, then
    // route them through finalize_weights_buf so the Σw=n + bounds_mode contract holds on
    // ALL exit paths (mirror greenkhorn kxna.1 / oris_finalize). Without this, logit's cell
    // iterate was returned unscaled → Σw≠n on BUDGET/STALL exit (measured 802/800), and the
    // reported cell best_error diverged from the residual the RETURNED weights realize.
    res.base.best_weights.resize(st.n);
    for (int i = 0; i < st.n; i++) {
        int c = ct.cell_of[i];
        res.base.best_weights[i] = (X_init[c] > 0.0)
            ? st.weights[i] * w_best[c] / X_init[c]
            : st.weights[i];
    }
    lbw::finalize_weights_buf(res.base.best_weights.data(), st.n, st, ct,
                              res.n_bounds_violated, res.n_bounds_clamped);

    // Final metrics — computed from the RETURNED weights (finalized obs weights remassed to
    // cells), NOT the pre-finalize cell iterate, so returned≡reported (mxcl.5) on BUDGET too.
    std::vector<double> w_ret(M, 0.0);
    for (int i = 0; i < st.n; i++) w_ret[ct.cell_of[i]] += res.base.best_weights[i];
    double W_best = 0.0;
    for (int c = 0; c < M; c++) W_best += w_ret[c];
    lbw::CellMetrics m_best = lbw::compute_cell_metrics(st, ct, w_ret, W_best, bucket_scratch);
    // eb79.16: recompute the RESIDUAL-class reported fields in ABSOLUTE-count space
    // (pop = targets[k][j]*n, NOT *W). compute_cell_metrics normalizes by the current
    // total W, so on a rank-deficient stall it reports max_error=0 despite a large true
    // margin violation. Mirror EXACTLY the per-metric aggregation compute_cell_metrics
    // uses, only against the true target total n:
    //   errRp     -> MAX over all (k,j) of r_kj                       (max_error/best_error)
    //   mean_err  -> MEAN over K margins of (MAX over j within margin) of r_kj
    //   grake     -> MAX over all (k,j) of r_kj / (1 + targets*n)
    // where r_kj = |targets[k][j]*n - Σ_{c in bucket} w_ret[c]| (RETURNED cell masses).
    // kl and chi2 are proportion-space divergence metrics (compare probability mass) —
    // legitimately scale-relative — so they are LEFT UNCHANGED from compute_cell_metrics.
    const double n_d = static_cast<double>(st.n);
    double abs_max      = 0.0;   // errRp analogue: global MAX of r_kj
    double abs_mean_sum = 0.0;   // mean_err analogue: sum over margins of per-margin MAX
    double abs_grake    = 0.0;   // grake_norm analogue: global MAX of normalized r_kj
    for (int k = 0; k < K; k++) {
        double max_k = 0.0;      // per-margin MAX over categories (mirror compute_cell_metrics)
        for (int j = 0; j < st.cat_counts[k]; j++) {
            double pop = st.targets[k][j] * n_d;
            double S_kj = 0.0;
            for (int c : cells_per_cat[k][j]) S_kj += w_ret[c];
            double r_kj = std::fabs(pop - S_kj);
            if (r_kj > max_k)   max_k = r_kj;
            if (r_kj > abs_max) abs_max = r_kj;
            double nm = r_kj / (1.0 + std::fabs(pop));
            if (nm > abs_grake) abs_grake = nm;
        }
        abs_mean_sum += max_k;
    }
    res.base.max_error  = abs_max / n_d;
    res.base.best_error = abs_max / n_d;
    res.base.mean_error = (K > 0) ? (abs_mean_sum / static_cast<double>(K)) / n_d : 0.0;
    res.base.kl         = m_best.kl;
    res.base.chi2       = m_best.chi2;
    res.base.grake_norm = abs_grake;

    // eb79.20: convergence_solver_objective = Σ_c D_L(w_ret[c]), the midpoint-referenced
    // Fermi-Dirac logit distance the solver minimizes, evaluated on the RETURNED cell
    // masses. Its gradient log((w−L)/(U−w)) = z matches the forward map w = L + range·σ(z)
    // (main loop above), so the reference is the MIDPOINT (σ(0)=½ ⇒ w=(L+U)/2), NOT X_init —
    // the sigmoid carries no z_target offset (lambda_0=z_target is a warm start only). Closed
    // form: D_L(w) = range·[p·log p + (1−p)·log(1−p) + log2],  p = (w−L)/range ∈ [0,1]
    // = range·(log2 − H(p)) ∈ [0, range·log2]; convex, 0 at the midpoint. Derivation
    // independently verified (agy). select_solver_objective (calib_dispatch.hpp) does NOT
    // cover logit (default branch returns NaN) — computed inline.
    double solver_obj = 0.0;
    for (int c = 0; c < M; c++) {
        double range = U_cell[c] - L_cell[c];
        if (range < 1e-12) continue;                  // degenerate L==U cell: contributes 0
        double p = std::clamp((w_ret[c] - L_cell[c]) / range, 0.0, 1.0);
        double t = std::log(2.0);                     // log2 − H(½) = 0 at the midpoint
        if (p > 1e-300)         t += p * std::log(p);           // x·log x → 0 as x → 0
        if (1.0 - p > 1e-300)   t += (1.0 - p) * std::log(1.0 - p);
        solver_obj += range * t;
    }
    res.base.convergence_solver_objective = solver_obj;

    // eb79.20: l1_weight_change = Σ_i|Δw|/n (leafblower.h:133), obs-level calibrated −
    // input design weight. Both vectors are length st.n (obs-level) regardless of
    // bounds_mode, so the n denominator is st.n unambiguously.
    double l1_wc = 0.0;
    for (int i = 0; i < st.n; i++)
        l1_wc += std::fabs(res.base.best_weights[i] - st.weights[i]);
    res.base.l1_weight_change = l1_wc / static_cast<double>(st.n);

    // kxna.20 (CR-C7c): in bounds_mode="unit" the per-obs water-fill above can leave a
    // fully-pinned cell off-target after status was set RK_OK on the pre-finalize cell
    // iterate. res.base.max_error is the RETURNED margin error (abs/n_d, proportion,
    // post-finalize at L553), so re-gate on it directly (no second recompute).
    lbw::regate_unit_status(res.base, st, res.base.max_error);

    return res;
}

} // namespace lbw
