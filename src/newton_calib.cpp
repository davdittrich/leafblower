#include "newton_calib.hpp"
#include "leafblower.h"
#include "calib_dispatch.hpp"
#include "calib_linalg.hpp"
#include "cell_table.hpp"
#include "calib_validate.hpp"
#ifndef LBW_NO_R
#  include <R_ext/Lapack.h>
#  include <R_ext/Print.h>  // Rprintf (Epic-H WH-e diagnostics)
#  include <R_ext/RS.h>     // F77_CALL
#  ifndef FCONE
#    define FCONE
#  endif
#else
#  include <cstdio>
// System LAPACK declarations for Python build
extern "C" {
    void dsyevd_(const char* jobz, const char* uplo, int* n, double* A, int* lda,
                 double* w, double* work, int* lwork, int* iwork,
                 int* liwork, int* info);
}
#  define F77_CALL(x) x ## _
#  define Rprintf(...) fprintf(stderr, __VA_ARGS__)
#  define FCONE
#endif
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
    res.base.convergence_solver_objective = std::numeric_limits<double>::infinity();

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
    // CR-C13 (kxna.13): the old `if (n_lam > 4096) BADARG` guard was unreachable and
    // is removed. solver_setup_ct (above, return-checked) runs calib_validate_preentry,
    // which rejects n_cats_total > kNCatsTotalMax=2048 (calib_validate.cpp) before we get
    // here. n_lam = Σ_k max(0, cat_counts[k]−1) ≤ Σ_k cat_counts[k] = n_cats_total ≤ 2048
    // < 4096, so the branch could never fire (bound holds regardless of per-margin counts).

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

    // WL-3 trust-region radius (Steihaug-CG step). Persistent across inner-Newton
    // iters within one solve (init at solve start); adaptive via Marquardt
    // ρ-ratio. NO upper cap on growth (plan rev-2 fix H). Composes additively
    // with LM lm_mu adaptation — both use the same actual/predicted ρ but each
    // updates its own state.
    double delta_radius = 1.0;

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
    // CR-C17 (kxna.17): u_max is the LSE stabilization shift and must be the max u
    // over rows that CONTRIBUTE to Z_stable = Σ st.weights[i]·exp(u_i−u_max) — i.e. the
    // positive-design-weight rows. A zero-weight (d_i=0) row adds 0 to Z_stable, so if
    // it held the max u it would only inflate the shift; and if its u exceeded every
    // positive row by >745, Z_stable would underflow to exactly 0 (log(0)=−inf;
    // scale=n/0=inf; w_i=0·inf=NaN in recovery). Filtering d_i>0 here makes the shift
    // correct BY CONSTRUCTION. Zero-weight rows are excluded from the ENTIRE LSE — this
    // scan, the Z_stable sums (eval_dual + step 5b) and weight recovery — so they have
    // zero effect: no underflow (positive rows define the shift) and no 0·inf overflow
    // (their terms are never evaluated). Byte-identical whenever a positive row already
    // holds the max (the norm). Σd_i>0 is validated upstream, so ≥1 row survives.
    auto compute_u_max = [&](const std::vector<double>& lam_v) -> double {
        double u_max = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < n; ++i) {
            if (st.weights[i] <= 0.0) continue;   // only positive-weight rows enter Z_stable
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
            if (st.weights[i] <= 0.0) continue;   // CR-C17: d=0 row adds 0; skip (no 0*inf)
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
    // Steihaug-CG on a DIAGONAL PSD operator H_diag (eigenbasis). Solves
    // H_diag·x ≈ g_kp under the trust constraint ||x||₂ ≤ Δ. Returns the boundary
    // intersection point if the iterate would exit the trust region, otherwise
    // returns the converged interior solution.
    //
    // Mathematics (no negative-curvature branch — Λ_damped > 0 by construction):
    //   x₀ = 0, r₀ = g_kp − H·x₀ = g_kp, p₀ = r₀
    //   for k = 0..max_iter:
    //       Hp = H_diag .* p
    //       α  = (rᵀr) / (pᵀHp)
    //       if ||x + α·p||₂ ≥ Δ: solve quadratic for τ ∈ (0,α], return x + τ·p
    //       x ← x + α·p
    //       r ← r − α·Hp
    //       if ||r||₂ < tol_cg: return x
    //       β  = (r_newᵀr_new) / (rᵀr); p ← r + β·p
    //
    // PSD ⇒ pᵀHp > 0 always (no division-by-zero). max_iter = 2·n_kp safety cap;
    // exact arithmetic CG converges in ≤n_kp iters but FP roundoff may cost a
    // small extra factor.
    auto steihaug_cg = [](const std::vector<double>& g_kp,
                          const std::vector<double>& H_diag,
                          double Delta,
                          double tol_cg,
                          int max_iter) -> std::vector<double> {
        const int n_kp = static_cast<int>(g_kp.size());
        std::vector<double> x(n_kp, 0.0);
        if (n_kp == 0) return x;
        std::vector<double> r(g_kp);            // r₀ = g_kp (since x₀ = 0)
        std::vector<double> p(g_kp);            // p₀ = r₀
        double r_dot_r = 0.0;
        for (int i = 0; i < n_kp; ++i) r_dot_r += r[i] * r[i];
        std::vector<double> Hp(n_kp);
        for (int k = 0; k < max_iter; ++k) {
            double p_dot_Hp = 0.0;
            for (int i = 0; i < n_kp; ++i) {
                Hp[i] = H_diag[i] * p[i];
                p_dot_Hp += p[i] * Hp[i];
            }
            // PSD ⇒ p_dot_Hp ≥ 0; Λ_damped > 0 strictly (lm_mu > 0 + retained
            // eigenvalues ≥ ratio·λ_max), so p_dot_Hp > 0 unless p = 0.
            if (p_dot_Hp <= 0.0) return x;  // defensive; should not occur
            const double alpha = r_dot_r / p_dot_Hp;
            // Trust-boundary check: if ||x + α·p||₂ ≥ Δ, solve quadratic.
            double xx = 0.0, xp = 0.0, pp = 0.0;
            for (int i = 0; i < n_kp; ++i) {
                xx += x[i] * x[i];
                xp += x[i] * p[i];
                pp += p[i] * p[i];
            }
            // ||x + α·p||² = pp·α² + 2·xp·α + xx
            const double xn2 = pp * alpha * alpha + 2.0 * xp * alpha + xx;
            if (xn2 >= Delta * Delta) {
                // Solve pp·τ² + 2·xp·τ + (xx − Δ²) = 0 for τ ∈ (0, α].
                // Take the larger root (gives the further intersection along p).
                const double aq = pp;
                const double bq = 2.0 * xp;
                const double cq = xx - Delta * Delta;
                const double disc = std::max(0.0, bq * bq - 4.0 * aq * cq);
                // CR-C14 (kxna.14): numerically stable positive root. The naive
                // (−bq+√disc)/(2aq) loses precision to catastrophic cancellation when
                // bq>0 and the CG iterate already sits near the trust boundary
                // (cq=xx−Δ²→0⁻ ⇒ √disc→|bq|). Citardauq form (Numerical Recipes §5.6):
                // q = −½(bq + sign(bq)·√disc); the two roots are q/aq and cq/q. Here
                // aq=pp>0 and cq=xx−Δ²<0 (the iterate is strictly interior on entry), so
                // the product cq/aq<0 ⇒ exactly one root is positive; pick it by sign(q),
                // which never subtracts near-equal magnitudes.
                const double q = -0.5 * (bq + std::copysign(std::sqrt(disc), bq));
                const double tau = (q > 0.0) ? (q / aq)
                                 : (q < 0.0) ? (cq / q)
                                 : (-bq / (2.0 * aq));  // q==0 ⇔ bq==cq==0 (both roots 0)
                for (int i = 0; i < n_kp; ++i) x[i] += tau * p[i];
                return x;
            }
            // Interior step: accept full α.
            for (int i = 0; i < n_kp; ++i) {
                x[i] += alpha * p[i];
                r[i] -= alpha * Hp[i];
            }
            double r_dot_r_new = 0.0;
            for (int i = 0; i < n_kp; ++i) r_dot_r_new += r[i] * r[i];
            if (std::sqrt(r_dot_r_new) < tol_cg) return x;
            const double beta = r_dot_r_new / r_dot_r;
            for (int i = 0; i < n_kp; ++i) p[i] = r[i] + beta * p[i];
            r_dot_r = r_dot_r_new;
        }
        return x;
    };

    auto run_newton_inner = [&](const std::vector<double>& T_eps,
                                int max_iter_inner) -> bool {
        int consecutive_failed = 0;  // consecutive line-search failures; cap at 3
        // Best-iterate tracking: rank-deficient near-optimum can drift λ to ∞ even
        // while Armijo accepts (model under-predicts decrease). Save λ at lowest
        // gap seen; restore before recovery if the loop exits in a worse state.
        // lam_best/best_Z/best_u_max: dual-space state for λ recovery (dual
        //   iterate is separate from weight-space output). best_gap is the
        //   lowest dual_gap seen; it also feeds convergence_solver_objective
        //   directly (CR-H7 removed the redundant BestIterTracker copy — it
        //   recorded dual_gap==best_gap on every iter, so its best_objective
        //   was provably identical to best_gap).
        std::vector<double> lam_best(lam);
        double best_gap     = std::numeric_limits<double>::infinity();
        double best_Z       = Z_curr;
        double best_u_max   = u_max_curr;
        int    best_iter_id = 0;

        // dk9l.2: hoist dsyevd workspace query + allocation outside iter loop.
        // Workspace sizes depend only on n_lam (constant per solve), so the
        // lwork=-1 query and work/iwork vector allocations are invariant.
        // The actual eigendecomposition still runs every iter (eigvals/vecs
        // change with H — must recompute). dsyevd writes work[0]/iwork[0] on
        // query. CR-H8: the lwork=-1 query path never dereferences A or W
        // (netlib dsyevd.f: arg-check lda>=n, compute LWMIN/LIWMIN, set
        // WORK(1)/IWORK(1), RETURN — all before any A/W access), so a
        // single-element scratch with lda=n_lam suffices; avoids a transient
        // n_lam² (up to ~33MB) allocation that would otherwise never be read.
        int dsy_lwork_cache = -1, dsy_liwork_cache = -1;
        std::vector<double> dsy_work;
        std::vector<int>    dsy_iwork;
        {
            int dsy_n_q = n_lam, dsy_lda_q = n_lam, dsy_info_q = 0;
            int q_lwork = -1, q_liwork = -1;
            double work_query = 0.0;
            int    iwork_query = 0;
            double H_query_scratch = 0.0;      // query never reads A (lda=n_lam)
            double eigvals_query_scratch = 0.0; // query never reads W
            F77_CALL(dsyevd)("V", "U", &dsy_n_q, &H_query_scratch, &dsy_lda_q,
                             &eigvals_query_scratch,
                             &work_query, &q_lwork, &iwork_query, &q_liwork,
                             &dsy_info_q FCONE FCONE);
            if (dsy_info_q == 0) {
                dsy_lwork_cache  = static_cast<int>(work_query);
                dsy_liwork_cache = iwork_query;
            } else {
                // Workspace query failed — pick safe upper bounds per LAPACK doc:
                //   lwork  >= 1 + 6N + 2N²   (jobz='V')
                //   liwork >= 3 + 5N
                dsy_lwork_cache  = 1 + 6 * n_lam + 2 * n_lam * n_lam;
                dsy_liwork_cache = 3 + 5 * n_lam;
            }
            dsy_work.resize(std::max(1, dsy_lwork_cache));
            dsy_iwork.resize(std::max(1, dsy_liwork_cache));
        }

        // lqex.1: hoist per-iter vector allocations outside Newton loop.
        // H_pre, H_eigvecs, eigvals, active_a, keep_idx, g_keep,
        // lambda_damped, delta_keep_pinv, delta_keep, lam_trial declared
        // inside the loop — reallocation every iter at O(n_lam^2) cost.
        // Declared here; resize/assign at iter start. Results bit-identical.
        std::vector<double> H_pre(static_cast<size_t>(n_lam) * n_lam);
        std::vector<double> H_eigvecs(static_cast<size_t>(n_lam) * n_lam);
        std::vector<double> eigvals(n_lam);
        std::vector<int>    active_a;    active_a.reserve(K);
        std::vector<int>    keep_idx;    keep_idx.reserve(n_lam);
        std::vector<double> g_keep;      g_keep.reserve(n_lam);
        std::vector<double> lambda_damped; lambda_damped.reserve(n_lam);
        std::vector<double> delta_keep_pinv; delta_keep_pinv.reserve(n_lam);
        std::vector<double> delta_keep;  delta_keep.reserve(n_lam);
        std::vector<double> lam_trial(n_lam);

        int iter = 0;
        for (iter = 0; iter < max_iter_inner; ++iter) {

            // WL-3 per-iter trust-region carryover (reset each iter).
            // m_pred_dec_for_trust: model-predicted decrease in eigenbasis,
            //   used as ρ denominator for adaptive Δ update.
            // trust_step_taken: true only when TSVD path computed δ_keep
            //   (fast path or Steihaug-CG); false on LDLT-fallback paths
            //   where Δ is left unchanged (no eigenbasis predicted-decrease).
            double m_pred_dec_for_trust = 0.0;
            bool   trust_step_taken     = false;

            // 5a. Zero accumulators for this step.
            std::fill(G.begin(), G.end(), 0.0);
            std::fill(H.begin(), H.end(), 0.0);
            double Z = 0.0;

            // 5b. LSE prep: scan u_max once, then accumulate using exp(u_i − u_max).
            // u_max factor cancels in all ratios (G/Z, H/Z) so this is a pure
            // numerical-stability shift. Cost: one extra O(n·K) pass per Newton iter.
            //
            // dk9l.1: hoist per-obs quantities out of inner K×K Hessian loop.
            //   - active_a[]: per-obs list of active dual-variable indices a =
            //     lam_off[k] + j - 1 for margins where j > 0 (j == 0 is the
            //     reference category, λ fixed to 0).
            //   - Built once per obs in increasing-k order. Since lam_off is
            //     strictly increasing across margins (disjoint index ranges),
            //     active_a is sorted ⇒ p1 ≤ p2 ⇒ active_a[p1] < active_a[p2],
            //     so no min/max needed for upper-triangle storage.
            //   - Reused for both gradient and Hessian accumulation.
            //   - Eliminates redundant group_ids[k2][i] re-fetch and skips
            //     inactive (j == 0) margins entirely. Drops the inner-loop
            //     work from O(n·K²) to O(n·K_active²).
            //   - u (sum of active λ entries) is accumulated in the same pass
            //     instead of via compute_u() (which iterates k=0..K again).
            const double u_max_step = compute_u_max(lam);
            // active_a hoisted pre-loop (lqex.1); capacity already reserved.
            for (int i = 0; i < n; ++i) {
                if (st.weights[i] <= 0.0) continue;   // CR-C17: d=0 row contributes 0 to Z/G/H
                // Per-obs hoist: active dual-index list and u in one pass.
                double u = 0.0;
                active_a.clear();
                for (int k = 0; k < K; ++k) {
                    const int j = st.group_ids[k][i];
                    if (j > 0) {
                        const int a_k = lam_off[k] + j - 1;
                        u += lam[a_k];
                        active_a.push_back(a_k);
                    }
                }
                const double f_i = st.weights[i] * std::exp(u - u_max_step);
                Z += f_i;

                const int n_active = static_cast<int>(active_a.size());

                // Gradient accumulation (active margins only).
                for (int p = 0; p < n_active; ++p) {
                    G[active_a[p]] += f_i;
                }

                // Hessian accumulation (upper triangle, active×active only).
                // active_a is strictly increasing in p ⇒ a1 < a2 for p1 < p2.
                for (int p1 = 0; p1 < n_active; ++p1) {
                    const int a1 = active_a[p1];
                    H[a1 * n_lam + a1] += f_i;  // diagonal
                    const size_t row_off = static_cast<size_t>(a1) * n_lam;
                    for (int p2 = p1 + 1; p2 < n_active; ++p2) {
                        H[row_off + active_a[p2]] += f_i;
                    }
                }
            }
            Z_curr     = Z;            // = Σ d_i exp(u_i − u_max_step)
            u_max_curr = u_max_step;   // shift used during this step's accumulation

            // 5c. Normalize G and H by Z.
            const double inv_Z = 1.0 / Z;
            for (int a = 0; a < n_lam; ++a) G[a] *= inv_Z;
            for (size_t idx = 0; idx < H.size(); ++idx) H[idx] *= inv_Z;

            // 5d. Form the Schur-complement Hessian H = M2 − G⊗G (covariance, upper tri).
            // CR-C18 (kxna.18): the DIAGONAL variance p_a − p_a² is a catastrophic
            // cancellation as p_a→1 (near-equal-magnitude subtraction). Compute it
            // cancellation-free as p_a·(1−p_a) — one rounding (~0.5 ulp) instead of the
            // subtraction's O(1/(1−p)) relative blow-up. Pre-subtraction H[a][a] holds
            // M2[a][a]=p_a=G[a], so the variance is exactly G[a]·(1−G[a]). Off-diagonals
            // p_ab − p_a·p_b are a genuine difference of two INDEPENDENT quantities (not a
            // p→1 self-cancellation), with no cheap algebraic reformulation; their residual
            // noise (~eps·O(1)) stays far below the TSVD threshold ratio_tsvd·λ_max
            // (ratio_tsvd=1e-8, :480/:510) and is left as the plain subtraction.
            for (int a = 0; a < n_lam; ++a) {
                const size_t row_off = static_cast<size_t>(a) * n_lam;
                H[row_off + a] = G[a] * (1.0 - G[a]);        // cancellation-free variance
                for (int b = a + 1; b < n_lam; ++b)
                    H[row_off + b] -= G[a] * G[b];            // off-diagonal covariance
            }

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
            // best_gap doubles as the solver objective for output; actual obs
            // weights are computed post-loop from lam.
            if (dual_gap < best_gap) {
                best_gap     = dual_gap;
                best_Z       = Z;
                best_u_max   = u_max_step;
                best_iter_id = iter;
                lam_best     = lam;
            }

            if (dual_gap < tol_abs) {
                res.lm_mu_final = lm_mu;
                mark_converged(res, st.convergence_cfg, iter, tol_abs);  // CR-C10b: local tol_abs is the governing fallback (1e-6 floor)
                break;
            }

            // 5h. LM damped Newton direction via truncated-SVD pseudoinverse.
            // Mechanism (Epic-Dβ WL-1):
            //   1. Eigendecompose pre-damp PSD H = V Λ Vᵀ via dsyevd.
            //   2. Threshold-truncate retained set R = {i : Λᵢ ≥ ratio · λ_max}.
            //   3. Damp in eigenbasis: Λ_dampᵢ = Λᵢ·(1+μ) + μ·d_floor_retained.
            //   4. Project G into eigenbasis (only retained columns), divide,
            //      back-project: δ = Σ_{i∈R} V[:,i] · (Vᵀ[:,i]·G) / Λ_dampᵢ.
            //   5. Null-space components (i∉R) contribute zero ⇒ λ unchanged
            //      in null space (well-defined behavior under rank deficiency).
            // Edge cases:
            //   • λ_max ≤ 0 (degenerate / all-zero H): δ ← 0, n_projected = n_λ.
            //   • n_keep == 0 (spectrum entirely below threshold): fall back to
            //     LM-LDLT path; n_projected_dims = n_λ.
            //   • dsyevd info ≠ 0: same fallback.
            //   • n_keep == n_λ: equivalent to plain LDLT (no truncation).
            // H_pre hoisted pre-loop (lqex.1); assign from H each iter.
            H_pre.assign(H.begin(), H.end());  // n_lam×n_lam copy
            // Tikhonov ridge: add τI to H_pre before LM damping and dsyevd.
            if (st.ridge_lambda > 0.0) {
                for (int k = 0; k < n_lam; k++)
                    H_pre[k * n_lam + k] += st.ridge_lambda;
            }

            // Epic-H WH-e: user-tunable TSVD truncation ratio; <=0 falls back to internal default 1e-8.
            const double ratio_tsvd = st.newton_tsvd_ratio > 0.0 ? st.newton_tsvd_ratio : 1e-8;
            // H_eigvecs, eigvals hoisted pre-loop (lqex.1); reinit each iter.
            H_eigvecs.assign(H_pre.begin(), H_pre.end());  // dsyevd overwrites in place with V (column-major)
            eigvals.assign(n_lam, 0.0);

            int dsy_n = n_lam, dsy_lda = n_lam, dsy_info = 0;
            // dk9l.2: workspace hoisted out of iter loop — query + allocation
            // ran every iter despite depending only on n_lam. Reuse cached
            // dsy_lwork_cache/dsy_liwork_cache and dsy_work/dsy_iwork buffers.
            int dsy_lwork  = dsy_lwork_cache;
            int dsy_liwork = dsy_liwork_cache;
            F77_CALL(dsyevd)("V", "U", &dsy_n, H_eigvecs.data(), &dsy_lda, eigvals.data(),
                             dsy_work.data(), &dsy_lwork, dsy_iwork.data(), &dsy_liwork,
                             &dsy_info FCONE FCONE);

            bool tsvd_used = false;   // if false → fall back to LDLT path below
            int  n_keep    = 0;

            if (dsy_info == 0) {
                // dsyevd returns eigvals in ascending order; λ_max = eigvals[n_lam-1].
                double lam_max = eigvals[n_lam - 1];
                if (lam_max <= 0.0) {
                    // Epic-H WH-e I2: surface degenerate spectrum on stderr-equivalent (Rprintf).
                    Rprintf("[newton_kl] degenerate lambda_max=%.3e <= 0; skipping TSVD\n", lam_max);
                    // Degenerate spectrum (all-zero or all-negative due to FP noise).
                    std::fill(delta.begin(), delta.end(), 0.0);
                    res.n_projected_dims = n_lam;
                    tsvd_used = true;
                    n_keep    = 0;
                } else {
                    const double thresh = ratio_tsvd * lam_max;
                    // Retained set: R = {i : eigvals[i] >= thresh}.
                    // Compute mean over retained eigenvalues for d_floor_retained.
                    double sum_keep = 0.0;
                    for (int i = 0; i < n_lam; ++i)
                        if (eigvals[i] >= thresh) { ++n_keep; sum_keep += eigvals[i]; }

                    if (n_keep == 0) {
                        // Spectrum entirely below threshold — Newton step undefined under
                        // truncation. Fall back to LM-LDLT path on H_pre.
                        res.n_projected_dims = n_lam;
                        tsvd_used = false;
                    } else {
                        const double d_floor_retained = sum_keep / static_cast<double>(n_keep);
                        // dsyevd column-major: V[:,i] occupies H_eigvecs[i*n_lam .. i*n_lam+n_lam-1].
                        //
                        // WL-3: build eigenbasis state (g_keep, lambda_damped) over retained
                        // dims, compute δ_keep_pinv (pseudoinverse step in eigenbasis),
                        // branch on ||δ_keep_pinv||₂ vs Δ:
                        //   ||·||₂ ≤ Δ  → fast path: δ_keep = δ_keep_pinv
                        //   ||·||₂ >  Δ → Steihaug-CG on diagonal Λ_damped to trust boundary
                        // Then back-project δ_keep into ambient space.
                        // keep_idx, g_keep, lambda_damped, delta_keep_pinv hoisted
                        // pre-loop (lqex.1); resize to actual n_keep each iter.
                        keep_idx.clear();
                        g_keep.resize(n_keep);
                        lambda_damped.resize(n_keep);
                        delta_keep_pinv.resize(n_keep);
                        int kp = 0;
                        for (int i = 0; i < n_lam; ++i) {
                            if (eigvals[i] < thresh) continue;
                            keep_idx.push_back(i);
                            // g_keep_i = V[:,i] · G  (project gradient into eigenbasis)
                            double gi = 0.0;
                            const double* vi = H_eigvecs.data() + static_cast<size_t>(i) * n_lam;
                            for (int j = 0; j < n_lam; ++j) gi += vi[j] * G[j];
                            const double lam_damp = eigvals[i] * (1.0 + lm_mu) + lm_mu * d_floor_retained;
                            g_keep[kp]          = gi;
                            lambda_damped[kp]   = lam_damp;
                            delta_keep_pinv[kp] = gi / lam_damp;  // pseudoinverse step
                            ++kp;
                        }
                        // Trust-region branch on pinv-step norm.
                        double pinv_norm_sq = 0.0;
                        for (int i = 0; i < n_keep; ++i)
                            pinv_norm_sq += delta_keep_pinv[i] * delta_keep_pinv[i];
                        const double pinv_norm = std::sqrt(pinv_norm_sq);

                        // delta_keep hoisted pre-loop (lqex.1); resize to n_keep.
                        delta_keep.resize(n_keep);
                        if (pinv_norm <= delta_radius) {
                            // Fast path: pinv step inside trust region — use it directly.
                            delta_keep = delta_keep_pinv;
                        } else {
                            // Steihaug-CG on diagonal Λ_damped to trust boundary Δ.
                            delta_keep = steihaug_cg(g_keep, lambda_damped, delta_radius,
                                                     /*tol_cg=*/1e-10,
                                                     /*max_iter=*/2 * n_keep);
                        }
                        // Back-project: δ[j] = Σ_{i∈R} V[j,i] · δ_keep[i].
                        std::fill(delta.begin(), delta.end(), 0.0);
                        for (int kp2 = 0; kp2 < n_keep; ++kp2) {
                            const int i = keep_idx[kp2];
                            const double* vi = H_eigvecs.data() + static_cast<size_t>(i) * n_lam;
                            const double dki = delta_keep[kp2];
                            for (int j = 0; j < n_lam; ++j) delta[j] += vi[j] * dki;
                        }
                        // WL-3: predicted decrease in eigenbasis (H is diagonal there):
                        //   m(δ_keep) = -g_keepᵀ·δ_keep - 0.5·δ_keepᵀ·(Λ_damped .* δ_keep)
                        // Used in ρ-ratio for adaptive Δ update post-Armijo accept.
                        double pred_dec_eig = 0.0;
                        for (int i = 0; i < n_keep; ++i) {
                            pred_dec_eig += -g_keep[i] * delta_keep[i]
                                            - 0.5 * delta_keep[i] * lambda_damped[i] * delta_keep[i];
                        }
                        // Stash for ρ computation post-Armijo.
                        m_pred_dec_for_trust = pred_dec_eig;
                        trust_step_taken     = true;
                        res.n_projected_dims = n_lam - n_keep;
                        tsvd_used = true;
                    }
                }
            } else {
                // Epic-H WH-e I1: surface dsyevd failure before LDLT fallback.
                Rprintf("[newton_kl] dsy_info=%d (LAPACK dsyevd failure)\n", dsy_info);
                // dsyevd failed — fall back to LDLT.
                res.n_projected_dims = n_lam;
                tsvd_used = false;
            }

            if (!tsvd_used) {
                // LDLT fallback path (preserves pre-WL-1 behavior for n_keep==0
                // and dsyevd-failure edge cases).
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
                    if (cholesky_factor_inplace(H.data(), static_cast<size_t>(n_lam), 1e-12) != RK_OK) {
                        lm_mu = std::min(lm_mu_max, lm_mu * 10.0);
                        continue;
                    }
                    cholesky_solve(H.data(), static_cast<size_t>(n_lam), delta.data());
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
                    res.base.best_error = best_gap;
                    res.base.best_iter  = best_iter_id;
                    res.base.convergence_solver_objective = best_gap;
                    return false;
                }
            }

            // δᵀ H_pre δ using pre-damp (PSD) Hessian for the gain-ratio denominator.
            // PSD is the mathematical truth; under severe rank deficiency (severe-skew
            // K=20, many empty cells), Schur-complement subtraction can produce
            // tiny-negative results from FP rounding — those we clamp to 0. But a
            // genuinely negative δᵀHδ (D3) means the model predicts ascent at this δ;
            // reject the step and shrink the trust region instead of pretending it's 0.
            //
            // lqex.2: O(n_keep) shortcut when TSVD path was used.
            // δ = Σ_{kp} V[:,keep_idx[kp]] · delta_keep[kp] (back-projected), so:
            //   δᵀ H_pre δ = Σ_{kp} eigvals[keep_idx[kp]] · delta_keep[kp]²
            // (null-space components carry zero eigenvalue contribution). O(n_keep)
            // vs O(n_lam²) dense matvec. LDLT fallback path uses dense form unchanged.
            double delta_H_delta = 0.0;
            if (tsvd_used && n_keep > 0) {
                for (int kp = 0; kp < n_keep; ++kp)
                    delta_H_delta += eigvals[keep_idx[kp]] * delta_keep[kp] * delta_keep[kp];
            } else {
                for (int a = 0; a < n_lam; ++a)
                    for (int b = 0; b < n_lam; ++b)
                        delta_H_delta += delta[a] * H_pre[a * n_lam + b] * delta[b];
            }

            // G·δ = G^T H_damp^{-1} G ≥ 0 (descent curvature).
            double g_dot_d = 0.0;
            for (int a = 0; a < n_lam; ++a) g_dot_d += G[a] * delta[a];

            // D3: distinguish FP rounding noise (~|g_dot_d|·eps) from true ascent.
            // If genuinely negative, skip the step and shrink trust region — do NOT
            // mask the issue by clamping to zero.
            if (delta_H_delta < -1e-10 * std::fabs(g_dot_d)) {
                delta_radius *= 0.25;
                lm_mu = std::min(lm_mu_max, lm_mu * 10.0);
                if (++consecutive_failed >= 3) {
                    res.base.status = RK_ERR_NOCONV;
                    res.lm_mu_final = lm_mu;
                    break;
                }
                continue;  // re-enter iter without updating λ
            }
            // Tiny-negative-from-rounding only: clamp to 0 for ρ denominator.
            if (delta_H_delta < 0.0) delta_H_delta = 0.0;

            // D4: Armijo descent guard. δ from PSD H_damp⁻¹ G satisfies G·δ ≥ 0
            // (descent for λ ← λ − α·δ); but TSVD truncation, Steihaug early-stop,
            // or null-space projection (n_keep == 0 → δ = 0) can yield a direction
            // with G·δ ≤ 0, in which case the Armijo test g_trial ≤ g_curr − c1·α·G·δ
            // degenerates (RHS ≥ g_curr). Fall back to steepest descent: δ = G,
            // G·δ = ||G||² > 0 (gap-test above already returned for ||G||_∞ < tol).
            const double g_norm_sq = [&]{
                double s = 0.0;
                for (int a = 0; a < n_lam; ++a) s += G[a] * G[a];
                return s;
            }();
            double delta_norm = 0.0;
            for (int a = 0; a < n_lam; ++a) delta_norm += delta[a] * delta[a];
            delta_norm = std::sqrt(delta_norm);
            const double g_norm = std::sqrt(g_norm_sq);
            // dqs3 audit (non-bug): delta_norm is a per-iteration local whose ONLY
            // consumer is this guard test, computed from the pre-swap delta — which is
            // correct, since the test decides whether to swap. After the swap below it
            // is never read again (no Armijo line-search or trust-region consumer; it is
            // re-derived fresh next iteration), so it is deliberately NOT recomputed for
            // the steepest-descent direction. Verified: sole uses are L671-675.
            if (g_dot_d <= 1e-12 * g_norm * delta_norm) {
                std::copy(G.begin(), G.end(), delta.begin());
                g_dot_d = g_norm_sq;
                // Recompute δᵀHδ for the steepest-descent direction so ρ is honest.
                delta_H_delta = 0.0;
                for (int a = 0; a < n_lam; ++a)
                    for (int b = 0; b < n_lam; ++b)
                        delta_H_delta += delta[a] * H_pre[a * n_lam + b] * delta[b];
                if (delta_H_delta < 0.0) delta_H_delta = 0.0;
                // Steepest-descent direction did not come from the trust path —
                // ρ_trust below would compare against a stale m_pred_dec_for_trust;
                // disable trust-region update for this iteration.
                trust_step_taken = false;
            }

            // 5i. Armijo backtracking line search.
            // Sufficient decrease: g(λ - α·δ) ≤ g(λ) - c1·α·(G·δ).
            double alpha = 1.0;
            // lam_trial hoisted pre-loop (lqex.1); no resize needed (n_lam fixed).
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
            // CXX.3: a non-positive predicted decrease means the quadratic model
            // does not predict improvement; clamping the denominator to 1e-300
            // (the old guard) inflated rho into the good-model branch and wrongly
            // tightened damping. Treat predicted<=0 as rho=0 (poor model ⇒ back off).
            const double predicted = alpha * g_dot_d
                                     - 0.5 * alpha * alpha * delta_H_delta;
            const double rho = (predicted <= 0.0)
                               ? 0.0
                               : (g_curr - g_trial) / predicted;

            if (alpha >= 0.999 && rho > 0.75) {
                // FULL_ACCEPT with good model quality — tighten damping.
                lm_mu = std::max(lm_mu_min, lm_mu / 3.0);
            } else if (rho < 0.25) {
                // Poor quadratic model (accepted Armijo but low gain) — back off.
                lm_mu = std::min(lm_mu_max, lm_mu * 10.0);
            }
            // else: BACKTRACKED or moderate ρ — keep lm_mu unchanged.

            // WL-3: adaptive trust radius Δ via Marquardt ρ-ratio.
            // Uses eigenbasis predicted-decrease m_pred_dec_for_trust (computed
            // when TSVD path took the step). On LDLT-fallback iters trust path
            // is inactive — leave Δ unchanged.
            //   ρ_trust > 0.75 ⇒ Δ ← 2·Δ   (NO upper cap, plan rev-2 fix H)
            //   ρ_trust < 0.25 ⇒ Δ ← Δ / 4
            //   otherwise      ⇒ Δ unchanged.
            // Composes additively with lm_mu adaptation above.
            if (trust_step_taken) {
                const double actual_dec  = g_curr - g_trial;
                const double rho_trust   = actual_dec / std::max(m_pred_dec_for_trust, 1e-300);
                if (rho_trust > 0.75) {
                    delta_radius *= 2.0;          // grow trust region (no upper cap)
                } else if (rho_trust < 0.25) {
                    delta_radius *= 0.25;         // shrink (Δ /= 4)
                }
            }

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
                mark_converged(res, st.convergence_cfg, iter, tol_abs);  // CR-C10b: local tol_abs is the governing fallback (1e-6 floor)
                break;
            }
        }
        // CR-C12 (kxna.12): on budget exhaustion the loop `for(iter=0; iter<max_iter_inner)`
        // exits with iter==max_iter_inner (the ++iter past the last EXECUTED step 0..max−1),
        // so iter+1 = max_iter_inner+1 over-counts the actual max_iter_inner steps. On a break
        // (iter<max_iter_inner) iter+1 is the true count. Cap at max_iter_inner to fix the
        // no-break case without disturbing the break paths.
        res.base.iterations = std::min(iter + 1, max_iter_inner);

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
        // best_error = best dual gap seen (Newton's objective); best_iter = that iteration.
        res.base.best_error = best_gap;
        res.base.best_iter  = best_iter_id;
        res.base.convergence_solver_objective = best_gap;
        return res.base.status == RK_OK;
    };

    bool ok = run_newton_inner(T, max_iter);
    if (!ok && res.base.status == RK_ERR_BADARG) {
        return res;  // do not proceed to weight recovery on bad state
    }

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
        if (st.weights[i] <= 0.0) { w[i] = 0.0; continue; }  // CR-C17: d=0 obs → w=0 (not in sample)
        const double u = compute_u(lam, i);
        const double w_i = st.weights[i] * std::exp(u - u_max_curr) * scale;
        w[i] = w_i;
        if (w_i > st.max_weight || w_i < st.min_weight) ++n_violated;
    }

    // Always surface violation count; 5% threshold gates status change only.
    res.n_bounds_violated = n_violated;
    const double frac_violated = static_cast<double>(n_violated) / n;
    if (frac_violated > 0.05) {
        res.base.status = RK_ERR_NOCONV;
        std::snprintf(res.message, sizeof(res.message),
            "newton_kl: %.1f%% of obs violate [min,max] weight bounds",
            frac_violated * 100.0);
    }

    // ── 7. Write calibrated weights and compute max_error ───────────────────
    // Guard: only write to caller's st.weights when no threshold violation.
    // When frac_violated > 0.05 (T4 contract), st.weights retains the
    // pre-violation iterate; res.base.best_weights is left empty/default.
    // max_error is still computed from w[] so callers get a finite diagnostic.
    const bool violation_threshold_crossed = (frac_violated > 0.05);
    if (!violation_threshold_crossed) {
        for (int i = 0; i < n; ++i) st.weights[i] = w[i];
    }

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
    res.base.max_error = max_err;
    if (!violation_threshold_crossed) {
        res.base.best_weights = w;
    }
    res.lm_mu_final = lm_mu;

    return res;
}

} // namespace lbw
