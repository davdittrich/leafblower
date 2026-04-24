#include "lbw_config.h"
#include "lbfgsb_solver.hpp"
#include "leafblower.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <vector>
#include <deque>
#include "lbw_math.hpp"

namespace lbw {

// Compute T_kj = targets[k][j] * W where W = sum(weights).
// Returns total W.
static double compute_targets_abs(const CalibState& st,
                                   std::vector<double>& T) {
    double W = 0.0;
    for (int i = 0; i < st.n; i++) W += st.weights[i];
    int off = 0;
    for (int k = 0; k < st.K; k++) {
        for (int j = 0; j < st.cat_counts[k]; j++) {
            T[off + j] = st.targets[k][j] * W;
        }
        off += st.cat_counts[k];
    }
    return W;
}

// Build offset array: offset[k] = index of first dual var for margin k.
static std::vector<int> build_offsets(const CalibState& st) {
    std::vector<int> off(st.K + 1, 0);
    for (int k = 0; k < st.K; k++) off[k+1] = off[k] + st.cat_counts[k];
    return off;
}

// Compute u_i = sum_k lambda[offset[k] + group_ids[k][i]]  (skip NA: g==-1).
static void compute_u(const CalibState& st, const std::vector<int>& off,
                       const std::vector<double>& lam, std::vector<double>& u) {
    std::fill(u.begin(), u.end(), 0.0);
    for (int k = 0; k < st.K; k++) {
        for (int i = 0; i < st.n; i++) {
            int g = st.group_ids[k][i];
            if (g >= 0) u[i] += lam[off[k] + g];
        }
    }
}

// phi_from_u: evaluate objective and gradient given pre-computed u.
// Avoids re-running compute_u; gradient in lambda-space still requires K*n pass.
static double phi_from_u(const CalibState& st,
                          const LinkFn& fn,
                          const std::vector<int>& off,
                          const std::vector<double>& lam,
                          const std::vector<double>& T,
                          const std::vector<double>& d,
                          std::vector<double>& grad,
                          const std::vector<double>& u) {
    int total = off[st.K];

    double obj = 0.0;
    for (int idx = 0; idx < total; idx++) obj += T[idx] * lam[idx];
    for (int i = 0; i < st.n; i++) obj -= d[i] * fn.H(u[i]);

    std::fill(grad.begin(), grad.end(), 0.0);
    for (int idx = 0; idx < total; idx++) grad[idx] = T[idx];
    for (int k = 0; k < st.K; k++) {
        for (int i = 0; i < st.n; i++) {
            int g = st.group_ids[k][i];
            if (g >= 0) grad[off[k] + g] -= d[i] * fn.F(u[i]);
        }
    }

    if (st.alm_mu > 0.0) {
        double sum_w = 0.0;
        for (int i = 0; i < st.n; i++) sum_w += d[i] * fn.F(u[i]);
        double residual = sum_w - static_cast<double>(st.n);
        double alm_scale = st.alm_lambda + st.alm_mu * residual;

        for (int k = 0; k < st.K; k++) {
            for (int i = 0; i < st.n; i++) {
                int g = st.group_ids[k][i];
                if (g >= 0)
                    grad[off[k] + g] += alm_scale * d[i] * fn.dF(u[i]);
            }
        }
        obj += st.alm_lambda * residual + (st.alm_mu / 2.0) * residual * residual;
    }

    return obj;
}

// phi(lambda) = sum_kj T_kj*lam_kj - sum_i d_i*H(u_i)
// grad[off_k+j] = T_kj - S_kj where S_kj = sum_{i:g_k(i)==j} d_i*F(u_i)
// Computes u from lam, then delegates to phi_from_u.
static double phi_and_grad(const CalibState& st,
                            const LinkFn& fn,
                            const std::vector<int>& off,
                            const std::vector<double>& lam,
                            const std::vector<double>& T,
                            const std::vector<double>& d,
                            std::vector<double>& grad,
                            std::vector<double>& u) {
    compute_u(st, off, lam, u);
    return phi_from_u(st, fn, off, lam, T, d, grad, u);
}

// Precompute per-observation directional derivative: du[i] = sum_k dir[off[k]+g_k(i)]
// Used to update u incrementally in Wolfe search: u_new = u_base + alpha * du
static void compute_du(const CalibState& st,
                        const std::vector<int>& off,
                        const std::vector<double>& dir,
                        std::vector<double>& du) {
    std::fill(du.begin(), du.end(), 0.0);
    for (int k = 0; k < st.K; k++) {
        for (int i = 0; i < st.n; i++) {
            int g = st.group_ids[k][i];
            if (g >= 0) du[i] += dir[off[k] + g];
        }
    }
}

static double maxAbs(const std::vector<double>& v) {
    double mx = 0.0;
    for (double x : v) mx = std::max(mx, std::fabs(x));
    return mx;
}

static double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    int n = (int)a.size();
    for (int i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

// L-BFGS 2-loop: given gradient g, produce ascent direction H*g.
static void lbfgs_direction(const std::deque<std::vector<double>>& svec,
                              const std::deque<std::vector<double>>& yvec,
                              const std::deque<double>& rho,
                              double gamma,
                              const std::vector<double>& g,
                              std::vector<double>& dir) {
    int m = (int)svec.size();
    int D = (int)g.size();
    dir = g;

    // Per-thread scratch, grown once per thread. m ≤ st.lbfgs_m (default 10,
    // runtime-configurable via leafblower.h:38). thread_local avoids a heap
    // alloc every outer iteration without capping m at a hard-coded bound.
    thread_local std::vector<double> alpha;
    if ((int)alpha.size() < m) alpha.resize(m);
    for (int i = m - 1; i >= 0; i--) {
        alpha[i] = rho[i] * dot(svec[i], dir);
        for (int j = 0; j < D; j++) dir[j] -= alpha[i] * yvec[i][j];
    }
    for (int j = 0; j < D; j++) dir[j] *= gamma;
    for (int i = 0; i < m; i++) {
        double beta = rho[i] * dot(yvec[i], dir);
        for (int j = 0; j < D; j++) dir[j] += (alpha[i] - beta) * svec[i][j];
    }
}

// Compute final calibrated weights and maximum calibration error.
static LBFGSResult compute_final_weights_and_error(
        CalibState& st, const LinkFn& fn,
        const std::vector<double>& d, const std::vector<int>& off,
        std::vector<double>& u, int iterations) {
    LBFGSResult res;
    res.iterations = iterations;

    for (int i = 0; i < st.n; i++) {
        double wi = d[i] * fn.F(u[i]);
        st.weights[i] = wi;
    }
    // Solver self-normalizes at exit: the outer lbfgsb_solve() wraps this call
    // and applies Σ w = n before returning (moved from wrapper 2026-04-24).
    // max_err below is computed as max_{k,j} |S_kj/Wn - target_kj|; S_kj and Wn
    // scale identically under uniform normalization, so max_err is
    // scale-invariant and remains valid after the outer normalize.

    double Wn = 0.0;
    for (int i = 0; i < st.n; i++) Wn += st.weights[i];
    double max_err = 0.0;
    for (int k = 0; k < st.K; k++) {
        std::vector<double> S(st.cat_counts[k], 0.0);
        for (int i = 0; i < st.n; i++) {
            int g = st.group_ids[k][i];
            if (g >= 0) S[g] += st.weights[i];
        }
        for (int j = 0; j < st.cat_counts[k]; j++) {
            max_err = std::max(max_err, std::fabs(S[j] / Wn - st.targets[k][j]));
        }
    }
    // Saturation diagnostic: when NOCONV and >50% of weights are pinned at
    // [L,U], the problem is likely infeasible given min/max_weight. Emit via
    // st.log (verbose-gated). Heuristic limitation: proxy tests d[i]·F(u[i])
    // against L, so saturated observations with d[i] ≫ 1 are false-negatives.
    if (max_err >= st.tol_abs && !fn.exponential) {
        int n_sat_low = 0, n_sat_high = 0;
        const double kSatTol = 1e-9;
        const double bw = fn.U - fn.L;
        for (int i = 0; i < st.n; i++) {
            if      (st.weights[i] <= fn.L + kSatTol * bw) n_sat_low++;
            else if (st.weights[i] >= fn.U - kSatTol * bw) n_sat_high++;
        }
        double sat_frac = static_cast<double>(n_sat_low + n_sat_high) / st.n;
        if (sat_frac > 0.5) {
            char msg[256];
            std::snprintf(msg, 256,
                "L-BFGS-B: %.1f%% of weights pinned at [L,U] bounds "
                "(low=%d high=%d); targets may be infeasible given "
                "min/max_weight.",
                100.0 * sat_frac, n_sat_low, n_sat_high);
            st.log(msg);
        }
    }
    res.max_error = max_err;

    // WU-B: compute pct_change (start weights d[i] vs. final st.weights[i]).
    double pct_change = 0.0;
    for (int i = 0; i < st.n; i++) {
        double rel = std::fabs(st.weights[i] - d[i]) / std::max(d[i], 1e-12);
        if (rel > pct_change) pct_change = rel;
    }

    // WU-B: alternative metrics from final weights.
    constexpr double kMetricEps = 1e-10;
    constexpr double kChi2Floor = 1.0;
    double mean_err_sum = 0.0;
    double kl_max       = 0.0;
    double chi2_total   = 0.0;
    if (Wn > 0.0) {
        for (int k = 0; k < st.K; k++) {
            std::vector<double> S2(st.cat_counts[k], 0.0);
            for (int i = 0; i < st.n; i++) {
                int g = st.group_ids[k][i];
                if (g >= 0) S2[g] += st.weights[i];
            }
            double max_k = 0.0;
            double kl_k  = 0.0;
            for (int j = 0; j < st.cat_counts[k]; j++) {
                double S_p = S2[j] / Wn;
                double T   = st.targets[k][j];
                double err = std::fabs(S_p - T);
                if (err > max_k) max_k = err;
                if (T > 0.0) {
                    kl_k += T * std::log((T + kMetricEps) / (S_p + kMetricEps));
                }
                double obs = S2[j];
                double exp_val = T * Wn;
                chi2_total += (obs - exp_val) * (obs - exp_val) / (exp_val + kChi2Floor);
            }
            mean_err_sum += max_k;
            if (kl_k > kl_max) kl_max = kl_k;
        }
    }
    double mean_err = (st.K > 0) ? (mean_err_sum / static_cast<double>(st.K)) : 0.0;

    res.pct_change = pct_change;
    res.mean_error = mean_err;
    res.kl         = kl_max;
    res.chi2       = chi2_total;

    // WU-B: active-criterion dispatch for status.
    // For L-BFGS-B, PCT criterion is applied post-hoc: pct_change measures
    // total shift (start → final), which scales with problem difficulty and is
    // not suitable as a standalone convergence gate. Require max_err < tol_abs
    // as a quality floor for PCT-only convergence.
    {
        const auto& cfg = st.convergence_cfg;
        double active_val = 0.0;
        if (cfg.absolute_tol > 0.0) {
            switch (cfg.criterion) {
                case lbw::CalibCriterion::MAX_ERR:  active_val = max_err;    break;
                case lbw::CalibCriterion::MEAN_ERR: active_val = mean_err;   break;
                case lbw::CalibCriterion::KL:       active_val = kl_max;     break;
                case lbw::CalibCriterion::CHI2:     active_val = chi2_total; break;
                case lbw::CalibCriterion::PCT:      active_val = pct_change; break;
            }
        }
        bool converged_abs = (cfg.absolute_tol > 0.0) && (active_val < cfg.absolute_tol);
        // Spec §1: PCT-only convergence is pct_change < pct_tol, no errRp floor.
        bool converged_pct = (cfg.pct_tol > 0.0) && (pct_change < cfg.pct_tol);

        bool have_pct = (cfg.pct_tol > 0.0);
        bool have_abs = (cfg.absolute_tol > 0.0);
        bool converged = false;
        if (have_pct && have_abs) {
            converged = (cfg.stop_when == lbw::CalibStopWhen::ALL)
                        ? (converged_pct && converged_abs)
                        : (converged_pct || converged_abs);
        } else if (have_pct) {
            converged = converged_pct;
        } else if (have_abs) {
            converged = converged_abs;
        }
        res.status = converged ? RK_OK : RK_ERR_NOCONV;
    }
    return res;
}

// Zoom phase (Nocedal & Wright Alg 3.6): bisect until strong Wolfe satisfied.
// Phi and slope computed in O(n) using precomputed du[]; full O(K*n) gradient
// computed once at the accepted alpha via phi_from_u.
static double wolfe_zoom(
        const CalibState& st, const LinkFn& fn,
        const std::vector<int>& off, const std::vector<double>& T,
        const std::vector<double>& d, double phi_0, double slope_0,
        double alpha_lo, double phi_lo, double alpha_hi,
        const std::vector<double>& u_base, const std::vector<double>& du,
        const std::vector<double>& lam,
        const std::vector<double>& dir,
        std::vector<double>& u_work,
        std::vector<double>& e_vec,
        std::vector<double>& lam_new, std::vector<double>& grad_new,
        double& phi_new) {
    constexpr double kC1 = 1e-4;
    constexpr double kC2 = 0.9;
    const int total = (int)lam.size();

    // Precompute T*dir and T*lam (both constant across all bisection trials)
    double Tdir = 0.0;
    double Tlam = 0.0;
    for (int idx = 0; idx < total; idx++) { Tdir += T[idx] * dir[idx]; Tlam += T[idx] * lam[idx]; }

    double alpha_accepted = -1.0;
    for (int j = 0; j < 20; j++) {
        double alpha = 0.5 * (alpha_lo + alpha_hi);
        {
            const double a = alpha;
            const double* __restrict__ ub = u_base.data();
            const double* __restrict__ dv = du.data();
            double*       __restrict__ uw = u_work.data();
#if LBW_HAS_OMP_SIMD
#pragma omp simd
#endif
            for (int i = 0; i < st.n; i++) uw[i] = ub[i] + a * dv[i];
        }
        if (!fn.exponential)
            lbw::bulk_scaled_exp(fn.logit_scale, u_work.data(), e_vec.data(), st.n);
        double phi_trial = Tlam + alpha * Tdir;
        double slope = Tdir;
        // WU-A3: branch-hoisted SIMD. reduction(+:...) avoids the OpenMP
        // reduction(-:x) combiner bug (would yield orig+Σ, not orig-Σ).
        if (st.alm_mu > 0.0) {
            // ALM scalar fallback: dead at runtime (alm_mu=0.0 forced) but
            // kept correct for future reactivation.
            double sum_w = 0.0, sum_dw = 0.0;
            for (int i = 0; i < st.n; i++) {
                double Fi, Hi;
                if (fn.exponential) { auto fh = fn.FH(u_work[i]); Fi = fh.F; Hi = fh.H; }
                else                { Fi = fn.F_from_e(e_vec[i]); Hi = fn.H_from_e(e_vec[i], u_work[i]); }
                phi_trial -= d[i] * Hi;
                slope     -= d[i] * Fi * du[i];
                sum_w     += d[i] * Fi;
                sum_dw    += d[i] * fn.dF(u_work[i]) * du[i];
            }
            double residual = sum_w - static_cast<double>(st.n);
            double alm_scale = st.alm_lambda + st.alm_mu * residual;
            phi_trial += st.alm_lambda * residual + (st.alm_mu / 2.0) * residual * residual;
            slope     += alm_scale * sum_dw;
        } else {
            double phi_acc = 0.0, slope_acc = 0.0;
            if (fn.exponential) {
#if LBW_HAS_OMP_SIMD
#pragma omp simd reduction(+:phi_acc) reduction(+:slope_acc)
#endif
                for (int i = 0; i < st.n; i++) {
                    double e = std::exp(std::min(u_work[i], 700.0));
                    phi_acc   += d[i] * e;
                    slope_acc += d[i] * e * du[i];
                }
            } else {
                const double L = fn.L, U = fn.U, ls = fn.logit_scale;
                const double A = L * (U - 1.0), B = U * (1.0 - L);
                const double P = (U - 1.0),     Q = (1.0 - L);
                const double R = (U - L) / ls;
                const double UmL = U - L;
#if LBW_HAS_OMP_SIMD
#pragma omp simd reduction(+:phi_acc) reduction(+:slope_acc)
#endif
                for (int i = 0; i < st.n; i++) {
                    double ei = e_vec[i];
                    double denom = P + Q * ei;
                    double Fi = (A + B * ei) / denom;
                    double Hi = L * u_work[i] + R * std::log(denom / UmL);
                    phi_acc   += d[i] * Hi;
                    slope_acc += d[i] * Fi * du[i];
                }
            }
            phi_trial -= phi_acc;
            slope     -= slope_acc;
        }

        // Armijo test (maximization: slope_0 > 0). Fail ⟹ phi_trial did NOT
        // exceed phi_0 by at least c1*α*slope_0 ⟹ step too long ⟹ narrow bracket.
        if (phi_trial < phi_0 + kC1 * alpha * slope_0 || phi_trial <= phi_lo) {
            alpha_hi = alpha;
        } else {
            if (std::fabs(slope) <= kC2 * std::fabs(slope_0)) {
                alpha_accepted = alpha;
                phi_new = phi_trial;
                break;
            }
            if (slope * (alpha_hi - alpha_lo) >= 0.0) alpha_hi = alpha_lo;
            alpha_lo = alpha; phi_lo = phi_trial;
        }
    }
    if (alpha_accepted < 0.0) alpha_accepted = 0.5 * (alpha_lo + alpha_hi);
    // Compute full gradient at accepted point (O(K*n) once per outer step)
    for (int idx = 0; idx < total; idx++) lam_new[idx] = lam[idx] + alpha_accepted * dir[idx];
    {
        const double a = alpha_accepted;
        const double* __restrict__ ub = u_base.data();
        const double* __restrict__ dv = du.data();
        double*       __restrict__ uw = u_work.data();
#if LBW_HAS_OMP_SIMD
#pragma omp simd
#endif
        for (int i = 0; i < st.n; i++) uw[i] = ub[i] + a * dv[i];
    }
    phi_new = phi_from_u(st, fn, off, lam_new, T, d, grad_new, u_work);
    return alpha_accepted;
}

// Strong Wolfe line search (Nocedal & Wright Alg 3.5/3.6).
// c1=1e-4 (Armijo), c2=0.9 (curvature). Guarantees sy>0 for L-BFGS update.
// Uses O(n) phi+slope evals during bracketing; O(K*n) only at accepted step.
static double wolfe_line_search(
        const CalibState& st, const LinkFn& fn,
        const std::vector<int>& off, const std::vector<double>& T,
        const std::vector<double>& d,
        const std::vector<double>& lam, double phi_0, double slope_0,
        const std::vector<double>& u_base, const std::vector<double>& du,
        std::vector<double>& u_work,
        std::vector<double>& e_vec,
        const std::vector<double>& dir,
        std::vector<double>& lam_new, std::vector<double>& grad_new,
        double& phi_new) {
    constexpr double kC1 = 1e-4;
    constexpr double kC2 = 0.9;
    const int total = (int)lam.size();

    // Precompute T*dir and T*lam (both constant per Wolfe search)
    double Tdir = 0.0;
    double Tlam = 0.0;
    for (int idx = 0; idx < total; idx++) {
        Tdir += T[idx] * dir[idx];
        Tlam += T[idx] * lam[idx];
    }

    double alpha_prev = 0.0, phi_prev = phi_0;
    double alpha = 1.0;

    for (int i = 0; i < 20; i++) {
        {
            const double a = alpha;
            const double* __restrict__ ub = u_base.data();
            const double* __restrict__ dv = du.data();
            double*       __restrict__ uw = u_work.data();
#if LBW_HAS_OMP_SIMD
#pragma omp simd
#endif
            for (int j = 0; j < st.n; j++) uw[j] = ub[j] + a * dv[j];
        }
        if (!fn.exponential)
            lbw::bulk_scaled_exp(fn.logit_scale, u_work.data(), e_vec.data(), st.n);
        double phi_trial = Tlam + Tdir * alpha;
        double slope = Tdir;
        // WU-A3: branch-hoisted SIMD. reduction(+:...) avoids combiner bug.
        if (st.alm_mu > 0.0) {
            // ALM scalar fallback (dead at runtime; kept correct).
            double sum_w = 0.0, sum_dw = 0.0;
            for (int j = 0; j < st.n; j++) {
                double Fj, Hj;
                if (fn.exponential) { auto fh = fn.FH(u_work[j]); Fj = fh.F; Hj = fh.H; }
                else                { Fj = fn.F_from_e(e_vec[j]); Hj = fn.H_from_e(e_vec[j], u_work[j]); }
                phi_trial -= d[j] * Hj;
                slope     -= d[j] * Fj * du[j];
                sum_w     += d[j] * Fj;
                sum_dw    += d[j] * fn.dF(u_work[j]) * du[j];
            }
            double residual = sum_w - static_cast<double>(st.n);
            double alm_scale = st.alm_lambda + st.alm_mu * residual;
            phi_trial += st.alm_lambda * residual + (st.alm_mu / 2.0) * residual * residual;
            slope     += alm_scale * sum_dw;
        } else {
            double phi_acc = 0.0, slope_acc = 0.0;
            if (fn.exponential) {
#if LBW_HAS_OMP_SIMD
#pragma omp simd reduction(+:phi_acc) reduction(+:slope_acc)
#endif
                for (int j = 0; j < st.n; j++) {
                    double e = std::exp(std::min(u_work[j], 700.0));
                    phi_acc   += d[j] * e;
                    slope_acc += d[j] * e * du[j];
                }
            } else {
                const double L = fn.L, U = fn.U, ls = fn.logit_scale;
                const double A = L * (U - 1.0), B = U * (1.0 - L);
                const double P = (U - 1.0),     Q = (1.0 - L);
                const double R = (U - L) / ls;
                const double UmL = U - L;
#if LBW_HAS_OMP_SIMD
#pragma omp simd reduction(+:phi_acc) reduction(+:slope_acc)
#endif
                for (int j = 0; j < st.n; j++) {
                    double ej = e_vec[j];
                    double denom = P + Q * ej;
                    double Fj = (A + B * ej) / denom;
                    double Hj = L * u_work[j] + R * std::log(denom / UmL);
                    phi_acc   += d[j] * Hj;
                    slope_acc += d[j] * Fj * du[j];
                }
            }
            phi_trial -= phi_acc;
            slope     -= slope_acc;
        }

        // Armijo test (maximization: slope_0 > 0). Fail ⟹ phi_trial did NOT
        // exceed phi_0 by at least c1*α*slope_0 ⟹ step too long ⟹ enter zoom.
        if (phi_trial < phi_0 + kC1 * alpha * slope_0 || (i > 0 && phi_trial <= phi_prev)) {
            return wolfe_zoom(st, fn, off, T, d, phi_0, slope_0,
                              alpha_prev, phi_prev, alpha,
                              u_base, du, lam, dir, u_work, e_vec,
                              lam_new, grad_new, phi_new);
        }
        if (std::fabs(slope) <= kC2 * std::fabs(slope_0)) {
            // Accepted: compute full gradient
            for (int j = 0; j < total; j++) lam_new[j] = lam[j] + alpha * dir[j];
            phi_new = phi_from_u(st, fn, off, lam_new, T, d, grad_new, u_work);
            return alpha;
        }
        if (slope <= 0) {
            return wolfe_zoom(st, fn, off, T, d, phi_0, slope_0,
                              alpha, phi_trial, alpha_prev,
                              u_base, du, lam, dir, u_work, e_vec,
                              lam_new, grad_new, phi_new);
        }
        alpha_prev = alpha; phi_prev = phi_trial;
        alpha = std::min(2.0 * alpha, 8.0);
    }
    st.log("L-BFGS-B: Wolfe bracket did not converge; using last step");
    for (int j = 0; j < total; j++) lam_new[j] = lam[j] + alpha * dir[j];
    {
        const double a = alpha;
        const double* __restrict__ ub = u_base.data();
        const double* __restrict__ dv = du.data();
        double*       __restrict__ uw = u_work.data();
#if LBW_HAS_OMP_SIMD
#pragma omp simd
#endif
        for (int j = 0; j < st.n; j++) uw[j] = ub[j] + a * dv[j];
    }
    phi_new = phi_from_u(st, fn, off, lam_new, T, d, grad_new, u_work);
    return alpha;
}

// Extracted from lbfgsb_solve to isolate the L-BFGS-B iteration kernel.
// If an ALM outer loop is added, lbfgsb_solve calls this per outer iteration.
static LBFGSResult lbfgsb_solve_inner(CalibState& st,
                                      const std::vector<int>& off,
                                      const std::vector<double>& T,
                                      const std::vector<double>& d,
                                      double W_sum) {
    LinkFn fn(st.min_weight, st.max_weight);
    int total = off[st.K];

    std::vector<double> lam(total, 0.0);
    std::vector<double> grad(total), grad_new(total);
    std::vector<double> u(st.n);       // u at current lam
    std::vector<double> du(st.n);      // per-obs directional derivative for Wolfe
    std::vector<double> u_work(st.n);  // scratch for Wolfe trial u
    std::vector<double> e_vec(st.n);   // scratch: exp(logit_scale * u_work[i]) per trial
    std::vector<double> dir(total);

    std::deque<std::vector<double>> svec, yvec;
    std::deque<double> rho_hist;
    double gamma = 1.0;

    double phi_curr = phi_and_grad(st, fn, off, lam, T, d, grad, u);

    int max_iter = st.outer_max_iter;
    int final_iter = 0;
    std::vector<double> lam_new(total), s_new(total), y_new(total);
    for (int iter = 0; iter < max_iter; iter++) {
        final_iter = iter + 1;

        double gn = maxAbs(grad) / W_sum;
        if (gn < st.tol_abs) break;

        if (svec.empty()) {
            dir = grad;
        } else {
            lbfgs_direction(svec, yvec, rho_hist, gamma, grad, dir);
        }

        double slope_0 = dot(grad, dir);
        // Precompute per-observation directional derivative O(K*n) once per step.
        // Wolfe evals then cost O(n) each instead of O(K*n).
        compute_du(st, off, dir, du);

        double phi_new = phi_curr;

        wolfe_line_search(st, fn, off, T, d, lam, phi_curr, slope_0,
                          u, du, u_work, e_vec, dir, lam_new, grad_new, phi_new);

        for (int i = 0; i < total; i++) {
            s_new[i] = lam_new[i] - lam[i];
            y_new[i] = grad_new[i] - grad[i];
        }
        double sy = dot(s_new, y_new);
        double yy = dot(y_new, y_new);
        // Relative curvature gate (Liu & Nocedal 1989 §3; Nocedal-Wright 2e §7.2):
        // sy > ε · ‖s‖ · ‖y‖. Squared form avoids sqrt: sy² > ε²·‖s‖²·‖y‖².
        // On rejection, fall through to steepest-ascent on this iteration
        // (L-BFGS history unchanged). Concave objective ⟹ still converges.
        constexpr double kCurvRel = 1e-8;
        double s_norm2 = dot(s_new, s_new);
        bool curv_ok = (sy > 0.0) && (sy * sy > kCurvRel * kCurvRel * s_norm2 * yy);
        if (curv_ok) {
            if ((int)svec.size() >= st.lbfgs_m) {
                svec.pop_front(); yvec.pop_front(); rho_hist.pop_front();
            }
            svec.push_back(s_new);
            yvec.push_back(y_new);
            rho_hist.push_back(1.0 / sy);
            gamma = sy / yy;
        }

        // O(1) pointer swap. lam/grad take the new values; lam_new/grad_new
        // retain the old buffers and are fully overwritten on the next Wolfe
        // call (lam_new at 257, 327, 341; grad_new via phi_from_u at 259,
        // 328, 343) before any further read.
        std::swap(lam, lam_new);
        std::swap(grad, grad_new);
        // Recompute u from scratch rather than using u_work from the Wolfe search.
        // Incremental du accumulation (u = u_base + alpha*du) can drift over many
        // outer steps; recomputing ensures u stays exactly consistent with lam.
        // Cost: one O(K*n) pass per outer iteration — same order as gradient computation.
        compute_u(st, off, lam, u);
        phi_curr = phi_new;
    }

    return compute_final_weights_and_error(st, fn, d, off, u, final_iter);
}

LBFGSResult lbfgsb_solve(CalibState& st) {
    auto off = build_offsets(st);

    std::vector<double> T(off[st.K]);
    double W_sum = compute_targets_abs(st, T);

    std::vector<double> d(st.n);
    for (int i = 0; i < st.n; i++) d[i] = st.weights[i];

    // ALM inactive. Solver-level invariant: at dual convergence ∇phi=0 ⟹
    // for each margin k, Σ_j S_kj = Σ_j T_kj = W (given targets sum to 1),
    // and summing observations once gives sum(w) = W = Σ d_i. Input contract
    // (R/harvest.R:218, python/_harvest.py): start_weights are pre-normalized
    // so Σ d_i = n. Output contract (moved in-solver 2026-04-24 per user
    // directive): solver self-normalizes to sum(w) = n before return.
    st.alm_lambda = 0.0;
    st.alm_mu     = 0.0;
    LBFGSResult res = lbfgsb_solve_inner(st, off, T, d, W_sum);

    // Post-solve normalization. Placement after the inner solve is safe because
    // res.max_error is computed inside compute_final_weights_and_error as
    // max_{k,j} |S_kj/Wn - target_kj| — scale-invariant under uniform weight
    // rescaling. Therefore normalizing st.weights here does NOT invalidate
    // res.max_error. Contract: total_w == 0 leaves weights untouched.
    double total_w = 0.0;
    for (int i = 0; i < st.n; i++) total_w += st.weights[i];
    if (std::isfinite(total_w) && total_w > 0.0) {
        const double norm = static_cast<double>(st.n) / total_w;
        for (int i = 0; i < st.n; i++) st.weights[i] *= norm;
    }
    return res;
}

} // namespace lbw