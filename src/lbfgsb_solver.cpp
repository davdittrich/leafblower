#include "lbfgsb_solver.hpp"
#include "leafblower.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <deque>

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

static double linf(const std::vector<double>& v) {
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

    std::vector<double> alpha(m);
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
        wi = std::max(st.min_weight, std::min(st.max_weight, wi));
        st.weights[i] = wi;
    }
    // Do NOT normalize: bridge normalizes start_weights to mean=1 before
    // rk_calibrate(); re-normalizing after clamping invalidates calibration constraints.

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
    res.max_error = max_err;
    res.status = (max_err < st.tol_abs) ? RK_OK : RK_ERR_NOCONV;
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
        for (int i = 0; i < st.n; i++) u_work[i] = u_base[i] + alpha * du[i];
        double phi_trial = Tlam + alpha * Tdir;
        double slope = Tdir;
        for (int i = 0; i < st.n; i++) {
            double Fi = fn.F(u_work[i]);
            phi_trial -= d[i] * fn.H(u_work[i]);
            slope -= d[i] * Fi * du[i];
        }

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
    for (int i = 0; i < st.n; i++) u_work[i] = u_base[i] + alpha_accepted * du[i];
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
        const std::vector<double>& dir,
        std::vector<double>& lam_new, std::vector<double>& grad_new,
        double& phi_new) {
    constexpr double kC1 = 1e-4;
    constexpr double kC2 = 0.9;
    const int total = (int)lam.size();

    // Precompute T*dir (constant per Wolfe search)
    double Tdir = 0.0;
    for (int idx = 0; idx < total; idx++) Tdir += T[idx] * dir[idx];

    double alpha_prev = 0.0, phi_prev = phi_0;
    double alpha = 1.0;

    for (int i = 0; i < 20; i++) {
        for (int j = 0; j < st.n; j++) u_work[j] = u_base[j] + alpha * du[j];
        double phi_trial = Tdir * alpha;
        for (int idx = 0; idx < total; idx++) phi_trial += T[idx] * lam[idx];
        double slope = Tdir;
        for (int j = 0; j < st.n; j++) {
            double Fj = fn.F(u_work[j]);
            phi_trial -= d[j] * fn.H(u_work[j]);
            slope -= d[j] * Fj * du[j];
        }

        if (phi_trial < phi_0 + kC1 * alpha * slope_0 || (i > 0 && phi_trial <= phi_prev)) {
            return wolfe_zoom(st, fn, off, T, d, phi_0, slope_0,
                              alpha_prev, phi_prev, alpha,
                              u_base, du, lam, dir, u_work,
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
                              u_base, du, lam, dir, u_work,
                              lam_new, grad_new, phi_new);
        }
        alpha_prev = alpha; phi_prev = phi_trial;
        alpha = std::min(2.0 * alpha, 8.0);
    }
    st.log("L-BFGS-B: Wolfe bracket did not converge; using last step");
    for (int j = 0; j < total; j++) lam_new[j] = lam[j] + alpha * dir[j];
    for (int j = 0; j < st.n; j++) u_work[j] = u_base[j] + alpha * du[j];
    phi_new = phi_from_u(st, fn, off, lam_new, T, d, grad_new, u_work);
    return alpha;
}

LBFGSResult lbfgsb_solve(CalibState& st) {
    LinkFn fn(st.min_weight, st.max_weight);
    auto off = build_offsets(st);
    int total = off[st.K];

    std::vector<double> d(st.n);
    for (int i = 0; i < st.n; i++) d[i] = st.weights[i];

    std::vector<double> lam(total, 0.0);
    std::vector<double> T(total);
    std::vector<double> grad(total), grad_new(total);
    std::vector<double> u(st.n);       // u at current lam
    std::vector<double> du(st.n);      // per-obs directional derivative for Wolfe
    std::vector<double> u_work(st.n);  // scratch for Wolfe trial u
    std::vector<double> dir(total);

    std::deque<std::vector<double>> svec, yvec;
    std::deque<double> rho_hist;
    double gamma = 1.0;

    double W_sum = compute_targets_abs(st, T);
    double phi_curr = phi_and_grad(st, fn, off, lam, T, d, grad, u);

    int max_iter = st.outer_max_iter;
    int final_iter = 0;
    for (int iter = 0; iter < max_iter; iter++) {
        final_iter = iter + 1;

        double gn = linf(grad) / W_sum;
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

        std::vector<double> lam_new(total);
        double phi_new = phi_curr;

        wolfe_line_search(st, fn, off, T, d, lam, phi_curr, slope_0,
                          u, du, u_work, dir, lam_new, grad_new, phi_new);

        std::vector<double> s_new(total), y_new(total);
        for (int i = 0; i < total; i++) {
            s_new[i] = lam_new[i] - lam[i];
            y_new[i] = grad_new[i] - grad[i];
        }
        double sy = dot(s_new, y_new);
        double yy = dot(y_new, y_new);
        constexpr double kCurvMin = 1e-20;
        if (sy > kCurvMin && yy > kCurvMin) {
            if ((int)svec.size() >= st.lbfgs_m) {
                svec.pop_front(); yvec.pop_front(); rho_hist.pop_front();
            }
            svec.push_back(s_new);
            yvec.push_back(y_new);
            rho_hist.push_back(1.0 / sy);
            gamma = sy / yy;
        }

        lam = lam_new;
        grad = grad_new;
        // Recompute u from scratch rather than using u_work from the Wolfe search.
        // Incremental du accumulation (u = u_base + alpha*du) can drift over many
        // outer steps; recomputing ensures u stays exactly consistent with lam.
        // Cost: one O(K*n) pass per outer iteration — same order as gradient computation.
        compute_u(st, off, lam, u);
        phi_curr = phi_new;
    }

    return compute_final_weights_and_error(st, fn, d, off, u, final_iter);
}

} // namespace lbw
