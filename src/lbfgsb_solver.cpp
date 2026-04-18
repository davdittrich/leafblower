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

// phi(lambda) = sum_kj T_kj*lam_kj - sum_i d_i*H(u_i)
// grad[off_k+j] = T_kj - S_kj where S_kj = sum_{i:g_k(i)==j} d_i*F(u_i)
static double phi_and_grad(const CalibState& st,
                            const LinkFn& fn,
                            const std::vector<int>& off,
                            const std::vector<double>& lam,
                            const std::vector<double>& T,
                            const std::vector<double>& d,
                            std::vector<double>& grad,
                            std::vector<double>& u) {
    int total = off[st.K];
    compute_u(st, off, lam, u);

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
static double wolfe_zoom(
        const CalibState& st, const LinkFn& fn,
        const std::vector<int>& off, const std::vector<double>& T,
        const std::vector<double>& d, double phi_0, double slope_0,
        double alpha_lo, double phi_lo, double alpha_hi,
        std::vector<double>& u, const std::vector<double>& lam,
        const std::vector<double>& dir,
        std::vector<double>& lam_new, std::vector<double>& grad_new,
        double& phi_new) {
    constexpr double kC1 = 1e-4;
    constexpr double kC2 = 0.9;
    const int total = (int)lam.size();

    for (int j = 0; j < 20; j++) {
        double alpha = 0.5 * (alpha_lo + alpha_hi);
        for (int i = 0; i < total; i++) lam_new[i] = lam[i] + alpha * dir[i];
        phi_new = phi_and_grad(st, fn, off, lam_new, T, d, grad_new, u);
        double slope = dot(grad_new, dir);

        if (phi_new < phi_0 + kC1 * alpha * slope_0 || phi_new <= phi_lo) {
            alpha_hi = alpha;
        } else {
            if (std::fabs(slope) <= kC2 * std::fabs(slope_0)) return alpha;
            if (slope * (alpha_hi - alpha_lo) >= 0.0) alpha_hi = alpha_lo;
            alpha_lo = alpha; phi_lo = phi_new;
        }
    }
    return 0.5 * (alpha_lo + alpha_hi);
}

// Strong Wolfe line search (Nocedal & Wright Alg 3.5/3.6).
// c1=1e-4 (Armijo), c2=0.9 (curvature). Guarantees sy>0 for L-BFGS update.
static double wolfe_line_search(
        const CalibState& st, const LinkFn& fn,
        const std::vector<int>& off, const std::vector<double>& T,
        const std::vector<double>& d,
        const std::vector<double>& lam, double phi_0, double slope_0,
        std::vector<double>& u, const std::vector<double>& dir,
        std::vector<double>& lam_new, std::vector<double>& grad_new,
        double& phi_new) {
    constexpr double kC1 = 1e-4;
    constexpr double kC2 = 0.9;
    const int total = (int)lam.size();

    double alpha_prev = 0.0, phi_prev = phi_0;
    double alpha = 1.0;

    for (int i = 0; i < 20; i++) {
        for (int j = 0; j < total; j++) lam_new[j] = lam[j] + alpha * dir[j];
        phi_new = phi_and_grad(st, fn, off, lam_new, T, d, grad_new, u);
        double slope = dot(grad_new, dir);

        if (phi_new < phi_0 + kC1 * alpha * slope_0 || (i > 0 && phi_new <= phi_prev)) {
            return wolfe_zoom(st, fn, off, T, d, phi_0, slope_0,
                              alpha_prev, phi_prev, alpha, u, lam, dir,
                              lam_new, grad_new, phi_new);
        }
        if (std::fabs(slope) <= kC2 * std::fabs(slope_0)) return alpha;
        if (slope <= 0) {
            return wolfe_zoom(st, fn, off, T, d, phi_0, slope_0,
                              alpha, phi_new, alpha_prev, u, lam, dir,
                              lam_new, grad_new, phi_new);
        }
        alpha_prev = alpha; phi_prev = phi_new;
        alpha = std::min(2.0 * alpha, 8.0);
    }
    st.log("L-BFGS-B: Wolfe bracket did not converge; using last step");
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
    std::vector<double> u(st.n);
    std::vector<double> dir(total);

    std::deque<std::vector<double>> svec, yvec;
    std::deque<double> rho_hist;
    double gamma = 1.0;

    compute_targets_abs(st, T);
    double W_sum = 0.0;
    for (int i = 0; i < st.n; i++) W_sum += d[i];
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
        std::vector<double> lam_new(total);
        double phi_new = phi_curr;

        wolfe_line_search(st, fn, off, T, d, lam, phi_curr, slope_0,
                          u, dir, lam_new, grad_new, phi_new);

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
        phi_curr = phi_new;
    }

    compute_u(st, off, lam, u);
    return compute_final_weights_and_error(st, fn, d, off, u, final_iter);
}

} // namespace lbw
