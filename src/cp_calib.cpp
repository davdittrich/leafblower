// Single-threaded solver; cell-table build inherits ieppa OpenMP behavior unchanged.
//
// Epic-K K-1: port of research/cp_calib (Algorithm 1 obs-level, Chambolle-Pock
// primal-dual on the box-constrained KL calibration problem). The CSR design
// matrix A is constructed on entry from CalibState (one nonzero per (obs i,
// margin k) at column cat_offset[k] + group_ids[k][i], value 1.0). Algorithm 2
// (accelerated PDHG) and cell compression land in K-3 / K-4 respectively;
// accelerate=TRUE here transiently falls back to Algorithm 1 with
// fell_back_to_pdhg=true.

#include "cp_calib.hpp"
#include "leafblower.h"
#include "calib_dispatch.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace lbw {

CpCalibResult cp_calibrate(CalibState& st) {
    CpCalibResult res;
    res.base.status         = RK_ERR_NOCONV;
    res.algorithm_requested = st.accelerate ? "accelerated_pdhg" : "pdhg";
    // K-1: Algorithm 2 not yet implemented — transient fallback to Algorithm 1.
    // Full Alg 2 (gamma-strong-convexity dispatch + adaptive step sizes) lands in K-3.
    res.algorithm_used      = "pdhg";
    res.fell_back_to_pdhg   = st.accelerate;

    auto t0 = std::chrono::steady_clock::now();
    auto stamp_walltime = [&]() {
        res.wall_time_ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - t0).count();
    };

    const int n = st.n;
    const int K = st.K;
    res.n_cells = n;  // K-1 = obs-level only; cell-compression lands in K-4

    // Build cat_offset and CSR matrix A (n_row = n, n_col = ΣJ_k, K nnz/row).
    std::vector<int> cat_offset;
    const int n_col = build_cat_offset(K, st.cat_counts, cat_offset);

    // Per-row design weight d_i (= input st.weights[i]) and bounds lo[i], hi[i].
    const double hi_scalar = resolve_hi(st);
    const double lo_scalar = st.min_weight;
    std::vector<double> d_vec(n), lo_vec(n), hi_vec(n);
    for (int i = 0; i < n; ++i) {
        d_vec[i]  = st.weights[i];
        lo_vec[i] = lo_scalar;
        hi_vec[i] = hi_scalar;
    }
    const double* d  = d_vec.data();
    const double* lo = lo_vec.data();
    const double* hi = hi_vec.data();

    // CSR: p[n+1], j[n*K], x[n*K] = 1.0. Each row has up to K nonzeros — those
    // margins where group_ids[k][i] >= 0. Rows with NA gids contribute fewer.
    std::vector<int>    p_csr(n + 1, 0);
    std::vector<int>    j_csr;
    std::vector<double> x_csr;
    j_csr.reserve(static_cast<size_t>(n) * K);
    x_csr.reserve(static_cast<size_t>(n) * K);
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < K; ++k) {
            const int g = st.group_ids[k][i];
            if (g < 0 || g >= st.cat_counts[k]) continue;
            j_csr.push_back(cat_offset[k] + g);
            x_csr.push_back(1.0);
        }
        p_csr[i + 1] = static_cast<int>(j_csr.size());
    }
    const int* p = p_csr.data();
    const int* j = j_csr.data();
    const double* x = x_csr.data();

    // Targets b[c] = targets[k][j] * Z, where Z = sum(d). Per spec Sec 3.1.
    double Z = 0.0;
    for (int i = 0; i < n; ++i) Z += d[i];
    std::vector<double> b_vec(n_col, 0.0);
    if (Z <= 0.0) {
        res.base.status      = 2;
        res.base.iterations  = 0;
        res.base.max_error   = 1.0;
        res.base.best_weights.assign(n, 1.0);
        for (int i = 0; i < n; ++i) st.weights[i] = 1.0;
        stamp_walltime();
        return res;
    }
    for (int k = 0; k < K; ++k)
        for (int jj = 0; jj < st.cat_counts[k]; ++jj)
            b_vec[cat_offset[k] + jj] = st.targets[k][jj] * Z;
    const double* b = b_vec.data();

    // Infeasibility check: lo[i] >= hi[i] - 2*delta. Returns init-clamp weights
    // per spec Sec 2 status_code=3 contract (NOT all-1.0).
    const double delta = 1e-8;
    bool infeasible = false;
    for (int i = 0; i < n; ++i) {
        if (lo[i] >= hi[i] - 2.0 * delta) { infeasible = true; break; }
    }
    if (infeasible) {
        std::vector<double> w_init(n);
        for (int i = 0; i < n; ++i)
            w_init[i] = std::clamp(d[i], lo[i] + delta, hi[i] - delta);
        res.base.status       = 3;          // infeasible bounds
        res.base.iterations   = 0;
        res.base.max_error    = 1.0;
        res.base.best_weights = w_init;
        for (int i = 0; i < n; ++i) st.weights[i] = w_init[i];
        stamp_walltime();
        return res;
    }

    const int max_iter = (st.outer_max_iter > 0) ? st.outer_max_iter : 500;

    // ── Power iteration for ||A|| (50-cap; rel-delta < 1e-6 for converged). ──
    std::vector<double> v(n_col, 1.0 / std::sqrt(static_cast<double>(n_col)));
    double lambda_old = 0.0;
    double rel_delta  = 1.0;
    bool   power_converged = false;
    int    n_power_iter    = 0;
    for (int iter = 0; iter < 50; ++iter) {
        n_power_iter = iter + 1;
        std::vector<double> u(n, 0.0);
        for (int i = 0; i < n; ++i) {
            double s = 0.0;
            for (int idx = p[i]; idx < p[i + 1]; ++idx) s += x[idx] * v[j[idx]];
            u[i] = s;
        }
        std::vector<double> v_new(n_col, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int idx = p[i]; idx < p[i + 1]; ++idx)
                v_new[j[idx]] += x[idx] * u[i];
        }
        double lambda_new = 0.0;
        for (int c = 0; c < n_col; ++c) lambda_new += v[c] * v_new[c];
        if (lambda_old > 0.0) {
            rel_delta = std::abs(lambda_new - lambda_old) / lambda_old;
            if (rel_delta < 1e-6) { power_converged = true; lambda_old = lambda_new; break; }
        }
        lambda_old = lambda_new;
        double norm_v_new = 0.0;
        for (int c = 0; c < n_col; ++c) norm_v_new += v_new[c] * v_new[c];
        norm_v_new = std::sqrt(norm_v_new);
        if (norm_v_new <= 0.0) break;
        for (int c = 0; c < n_col; ++c) v[c] = v_new[c] / norm_v_new;
    }
    res.n_power_iter = n_power_iter;

    if (!power_converged && rel_delta > 1e-3) {
        res.base.status      = 4;          // power-iter divergence
        res.base.iterations  = 0;
        res.base.max_error   = 1.0;
        std::vector<double> w_init(n);
        for (int i = 0; i < n; ++i)
            w_init[i] = std::clamp(d[i], lo[i] + delta, hi[i] - delta);
        res.base.best_weights = w_init;
        for (int i = 0; i < n; ++i) st.weights[i] = w_init[i];
        stamp_walltime();
        return res;
    }

    const double norm_A      = std::sqrt(std::max(lambda_old, 0.0));
    const double norm_A_used = norm_A * 1.05;
    res.A_norm_estimate = norm_A;
    if (norm_A_used <= 0.0) {
        res.base.status      = 2;
        res.base.iterations  = 0;
        res.base.max_error   = 1.0;
        std::vector<double> w_init(n);
        for (int i = 0; i < n; ++i)
            w_init[i] = std::clamp(d[i], lo[i] + delta, hi[i] - delta);
        res.base.best_weights = w_init;
        for (int i = 0; i < n; ++i) st.weights[i] = w_init[i];
        stamp_walltime();
        return res;
    }
    double sigma = 1.0 / norm_A_used;
    double tau   = 1.0 / norm_A_used;

    // ── Init w, w_bar, y. ──
    std::vector<double> w(n, 0.0);
    for (int i = 0; i < n; ++i)
        w[i] = std::clamp(d[i], lo[i] + delta, hi[i] - delta);
    std::vector<double> w_bar = w;
    std::vector<double> y(n_col, 0.0);

    int    final_iter = max_iter;
    int    status     = 1;          // max_iter exhausted

    for (int k = 0; k < max_iter; ++k) {
        // y ← y + sigma (A^T w_bar − b)
        std::vector<double> AT_w_bar(n_col, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int idx = p[i]; idx < p[i + 1]; ++idx)
                AT_w_bar[j[idx]] += x[idx] * w_bar[i];
        }
        for (int c = 0; c < n_col; ++c) y[c] += sigma * (AT_w_bar[c] - b[c]);

        // Ay
        std::vector<double> Ay(n, 0.0);
        for (int i = 0; i < n; ++i) {
            double s = 0.0;
            for (int idx = p[i]; idx < p[i + 1]; ++idx) s += x[idx] * y[j[idx]];
            Ay[i] = s;
        }

        bool has_nan_inf = false;
        for (int i = 0; i < n; ++i) {
            const double z_i      = w[i] - tau * Ay[i];
            const double w_old_i  = w[i];
            const double tau_d_i  = tau * d[i];
            double w_n = std::clamp(d[i], lo[i] + delta, hi[i] - delta);

            if (tau_d_i <= 0.0) {
                w_n = z_i;
            } else if (std::abs(z_i / tau_d_i - 1.0) > 700.0) {
                w_n = z_i;
            } else {
                for (int p_iter = 0; p_iter < 20; ++p_iter) {
                    const double f_prime = tau_d_i * std::log(w_n / d[i]) + (w_n - z_i);
                    if (std::abs(f_prime) < 1e-12) break;
                    const double f_double_prime = tau_d_i / w_n + 1.0;
                    w_n = w_n - f_prime / f_double_prime;
                    if (w_n <= 0.0) w_n = 1e-12;
                }
            }
            w_n = std::clamp(w_n, lo[i] + delta, hi[i] - delta);
            w[i]     = w_n;
            w_bar[i] = 2.0 * w_n - w_old_i;
        }

        for (int c = 0; c < n_col; ++c) if (!std::isfinite(y[c])) { has_nan_inf = true; break; }
        if (!has_nan_inf)
            for (int i = 0; i < n; ++i) if (!std::isfinite(w[i])) { has_nan_inf = true; break; }

        if (has_nan_inf) {
            status     = 2;
            final_iter = k + 1;
            break;
        }

        // Convergence: ||A^T w - b||_inf / max(Z, 1) < 1e-7  → status=0.
        std::vector<double> AT_w(n_col, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int idx = p[i]; idx < p[i + 1]; ++idx)
                AT_w[j[idx]] += x[idx] * w[i];
        }
        double primal_resid = 0.0;
        for (int c = 0; c < n_col; ++c) {
            const double r = std::abs(AT_w[c] - b[c]);
            if (r > primal_resid) primal_resid = r;
        }
        final_iter = k + 1;
        if (primal_resid / std::max(Z, 1.0) < 1e-7) {
            status = 0;
            break;
        }
    }

    // ── Finalize: write calibrated weights and compute marginal max_error. ──
    for (int i = 0; i < n; ++i) st.weights[i] = w[i];

    double total = 0.0;
    for (int i = 0; i < n; ++i) total += w[i];
    double max_err = 0.0;
    if (total > 0.0) {
        const double inv_tot = 1.0 / total;
        for (int kk = 0; kk < K; ++kk) {
            std::vector<double> achieved(st.cat_counts[kk], 0.0);
            for (int i = 0; i < n; ++i) {
                const int g = st.group_ids[kk][i];
                if (g >= 0 && g < st.cat_counts[kk]) achieved[g] += w[i];
            }
            for (int jj = 0; jj < st.cat_counts[kk]; ++jj) {
                const double e = std::abs(achieved[jj] * inv_tot - st.targets[kk][jj]);
                if (e > max_err) max_err = e;
            }
        }
    } else {
        max_err = 1.0;
    }

    res.base.status       = status;
    res.base.iterations   = final_iter;
    res.base.max_error    = max_err;
    res.base.best_weights = w;
    res.base.best_error   = max_err;
    res.base.best_iter    = final_iter;
    if (status == 0) {
        res.base.convergence_metric = static_cast<int>(CalibMetric::MAX_ERR);
        res.base.convergence_rule   = static_cast<int>(CalibRule::THRESHOLD);
        res.base.convergence_tol    = 1e-7;
        res.base.convergence_iter   = final_iter;
    }

    res.final_theta = std::numeric_limits<double>::quiet_NaN();
    res.final_tau   = tau;
    res.final_sigma = sigma;

    stamp_walltime();
    return res;
}

} // namespace lbw
