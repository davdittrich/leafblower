// Epic-J WU-5: Interior-Point Newton (IPM) implementation per spec Sec 3.2.
//
// Faithful textbook vanilla central path (NO Mehrotra predictor-corrector,
// NO momentum, NO warm-start). All state in double. Diagonal Hessian + dense
// Schur S = A^T H_w^{-1} A solved via TSVD (LAPACK dsyevd, ratio 1e-8) — mirrors
// Epic-Dβ WL-1 in src/newton_calib.cpp.

#include "ipm_calib.hpp"

#include <R_ext/Lapack.h>
#include <R_ext/RS.h>   // F77_CALL
#ifndef FCONE
# define FCONE
#endif

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <chrono>
#include <limits>
#include <vector>

IPMResult ipm_calibrate(
    const int n_row,
    const int n_col,
    const int* p,
    const int* j,
    const double* x,
    const double* b,
    const double* d,
    const double* lo,
    const double* hi,
    const int max_iterations,
    const bool capture_trace
) {
    IPMResult res;
    res.weights.assign(n_row, 1.0);
    res.status_code = 1;  // outer/inner cap exceeded by default
    res.status_msg = "max iterations reached";
    res.iterations = 0;
    res.wall_time_ms = 0.0;

    auto start_time = std::chrono::steady_clock::now();

    // Init w^0 = clamp(d_i, lo + delta, hi - delta). Strict interior.
    const double delta = 1e-8;
    std::vector<double> w(n_row, 0.0);
    for (int i = 0; i < n_row; i++) {
        w[i] = std::clamp(d[i], lo[i] + delta, hi[i] - delta);
    }
    std::vector<double> lam(n_col, 0.0);  // dual multipliers

    // Outer barrier loop schedule.
    const int    outer_cap   = 50;
    const int    inner_cap   = 5;
    const double mu0         = 1.0;
    const double mu_floor    = 1e-9;
    const double tsvd_ratio  = 1e-8;
    const double final_tol   = 1e-7;
    // Trace stride: cap at 1000 rows. Outer cap 50, so stride 1 always fits.
    const int trace_stride = std::max(1, (outer_cap + 999) / 1000);

    int total_inner = 0;
    bool converged = false;
    double mu = mu0;

    // Workspace for dense Schur S (m × m) using Eigen column-major.
    Eigen::MatrixXd S(n_col, n_col);
    std::vector<double> Hinv(n_row, 0.0);    // 1 / H_w diagonal
    std::vector<double> r_w(n_row, 0.0);
    std::vector<double> r_lam(n_col, 0.0);
    std::vector<double> rhs(n_col, 0.0);     // A^T H_w^{-1} r_w - r_lam
    std::vector<double> dlam(n_col, 0.0);
    std::vector<double> dw(n_row, 0.0);
    std::vector<double> Adlam(n_row, 0.0);

    // dsyevd workspace re-used across inner iterations; allocate after LWORK query.
    std::vector<double> S_eig;        // copy of S that dsyevd overwrites with V
    std::vector<double> eigvals(n_col, 0.0);
    std::vector<double> work;
    std::vector<int>    iwork;
    int dsy_lwork  = 0;
    int dsy_liwork = 0;

    for (int outer = 0; outer < outer_cap; outer++) {
        bool inner_converged = false;
        int  inner_count     = 0;
        double last_kkt      = std::numeric_limits<double>::infinity();
        int    last_n_proj   = 0;
        double last_alpha    = 0.0;

        for (int inner = 0; inner < inner_cap; inner++) {
            inner_count = inner + 1;
            total_inner++;

            // ----- Compute residuals r_w, r_lam at current (w, lam, mu).
            // r_w = log(w/d) + mu (1/(u-w) - 1/(w-l)) + A lam
            // r_lam = A^T w - b
            for (int c = 0; c < n_col; c++) r_lam[c] = -b[c];
            for (int i = 0; i < n_row; i++) {
                for (int idx = p[i]; idx < p[i+1]; idx++) {
                    r_lam[j[idx]] += x[idx] * w[i];
                }
            }

            for (int i = 0; i < n_row; i++) {
                double ulo = w[i] - lo[i];
                double uhi = hi[i] - w[i];
                // A * lam component for row i: sum over nz in row i.
                double A_lam_i = 0.0;
                for (int idx = p[i]; idx < p[i+1]; idx++) {
                    A_lam_i += x[idx] * lam[j[idx]];
                }
                r_w[i] = std::log(w[i] / d[i]) + mu * (1.0 / uhi - 1.0 / ulo) + A_lam_i;
            }

            // KKT residual ||r_w||_inf + ||r_lam||_inf.
            double r_w_inf = 0.0, r_lam_inf = 0.0;
            bool has_nan_inf = false;
            for (int i = 0; i < n_row; i++) {
                if (!std::isfinite(r_w[i]) || !std::isfinite(w[i])) { has_nan_inf = true; break; }
                double a = std::abs(r_w[i]);
                if (a > r_w_inf) r_w_inf = a;
            }
            if (!has_nan_inf) {
                for (int c = 0; c < n_col; c++) {
                    if (!std::isfinite(r_lam[c]) || !std::isfinite(lam[c])) { has_nan_inf = true; break; }
                    double a = std::abs(r_lam[c]);
                    if (a > r_lam_inf) r_lam_inf = a;
                }
            }
            if (has_nan_inf) {
                res.status_code = 2;
                res.status_msg  = "NaN/Inf detected in residuals";
                res.iterations  = total_inner;
                goto finalize;
            }
            last_kkt = r_w_inf + r_lam_inf;

            // Inner stop: ||r_w||_inf + ||r_lam||_inf < 10*mu.
            if (last_kkt < 10.0 * mu) {
                inner_converged = true;
                break;
            }

            // ----- Build H_w (diagonal) and Hinv = 1/H_w.
            // H_w_i = 1/w_i + mu/(w_i-l_i)^2 + mu/(u_i-w_i)^2
            for (int i = 0; i < n_row; i++) {
                double ulo = w[i] - lo[i];
                double uhi = hi[i] - w[i];
                double H_i = 1.0 / w[i] + mu / (ulo * ulo) + mu / (uhi * uhi);
                if (!std::isfinite(H_i) || H_i <= 0.0) {
                    res.status_code = 2;
                    res.status_msg  = "non-positive or non-finite H_w";
                    res.iterations  = total_inner;
                    goto finalize;
                }
                Hinv[i] = 1.0 / H_i;
            }

            // ----- Build dense Schur S = A^T H_w^{-1} A (m × m).
            // Block-incidence A has K=p[i+1]-p[i] nonzeros per row, all 1.0 (or x[idx]).
            // S[c1, c2] = sum_{i: A[i,c1]≠0 and A[i,c2]≠0} x[i,c1] * x[i,c2] / H_w_i.
            // Expand row-by-row outer products: S += Hinv[i] * (a_i a_i^T) where a_i is row i of A.
            S.setZero();
            for (int i = 0; i < n_row; i++) {
                double hi_inv = Hinv[i];
                for (int idx1 = p[i]; idx1 < p[i+1]; idx1++) {
                    int c1 = j[idx1];
                    double xv1 = x[idx1] * hi_inv;
                    for (int idx2 = p[i]; idx2 < p[i+1]; idx2++) {
                        int c2 = j[idx2];
                        S(c1, c2) += xv1 * x[idx2];
                    }
                }
            }

            // ----- Build RHS: rhs = A^T (H_w^{-1} r_w) - r_lam
            for (int c = 0; c < n_col; c++) rhs[c] = -r_lam[c];
            for (int i = 0; i < n_row; i++) {
                double v = Hinv[i] * r_w[i];
                for (int idx = p[i]; idx < p[i+1]; idx++) {
                    rhs[j[idx]] += x[idx] * v;
                }
            }

            // ----- TSVD on S via LAPACK dsyevd. Mirror Epic-Dβ WL-1.
            int dsy_n = n_col, dsy_lda = n_col, dsy_info = 0;
            S_eig.assign(S.data(), S.data() + (size_t)n_col * n_col);

            // LWORK query (only on first call; cache afterward).
            if (dsy_lwork == 0) {
                int qlw = -1, qliw = -1;
                double work_q = 0.0;
                int    iwork_q = 0;
                F77_CALL(dsyevd)("V", "U", &dsy_n, S_eig.data(), &dsy_lda, eigvals.data(),
                                 &work_q, &qlw, &iwork_q, &qliw, &dsy_info FCONE FCONE);
                if (dsy_info != 0) {
                    res.status_code = 2;
                    res.status_msg  = "dsyevd LWORK query failed";
                    res.iterations  = total_inner;
                    goto finalize;
                }
                dsy_lwork  = static_cast<int>(work_q);
                dsy_liwork = iwork_q;
                work.assign(std::max(1, dsy_lwork), 0.0);
                iwork.assign(std::max(1, dsy_liwork), 0);
                // S_eig was overwritten by query (no, dsyevd query does not overwrite);
                // restore from S anyway to be safe.
                S_eig.assign(S.data(), S.data() + (size_t)n_col * n_col);
            }

            F77_CALL(dsyevd)("V", "U", &dsy_n, S_eig.data(), &dsy_lda, eigvals.data(),
                             work.data(), &dsy_lwork, iwork.data(), &dsy_liwork, &dsy_info FCONE FCONE);
            if (dsy_info != 0) {
                res.status_code = 2;
                res.status_msg  = "dsyevd factorization failed";
                res.iterations  = total_inner;
                goto finalize;
            }

            // dsyevd returns eigenvalues ascending. lam_max = eigvals[n_col - 1].
            double lam_max = eigvals[n_col - 1];
            int n_keep = 0;
            int n_proj = 0;
            // Solve S^{-1} rhs in eigenbasis, truncating dims with eig < ratio * lam_max.
            // dlam = sum_{i kept} (V[:,i]^T rhs / eig_i) * V[:,i]
            for (int c = 0; c < n_col; c++) dlam[c] = 0.0;
            if (lam_max > 0.0) {
                double thresh = tsvd_ratio * lam_max;
                for (int e = 0; e < n_col; e++) {
                    double ev = eigvals[e];
                    if (ev < thresh) {
                        n_proj++;
                        continue;
                    }
                    n_keep++;
                    // V[:,e] is column e of S_eig (column-major), occupies S_eig[e*n_col .. e*n_col + n_col - 1].
                    const double* ve = S_eig.data() + (size_t)e * n_col;
                    double proj = 0.0;
                    for (int c = 0; c < n_col; c++) proj += ve[c] * rhs[c];
                    double coef = proj / ev;
                    for (int c = 0; c < n_col; c++) dlam[c] += coef * ve[c];
                }
            } else {
                // Degenerate spectrum: dlam = 0; all dims projected.
                n_proj = n_col;
            }
            last_n_proj = n_proj;

            // ----- Δw = -H_w^{-1} (r_w + A Δλ)
            // Adlam_i = sum_{idx in row i} x[idx] * dlam[j[idx]]
            for (int i = 0; i < n_row; i++) {
                double s = 0.0;
                for (int idx = p[i]; idx < p[i+1]; idx++) {
                    s += x[idx] * dlam[j[idx]];
                }
                Adlam[i] = s;
            }
            for (int i = 0; i < n_row; i++) {
                dw[i] = -Hinv[i] * (r_w[i] + Adlam[i]);
            }

            // NaN/Inf check on steps.
            bool step_bad = false;
            for (int i = 0; i < n_row; i++) {
                if (!std::isfinite(dw[i])) { step_bad = true; break; }
            }
            if (!step_bad) {
                for (int c = 0; c < n_col; c++) {
                    if (!std::isfinite(dlam[c])) { step_bad = true; break; }
                }
            }
            if (step_bad) {
                res.status_code = 2;
                res.status_msg  = "NaN/Inf detected in Newton step";
                res.iterations  = total_inner;
                goto finalize;
            }

            // ----- Fraction-to-boundary: alpha_max = max alpha s.t. l + delta <= w + alpha dw <= u - delta.
            double alpha_max = 1.0;
            for (int i = 0; i < n_row; i++) {
                if (dw[i] < 0.0) {
                    double slack = (w[i] - (lo[i] + delta));
                    if (slack <= 0.0) { alpha_max = 0.0; break; }
                    double a = -slack / dw[i];  // positive since dw[i] < 0
                    if (a < alpha_max) alpha_max = a;
                } else if (dw[i] > 0.0) {
                    double slack = ((hi[i] - delta) - w[i]);
                    if (slack <= 0.0) { alpha_max = 0.0; break; }
                    double a = slack / dw[i];
                    if (a < alpha_max) alpha_max = a;
                }
            }
            if (alpha_max <= 0.0) {
                // Already at boundary in step direction; cannot move. Treat as failure.
                res.status_code = 2;
                res.status_msg  = "fraction-to-boundary alpha_max = 0";
                res.iterations  = total_inner;
                goto finalize;
            }
            double alpha = 0.99 * alpha_max;
            last_alpha = alpha;

            // ----- Update w, lam.
            for (int i = 0; i < n_row; i++) w[i]   += alpha * dw[i];
            for (int c = 0; c < n_col; c++) lam[c] += alpha * dlam[c];
        } // inner

        // ----- Trace capture per outer iter (record n_projected_dims even when zero).
        bool record_trace = capture_trace && (outer % trace_stride == 0);
        if (record_trace) {
            // Recompute max_err on current weights: |A^T w / Z - b / Z|_inf where Z=sum(d).
            // (Use absolute error in margin frequency space.)
            std::vector<double> AT_w(n_col, 0.0);
            for (int i = 0; i < n_row; i++) {
                for (int idx = p[i]; idx < p[i+1]; idx++) {
                    AT_w[j[idx]] += x[idx] * w[i];
                }
            }
            double Z = 0.0;
            for (int i = 0; i < n_row; i++) Z += d[i];
            if (Z <= 0.0) Z = 1.0;
            double W_total = 0.0;
            for (int i = 0; i < n_row; i++) W_total += w[i];
            if (W_total <= 0.0) W_total = 1.0;
            double max_err = 0.0;
            for (int c = 0; c < n_col; c++) {
                double e = std::abs(AT_w[c] / W_total - b[c] / Z);
                if (e > max_err) max_err = e;
            }

            auto now = std::chrono::steady_clock::now();
            double time_ms = std::chrono::duration<double, std::milli>(now - start_time).count();

            res.trace_data.push_back((double)outer);
            res.trace_data.push_back((double)inner_count);
            res.trace_data.push_back(time_ms);
            res.trace_data.push_back(mu);
            res.trace_data.push_back(max_err);
            res.trace_data.push_back(last_kkt);
            res.trace_data.push_back((double)last_n_proj);
            res.trace_data.push_back(last_alpha);
        }

        res.iterations = total_inner;

        if (!inner_converged) {
            // Inner cap exceeded for this mu. status_code=1, halt.
            res.status_code = 1;
            res.status_msg  = "inner Newton cap exceeded";
            goto finalize;
        }

        // ----- Outer stop: mu < mu_floor AND KKT < 1e-7.
        if (mu < mu_floor && last_kkt < final_tol) {
            converged = true;
            res.status_code = 0;
            res.status_msg  = "converged";
            goto finalize;
        }

        // Decay barrier.
        mu *= 0.5;
    } // outer

    // Outer cap exceeded.
    if (!converged) {
        // If KKT meets final_tol already, treat as converged even though mu hasn't reached floor.
        // (Spec: status_code=0 requires both. Otherwise code=1.)
        res.status_code = 1;
        res.status_msg  = "outer barrier cap exceeded";
    }

finalize:
    res.weights = w;
    auto final_time = std::chrono::steady_clock::now();
    res.wall_time_ms = std::chrono::duration<double, std::milli>(final_time - start_time).count();
    return res;
}
