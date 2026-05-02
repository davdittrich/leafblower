#include "cp_calib.hpp"
#include <cmath>
#include <algorithm>
#include <chrono>
#include <limits>

CPResult cp_calibrate(
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
    CPResult res;
    res.weights.assign(n_row, 1.0);
    res.status_code = 1; // max_iter exhausted
    res.status_msg = "max iterations reached";
    res.iterations = max_iterations;
    res.wall_time_ms = 0.0;

    auto start_time = std::chrono::steady_clock::now();

    double Z = 0.0;
    for (int i = 0; i < n_row; i++) {
        Z += d[i];
    }
    if (Z <= 0.0) {
        res.status_code = 2;
        res.status_msg = "Z = sum(d) <= 0";
        return res;
    }

    // Power iteration for ||A||
    std::vector<double> v(n_col, 1.0 / std::sqrt((double)n_col));
    double lambda_old = 0.0;
    double rel_delta = 1.0;
    bool power_converged = false;
    for (int iter = 0; iter < 50; iter++) {
        std::vector<double> u(n_row, 0.0);
        for (int i = 0; i < n_row; i++) {
            double sum = 0.0;
            for (int idx = p[i]; idx < p[i+1]; idx++) {
                sum += x[idx] * v[j[idx]];
            }
            u[i] = sum;
        }
        std::vector<double> v_new(n_col, 0.0);
        for (int i = 0; i < n_row; i++) {
            for (int idx = p[i]; idx < p[i+1]; idx++) {
                v_new[j[idx]] += x[idx] * u[i];
            }
        }
        double lambda_new = 0.0;
        for (int c = 0; c < n_col; c++) {
            lambda_new += v[c] * v_new[c];
        }
        if (lambda_old > 0.0) {
            rel_delta = std::abs(lambda_new - lambda_old) / lambda_old;
            if (rel_delta < 1e-6) {
                power_converged = true;
                lambda_old = lambda_new;
                break;
            }
        }
        lambda_old = lambda_new;
        double norm_v_new = 0.0;
        for (int c = 0; c < n_col; c++) {
            norm_v_new += v_new[c] * v_new[c];
        }
        norm_v_new = std::sqrt(norm_v_new);
        for (int c = 0; c < n_col; c++) {
            v[c] = v_new[c] / norm_v_new;
        }
    }

    if (!power_converged && rel_delta > 1e-3) {
        res.status_code = 4;
        res.status_msg = "power iter divergence";
        return res;
    }

    double norm_A = std::sqrt(lambda_old);
    double norm_A_used = norm_A * 1.05;
    if (norm_A_used <= 0.0) {
        res.status_code = 2;
        res.status_msg = "||A|| <= 0";
        return res;
    }
    double sigma = 1.0 / norm_A_used;
    double tau = 1.0 / norm_A_used;

    const double delta = 1e-8;
    std::vector<double> w(n_row, 0.0);
    for (int i = 0; i < n_row; i++) {
        w[i] = std::clamp(d[i], lo[i] + delta, hi[i] - delta);
    }
    std::vector<double> w_bar = w;
    std::vector<double> y(n_col, 0.0);
    std::vector<double> w_ergodic = w;

    int trace_stride = std::max(1, (max_iterations + 999) / 1000);

    for (int k = 0; k < max_iterations; k++) {
        std::vector<double> AT_w_bar(n_col, 0.0);
        for (int i = 0; i < n_row; i++) {
            for (int idx = p[i]; idx < p[i+1]; idx++) {
                AT_w_bar[j[idx]] += x[idx] * w_bar[i];
            }
        }
        for (int c = 0; c < n_col; c++) {
            y[c] += sigma * (AT_w_bar[c] - b[c]);
        }

        std::vector<double> Ay(n_row, 0.0);
        for (int i = 0; i < n_row; i++) {
            double sum = 0.0;
            for (int idx = p[i]; idx < p[i+1]; idx++) {
                sum += x[idx] * y[j[idx]];
            }
            Ay[i] = sum;
        }

        double primal_stationarity_max = 0.0;
        bool has_nan_inf = false;

        for (int i = 0; i < n_row; i++) {
            double z_i = w[i] - tau * Ay[i];
            double w_old_i = w[i];
            double tau_d_i = tau * d[i];
            double w_n = std::clamp(d[i], lo[i] + delta, hi[i] - delta);

            if (tau_d_i <= 0.0) {
                w_n = z_i;
            } else if (std::abs(z_i / tau_d_i - 1.0) > 700.0) {
                w_n = z_i;
            } else {
                for (int p_iter = 0; p_iter < 20; p_iter++) {
                    double f_prime = tau_d_i * std::log(w_n / d[i]) + (w_n - z_i);
                    if (std::abs(f_prime) < 1e-12) break;
                    double f_double_prime = tau_d_i / w_n + 1.0;
                    w_n = w_n - f_prime / f_double_prime;
                    if (w_n <= 0.0) w_n = 1e-12;
                }
            }
            w_n = std::clamp(w_n, lo[i] + delta, hi[i] - delta);
            w[i] = w_n;
            w_bar[i] = 2.0 * w_n - w_old_i;
            w_ergodic[i] = w_ergodic[i] + (w_n - w_ergodic[i]) / (k + 2.0);

            double diff = std::abs(w[i] - w_old_i);
            if (diff > primal_stationarity_max) {
                primal_stationarity_max = diff;
            }
        }
        
        for (int c = 0; c < n_col; c++) {
            if (!std::isfinite(y[c])) has_nan_inf = true;
        }
        for (int i = 0; i < n_row; i++) {
            if (!std::isfinite(w[i])) has_nan_inf = true;
        }

        if (has_nan_inf) {
            res.status_code = 2;
            res.status_msg = "NaN/Inf detected";
            res.iterations = k + 1;
            break;
        }

        double primal_stationarity_proxy = primal_stationarity_max / tau;

        std::vector<double> AT_w(n_col, 0.0);
        for (int i = 0; i < n_row; i++) {
            for (int idx = p[i]; idx < p[i+1]; idx++) {
                AT_w[j[idx]] += x[idx] * w[i];
            }
        }

        double primal_resid = 0.0;
        for (int c = 0; c < n_col; c++) {
            double resid = std::abs(AT_w[c] - b[c]);
            if (resid > primal_resid) primal_resid = resid;
        }
        
        bool record_trace = capture_trace && (k % trace_stride == 0);
        if (record_trace) {
            std::vector<double> AT_w_ergodic(n_col, 0.0);
            for (int i = 0; i < n_row; i++) {
                for (int idx = p[i]; idx < p[i+1]; idx++) {
                    AT_w_ergodic[j[idx]] += x[idx] * w_ergodic[i];
                }
            }
            double W_total_last = 0.0;
            double W_total_ergodic = 0.0;
            for (int i = 0; i < n_row; i++) {
                W_total_last += w[i];
                W_total_ergodic += w_ergodic[i];
            }
            double max_err_last = 0.0;
            double max_err_ergodic = 0.0;
            for (int c = 0; c < n_col; c++) {
                double err_last = std::abs(AT_w[c] / std::max(W_total_last, 1.0) - b[c] / Z);
                if (err_last > max_err_last) max_err_last = err_last;
                double err_ergodic = std::abs(AT_w_ergodic[c] / std::max(W_total_ergodic, 1.0) - b[c] / Z);
                if (err_ergodic > max_err_ergodic) max_err_ergodic = err_ergodic;
            }

            auto now = std::chrono::steady_clock::now();
            double time_ms = std::chrono::duration<double, std::milli>(now - start_time).count();
            
            res.trace_data.push_back((double)k);
            res.trace_data.push_back(time_ms);
            res.trace_data.push_back(max_err_last);
            res.trace_data.push_back(max_err_ergodic);
            res.trace_data.push_back(primal_resid);
            res.trace_data.push_back(primal_stationarity_proxy);
        }

        if (primal_resid / std::max(Z, 1.0) < 1e-7) {
            res.status_code = 0;
            res.status_msg = "converged";
            res.iterations = k + 1;
            break;
        }

        res.iterations = k + 1;
    }

    res.weights = w;
    auto final_time = std::chrono::steady_clock::now();
    res.wall_time_ms = std::chrono::duration<double, std::milli>(final_time - start_time).count();

    return res;
}
