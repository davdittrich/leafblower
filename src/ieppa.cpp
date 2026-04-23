#include "lbw_config.h"
#include "ieppa.hpp"
#include "cell_table.hpp"
#include "leafblower.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <set>
#include <vector>

namespace lbw {

// Log-space algBCD at C=0 — see design spec §2.3, §8 for math.
// lf[k][j]: log Sinkhorn factor (per margin k, category j)
// W[c]:    capacity multiplier per cell (linear-space; bounded in [L_c/X_tilde, U_c/X_tilde])
// X_tilde[c] = X_init[c] * exp(sum_k lf[k][g_k(c)])
// X[c] = clamp(X_tilde[c] * W[c], L_c, U_c) — one capacity BCD block per outer iter
//
// Log-sum-exp stabilization on S_kj prevents overflow when partial log-sums
// approach log(DBL_MAX) ≈ 709.
IEPPAResult ieppa_solve(CalibState& st) {
    constexpr int    kErrCheckInterval = 10;
    constexpr double kEmptyBucketThreshold = 1e-15;
    constexpr double kLogClip = 700.0;  // exp(700) < DBL_MAX

    IEPPAResult res;
    res.status = RK_ERR_NOCONV;
    res.iterations = 0;
    res.max_error = 1.0;
    res.M_cell = 0;
    res.n_cap_active = 0;

    // Build cell table.
    CellTable ct;
    {
        int rc = build_cell_table(st.n, st.K, st.group_ids,
                                  st.cat_counts, st.weights, ct);
        if (rc != 0) {
            res.status = RK_ERR_BADARG;
            return res;
        }
    }
    res.M_cell = ct.M_cell;

    // Cell-aggregate initial weights.
    std::vector<double> X_init(ct.M_cell, 0.0);
    for (int i = 0; i < st.n; i++) {
        X_init[ct.cell_of[i]] += st.weights[i];
    }

    // Precompute log(X_init[c]) once; reused in margin sweep + X_tilde.
    std::vector<double> log_X_init(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) {
        log_X_init[c] = (X_init[c] > 0.0) ? std::log(X_init[c]) : -std::numeric_limits<double>::infinity();
    }

    // Per-cell bounds.
    std::vector<double> L_cell(ct.M_cell), U_cell(ct.M_cell);
    double hi = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
    double lo = st.min_weight;
    for (int c = 0; c < ct.M_cell; c++) {
        L_cell[c] = lo * ct.n_per_cell[c];
        U_cell[c] = hi * ct.n_per_cell[c];
    }

    // Log-space Sinkhorn factors per margin-category.
    int total_cats = 0;
    std::vector<int> cat_offset(st.K + 1, 0);
    for (int k = 0; k < st.K; k++) {
        cat_offset[k + 1] = cat_offset[k] + st.cat_counts[k] + 1;  // +1 for NA bucket
    }
    total_cats = cat_offset[st.K];
    std::vector<double> lf(total_cats, 0.0);  // lf[cat_offset[k] + j]

    // Per-cell capacity multiplier (linear-space).
    std::vector<double> W(ct.M_cell, 1.0);
    std::vector<double> X_tilde(ct.M_cell);
    std::vector<double> X(ct.M_cell);

    // Scratch for margin sweep.
    std::vector<std::vector<int>> cells_by_margin_cat(total_cats);
    for (int k = 0; k < st.K; k++) {
        for (int c = 0; c < ct.M_cell; c++) {
            int j = ct.g_per_cell[k][c];  // j in [0, cat_counts[k]] (NA → cat_counts[k])
            cells_by_margin_cat[cat_offset[k] + j].push_back(c);
        }
    }

    // Targets (user-provided, positive marginals per (k, j)).
    // Note: st.targets[k] has cat_counts[k] entries (no NA slot).
    // Paper: margin sum τ_{k,j} × W_total should match S_kj for j in [0, cat_counts[k]).
    // For NA bucket (j = cat_counts[k]), no constraint; f remains 1.0.

    bool is_infeasible = false;
    std::set<std::pair<int,int>> infeasible_pairs;  // dedup via set ordering

    if (st.verbose >= 1) {
        char msg[256];
        // Caller (c_api.cpp) sets st.ieppa_auto_selected=true when routing came
        // via AUTO; solver prepends [AUTO->iEPPA] marker. Otherwise plain entry.
        const char* prefix = (st.ieppa_auto_selected ? "[AUTO->iEPPA] " : "");
        std::snprintf(msg, sizeof(msg),
                      "%siEPPA: n=%d K=%d M_cell=%d compression=%.1fx",
                      prefix, st.n, st.K, ct.M_cell,
                      (double)st.n / (double)std::max(ct.M_cell, 1));
        st.log(msg);
    }

    std::vector<double> lv;
    lv.reserve(ct.M_cell);

    for (int iter = 1; iter <= st.inner_max_iter; iter++) {
        res.iterations = iter;

        // Margin sweep: one block per margin k.
        for (int k = 0; k < st.K; k++) {
            for (int j = 0; j < st.cat_counts[k]; j++) {
                const auto& cells = cells_by_margin_cat[cat_offset[k] + j];
                if (cells.empty()) {
                    if (st.targets[k][j] > 0.0) {
                        is_infeasible = true;
                        infeasible_pairs.emplace(k, j);
                    }
                    continue;
                }
                // log-sum-exp stabilization: compute lv_c = log X_init[c] + sum_{m!=k} lf[m]
                //                                       + log W[c]
                lv.assign(cells.size(), -std::numeric_limits<double>::infinity());
                double lv_max = -std::numeric_limits<double>::infinity();
                for (size_t r = 0; r < cells.size(); r++) {
                    int c = cells[r];
                    if (X_init[c] <= 0.0 || W[c] <= 0.0) {
                        continue;
                    }
                    double s = log_X_init[c] + std::log(W[c]);
                    for (int m = 0; m < st.K; m++) {
                        if (m == k) continue;
                        int gm = ct.g_per_cell[m][c];
                        s += lf[cat_offset[m] + gm];
                    }
                    lv[r] = s;
                    if (s > lv_max) lv_max = s;
                }
                if (!std::isfinite(lv_max)) {
                    // All cells are degenerate for this (k,j).
                    if (st.targets[k][j] > 0.0) {
                        is_infeasible = true;
                        infeasible_pairs.emplace(k, j);
                    }
                    continue;
                }
                double sum = 0.0;
                for (size_t r = 0; r < lv.size(); r++) {
                    if (std::isfinite(lv[r])) sum += std::exp(lv[r] - lv_max);
                }
                double log_S_kj = lv_max + std::log(sum);
                // Compare in log-space; exp(lv_max) * sum defeats LSE stabilization when lv_max → 700.
                double log_threshold = std::log(kEmptyBucketThreshold * ct.W_input);
                if (!std::isfinite(log_S_kj) || log_S_kj < log_threshold) {
                    if (st.targets[k][j] > 0.0) {
                        is_infeasible = true;
                        infeasible_pairs.emplace(k, j);
                    }
                    continue;
                }
                double log_target = std::log(st.targets[k][j] * ct.W_input);
                lf[cat_offset[k] + j] = log_target - log_S_kj;
            }
        }

        // Compute X_tilde via clip-before-exp. Detect overflow on uncapped cells.
        bool overflow_detected = false;
        double max_log_X_tilde = -std::numeric_limits<double>::infinity();
        for (int c = 0; c < ct.M_cell; c++) {
            if (X_init[c] <= 0.0) { X_tilde[c] = 0.0; continue; }
            double s = log_X_init[c];
            for (int m = 0; m < st.K; m++) {
                int gm = ct.g_per_cell[m][c];
                s += lf[cat_offset[m] + gm];
            }
            if (s > max_log_X_tilde) max_log_X_tilde = s;
            double s_clip = (s > kLogClip) ? kLogClip : s;
            if (s > kLogClip && U_cell[c] >= 1e299) {
                // Uncapped cell would overflow: log-factor drift beyond double precision.
                overflow_detected = true;
            }
            X_tilde[c] = std::exp(s_clip);
        }
        if (overflow_detected) {
            res.status = RK_ERR_NOCONV;
            res.max_error = std::numeric_limits<double>::infinity();
            if (st.verbose >= 2) {
                // Per design §8b: log-factor overflow event log lives in verbose=2.
                char msg[256];
                std::snprintf(msg, sizeof(msg),
                              "iEPPA: log-factor overflow (max_log_X_tilde=%.1f > 700) "
                              "indicates ill-conditioning; try looser max_weight or "
                              "tighter tol_abs, or method=raking.", max_log_X_tilde);
                st.log(msg);
            }
            break;
        }

        // Capacity block: X[c] = clamp(X_tilde[c], L_c, U_c); W[c] updated for next iter.
        int n_cap = 0;
        for (int c = 0; c < ct.M_cell; c++) {
            double xc = std::clamp(X_tilde[c], L_cell[c], U_cell[c]);
            X[c] = xc;
            if (X_tilde[c] > 0.0) {
                W[c] = xc / X_tilde[c];
            } else {
                W[c] = 1.0;
            }
            if (W[c] != 1.0) n_cap++;
        }
        res.n_cap_active = n_cap;

        // Convergence check.
        if (iter == 1 || iter % kErrCheckInterval == 0 || iter == st.inner_max_iter) {
            double W_total = 0.0;
            for (int c = 0; c < ct.M_cell; c++) W_total += X[c];
            double errRp = 0.0;
            for (int k = 0; k < st.K; k++) {
                for (int j = 0; j < st.cat_counts[k]; j++) {
                    const auto& cells = cells_by_margin_cat[cat_offset[k] + j];
                    double Skj = 0.0;
                    for (int c : cells) Skj += X[c];
                    double e = std::fabs(Skj / W_total - st.targets[k][j]);
                    if (e > errRp) errRp = e;
                }
            }
            res.max_error = errRp;
            if (st.verbose >= 1) {
                char msg[256];
                // Per design §8b: verbose=1 reports only errRp; n_cap lives in verbose=2.
                std::snprintf(msg, sizeof(msg),
                              "iEPPA iter %d: errRp=%.3e", iter, errRp);
                st.log(msg);
            }
            if (st.verbose >= 2) {
                // verbose=2: n_cap_active + log10 max/min of f[k][j] per margin for ill-conditioning debug.
                char msg[256];
                std::snprintf(msg, sizeof(msg),
                              "  n_cap_active=%d", n_cap);
                st.log(msg);
                for (int k = 0; k < st.K; k++) {
                    double lf_max = -std::numeric_limits<double>::infinity();
                    double lf_min =  std::numeric_limits<double>::infinity();
                    for (int j = 0; j < st.cat_counts[k]; j++) {
                        double v = lf[cat_offset[k] + j];
                        if (v > lf_max) lf_max = v;
                        if (v < lf_min) lf_min = v;
                    }
                    // log10(exp(v)) = v / ln(10)
                    std::snprintf(msg, sizeof(msg),
                                  "  margin=%d: log10(f) range [%.2f, %.2f]",
                                  k + 1,
                                  lf_min / 2.302585,
                                  lf_max / 2.302585);
                    st.log(msg);
                }
            }
            if (errRp < st.tol_abs) {
                res.status = is_infeasible ? RK_ERR_INFEAS : RK_OK;
                break;
            }
        }
    }

    // Expand to obs weights: w[i] = d[i] * X[c] / X_init[c]
    //
    // NOTE: No per-observation clamp here. The iEPPA solver enforces the
    // per-cell aggregate bound X[c] in [min_weight, max_weight] * |cell|.
    // Clamping per-observation weights (w[i] = d[i] * X[c] / X_init[c]) to
    // [min_weight, max_weight] would break the cell aggregate invariant
    // sum_{i in c} w[i] == X[c] whenever d[i] is non-uniform: the ratio
    // d_max / d_mean can legitimately push individual weights outside the
    // per-cell bounds even when the cell total is in range. Clamping silently
    // violated sum w == target marginals in that regime. See test
    // tests/testthat/test-ieppa-nonuniform-d.R.
    std::vector<double> mult(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) {
        mult[c] = (X_init[c] > 0.0) ? X[c] / X_init[c] : 0.0;
    }
    for (int i = 0; i < st.n; i++) {
        st.weights[i] = st.weights[i] * mult[ct.cell_of[i]];
    }

    if (is_infeasible && res.status == RK_ERR_NOCONV) {
        res.status = RK_ERR_INFEAS;
    }

    if (st.verbose >= 1) {
        const char* status_label =
            (res.status == RK_OK) ? "converged" :
            (res.status == RK_ERR_NOCONV) ? "max_iter exhausted (NOCONV)" :
            (res.status == RK_ERR_INFEAS) ? "infeasible" : "error";
        char msg[256];
        std::snprintf(msg, sizeof(msg),
                      "iEPPA %s in %d iters, errRp=%.3e",
                      status_label, res.iterations, res.max_error);
        st.log(msg);
    }

    if (st.verbose >= 1 && is_infeasible) {
        char msg[256];
        size_t off = 0;
        off += std::snprintf(msg + off, sizeof(msg) - off,
                             "iEPPA infeasible cells: ");
        size_t idx = 0;
        const size_t total = infeasible_pairs.size();
        for (auto it = infeasible_pairs.begin();
             it != infeasible_pairs.end() && off < sizeof(msg) - 32;
             ++it, ++idx) {
            off += std::snprintf(msg + off, sizeof(msg) - off,
                                 "margin=%d cat=%d%s",
                                 it->first + 1,
                                 it->second + 1,
                                 (idx + 1 < total) ? ", " : "");
        }
        st.log(msg);
    }

    return res;
}

} // namespace lbw
