#include "lbw_config.h"
#include "ieppa.hpp"
#include "cell_table.hpp"
#include "leafblower.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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

    // WU-2 dispatch: linear-space path when M_cell/n > 0.5 (i.e., compression <= 2x).
    // Env var LBW_IEPPA_FORCE_PATH in {"linear", "log"} overrides for tests; always
    // compiled (no #ifdef) -- getenv cost is microseconds, amortized over the solve.
    constexpr double kLinearSpaceThreshold = 2.0;  // compression ratio cutoff
    bool use_linear = (static_cast<double>(st.n) /
                       static_cast<double>(std::max(ct.M_cell, 1))) <
                      kLinearSpaceThreshold;
    if (const char* force = std::getenv("LBW_IEPPA_FORCE_PATH")) {
        if (std::strcmp(force, "linear") == 0) use_linear = true;
        else if (std::strcmp(force, "log") == 0) use_linear = false;
    }

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

    // Structural vs transient infeasibility split (post-WU-1 rev; stepstone finding):
    //   structural_infeas_pairs — bucket has cells.empty() AND target > 0. Static;
    //     a bucket with zero observations can never be satisfied. Latch immediately.
    //   infeas_streak — counts consecutive transient empties (log_S_kj < threshold).
    //     Drives WU-3 damping auto-trigger. Does NOT latch into persistent set;
    //     transient near-zeros arise from co-margin log-factor drift and often
    //     recover as the Sinkhorn iterates settle. Stepstone probe (commit state:
    //     mw=5, 500 iters, persistence disabled) confirmed iEPPA reaches
    //     errRp=2.24e-3 (best of 3 solvers) despite a bucket staying below
    //     threshold for 250+ iters — this is slow settling, not infeasibility.
    constexpr int kInfeasPersistence = 5;  // streak count that engages damping
    std::vector<int> infeas_streak(total_cats, 0);
    std::set<std::pair<int,int>> structural_infeas_pairs;

    auto record_empty = [&](int k, int j) {
        if (st.targets[k][j] <= 0.0) return;
        infeas_streak[cat_offset[k] + j]++;
    };
    auto record_nonempty = [&](int k, int j) {
        if (st.targets[k][j] <= 0.0) return;
        infeas_streak[cat_offset[k] + j] = 0;
    };

    // Pre-loop structural-infeas scan: buckets with cells.empty() AND target > 0
    // cannot be satisfied regardless of iteration. Latched once; static condition.
    for (int k = 0; k < st.K; k++) {
        for (int j = 0; j < st.cat_counts[k]; j++) {
            if (st.targets[k][j] > 0.0 &&
                cells_by_margin_cat[cat_offset[k] + j].empty()) {
                structural_infeas_pairs.emplace(k, j);
            }
        }
    }

    if (st.verbose >= 1) {
        char msg[256];
        // Caller (c_api.cpp) sets st.ieppa_auto_selected=true when routing came
        // via AUTO; solver prepends [AUTO->iEPPA] marker. Otherwise plain entry.
        const char* prefix = (st.ieppa_auto_selected ? "[AUTO->iEPPA] " : "");
        std::snprintf(msg, sizeof(msg),
                      "%siEPPA: n=%d K=%d M_cell=%d compression=%.1fx path=%s",
                      prefix, st.n, st.K, ct.M_cell,
                      (double)st.n / (double)std::max(ct.M_cell, 1),
                      use_linear ? "linear" : "log");
        st.log(msg);
    }

    std::vector<double> lv;
    lv.reserve(ct.M_cell);

    // Runtime trip: factor = X_init[c] * W[c] * prod_m f[m] <= max_X_init * trip^K < DBL_MAX/2.
    // Accounts for both K-way product accumulation AND X_init prefactor (Security B2).
    double max_X_init_val = 1.0;
    for (int c = 0; c < ct.M_cell; c++) {
        if (X_init[c] > max_X_init_val) max_X_init_val = X_init[c];
    }
    const double kLinearOverflowTrip = std::pow(
        std::numeric_limits<double>::max() / (2.0 * max_X_init_val),
        1.0 / static_cast<double>(st.K));

    // Linear-space Sinkhorn factors (mirror of lf[], but in linear domain).
    // Populated lazily -- only read by the linear-space sweep.
    std::vector<double> f_lin(total_cats, 1.0);
    // Prefactored per-cell accumulator X_cur[c] = X_init[c] * W[c] * prod_m f_lin[m][g_m(c)]
    // (spec rev 5 §5). Per-iter summation becomes O(|bucket|) instead of O(K*|bucket|),
    // cutting the linear-space per-iter cost from O(K^2*M_cell) to O(K*M_cell).
    std::vector<double> X_cur(ct.M_cell, 0.0);
    if (use_linear) {
        for (int c = 0; c < ct.M_cell; c++) {
            X_cur[c] = X_init[c] * W[c];  // f_lin == 1 at entry
        }
    }
    bool linear_fallback_used = false;

    // WU-3: adaptive damping. Default alpha=1.0 (stable mode, byte-identical to pre-WU-3).
    // Auto-switch to alpha=0.5 (damped mode) when any bucket's streak reaches
    // kInfeasPersistence/2 = 2. Latched per-solve; does not revert.
    double alpha = 1.0;
    bool damped_latched = false;
    // Test-only override (parallel to LBW_IEPPA_FORCE_PATH): "on"|"off"|unset.
    // Always compiled; microsecond getenv cost. Enables falsifiable
    // iter_damped > iter_stable assertion (spec §7, CTO B5).
    const char* force_damp = std::getenv("LBW_IEPPA_FORCE_DAMPING");
    bool force_damping_on  = (force_damp != nullptr && std::strcmp(force_damp, "on")  == 0);
    bool force_damping_off = (force_damp != nullptr && std::strcmp(force_damp, "off") == 0);
    if (force_damping_on) { alpha = 0.5; damped_latched = true; }
    auto maybe_engage_damping = [&]() {
        if (damped_latched || force_damping_off) return;
        for (int idx = 0; idx < total_cats; idx++) {
            if (infeas_streak[idx] >= kInfeasPersistence / 2) {
                alpha = 0.5;
                damped_latched = true;
                if (st.verbose >= 1) {
                    st.log("iEPPA: mid-streak detected; damping engaged (alpha=0.5).");
                }
                return;
            }
        }
    };

    // Scratch for linear-space sweep: per-category accumulators (max categories per margin).
    int max_cat = 0;
    for (int k = 0; k < st.K; k++) if (st.cat_counts[k] > max_cat) max_cat = st.cat_counts[k];
    std::vector<double> S_lin(max_cat, 0.0);
    std::vector<double> rescale_lin(max_cat, 1.0);
    std::vector<double> inv_f_old_lin(max_cat, 1.0);

    for (int iter = 1; iter <= st.inner_max_iter; iter++) {
        res.iterations = iter;
        maybe_engage_damping();

        // Margin sweep: branched by path.
        bool overflow_trip = false;

        if (use_linear) {
            // WU-2 prefactored linear-space sweep (spec rev 5 §5).
            // Sequential (raking-style) access to X_cur for cache efficiency.
            // Each margin k: 2 sequential passes over M_cell cells = O(K*M_cell) total.
            // Per bucket (k,j): S_kj = sum X_cur[c]/f_kj accumulated in pass 1;
            //   X_cur[c] rescaled by f_new/f_old in pass 2. Identical to bucket loop.

            for (int k = 0; k < st.K && !overflow_trip; k++) {
                const int nj = st.cat_counts[k];
                const int off = cat_offset[k];
                // Precompute reciprocals of old factors (one division per category, not per cell).
                std::fill(S_lin.begin(), S_lin.begin() + nj, 0.0);
                for (int j = 0; j < nj; j++) inv_f_old_lin[j] = 1.0 / f_lin[off + j];

                // Pass 1: accumulate S_kj via sequential X_cur scan.
                // g_per_cell[k][c] gives the category for cell c under margin k.
                const int* gk = ct.g_per_cell[k].data();
                for (int c = 0; c < ct.M_cell; c++) {
                    int j = gk[c];
                    if (j < 0 || j >= nj) continue;  // NA category
                    if (X_init[c] <= 0.0) continue;
                    S_lin[j] += X_cur[c] * inv_f_old_lin[j];
                }

                // Compute f_new and rescale factors per category.
                // Empty check uses S_lin threshold (covers both zero-cell buckets and
                // near-zero accumulations); cells_by_margin_cat not needed here.
                bool any_trip = false;
                for (int j = 0; j < nj; j++) {
                    if (!(S_lin[j] >= kEmptyBucketThreshold * ct.W_input) ||
                        !std::isfinite(S_lin[j])) {
                        record_empty(k, j);
                        rescale_lin[j] = 1.0;
                        continue;
                    }
                    record_nonempty(k, j);
                    double naive = (st.targets[k][j] * ct.W_input) / S_lin[j];
                    if (!std::isfinite(naive) || naive > kLinearOverflowTrip) {
                        overflow_trip = true;
                        any_trip = true;
                        break;
                    }
                    double new_f;
                    if (alpha == 1.0) {
                        new_f = naive;
                    } else {
                        double f_old = f_lin[off + j];
                        new_f = std::pow(f_old, 1.0 - alpha)
                              * std::pow(naive, alpha);
                    }
                    if (!std::isfinite(new_f) || new_f > kLinearOverflowTrip) {
                        overflow_trip = true;
                        any_trip = true;
                        break;
                    }
                    rescale_lin[j] = new_f * inv_f_old_lin[j];
                    f_lin[off + j] = new_f;
                }
                if (any_trip) break;

                // Pass 2: rescale X_cur via sequential scan.
                for (int c = 0; c < ct.M_cell; c++) {
                    int j = gk[c];
                    if (j < 0 || j >= nj) continue;
                    if (X_init[c] <= 0.0) continue;
                    X_cur[c] *= rescale_lin[j];
                }
            }
            if (overflow_trip && !linear_fallback_used) {
                // One-shot fallback: reset all solver state, switch to log-space,
                // restart outer loop from iter 0. State-clean list per spec rev 5 §5.
                linear_fallback_used = true;
                use_linear = false;
                std::fill(lf.begin(), lf.end(), 0.0);
                std::fill(f_lin.begin(), f_lin.end(), 1.0);
                std::fill(X_cur.begin(), X_cur.end(), 0.0);
                std::fill(W.begin(), W.end(), 1.0);
                for (int c = 0; c < ct.M_cell; c++) {
                    X_tilde[c] = X_init[c];
                    X[c] = X_init[c];
                }
                std::fill(infeas_streak.begin(), infeas_streak.end(), 0);
                // structural_infeas_pairs NOT cleared — static bucket topology
                alpha = 1.0;
                damped_latched = false;
                if (force_damping_on) { alpha = 0.5; damped_latched = true; }
                iter = 0;  // outer for-loop increments -> iter=1 next round
                if (st.verbose >= 1) {
                    st.log("iEPPA: linear-space overflow trip; fallback to log-space.");
                }
                continue;  // skip the post-sweep X_tilde / capacity / errRp blocks this round
            }
        } else {
            // Original log-space sweep (unchanged except record_empty/record_nonempty
            // hooks already installed by WU-1).
            for (int k = 0; k < st.K; k++) {
                for (int j = 0; j < st.cat_counts[k]; j++) {
                    const auto& cells = cells_by_margin_cat[cat_offset[k] + j];
                    if (cells.empty()) { record_empty(k, j); continue; }
                    lv.assign(cells.size(), -std::numeric_limits<double>::infinity());
                    double lv_max = -std::numeric_limits<double>::infinity();
                    for (size_t r = 0; r < cells.size(); r++) {
                        int c = cells[r];
                        if (X_init[c] <= 0.0 || W[c] <= 0.0) continue;
                        double s = log_X_init[c] + std::log(W[c]);
                        for (int m = 0; m < st.K; m++) {
                            if (m == k) continue;
                            int gm = ct.g_per_cell[m][c];
                            s += lf[cat_offset[m] + gm];
                        }
                        lv[r] = s;
                        if (s > lv_max) lv_max = s;
                    }
                    if (!std::isfinite(lv_max)) { record_empty(k, j); continue; }
                    double sum = 0.0;
                    for (size_t r = 0; r < lv.size(); r++) {
                        if (std::isfinite(lv[r])) sum += std::exp(lv[r] - lv_max);
                    }
                    double log_S_kj = lv_max + std::log(sum);
                    double log_threshold = std::log(kEmptyBucketThreshold * ct.W_input);
                    if (!std::isfinite(log_S_kj) || log_S_kj < log_threshold) {
                        record_empty(k, j);
                        continue;
                    }
                    record_nonempty(k, j);
                    double log_target = std::log(st.targets[k][j] * ct.W_input);
                    if (alpha == 1.0) {
                        lf[cat_offset[k] + j] = log_target - log_S_kj;
                    } else {
                        double lf_old = lf[cat_offset[k] + j];
                        lf[cat_offset[k] + j] =
                            (1.0 - alpha) * lf_old
                            + alpha * (log_target - log_S_kj);
                    }
                }
            }
        }

        // Compute X_tilde. Path-dependent:
        //   linear: X_cur[c] = X_init[c]*W[c]*prod_m f_lin; X_tilde = X_cur[c]/W[c] — O(M_cell).
        //   log:    clip-before-exp as before — O(K*M_cell).
        bool overflow_detected = false;
        double max_log_X_tilde = -std::numeric_limits<double>::infinity();
        for (int c = 0; c < ct.M_cell; c++) {
            if (X_init[c] <= 0.0) { X_tilde[c] = 0.0; continue; }
            if (use_linear) {
                // X_cur[c] = X_init[c] * W[c] * prod_m f_lin (maintained by sweep).
                // X_tilde = X_init[c] * prod_m f_lin = X_cur[c] / W[c].
                double v = (W[c] > 0.0) ? X_cur[c] / W[c] : 0.0;
                if (!std::isfinite(v)) overflow_detected = true;
                X_tilde[c] = v;
            } else {
                double s = log_X_init[c];
                for (int m = 0; m < st.K; m++) {
                    int gm = ct.g_per_cell[m][c];
                    s += lf[cat_offset[m] + gm];
                }
                if (s > max_log_X_tilde) max_log_X_tilde = s;
                double s_clip = (s > kLogClip) ? kLogClip : s;
                if (s > kLogClip && U_cell[c] >= 1e299) {
                    overflow_detected = true;
                }
                X_tilde[c] = std::exp(s_clip);
            }
        }
        if (overflow_detected) {
            if (use_linear && !linear_fallback_used) {
                // Treat same as overflow_trip: fall back once. State-clean per spec rev 5 §5.
                linear_fallback_used = true;
                use_linear = false;
                std::fill(lf.begin(), lf.end(), 0.0);
                std::fill(f_lin.begin(), f_lin.end(), 1.0);
                std::fill(X_cur.begin(), X_cur.end(), 0.0);
                std::fill(W.begin(), W.end(), 1.0);
                for (int c2 = 0; c2 < ct.M_cell; c2++) {
                    X_tilde[c2] = X_init[c2];
                    X[c2] = X_init[c2];
                }
                std::fill(infeas_streak.begin(), infeas_streak.end(), 0);
                // structural_infeas_pairs NOT cleared — static bucket topology
                alpha = 1.0;
                damped_latched = false;
                if (force_damping_on) { alpha = 0.5; damped_latched = true; }
                iter = 0;
                if (st.verbose >= 1) {
                    st.log("iEPPA: linear-space X_tilde overflow; fallback to log-space.");
                }
                continue;
            }
            res.status = RK_ERR_NOCONV;
            res.max_error = std::numeric_limits<double>::infinity();
            if (st.verbose >= 2) {
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

        // Update X_cur after capacity update (linear path only).
        // X_cur[c] = X_init[c] * W_new[c] * prod_m f_lin
        //          = X_tilde[c] * W_new[c]   (X_tilde already computed above)
        // O(M_cell) — no K-loop. Spec §5: "scale X_cur[c] by W_new/W_old".
        if (use_linear) {
            for (int c = 0; c < ct.M_cell; c++) {
                X_cur[c] = X_tilde[c] * W[c];  // W[c] is now W_new after capacity block
            }
        }

        // Convergence check.
        if (iter == 1 || iter % kErrCheckInterval == 0 || iter == st.inner_max_iter) {
            double W_total = 0.0;
            for (int c = 0; c < ct.M_cell; c++) W_total += X[c];
            double errRp = 0.0;
            if (use_linear) {
                // Sequential (raking-style) errRp: 2 sequential passes over M_cell per margin.
                // Avoids stride-5 cache pattern of bucket approach.
                for (int k = 0; k < st.K; k++) {
                    const int nj = st.cat_counts[k];
                    const int off = cat_offset[k];
                    std::fill(S_lin.begin(), S_lin.begin() + nj, 0.0);
                    const int* gk = ct.g_per_cell[k].data();
                    for (int c = 0; c < ct.M_cell; c++) {
                        int j = gk[c];
                        if (j >= 0 && j < nj) S_lin[j] += X[c];
                    }
                    for (int j = 0; j < nj; j++) {
                        double e = std::fabs(S_lin[j] / W_total - st.targets[k][j]);
                        if (e > errRp) errRp = e;
                    }
                }
            } else {
                for (int k = 0; k < st.K; k++) {
                    for (int j = 0; j < st.cat_counts[k]; j++) {
                        const auto& cells = cells_by_margin_cat[cat_offset[k] + j];
                        double Skj = 0.0;
                        for (int c : cells) Skj += X[c];
                        double e = std::fabs(Skj / W_total - st.targets[k][j]);
                        if (e > errRp) errRp = e;
                    }
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
                res.status = structural_infeas_pairs.empty() ? RK_OK : RK_ERR_INFEAS;
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

    if (!structural_infeas_pairs.empty() && res.status == RK_ERR_NOCONV) {
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

    if (st.verbose >= 1 && !structural_infeas_pairs.empty()) {
        char msg[256];
        size_t off = 0;
        off += std::snprintf(msg + off, sizeof(msg) - off,
                             "iEPPA persistent infeasible cells: ");
        size_t idx = 0;
        const size_t total = structural_infeas_pairs.size();
        for (auto it = structural_infeas_pairs.begin();
             it != structural_infeas_pairs.end() && off < sizeof(msg) - 32;
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
