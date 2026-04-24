#include "lbw_config.h"
#include "ieppa.hpp"
#include "cell_table.hpp"
#include "leafblower.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fstream>
#include <limits>
#include <set>
#include <vector>

namespace lbw {

static std::vector<int> parse_trajectory_iters() {
    const char* s = std::getenv("LBW_TRAJECTORY_AT");
    if (!s || !*s) return {};
    std::vector<int> out;
    std::string buf;
    for (const char* p = s;; ++p) {
        if (*p == ',' || *p == '\0') {
            if (!buf.empty()) {
                try {
                    out.push_back(std::stoi(buf));
                } catch (const std::exception&) {
                    // skip malformed token
                }
                buf.clear();
            }
            if (*p == '\0') break;
        } else buf.push_back(*p);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}


static void write_trajectory_csv(
    const std::vector<std::pair<int,double>>& samples)
{
    if (samples.empty()) return;
    const char* path = std::getenv("LBW_TRAJECTORY_OUT");
    if (!path || !*path) return;
    std::ofstream f(path);
    if (!f.is_open()) return;
    f << "iter,errRp\n";
    for (const auto& p : samples)
        f << p.first << "," << p.second << "\n";
}

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
    // Intermediate homotopy levels stop early once errRp drops below this loose
    // threshold; only the final level uses the user's tol_abs for full precision.
    constexpr double kHomotopyIntermediateTol = 1e-5;

    IEPPAResult res;
    res.status = RK_ERR_NOCONV;
    res.iterations = 0;
    res.max_error = 1.0;
    res.M_cell = 0;
    res.n_cap_active = 0;
    res.n_xcur_writes_per_iter_linear = 0;
    res.min_alpha_seen = 1.0;
    res.final_alpha = 1.0;
    res.n_bounds_violated = 0;
    res.n_bounds_clamped  = 0;

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

    // Per-cell bounds. U_cell is recomputed per homotopy level below from
    // current_max_weight; L_cell is min_weight-only and level-invariant.
    std::vector<double> L_cell(ct.M_cell), U_cell(ct.M_cell);
    double lo = st.min_weight;
    for (int c = 0; c < ct.M_cell; c++) {
        L_cell[c] = lo * ct.n_per_cell[c];
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

    // Structural vs transient infeasibility split (stepstone finding):
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

    // P2.1 adaptive damping schedule. alpha = 1 / (1 + β · stress); stress = max streak.
    // β = 0.5: at stress=2, alpha=0.5; at stress=10, alpha≈0.17. Unlatched — recovers as
    // streaks reset. Preserves Peyré-Cuturi §4.4 convergence (alpha ∈ (0, 1]).
    // β is mutable so eta_schedule can vary it per homotopy level.
    double beta  = 0.5;
    double alpha = 1.0;
    // Test-only override (parallel to LBW_IEPPA_FORCE_PATH): "on"|"off"|unset.
    // Always compiled; microsecond getenv cost. Enables falsifiable
    // min_alpha_seen assertion (spec §7, CTO B5).
    const char* force_damp = std::getenv("LBW_IEPPA_FORCE_DAMPING");
    bool force_damping_on  = (force_damp != nullptr && std::strcmp(force_damp, "on")  == 0);
    bool force_damping_off = (force_damp != nullptr && std::strcmp(force_damp, "off") == 0);
    auto compute_alpha = [&]() -> double {
        if (force_damping_on)  return 0.5;
        if (force_damping_off) return 1.0;
        int stress = 0;
        for (int idx = 0; idx < total_cats; idx++) {
            if (infeas_streak[idx] > stress) stress = infeas_streak[idx];
        }
        if (stress == 0) return 1.0;
        // Cap stress at kInfeasPersistence to prevent unbounded streak growth from
        // driving alpha to near-zero on inputs with long transient infeasibility.
        // At stress=kInfeasPersistence/2=2 (original trigger point), alpha=0.5;
        // at stress=kInfeasPersistence=5, alpha≈0.29. Bounded, recoverable.
        int capped = std::min(stress, kInfeasPersistence);
        return 1.0 / (1.0 + beta * static_cast<double>(capped));
    };

    // Scratch for linear-space sweep: per-category accumulators (max categories per margin).
    int max_cat = 0;
    for (int k = 0; k < st.K; k++) if (st.cat_counts[k] > max_cat) max_cat = st.cat_counts[k];
    std::vector<double> S_lin(max_cat, 0.0);
    std::vector<double> rescale_lin(max_cat, 1.0);
    std::vector<double> inv_f_old_lin(max_cat, 1.0);

    // Trajectory probe state (internal-only; driven by LBW_TRAJECTORY_AT env var).
    // Probe deque + samples are shared across homotopy levels so iter counts are global.
    const std::vector<int> probe_targets = parse_trajectory_iters();
    std::deque<int> probe_queue(probe_targets.begin(), probe_targets.end());
    std::vector<std::pair<int,double>> probe_samples;

    // WU-B: X_prev tracks X[c] at the last convergence check for pct_change computation.
    // Initialized from X_init (uniform W[c]=1 at entry, X[c] = X_init[c]).
    std::vector<double> X_prev(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) X_prev[c] = X_init[c];

    // ── Homotopy outer driver ──
    // N=1 (default) degenerates to single level at max_weight — bit-identical to
    // prior flat-loop behaviour (net-zero identity). N>1 progressively tightens
    // max_weight from start_factor*max_weight down to end_factor*max_weight over
    // N levels with Chizat-inspired budget split (∝ (lvl+1)^p).
    // Drive homotopy from n_levels alone; enabled is derived from n_levels in the
    // bridge and checking it here is tautological (and creates a latent inconsistency
    // if callers query enabled independently, e.g. eta_schedule N>1 branch).
    const int    N_levels = (st.homotopy.n_levels > 1)
                            ? st.homotopy.n_levels : 1;
    const double k_start  = st.homotopy.start_factor;
    const double k_end    = st.homotopy.end_factor;
    const double p_budget = st.homotopy.budget_split_p;
    double budget_weight_sum = 0.0;
    for (int lvl = 0; lvl < N_levels; lvl++) {
        budget_weight_sum += std::pow(lvl + 1.0, p_budget);
    }
    const double alm_mu_base = st.alm_mu;
    int total_iters = 0;
    bool homotopy_break = false;

    for (int lvl = 0; lvl < N_levels && !homotopy_break; lvl++) {
        const double frac   = (N_levels == 1) ? 0.0
            : static_cast<double>(lvl) / static_cast<double>(N_levels - 1);
        const double factor = (N_levels == 1) ? 1.0
            : k_start * std::pow(k_end / k_start, frac);
        const double current_max_weight = st.max_weight * factor;
        const int budget_lvl = (N_levels == 1) ? st.inner_max_iter
            : std::max(1, static_cast<int>(std::round(
                static_cast<double>(st.inner_max_iter) *
                std::pow(lvl + 1.0, p_budget) / budget_weight_sum)));
        // Loose intermediate tol until final level (enables warm-jump mid-level exit).
        const double tol_lvl = (lvl == N_levels - 1) ? st.tol_abs : kHomotopyIntermediateTol;

        // Tang 2024 (arxiv:2403.05054 §4) dynamic-eta schedule.
        // Geometric interpolation: eta_i = eta_start * (eta_end/eta_start)^(frac^power).
        // frac=0 at lvl=0 → eta_i=eta_start; frac=1 at lvl=N-1 → eta_i=eta_end.
        // alm_mu>0: scale ALM penalty. alm_mu=0 (default): scale damping β instead,
        // so early levels apply stronger per-step smoothing and final level uses β=0.5.
        if (st.eta_schedule.mode == EtaScheduleMode::TANG_DYNAMIC && N_levels > 1) {
            const double scaled_frac = std::pow(frac, st.eta_schedule.schedule_power);
            const double eta_i = st.eta_schedule.eta_start *
                std::pow(st.eta_schedule.eta_end / st.eta_schedule.eta_start, scaled_frac);
            res.eta_final = eta_i;
            if (alm_mu_base > 0.0) {
                st.alm_mu = eta_i * alm_mu_base;
            } else {
                beta = 0.5 * eta_i;
            }
        }

        // Recompute U_cell for this level's current_max_weight.
        double hi = std::isfinite(current_max_weight) ? current_max_weight : 1e300;
        for (int c = 0; c < ct.M_cell; c++) {
            U_cell[c] = hi * ct.n_per_cell[c];
        }

        res.homotopy_levels_used  = lvl + 1;
        res.homotopy_final_factor = factor;

        bool level_converged = false;

        // WU-B Fix 1: reset X_prev at the start of each homotopy level so that
        // pct_change measures iteration-to-iteration shift within a level, not
        // cross-level drift from the previous level's final X.
        for (int c = 0; c < ct.M_cell; c++) X_prev[c] = X[c];

    for (int iter_in_lvl = 1; iter_in_lvl <= budget_lvl; iter_in_lvl++) {
        const int iter = total_iters + iter_in_lvl;
        res.iterations = iter;
        alpha = compute_alpha();
        if (alpha < res.min_alpha_seen) res.min_alpha_seen = alpha;

        // Margin sweep: branched by path.
        bool overflow_trip = false;

        // Do NOT read st.max_weight here — homotopy levels pass
        // current_max_weight indirectly via the shared U_cell already built.
        auto apply_single_margin_linear = [&](int k) -> bool {
            const int nj = st.cat_counts[k];
            const int off = cat_offset[k];
            std::fill(S_lin.begin(), S_lin.begin() + nj, 0.0);
            for (int j = 0; j < nj; j++) inv_f_old_lin[j] = 1.0 / f_lin[off + j];
            const int* gk = ct.g_per_cell[k].data();
            for (int c = 0; c < ct.M_cell; c++) {
                int j = gk[c];
                if (j < 0 || j >= nj) continue;
                if (X_init[c] <= 0.0) continue;
                S_lin[j] += X_cur[c] * inv_f_old_lin[j];
            }
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
                    return true;
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
                    return true;
                }
                rescale_lin[j] = new_f * inv_f_old_lin[j];
                f_lin[off + j] = new_f;
            }
            for (int c = 0; c < ct.M_cell; c++) {
                int j = gk[c];
                if (j < 0 || j >= nj) continue;
                if (X_init[c] <= 0.0) continue;
                X_cur[c] *= rescale_lin[j];
            }
            return false;
        };

        auto apply_single_margin_log = [&](int k) -> bool {
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
            return false;  // log path does not trip overflow mid-sweep
        };

        // Per-margin residual errRp_k = max_j |S_lin[j]/W_total - targets[k][j]|.
        // Uses the current X_cur (linear path) or rebuilt S via cells_by_margin_cat (log).
        auto compute_margin_errRp_linear = [&](int k) -> double {
            const int nj = st.cat_counts[k];
            std::fill(S_lin.begin(), S_lin.begin() + nj, 0.0);
            const int* gk = ct.g_per_cell[k].data();
            double W_total = 0.0;
            for (int c = 0; c < ct.M_cell; c++) W_total += X_cur[c];
            if (W_total <= 0.0) return 0.0;
            for (int c = 0; c < ct.M_cell; c++) {
                int j = gk[c];
                if (j >= 0 && j < nj) S_lin[j] += X_cur[c];
            }
            double err = 0.0;
            for (int j = 0; j < nj; j++) {
                double e = std::fabs(S_lin[j] / W_total - st.targets[k][j]);
                if (e > err) err = e;
            }
            return err;
        };

        auto compute_margin_errRp_log = [&](int k) -> double {
            double W_total = 0.0;
            for (int c = 0; c < ct.M_cell; c++) W_total += X_tilde[c];
            if (W_total <= 0.0) return 0.0;
            double err = 0.0;
            for (int j = 0; j < st.cat_counts[k]; j++) {
                const auto& cells = cells_by_margin_cat[cat_offset[k] + j];
                double Skj = 0.0;
                for (int c : cells) Skj += X_tilde[c];
                double e = std::fabs(Skj / W_total - st.targets[k][j]);
                if (e > err) err = e;
            }
            return err;
        };

        const bool use_greedy = (st.scheduler.mode == SchedulerMode::GREEDY);

        if (use_linear) {
            // WU-2 prefactored linear-space sweep (spec rev 5 §5).
            // Sequential (raking-style) access to X_cur for cache efficiency.
            // Each margin k: 2 sequential passes over M_cell cells = O(K*M_cell) total.
            // Per bucket (k,j): S_kj = sum X_cur[c]/f_kj accumulated in pass 1;
            //   X_cur[c] rescaled by f_new/f_old in pass 2. Identical to bucket loop.

            if (use_greedy) {
                std::vector<double> per_margin_err(st.K);
                for (int k = 0; k < st.K; k++) {
                    per_margin_err[k] = compute_margin_errRp_linear(k);
                }
                double initial_sum = 0.0;
                for (int k = 0; k < st.K; k++) initial_sum += per_margin_err[k];
                const double stop_threshold =
                    st.scheduler.residual_recheck_fraction * initial_sum;
                for (int step = 0; step < st.K && !overflow_trip; step++) {
                    int k_star = 0;
                    double best = per_margin_err[0];
                    for (int k = 1; k < st.K; k++) {
                        if (per_margin_err[k] > best) { best = per_margin_err[k]; k_star = k; }
                    }
                    if (apply_single_margin_linear(k_star)) {
                        overflow_trip = true; break;
                    }
                    res.greedy_sweeps_taken++;
                    per_margin_err[k_star] = compute_margin_errRp_linear(k_star);
                    double total_err = 0.0;
                    for (int k = 0; k < st.K; k++) total_err += per_margin_err[k];
                    if (total_err < stop_threshold) break;
                }
            } else {
                for (int k = 0; k < st.K && !overflow_trip; k++) {
                    if (apply_single_margin_linear(k)) overflow_trip = true;
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
                // structural_infeas_pairs NOT cleared — static bucket topology.
                // alpha recomputed at top of next iter by compute_alpha().
                // Reset iter_in_lvl: for-loop increment makes iter_in_lvl=1 next round
                // (preserves prior iter=0 semantics; global iter counter may regress
                // by one level's worth of iters but that matches pre-refactor behaviour).
                iter_in_lvl = 0;
                // WU-B Fix 2: reset X_prev after fallback — X semantics changed (log-path).
                for (int c = 0; c < ct.M_cell; c++) X_prev[c] = X[c];
                if (st.verbose >= 1) {
                    st.log("iEPPA: linear-space overflow trip; fallback to log-space.");
                }
                continue;  // skip the post-sweep X_tilde / capacity / errRp blocks this round
            }
        } else {
            // Log-space sweep via helper. Greedy branch needs X_tilde for residual
            // probe; rebuild it cheaply (residual-only; no capacity side-effect).
            if (use_greedy) {
                // Rebuild X_tilde from current lf for residual scoring.
                for (int c = 0; c < ct.M_cell; c++) {
                    if (X_init[c] <= 0.0) { X_tilde[c] = 0.0; continue; }
                    double s = log_X_init[c];
                    for (int m = 0; m < st.K; m++) {
                        int gm = ct.g_per_cell[m][c];
                        s += lf[cat_offset[m] + gm];
                    }
                    double s_clip = (s > kLogClip) ? kLogClip : s;
                    X_tilde[c] = std::exp(s_clip);
                }
                std::vector<double> per_margin_err(st.K);
                for (int k = 0; k < st.K; k++)
                    per_margin_err[k] = compute_margin_errRp_log(k);
                double initial_sum = 0.0;
                for (int k = 0; k < st.K; k++) initial_sum += per_margin_err[k];
                const double stop_threshold =
                    st.scheduler.residual_recheck_fraction * initial_sum;
                for (int step = 0; step < st.K; step++) {
                    int k_star = 0;
                    double best = per_margin_err[0];
                    for (int k = 1; k < st.K; k++) {
                        if (per_margin_err[k] > best) { best = per_margin_err[k]; k_star = k; }
                    }
                    (void) apply_single_margin_log(k_star);
                    res.greedy_sweeps_taken++;
                    // Refresh X_tilde for the touched margin's effect; cheap:
                    // lf changed only at k_star, so X_tilde multiplies by exp(lf_new-lf_old)
                    // per cell — simpler to recompute X_tilde fully (O(K*M_cell) but correct).
                    for (int c = 0; c < ct.M_cell; c++) {
                        if (X_init[c] <= 0.0) { X_tilde[c] = 0.0; continue; }
                        double s = log_X_init[c];
                        for (int m = 0; m < st.K; m++) {
                            int gm = ct.g_per_cell[m][c];
                            s += lf[cat_offset[m] + gm];
                        }
                        double s_clip = (s > kLogClip) ? kLogClip : s;
                        X_tilde[c] = std::exp(s_clip);
                    }
                    per_margin_err[k_star] = compute_margin_errRp_log(k_star);
                    double total_err = 0.0;
                    for (int k = 0; k < st.K; k++) total_err += per_margin_err[k];
                    if (total_err < stop_threshold) break;
                }
            } else {
                for (int k = 0; k < st.K; k++) {
                    (void) apply_single_margin_log(k);
                }
            }
        }

        if (use_linear) {
            // P1.1 fused block: X_tilde derived inline as X_cur / W; capacity + X_cur rebuild fused.
            bool overflow_detected = false;
            int n_cap = 0;
            for (int c = 0; c < ct.M_cell; c++) {
                if (X_init[c] <= 0.0 || W[c] <= 0.0) {
                    X[c] = 0.0;
                    X_cur[c] = 0.0;
                    W[c] = 1.0;
                    continue;
                }
                double X_tilde_c = X_cur[c] / W[c];
                if (!std::isfinite(X_tilde_c) || X_tilde_c > kLinearOverflowTrip) {
                    overflow_detected = true;
                    break;
                }
                if (X_tilde_c <= 0.0) {
                    X[c] = 0.0;
                    X_cur[c] = 0.0;
                    W[c] = 1.0;
                    continue;
                }
                double xc = std::clamp(X_tilde_c, L_cell[c], U_cell[c]);
                X[c] = xc;
                W[c] = xc / X_tilde_c;
                X_cur[c] = xc;
                res.n_xcur_writes_per_iter_linear++;
                if (xc != X_tilde_c) n_cap++;
            }
            res.n_cap_active = n_cap;
            if (overflow_detected) {
                // Full state reset on mid-loop break; partial writes to W/X/X_cur undone.
                std::fill(X_cur.begin(),   X_cur.end(),   0.0);
                std::fill(W.begin(),       W.end(),       1.0);
                std::fill(X.begin(),       X.end(),       0.0);
                std::fill(X_tilde.begin(), X_tilde.end(), 0.0);
                std::fill(lf.begin(),      lf.end(),      0.0);
                std::fill(f_lin.begin(),   f_lin.end(),   1.0);
                std::fill(infeas_streak.begin(), infeas_streak.end(), 0);
                res.n_xcur_writes_per_iter_linear = 0;
                use_linear = false;
                linear_fallback_used = true;
                // WU-B Fix 2: reset X_prev after fallback — X semantics changed (log-path).
                for (int c = 0; c < ct.M_cell; c++) X_prev[c] = X[c];
                if (st.verbose >= 1) st.log("iEPPA: linear-space overflow trip; fallback to log-space.");
                continue;
            }
        } else {
            // Log-path: X_tilde + capacity + X_cur unchanged from current implementation.
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
                    overflow_detected = true;
                }
                X_tilde[c] = std::exp(s_clip);
            }
            if (overflow_detected) {
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
        }

        // Convergence check.
        if (iter == 1 || iter % kErrCheckInterval == 0 || iter_in_lvl == budget_lvl) {
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

            // WU-B: compute pct_change (max relative shift in cell mass since last check).
            double pct_change = 0.0;
            if (W_total > 0.0) {
                for (int c = 0; c < ct.M_cell; c++) {
                    double rel = std::fabs(X[c] - X_prev[c]) / std::max(X_prev[c], 1e-12);
                    if (rel > pct_change) pct_change = rel;
                }
            }

            // WU-B: alternative metrics — accumulate per-margin S_k, then compute
            // mean_err (mean of per-margin max), kl_max (max per-margin KL), chi2.
            constexpr double kMetricEps = 1e-10;
            constexpr double kChi2Floor = 1.0;
            double mean_err_sum = 0.0;
            double kl_max       = 0.0;
            double chi2_total   = 0.0;
            if (W_total > 0.0) {
                for (int k = 0; k < st.K; k++) {
                    const int nj = st.cat_counts[k];
                    std::fill(S_lin.begin(), S_lin.begin() + nj, 0.0);
                    const int* gk = ct.g_per_cell[k].data();
                    for (int c = 0; c < ct.M_cell; c++) {
                        int j = gk[c];
                        if (j >= 0 && j < nj) S_lin[j] += X[c];
                    }
                    double max_k = 0.0;
                    double kl_k  = 0.0;
                    for (int j = 0; j < nj; j++) {
                        double S_p = S_lin[j] / W_total;
                        double T   = st.targets[k][j];
                        double err = std::fabs(S_p - T);
                        if (err > max_k) max_k = err;
                        if (T > 0.0) {
                            kl_k += T * std::log((T + kMetricEps) / (S_p + kMetricEps));
                        }
                        double obs = S_lin[j];
                        double exp_val = T * W_total;
                        chi2_total += (obs - exp_val) * (obs - exp_val) / (exp_val + kChi2Floor);
                    }
                    mean_err_sum += max_k;
                    if (kl_k > kl_max) kl_max = kl_k;
                }
            }
            double mean_err = (st.K > 0) ? (mean_err_sum / static_cast<double>(st.K)) : 0.0;

            // Store metrics in result struct.
            res.pct_change  = pct_change;
            res.mean_error  = mean_err;
            res.kl          = kl_max;
            res.chi2        = chi2_total;

            // Update X_prev AFTER computing pct_change.
            for (int c = 0; c < ct.M_cell; c++) X_prev[c] = X[c];

            // Trajectory probe: capture at requested iterations (captured on the
            // nearest check-interval iter >= requested iter, since errRp is only
            // computed at kErrCheckInterval boundaries, iter==1, and last iter).
            if (!probe_queue.empty() && iter >= probe_queue.front()) {
                bool first_write = true;
                while (!probe_queue.empty() && iter >= probe_queue.front()) {
                    if (first_write) {
                        probe_samples.push_back({iter, errRp});
                        first_write = false;
                    }
                    probe_queue.pop_front();
                }
            }
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
                              "  n_cap_active=%d", res.n_cap_active);
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

            // WU-B: active-criterion dispatch + stop_when ANY/ALL logic.
            // Replaces the old single-criterion `errRp < tol_lvl` check.
            // PCT criterion uses pct_change; all others use absolute_tol.
            // tol_lvl is the homotopy-level-adjusted absolute tolerance; for
            // intermediate levels it is kHomotopyIntermediateTol (loose), and for
            // the final level it equals st.tol_abs. The PCT gate always uses
            // the final pct_tol (not tol_lvl) because pct_change is scale-invariant.
            {
                double active_val = 0.0;
                const auto& cfg = st.convergence_cfg;
                if (cfg.absolute_tol > 0.0) {
                    switch (cfg.criterion) {
                        case lbw::CalibCriterion::MAX_ERR:  active_val = errRp;      break;
                        case lbw::CalibCriterion::MEAN_ERR: active_val = mean_err;   break;
                        case lbw::CalibCriterion::KL:       active_val = kl_max;     break;
                        case lbw::CalibCriterion::CHI2:     active_val = chi2_total; break;
                        case lbw::CalibCriterion::PCT:      active_val = pct_change; break;
                    }
                }

                // Legacy tol_lvl: applies to MAX_ERR on all homotopy levels.
                // For non-MAX_ERR criteria on intermediate levels, use the same
                // kHomotopyIntermediateTol gate via errRp (conservative: still exit
                // intermediate levels when errRp is loose enough).
                bool converged_abs = (cfg.absolute_tol > 0.0) && (active_val < cfg.absolute_tol);
                // Spec §1: PCT-only convergence is pct_change < pct_tol, no errRp floor.
                bool converged_pct = (cfg.pct_tol > 0.0) && (pct_change < cfg.pct_tol);
                // Intermediate homotopy levels also exit early when errRp < tol_lvl
                // regardless of criterion (warm-jump semantics; tol_lvl is loose).
                bool converged_intermediate = (errRp < tol_lvl);

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
                // For intermediate homotopy levels: always allow warm-jump when errRp loose.
                if (!converged && lvl < N_levels - 1) {
                    converged = converged_intermediate;
                }

                if (converged) {
                    // Early exit from this homotopy level. If this is the final level,
                    // set terminal status; else warm-jump to the next (tighter) level.
                    level_converged = true;
                    if (lvl == N_levels - 1) {
                        res.status = structural_infeas_pairs.empty() ? RK_OK : RK_ERR_INFEAS;
                    }
                    res.final_alpha = alpha;
                    break;
                }
            }
        }
        res.final_alpha = alpha;
    }
        // Log-space overflow: ieppa inner loop sets res.status = RK_ERR_NOCONV
        // and res.max_error = +inf via `break`. Must not proceed to next level.
        if (!std::isfinite(res.max_error)) {
            homotopy_break = true;
        }
        total_iters = res.iterations;  // global iter counter; picks up partial level
        if (level_converged && lvl == N_levels - 1) {
            homotopy_break = true;
        }
    }  // end homotopy level loop

    // Expansion to observation weights.
    std::vector<double> mult(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) {
        mult[c] = (X_init[c] > 0.0) ? X[c] / X_init[c] : 0.0;
    }
    for (int i = 0; i < st.n; i++) {
        st.weights[i] = st.weights[i] * mult[ct.cell_of[i]];
    }

    // Solver-owned normalization (moved from wrapper 2026-04-24 per user directive).
    // Applied AFTER expansion and BEFORE bounds_mode post-processing so unit-mode
    // water-fill sees final-scale weights and can strictly enforce
    // [min_weight, max_weight]. Contract: if total_w == 0 (degenerate: all zero
    // X_init[c] or all-zero mult[c]), weights remain unchanged (all zero) and
    // solver status is left as set upstream. X[c] is NOT rescaled — it is dead
    // after expansion (water-fill uses (void)target_sum below and redistributes
    // within each cell, so cell-aggregate scale is irrelevant).
    double total_w = 0.0;
    for (int i = 0; i < st.n; i++) total_w += st.weights[i];
    if (std::isfinite(total_w) && total_w > 0.0) {
        const double norm = static_cast<double>(st.n) / total_w;
        for (int i = 0; i < st.n; i++) st.weights[i] *= norm;
    }

    if (st.bounds_mode == RK_BOUNDS_CELL) {
        // Diagnostic scan: count violations without action.
        int violations = 0;
        for (int i = 0; i < st.n; i++) {
            if (st.weights[i] > st.max_weight || st.weights[i] < st.min_weight) {
                violations++;
            }
        }
        res.n_bounds_violated = violations;
    } else {
        // Unit mode: per-cell water-filling.
        // Water-fill redistributes excess within each cell, preserving the
        // post-normalize cell sum X[c] exactly via three-way classification
        // in the scan (violator / pinned / free). Only strictly-free obs
        // enter free_sum, so factor = 1 + excess/free_sum_really_free
        // distributes excess in full.
        // Build cells_of_obs (list of obs indices per cell) in one pass.
        std::vector<std::vector<int>> cells_of_obs(ct.M_cell);
        for (int i = 0; i < st.n; i++) cells_of_obs[ct.cell_of[i]].push_back(i);

        constexpr int kWaterFillMaxIter = 50;
        int total_clamped = 0;
        for (int c = 0; c < ct.M_cell; c++) {
            const auto& idxs = cells_of_obs[c];
            if (idxs.empty()) continue;

            for (int it = 0; it < kWaterFillMaxIter; it++) {
                double excess = 0.0;
                double free_sum = 0.0;
                int    n_free = 0;
                bool   any_violation = false;
                // Running counter: increment total_clamped at each normal-path
                // clamp (strict > max / < min). Pinned weights stay at bound
                // exactly, so the next-iter scan (strict-inequality here) and
                // the redistribute guard (strict < max && > min at "free"
                // branch below) both exclude them — no double count.
                // Pathological re-clamps (n_free==0, budget exhausted) are
                // redundant with these increments and do not re-count.
                for (int i : idxs) {
                    if (st.weights[i] > st.max_weight) {
                        excess += st.weights[i] - st.max_weight;
                        st.weights[i] = st.max_weight;
                        any_violation = true;
                        total_clamped++;
                    } else if (st.weights[i] < st.min_weight) {
                        excess -= st.min_weight - st.weights[i];
                        st.weights[i] = st.min_weight;
                        any_violation = true;
                        total_clamped++;
                    } else if (st.weights[i] == st.max_weight || st.weights[i] == st.min_weight) {
                        // Pinned from prior iter (set exactly via direct
                        // assignment above). Excluded from free_sum so
                        // factor = 1 + excess/free_sum_really_free
                        // distributes excess in full; cell-sum conservation
                        // holds exactly. FP equality is safe here because
                        // pinned obs was assigned, not computed. See
                        // leafblower-6s1o for the pre-fix under-distribution.
                    } else {
                        free_sum += st.weights[i];
                        n_free++;
                    }
                }
                if (!any_violation) break;
                if (n_free == 0 || free_sum <= 0.0) {
                    // Pathological: no room to redistribute. All violators
                    // already pinned by the scan above (line 585/589); cell
                    // sum may deviate from target by accumulated excess.
                    break;
                }
                // Redistribute excess proportionally over free observations.
                double factor = 1.0 + excess / free_sum;
                for (int i : idxs) {
                    if (st.weights[i] > st.min_weight && st.weights[i] < st.max_weight) {
                        st.weights[i] *= factor;
                    }
                }
                // Budget-exhaustion case (it == kWaterFillMaxIter - 1): if
                // factor-redistribution newly pushes a free obs above max,
                // the scan in the NEXT iter would clamp it — but there is no
                // next iter. In practice iter count is always enough because
                // water-fill converges geometrically for feasible problems.
            }
        }
        res.n_bounds_clamped = total_clamped;
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

    write_trajectory_csv(probe_samples);
    return res;
}

} // namespace lbw
