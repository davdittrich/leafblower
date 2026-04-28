#include "lbw_config.h"
#include "ieppa.hpp"
#include "calib_dispatch.hpp"
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
    std::vector<double> X_tilde;  // deferred: allocated at first log-path/fallback use
    std::vector<double> X(ct.M_cell);
    // T1.B: per-cell log-product shadow. cell_lf[c] = Σ_k lf[k][g_k(c)].
    // Reuses lf[] (line 135) as log(f_lin) in the linear path.
    // cell_lf also used by T2.A in the log path.
    std::vector<double> cell_lf(ct.M_cell, 0.0);
    // High-water mark: max_c(log_X_init[c] + cell_lf[c]) ≈ max_c log(X_tilde[c]).
    // Monotone-nondecreasing between corrections — stale-high is intentional.
    double cell_lf_hwm = std::numeric_limits<double>::lowest();

    // ADMM dual variable for capacity constraint.
    // u[c] accumulates X_tilde - z violations; converges to 0 at fixed point.
    std::vector<double> u;  // allocated only for ieppa_soft (ADMM)
    if (st.use_admm_capacity) u.assign(ct.M_cell, 0.0);

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
    // Renorm threshold: halfway in log-space before the trip.
    // At this point, f_lin geometric mean = sqrt(trip); product f_lin^K = trip^(K/2),
    // leaving sqrt(trip) multiplicative headroom before the trip fires.
    const double kLinearOverflowThreshold = std::sqrt(kLinearOverflowTrip);
    const double kLogOverflowThreshold = std::log(kLinearOverflowThreshold);

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

    // WU-C: improvement/plateau rule — track active metric across kErrCheckInterval checks.
    // Initialized to +∞ so first check never triggers convergence.
    double prev_metric_for_rule = std::numeric_limits<double>::infinity();

    // WU-E: best-iterate tracking (cell-level W snapshot at min observed active metric).
    // Tracks the active convergence metric so best_weights minimizes the user's chosen objective.
    // W_best is initialized to all-zeros; best_metric_seen = ∞ ensures first check triggers copy.
    double best_metric_seen    = std::numeric_limits<double>::infinity();
    double best_objective_seen = 0.0;  // weight KL at best_iter (A1 fix)
    int    best_iter_val   = 0;
    std::vector<double> W_best(ct.M_cell, 0.0);

    // WU-D: SOR adaptive under-relaxation state (iEPPA-only).
    // Per-margin omega[k] starts at omega_init (1.0 = no damping = fast path).
    // Adaptation: sign-flip in per-margin errRp trajectory → omega *= 0.7 (floor: omega_min).
    // Monotone decrease → omega *= 1.05, capped at 1.0 (recovery).
    // Adaptation is suppressed for sor_burnin iterations so early transient oscillation
    // (driven by infeas-streak damping) does not prematurely reduce omega.
    const bool sor_active     = st.sor_cfg.enabled;
    const bool sor_auto       = st.sor_cfg.auto_adapt;
    const double omega_init_v = st.sor_cfg.omega_init;    // default 1.0
    const double omega_min_v  = st.sor_cfg.omega_min;     // default 0.3
    const double omega_fixed_v = st.sor_cfg.omega_fixed;  // -1.0 = use auto
    const int    sor_burnin_v  = st.sor_cfg.burnin;       // default 20
    std::vector<double> sor_omega(st.K, omega_init_v);
    std::vector<double> sor_prev_errRp(st.K, std::numeric_limits<double>::infinity());
    std::vector<bool>   sor_prev_decreasing(st.K, false);
    double sor_min_omega = 1.0;
    int    sor_n_damped  = 0;

    // Weight-space KL objective: Σ_c X[c]*log(X[c]/X_init[c])/n
    // Distinct from m.kl (marginal KL) — this is what ieppa actually minimizes.
    auto compute_weight_kl = [&]() -> double {
        double wkl = 0.0;
        const double inv_n = 1.0 / static_cast<double>(st.n);
        for (int c = 0; c < ct.M_cell; c++) {
            if (X_init[c] > 0.0 && X[c] > 0.0)
                wkl += X[c] * std::log(X[c] / X_init[c]) * inv_n;
        }
        return std::isfinite(wkl) ? wkl : 0.0;
    };

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
        // WU-C: reset improvement/plateau baseline at each homotopy level.
        prev_metric_for_rule = std::numeric_limits<double>::infinity();

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
            // WU-D: effective omega for this margin.
            // sor_active==false → eff_omega=1.0 (fast path, no pow()).
            // sor_active && !sor_auto && omega_fixed_v>0 → fixed omega.
            // sor_active && sor_auto → per-margin adaptive omega[k].
            double eff_omega;
            if (!sor_active) {
                eff_omega = 1.0;
            } else if (!sor_auto && omega_fixed_v > 0.0) {
                eff_omega = omega_fixed_v;
            } else {
                eff_omega = sor_omega[k];
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
                    // WU-D: SOR under-relaxation in linear space.
                    // ratio = T_kj*W_input / S_kj = naive (the full Sinkhorn step).
                    // Applying pow(ratio, eff_omega) is equivalent to a fractional step
                    // in log-space: lf_new = lf_old + eff_omega*(log(T)-log(S)).
                    // Fast path: eff_omega==1.0 skips pow() entirely (no regression for
                    // users who do not enable SOR, preserving the prior O(1) code path).
                    if (eff_omega == 1.0) {
                        new_f = naive;
                    } else {
                        double f_old = f_lin[off + j];
                        // new_f = f_old * (naive/f_old)^eff_omega = f_old^(1-eff_omega) * naive^eff_omega
                        new_f = std::pow(f_old, 1.0 - eff_omega)
                              * std::pow(naive, eff_omega);
                    }
                } else {
                    double f_old = f_lin[off + j];
                    // alpha damping (P2.1 infeas-streak) and SOR eff_omega compose:
                    // net exponent on the Sinkhorn ratio is alpha*eff_omega.
                    double net = alpha * eff_omega;
                    if (net == 1.0) {
                        new_f = naive;
                    } else {
                        new_f = std::pow(f_old, 1.0 - net)
                              * std::pow(naive, net);
                    }
                }
                if (!std::isfinite(new_f) || new_f > kLinearOverflowTrip) {
                    return true;
                }
                rescale_lin[j] = new_f * inv_f_old_lin[j];
                f_lin[off + j] = new_f;
                // T1.B: update lf shadow and propagate delta to cell_lf.
                // new_f > 0 guaranteed: empty buckets continue before reaching here,
                // and the overflow check above ensures isfinite(new_f).
                // Guard new_f > 0 (not 1e-300) to keep lf consistent even for subnormals —
                // log(subnormal) is finite and correct; skipping lf update would diverge
                // cell_lf from Σ lf[k][g_k(c)] after an over-correction.
                if (new_f > 0.0) {
                    double lf_new = std::log(new_f);
                    double delta = lf_new - lf[off + j];
                    lf[off + j] = lf_new;
                    if (std::fabs(delta) > 1e-12) {
                        if (delta > 0.0) {
                            // Positive delta: update cell_lf and high-water mark together.
                            for (int c : cells_by_margin_cat[off + j]) {
                                cell_lf[c] += delta;
                                double val = cell_lf[c] + log_X_init[c];
                                if (std::isfinite(val) && val > cell_lf_hwm)
                                    cell_lf_hwm = val;
                            }
                        } else {
                            for (int c : cells_by_margin_cat[off + j])
                                cell_lf[c] += delta;
                        }
                    }
                }
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
            // WU-D: effective omega for this margin (log-space path).
            double eff_omega_log;
            if (!sor_active) {
                eff_omega_log = 1.0;
            } else if (!sor_auto && omega_fixed_v > 0.0) {
                eff_omega_log = omega_fixed_v;
            } else {
                eff_omega_log = sor_omega[k];
            }
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
                // WU-D: in log-space, SOR applies as a fractional step:
                // lf_new = lf_old + net*(log_target - log_S_kj), where net = alpha * eff_omega_log.
                // alpha=1 && eff_omega_log=1 → lf_new = log_target - log_S_kj (fast path, no change).
                double net_log = alpha * eff_omega_log;
                {
                    double lf_old = lf[cat_offset[k] + j];
                    double lf_new = (net_log == 1.0)
                        ? (log_target - log_S_kj)
                        : ((1.0 - net_log) * lf_old + net_log * (log_target - log_S_kj));
                    lf[cat_offset[k] + j] = lf_new;
                    // T2.A: maintain cell_lf incrementally (eliminates K=20 DRAM streams)
                    double delta = lf_new - lf_old;
                    if (std::fabs(delta) > 1e-12) {
                        for (int c : cells_by_margin_cat[cat_offset[k] + j])
                            cell_lf[c] += delta;
                    }
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
            if (X_tilde.empty()) return std::numeric_limits<double>::infinity();
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

            // T1.B: correct before X_tilde product overflows.
            // Fires when max_c log(X_init[c] * Π_k f_lin[k][g_k(c)]) >= log(threshold).
            // shift > 0 guaranteed by >=; x_scale in (0,1].
            if (!overflow_trip && cell_lf_hwm >= kLogOverflowThreshold) {
                double shift = cell_lf_hwm - kLogOverflowThreshold;
                double lf_correction = -shift / static_cast<double>(st.K);
                double x_scale = std::exp(-shift);
                // j <= cat_counts[k]: include NA bucket (cat_offset has +1 per margin).
                // Without NA shift, cells NA for a margin would have cell_lf decremented
                // by full shift but lf[k][NA] unchanged — invariant violated.
                for (int k = 0; k < st.K; k++) {
                    for (int j = 0; j <= st.cat_counts[k]; j++) {
                        lf[cat_offset[k] + j] += lf_correction;
                        f_lin[cat_offset[k] + j] = std::exp(lf[cat_offset[k] + j]);
                    }
                }
                // X_cur *= exp(-shift): maintains X_cur = X_init × W × Π f_lin.
                for (int c = 0; c < ct.M_cell; c++) {
                    cell_lf[c] -= shift;
                    X_cur[c] *= x_scale;
                }
                // Reset to lowest so hwm repopulates honestly from next positive delta,
                // avoiding spurious re-trigger when true max is exactly at threshold.
                cell_lf_hwm = std::numeric_limits<double>::lowest();
                if (st.use_admm_capacity) std::fill(u.begin(), u.end(), 0.0);
                if (st.verbose >= 2) {
                    char msg[128];
                    std::snprintf(msg, sizeof(msg), "iEPPA T1.B renorm shift=%.2e", shift);
                    st.log(msg);
                }
            }

            if (overflow_trip && !linear_fallback_used) {
                // One-shot fallback: reset all solver state, switch to log-space,
                // restart outer loop from iter 0. State-clean list per spec rev 5 §5.
                linear_fallback_used = true;
                use_linear = false;
                std::fill(lf.begin(), lf.end(), 0.0);
                std::fill(cell_lf.begin(), cell_lf.end(), 0.0);
                cell_lf_hwm = std::numeric_limits<double>::lowest();
                if (st.use_admm_capacity) std::fill(u.begin(), u.end(), 0.0);
                std::fill(f_lin.begin(), f_lin.end(), 1.0);
                std::fill(X_cur.begin(), X_cur.end(), 0.0);
                std::fill(W.begin(), W.end(), 1.0);
                if (X_tilde.empty()) X_tilde.assign(ct.M_cell, 0.0);
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
                if (X_tilde.empty()) X_tilde.assign(ct.M_cell, 0.0);
                for (int c = 0; c < ct.M_cell; c++) {
                    if (X_init[c] <= 0.0) { X_tilde[c] = 0.0; continue; }
                    // T2.A: single-stream exp via cell_lf (was K=20 DRAM streams)
                    double s = log_X_init[c] + cell_lf[c];
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
                    // T2.A: lf delta is captured in cell_lf[c], single-stream exp
                    for (int c = 0; c < ct.M_cell; c++) {
                        if (X_init[c] <= 0.0) { X_tilde[c] = 0.0; continue; }
                        double s = log_X_init[c] + cell_lf[c];
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
                if (st.use_admm_capacity) {
                    // ADMM: z-update projects adjusted X_tilde onto capacity box.
                    // u[c] accumulates violations across iterations; converges to 0.
                    double z = std::clamp(X_tilde_c + u[c], L_cell[c], U_cell[c]);
                    u[c] += X_tilde_c - z;
                    X[c] = z; W[c] = z / X_tilde_c; X_cur[c] = z;
                } else {
                    double xc = std::clamp(X_tilde_c, L_cell[c], U_cell[c]);
                    X[c] = xc; W[c] = xc / X_tilde_c; X_cur[c] = xc;
                }
                res.n_xcur_writes_per_iter_linear++;
                if (X[c] != X_tilde_c) n_cap++;
            }
            res.n_cap_active = n_cap;
            if (overflow_detected) {
                // Full state reset on mid-loop break; partial writes to W/X/X_cur undone.
                std::fill(X_cur.begin(),   X_cur.end(),   0.0);
                std::fill(W.begin(),       W.end(),       1.0);
                std::fill(X.begin(),       X.end(),       0.0);
                std::fill(X_tilde.begin(), X_tilde.end(), 0.0);
                std::fill(lf.begin(),      lf.end(),      0.0);
                std::fill(cell_lf.begin(), cell_lf.end(), 0.0);
                cell_lf_hwm = std::numeric_limits<double>::lowest();
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
            if (X_tilde.empty()) X_tilde.assign(ct.M_cell, 0.0);
            for (int c = 0; c < ct.M_cell; c++) {
                if (X_init[c] <= 0.0) { X_tilde[c] = 0.0; continue; }
                // T2.A: single-stream exp via cell_lf (was K=20 DRAM streams)
                double s = log_X_init[c] + cell_lf[c];
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
            double marg_kl = 0.0;
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
                    // Marginal KL: Σ_k Σ_j t_kj log(t_kj / achieved_kj)
                    for (int j = 0; j < nj; j++) {
                        double tkj = st.targets[k][j];
                        double skj = (W_total > 0.0) ? S_lin[j] / W_total : 0.0;
                        if (tkj > 1e-300 && skj > 1e-300)
                            marg_kl += tkj * std::log(tkj / skj);
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
                        // Marginal KL contribution for this (k, j) cell
                        double tkj = st.targets[k][j];
                        double skj = (W_total > 0.0) ? Skj / W_total : 0.0;
                        if (tkj > 1e-300 && skj > 1e-300)
                            marg_kl += tkj * std::log(tkj / skj);
                    }
                }
            }
            res.marginal_kl_at_iter = marg_kl;
            res.max_error = errRp;

            // WU-E / g4oj: BLOCK 1 — MAX_ERR best-iterate (errRp always valid here,
            // outside need_extra_metrics gate). Tracks min errRp when MAX_ERR is active.
            if (st.convergence_cfg.metric == lbw::CalibMetric::MAX_ERR) {
                if (errRp < best_metric_seen) {
                    // === BEST-ITER UPDATE (metric, iter, objective MUST stay co-located) ===
                    best_metric_seen    = errRp;
                    best_iter_val       = iter;
                    best_objective_seen = compute_weight_kl();
                    for (int c = 0; c < ct.M_cell; c++)
                        W_best[c] = (X_init[c] > 0.0) ? X[c] / X_init[c] : 0.0;
                    // === END BEST-ITER UPDATE ===
                }
            }
            // BLOCK 1b — MARGINAL_KL best-iterate (marg_kl always valid alongside errRp).
            // Tracks min marginal KL when MARGINAL_KL is active.
            if (st.convergence_cfg.metric == lbw::CalibMetric::MARGINAL_KL) {
                if (res.marginal_kl_at_iter < best_metric_seen) {
                    // === BEST-ITER UPDATE ===
                    best_metric_seen    = res.marginal_kl_at_iter;
                    best_iter_val       = iter;
                    best_objective_seen = compute_weight_kl();
                    for (int c = 0; c < ct.M_cell; c++)
                        W_best[c] = (X_init[c] > 0.0) ? X[c] / X_init[c] : 0.0;
                    // === END BEST-ITER UPDATE ===
                }
            }

            // WU-D: per-margin omega adaptation (both linear and log paths; auto mode only;
            // suppressed during burnin to let the infeas-streak damping settle first).
            if (sor_active && sor_auto && iter >= sor_burnin_v) {
                if (W_total > 0.0) {
                    for (int k = 0; k < st.K; k++) {
                        const int nj_k = st.cat_counts[k];
                        // Compute per-margin errRp_k using X[c] (post-capacity, available
                        // in both paths). Reuse S_lin scratch (sized to max_cat >= nj_k).
                        std::fill(S_lin.begin(), S_lin.begin() + nj_k, 0.0);
                        const int* gk_s = ct.g_per_cell[k].data();
                        for (int c = 0; c < ct.M_cell; c++) {
                            int j = gk_s[c];
                            if (j >= 0 && j < nj_k) S_lin[j] += X[c];
                        }
                        double errRp_k = 0.0;
                        for (int j = 0; j < nj_k; j++) {
                            double e = std::fabs(S_lin[j] / W_total - st.targets[k][j]);
                            if (e > errRp_k) errRp_k = e;
                        }
                        bool decreasing = (errRp_k < sor_prev_errRp[k]);
                        bool sign_flip  = !decreasing && sor_prev_decreasing[k];
                        if (sign_flip) {
                            // Oscillation detected: damp omega by 0.7, clamp to floor.
                            sor_omega[k] = std::max(omega_min_v, sor_omega[k] * 0.7);
                            sor_n_damped++;
                        } else if (decreasing) {
                            // Monotone convergence: cautiously recover omega toward 1.0.
                            sor_omega[k] = std::min(1.0, sor_omega[k] * 1.05);
                        }
                        if (sor_omega[k] < sor_min_omega) sor_min_omega = sor_omega[k];
                        sor_prev_decreasing[k] = decreasing;
                        sor_prev_errRp[k]      = errRp_k;
                    }
                }
            }

            // WU-B: compute pct_change (max relative shift in cell mass since last check).
            double pct_change = 0.0;
            if (W_total > 0.0) {
                for (int c = 0; c < ct.M_cell; c++) {
                    double rel = std::fabs(X[c] - X_prev[c]) / std::max(X_prev[c], 1e-12);
                    if (rel > pct_change) pct_change = rel;
                }
            }

            // WU-C: l1_weight = Σ|ΔX| / W_input (normalized absolute shift in cell mass).
            double l1_weight = 0.0;
            for (int c = 0; c < ct.M_cell; c++)
                l1_weight += std::fabs(X[c] - X_prev[c]);
            if (ct.W_input > 0.0) l1_weight /= ct.W_input;

            // WU-B/C extra metrics gate: skip O(K*M_cell) passes when the active
            // stopping criterion does not need them. Always compute on the final
            // iteration (iter_in_lvl == budget_lvl) so result struct is fully
            // populated on exit. grake_norm moved inside this gate (was unconditional,
            // causing 9x slowdown on K=20/M_cell=1M with MAX_ERR+IMPROVEMENT default).
            const lbw::CalibMetric metric = st.convergence_cfg.metric;
            const auto& cfg_m = st.convergence_cfg;
            const bool need_extra_metrics =
                (metric == lbw::CalibMetric::MEAN_ERR    ||
                 metric == lbw::CalibMetric::KL          ||
                 metric == lbw::CalibMetric::CHI2        ||
                 metric == lbw::CalibMetric::GRAKE_NORM  ||
                 metric == lbw::CalibMetric::MARGINAL_KL ||
                 iter_in_lvl == budget_lvl);
            // Note: MARGINAL_KL is in need_extra_metrics to ensure grake_norm, kl, etc.
            // are populated in the result struct at convergence (marg_kl itself is already
            // computed in the errRp loop above at zero extra cost, but the full metrics
            // pass is needed to populate the other result fields for diagnostics).
            double grake_norm = 0.0;

            constexpr double kMetricEps = 1e-10;
            constexpr double kChi2Floor = 1.0;
            double mean_err_sum = 0.0;
            double kl_max       = 0.0;
            double chi2_total   = 0.0;
            if (need_extra_metrics && W_total > 0.0) {
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
                        // grake_norm computed here (reuses S_lin, no extra O(K*M_cell) pass)
                        double nm = std::fabs(obs - exp_val) / (1.0 + std::fabs(exp_val));
                        if (nm > grake_norm) grake_norm = nm;
                    }
                    mean_err_sum += max_k;
                    if (kl_k > kl_max) kl_max = kl_k;
                }
                // WU-E / g4oj: BLOCK 2 — best-iterate for non-MAX_ERR, non-MARGINAL_KL metrics.
                // All metric values (mean_err_sum, kl_max, chi2_total) are valid here.
                // MARGINAL_KL is handled by BLOCK 1b (marg_kl already computed in errRp loop).
                if (st.convergence_cfg.metric != lbw::CalibMetric::MAX_ERR &&
                    st.convergence_cfg.metric != lbw::CalibMetric::MARGINAL_KL) {
                    const double mean_err_blk2 = (st.K > 0)
                        ? (mean_err_sum / static_cast<double>(st.K)) : 0.0;
                    const double curr_best = lbw::select_metric(
                        st.convergence_cfg.metric,
                        errRp, mean_err_blk2, kl_max, chi2_total, grake_norm, l1_weight);
                    if (std::isfinite(curr_best) && curr_best < best_metric_seen) {
                        // === BEST-ITER UPDATE ===
                        best_metric_seen    = curr_best;
                        best_iter_val       = iter;
                        best_objective_seen = compute_weight_kl();
                        for (int c = 0; c < ct.M_cell; c++)
                            W_best[c] = (X_init[c] > 0.0) ? X[c] / X_init[c] : 0.0;
                        // === END BEST-ITER UPDATE ===
                    }
                }
            }
            double mean_err = (st.K > 0) ? (mean_err_sum / static_cast<double>(st.K)) : 0.0;

            // Store metrics in result struct (unconditional — intermediate checks
            // store 0 for gated metrics; final iter always populates all fields).
            res.l1_weight_change = l1_weight;    // WU-C: real Σ|ΔX|/W_input (replaces pct_change stub)
            res.grake_norm       = grake_norm;   // WU-C: max_kj normalized margin residual
            res.mean_error       = mean_err;
            res.kl               = kl_max;
            res.chi2             = chi2_total;

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

            // WU-C: metric + rule dispatch.
            // Select the active metric value, apply the stopping rule, then combine
            // with the optional secondary absolute threshold via stop_when semantics.
            // Intermediate homotopy levels additionally allow warm-jump when errRp < tol_lvl.
            {
                const auto& cfg = st.convergence_cfg;

                // Step 1: select active metric value (shared helper).
                // marg_kl is passed as the 7th arg (marginal_kl slot); valid for all metrics
                // because it is computed unconditionally in the errRp loop above.
                const double curr_metric = lbw::select_metric(
                    cfg.metric, errRp, mean_err, kl_max, chi2_total, grake_norm, l1_weight,
                    marg_kl);

                // Step 2: apply stopping rule. apply_rule(prev&) updates prev_metric_for_rule
                // in-place (sets prev=curr). This is the ONLY update inside the loop body.
                const bool converged_primary = lbw::apply_rule(
                    cfg.rule, curr_metric, prev_metric_for_rule, cfg.pct_tol);

                // Step 3: secondary (or sole) absolute threshold on the active metric.
                // absolute_tol applies to curr_metric (the selected metric), not always errRp.
                bool converged_abs = (cfg.absolute_tol > 0.0) && (curr_metric < cfg.absolute_tol);

                bool converged = false;
                if (cfg.absolute_tol > 0.0) {
                    converged = (cfg.stop_when == lbw::CalibStopWhen::ALL)
                                ? (converged_primary && converged_abs)
                                : (converged_primary || converged_abs);
                } else {
                    converged = converged_primary;
                }

                // Intermediate homotopy levels: allow warm-jump when errRp is loose.
                if (!converged && lvl < N_levels - 1) {
                    converged = (errRp < tol_lvl);
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

    // PCT stall detection: pct_change < pct_tol (PCT converged) but max_error >> pct_tol
    // signals infeasible problem. Threshold 10x: well-posed problems have
    // errRp/pct_change ratio 1-5x; infeasible stalls show 100x+; 10x separates them.
    // Warning only — status unchanged for backward compatibility.
    {
        const auto& cfg = st.convergence_cfg;
        if (cfg.pct_tol > 0.0 &&
            cfg.metric != CalibMetric::L1_WEIGHT &&
            res.max_error > 10.0 * cfg.pct_tol &&
            st.log_fn != nullptr) {
            char msg[256];
            std::snprintf(msg, sizeof(msg),
                "PCT convergence stall: pct_change < %.3g but max_error=%.3g "
                "(%.0fx pct_tol). Possible contradictory or infeasible targets.",
                cfg.pct_tol, res.max_error,
                res.max_error / cfg.pct_tol);
            st.log_fn(msg, st.log_ctx);
        }
    }

    // WU-C: populate convergence diagnostics at solver exit.
    res.convergence_metric             = static_cast<int>(st.convergence_cfg.metric);
    res.convergence_rule               = static_cast<int>(st.convergence_cfg.rule);
    res.convergence_tol                = st.convergence_cfg.pct_tol;
    res.convergence_iter               = (res.status == RK_OK) ? res.iterations : -1;
    res.convergence_solver_objective   = best_objective_seen;
    res.convergence_minimized_metric   = static_cast<int>(st.convergence_cfg.metric);

    // WU-E: expand W_best (cell-level snapshot) to obs-level best_weights.
    // Rule: scalar mult of initial obs weight by cell multiplier, then sum-normalize to n.
    // NO water-fill, NO bounds-clamping — this is a mid-loop snapshot.
    // If best_metric_seen == ∞ (solver exited before first check), best_weights is all zeros.
    res.best_error = best_metric_seen;
    res.best_iter  = best_iter_val;
    if (std::isfinite(best_metric_seen)) {
        std::vector<double> best_weights_obs(st.n);
        for (int i = 0; i < st.n; i++) {
            best_weights_obs[i] = st.weights[i] * W_best[ct.cell_of[i]];
        }
        // Sum-normalize to n (scalar mult + sum-normalize only, per spec).
        double s = 0.0;
        for (int i = 0; i < st.n; i++) s += best_weights_obs[i];
        if (s > 0.0) {
            const double scale = static_cast<double>(st.n) / s;
            for (int i = 0; i < st.n; i++) best_weights_obs[i] *= scale;
        }
        res.best_weights = std::move(best_weights_obs);
    } else {
        res.best_weights.assign(st.n, 0.0);
    }

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
        // Clamp post-normalization violations: normalization (w *= n/total_w)
        // can push weights above max_weight when cells were clamped (W_total<n).
        // Cell mode previously only counted — now clamps too, mirroring unit mode.
        // n_bounds_clamped = n_bounds_violated: every violation is clamped.
        int violations = 0;
        for (int i = 0; i < st.n; i++) {
            if (st.weights[i] > st.max_weight) {
                st.weights[i] = st.max_weight;
                violations++;
            } else if (st.weights[i] < st.min_weight) {
                st.weights[i] = st.min_weight;
                violations++;
            }
        }
        res.n_bounds_violated = violations;
        res.n_bounds_clamped  = violations;   // every violation was clamped
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

    // WU-D: store SOR diagnostics in result.
    res.sor_min_omega = sor_min_omega;
    res.sor_n_damped  = sor_n_damped;

    write_trajectory_csv(probe_samples);
    return res;
}

} // namespace lbw
