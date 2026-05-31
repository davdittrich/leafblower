#include "lbw_config.h"
#include "lbw_math.hpp"
#include "oris.hpp"
#include "calib_dispatch.hpp"
#include "cell_table.hpp"
#include "leafblower.h"
#include "sraa.hpp"
#include "oris_internal.hpp"
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

// parse_trajectory_iters() and write_trajectory_csv() moved to
// oris_trajectory.cpp (declared in oris_internal.hpp) — uu8r.1.

// Log-space algBCD at C=0 — see design spec §2.3, §8 for math.
// lf[k][j]: log Sinkhorn factor (per margin k, category j)
// W[c]:    capacity multiplier per cell (linear-space; bounded in [L_c/X_tilde, U_c/X_tilde])
// X_tilde[c] = X_init[c] * exp(sum_k lf[k][g_k(c)])
// X[c] = clamp(X_tilde[c] * W[c], L_c, U_c) — one capacity BCD block per outer iter
//
// Log-sum-exp stabilization on S_kj prevents overflow when partial log-sums
// approach log(DBL_MAX) ≈ 709.

// SRAA-m helpers — file-local. lf is the flat per-margin log-factor vector
// of size cat_offset[K] = Σ_k (cat_counts[k]+1). These helpers are O(M_cell)
// and pre-allocate nothing; all destination buffers are owned by oris_solve.
//
// pack_lf: copy lf -> dst (size n_cats_total_with_na). NA slots (j == cat_counts[k])
// are inert zeros; they participate in the SRAA linear system with ΔX=ΔR=0.
static inline void pack_lf(const std::vector<double>& lf,
                           std::vector<double>& dst) {
    std::copy(lf.begin(), lf.end(), dst.begin());
}

// unpack_lf: from a flat lf iterate, rebuild the derived state required by
// apply_single_margin_linear for the next sweep:
//   1) lf      <- src
//   2) f_lin[i] = exp(lf[i])     for i in 0..total_cats-1
//   3) cell_lf[c] = Σ_k lf[cat_offset[k] + g_per_cell[k][c]]   (g<0 skipped)
//   4) X_cur[c]   = X_init[c] * exp(cell_lf[c])
//
// inv_f_old_lin is NOT set here; apply_single_margin_linear recomputes it
// at the head of each margin call.
static inline void unpack_lf(const std::vector<double>& src,
                             std::vector<double>& lf,
                             std::vector<double>& f_lin,
                             std::vector<double>& cell_lf,
                             std::vector<double>& X_cur,
                             const lbw::CellTable& ct,
                             const std::vector<double>& X_init,
                             const std::vector<double>& log_X_init,
                             int K,
                             const std::vector<int>& cat_offset,
                             double& cell_lf_hwm) {
    const int total_cats = cat_offset[K];
    for (int i = 0; i < total_cats; i++) {
        lf[i]    = src[i];
        f_lin[i] = std::exp(src[i]);
    }
    const int M = ct.M_cell;
    std::fill(cell_lf.begin(), cell_lf.end(), 0.0);
    for (int k = 0; k < K; k++) {
        const int* gk  = ct.g_per_cell[k].data();
        const int  off = cat_offset[k];
        for (int c = 0; c < M; c++) {
            int g = gk[c];
            if (g < 0) continue;
            cell_lf[c] += src[off + g];
        }
    }
    double hwm = std::numeric_limits<double>::lowest();
    for (int c = 0; c < M; c++) {
        if (X_init[c] > 0.0) {
            X_cur[c] = X_init[c] * std::exp(cell_lf[c]);
            double v = log_X_init[c] + cell_lf[c];
            if (std::isfinite(v) && v > hwm) hwm = v;
        } else {
            X_cur[c] = 0.0;
        }
    }
    cell_lf_hwm = hwm;
}

// oris_finalize() moved to oris_finalize.cpp (declared in oris_internal.hpp) — uu8r.2.

ORISResult oris_solve(CalibState& st, std::vector<double>* lf_capture) {
    constexpr int    kErrCheckInterval = 10;
    constexpr double kEmptyBucketThreshold = 1e-15;
    constexpr double kLogClip = 700.0;  // exp(700) < DBL_MAX
    // Intermediate homotopy levels stop early once errRp drops below this loose
    // threshold; only the final level uses the user's tol_abs for full precision.
    constexpr double kHomotopyIntermediateTol = 1e-5;
    constexpr double kAlphaBeta          = 0.5;   // P2.1 stress→alpha mapping: alpha = 1/(1+β·stress)
    constexpr double kSorOscillationDamp = 0.7;   // SOR sign-flip: reduce omega by this factor
    constexpr double kSorRecoveryGrowth  = 1.05;  // SOR monotone: recover omega by this factor
    // kInfeasStallRatio defined once in oris_finalize() at line 148; this scope

    ORISResult res;
    res.base.status = RK_ERR_NOCONV;
    res.base.iterations = 0;
    res.base.max_error = 1.0;
    res.M_cell = 0;
    res.n_cap_active = 0;
    res.n_xcur_writes_per_iter_last = 0;
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
            res.base.status = RK_ERR_BADARG;
            return res;
        }
    }
    res.M_cell = ct.M_cell;

    // WU-2 dispatch: linear-space path when M_cell/n > 0.5 (i.e., compression <= 2x).
    // Env var LBW_ORIS_FORCE_PATH in {"linear", "log"} overrides for tests; always
    // compiled (no #ifdef) -- getenv cost is microseconds, amortized over the solve.
    constexpr double kLinearSpaceThreshold = 2.0;  // compression ratio cutoff
    bool use_linear = (static_cast<double>(st.n) /
                       static_cast<double>(std::max(ct.M_cell, 1))) <
                      kLinearSpaceThreshold;
    if (const char* force = std::getenv("LBW_ORIS_FORCE_PATH")) {
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

    // WI-1: best-iterate lf snapshot for optional lf_capture output.
    // Updated alongside W_best in the SRAA branch (mirrors lf_best = lf_flat).
    // Fallback: if best-iterate is never recorded (e.g. early-exit before any
    // improvement, or non-SRAA path that does not maintain a lf_best), the
    // RAII guard writes the final `lf` instead — caller documents this fallback.
    std::vector<double> lf_best_snap(total_cats, 0.0);
    bool                lf_best_snap_set = false;
    // RAII guard: writes the captured lf to *lf_capture on destruction so all
    // exit paths (returns, exceptions) are covered without manually wiring
    // every return statement in this long function.
    struct LfCaptureGuard {
        std::vector<double>*       out;
        const std::vector<double>* snap;
        const bool*                snap_set;
        const std::vector<double>* fallback;
        ~LfCaptureGuard() {
            if (out == nullptr) return;
            if (snap_set != nullptr && *snap_set && snap != nullptr) {
                *out = *snap;
            } else if (fallback != nullptr) {
                *out = *fallback;
            }
        }
    };
    LfCaptureGuard lf_capture_guard{lf_capture, &lf_best_snap,
                                    &lf_best_snap_set, &lf};

    // Per-cell capacity multiplier (linear-space).
    std::vector<double> W(ct.M_cell, 1.0);
    std::vector<double> log_W(ct.M_cell, 0.0);  // log_W[c] = log(W[c]); W init=1 → log=0
    std::vector<double> s_buf(ct.M_cell, 0.0);   // staging for bulk_exp_clipped
    std::vector<double> X_tilde;  // deferred: allocated at first log-path/fallback use
    std::vector<double> X(ct.M_cell);
    // T1.B: per-cell log-product shadow. cell_lf[c] = Σ_k lf[k][g_k(c)].
    // Reuses lf[] (line 135) as log(f_lin) in the linear path.
    // cell_lf also used by T2.A in the log path.
    std::vector<double> cell_lf(ct.M_cell, 0.0);
    // J2: Jacobi snapshot of cell_lf taken at outer-iter start (log path only).
    // Lambda reads from this when st.jacobi_log==true; else reads cell_lf directly.
    // Allocated only when jacobi_log is set (zero memory cost for default GS path).
    std::vector<double> cell_lf_frozen;
    if (st.jacobi_log) cell_lf_frozen.assign(ct.M_cell, 0.0);
    // High-water mark: max_c(log_X_init[c] + cell_lf[c]) ≈ max_c log(X_tilde[c]).
    // Monotone-nondecreasing between corrections — stale-high is intentional.
    double cell_lf_hwm = std::numeric_limits<double>::lowest();

    // ADMM dual variable for capacity constraint.
    // u[c] accumulates X_tilde - z violations; converges to 0 at fixed point.
    std::vector<double> u;  // allocated only for oris_soft (ADMM)
    if (st.use_admm_capacity) u.assign(ct.M_cell, 0.0);

    // ALM persistent state (oris_soft only). N_levels-dependent vars set below.
    const bool alm_active = st.use_admm_capacity && st.alm.capacity_mu > 0.0;
    const double capacity_mu_base    = st.alm.capacity_mu;
    double capacity_mu_adaptive      = capacity_mu_base;
    double eta_i_current             = 1.0;
    int    alm_violation_streak      = 0;
    std::vector<double> lambda_cell;
    if (alm_active) lambda_cell.assign(ct.M_cell, 0.0);

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
    //     mw=5, 500 iters, persistence disabled) confirmed ORIS reaches
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
        // Caller (c_api.cpp) sets st.oris_auto_selected=true when routing came
        // via AUTO; solver prepends [AUTO->ORIS] marker. Otherwise plain entry.
        const char* prefix = (st.oris_auto_selected ? "[AUTO->ORIS] " : "");
        std::snprintf(msg, sizeof(msg),
                      "%sORIS: n=%d K=%d M_cell=%d compression=%.1fx path=%s",
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
    const double kLinearOverflowTrip = std::min(
        std::pow(
            std::numeric_limits<double>::max() / (2.0 * max_X_init_val),
            1.0 / static_cast<double>(st.K)),
        1e15);
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
    double beta  = kAlphaBeta;
    double alpha = 1.0;
    // Test-only override (parallel to LBW_ORIS_FORCE_PATH): "on"|"off"|unset.
    // Always compiled; microsecond getenv cost. Enables falsifiable
    // min_alpha_seen assertion (spec §7, CTO B5).
    const char* force_damp = std::getenv("LBW_ORIS_FORCE_DAMPING");
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

    // X_prev tracks X[c] at the last convergence check for pct_change computation.
    // Initialized from X_init (uniform W[c]=1 at entry, X[c] = X_init[c]).
    std::vector<double> X_prev(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) X_prev[c] = X_init[c];

    // improvement/plateau rule — track active metric across kErrCheckInterval checks.
    // Initialized to +∞ so first check never triggers convergence.
    double prev_metric_for_rule = std::numeric_limits<double>::infinity();

    // WU-E / G8b: best-iterate tracking via BestIterTracker (replaces ad-hoc vars).
    // best.best_weights stores cell-level W ratio (X/X_init) at the best observed metric.
    BestIterTracker best;

    // SOR adaptive under-relaxation state (ORIS-only).
    // Per-margin omega[k] starts at omega_init (1.0 = no damping = fast path).
    // Adaptation: sign-flip in per-margin errRp trajectory → omega *= 0.7 (floor: omega_min).
    // Monotone decrease → omega *= 1.05, capped at 1.0 (recovery).
    // Adaptation is suppressed for sor_burnin iterations so early transient oscillation
    // (driven by infeas-streak damping) does not prematurely reduce omega.
    const bool sor_active     = st.sor_cfg.enabled;
    // When SRAA is active, SOR adaptation is allowed on plain steps
    // (aa_accepted=false) but suppressed on AA-accepted steps (their trajectory
    // is non-monotone from extrapolation; SOR would fight AA). The SRAA
    // while-loop updates sor_auto_v per-step via r.aa_accepted; non-SRAA
    // (accelerate=false) keeps the original masked semantics.
    const bool sor_base       = st.sor_cfg.auto_adapt;
    bool sor_auto_v           = sor_base && !st.accelerate;
    const double omega_init_v  = st.sor_cfg.omega_init;    // default 1.0
    const double omega_min_v   = st.sor_cfg.omega_min;     // default 0.3
    const double omega_max_v   = st.sor_cfg.omega_max;     // default 1.5; proven (0,2) (Thibault 2021)
    const double omega_fixed_v = st.sor_cfg.omega_fixed;   // -1.0 = use auto
    const int    sor_burnin_v  = st.sor_cfg.burnin;        // default 20
    // omega_mode_v: 0=heuristic(0.7/1.05), 1=fixed(omega_max), 2=spectral(Lehmann 2022)
    const int    omega_mode_v  = st.sor_cfg.omega_mode_id; // default 2

    // Spectral omega helpers — Lehmann et al. 2022 [lehmann2022overrelaxation].
    // theta2 = (||e_{k+1}|| / ||e_k||)^2 estimated from successive-residual ratio.
    // omega_opt = 2 / (1 + sqrt(1 - theta2)),  theta2 in [0,1).
    // Returns -1 from estimate_theta2 when non-informative (pre-burnin, NaN, ratio>=1).
    // omega_from_theta2 falls back to ceiling on non-informative theta2.
    // Both functions are call-overhead-free in the common case (stored residuals reused).
    //
    // Spectral mode uses kSorSpectralCeiling (1.99) — the Thibault 2021 strict-<2 boundary.
    // omega_max (default 1.5) is the ceiling for fixed (mode=1) and heuristic (mode=0) only.
    // These are independent: changing omega_max has NO effect on spectral mode.
    static constexpr double kSorSpectralCeiling = 1.99;  // Thibault 2021: convergence strict <2
    auto estimate_theta2 = [](double prev, double curr) -> double {
        if (prev <= 0.0 || curr <= 0.0 ||
            prev == std::numeric_limits<double>::infinity() ||
            !std::isfinite(prev) || !std::isfinite(curr))
            return -1.0;
        double ratio = curr / prev;
        if (ratio >= 1.0 || !std::isfinite(ratio)) return -1.0;
        return ratio * ratio;  // theta2 = (||e_{k+1}||/||e_k||)^2
    };
    // ceiling: kSorSpectralCeiling for spectral mode, omega_max_v for fixed/heuristic.
    auto omega_from_theta2 = [&](double theta2, double ceiling) -> double {
        if (theta2 < 0.0) return ceiling;  // non-informative -> fallback to ceiling
        double omega = 2.0 / (1.0 + std::sqrt(1.0 - theta2));
        return std::max(omega_min_v, std::min(ceiling, omega));
    };
    std::vector<double> sor_omega(st.K, omega_init_v);
    std::vector<double> sor_prev_errRp(st.K, std::numeric_limits<double>::infinity());
    std::vector<bool>   sor_prev_decreasing(st.K, false);
    double sor_min_omega = 1.0;
    int    sor_n_damped  = 0;

    // Scratch buffers for compute_weight_kl.
    std::vector<double> kl_ratio_buf(ct.M_cell);
    std::vector<double> kl_weight_buf(ct.M_cell);

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
    const bool tang_active = (st.eta_schedule.mode == EtaScheduleMode::TANG_DYNAMIC && N_levels > 1);
    const double k_start  = st.homotopy.start_factor;
    const double k_end    = st.homotopy.end_factor;
    const double p_budget = st.homotopy.budget_split_p;
    double budget_weight_sum = 0.0;
    for (int lvl = 0; lvl < N_levels; lvl++) {
        budget_weight_sum += std::pow(lvl + 1.0, p_budget);
    }
    const double alm_mu_base = st.alm.mu;
    int total_iters = 0;
    bool homotopy_break = false;
    bool absolute_tol_fired = false;  // set when absolute_tol triggers convergence
    // R3: kEmptyBucketThreshold and ct.W_input are loop-invariant; hoist products once.
    const double empty_threshold_abs  = kEmptyBucketThreshold * ct.W_input;
    const double log_empty_threshold  = std::log(empty_threshold_abs);

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
                st.alm.mu = eta_i * alm_mu_base;
            } else {
                beta = 0.5 * eta_i;
            }
        }

        // ALM (oris_soft) per-level state: scale capacity_mu, reset duals.
        if (alm_active) {
            if (tang_active) {
                const double scaled_frac = std::pow(frac, st.eta_schedule.schedule_power);
                eta_i_current = st.eta_schedule.eta_start *
                    std::pow(st.eta_schedule.eta_end / st.eta_schedule.eta_start, scaled_frac);
                res.eta_final = eta_i_current;
                st.alm.capacity_mu = eta_i_current * capacity_mu_adaptive;
            } else {
                eta_i_current = 1.0;
                st.alm.capacity_mu = capacity_mu_adaptive;
            }
            alm_violation_streak = 0;
            std::fill(lambda_cell.begin(), lambda_cell.end(), 0.0);
        }

        // Recompute U_cell for this level's current_max_weight.
        double hi = std::isfinite(current_max_weight) ? current_max_weight : kUnboundedSentinel;
        for (int c = 0; c < ct.M_cell; c++) {
            U_cell[c] = hi * ct.n_per_cell[c];
        }

        res.homotopy_levels_used  = lvl + 1;
        res.homotopy_final_factor = factor;

        bool level_converged = false;

        // reset X_prev at the start of each homotopy level so that
        // pct_change measures iteration-to-iteration shift within a level, not
        // cross-level drift from the previous level's final X.
        // B11: At lvl=0, X is still all-zeros (not yet populated); the X_prev
        // initialization from X_init above is the correct baseline. Only reset
        // for subsequent levels where X holds the warm-start from the prior level.
        if (lvl > 0) {
            for (int c = 0; c < ct.M_cell; c++) X_prev[c] = X[c];
        }
        // reset improvement/plateau baseline at each homotopy level.
        prev_metric_for_rule = std::numeric_limits<double>::infinity();

        // ────────────────────────────────────────────────────────────────────
        // SRAA-m linear-path acceleration — outer-loop replacement.
        // The lambdas below are hoisted out of the inner for-loop so f_eval_lf
        // can be invoked from the SRAA while-loop branch.
        // apply_single_margin_linear captures alpha (declared at line 319,
        // outside the for-loop) and writes through it via compute_alpha; it
        // does NOT reference iter_in_lvl or iter directly.
        // ────────────────────────────────────────────────────────────────────
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
            // effective omega for this margin.
            // sor_active==false → eff_omega=1.0 (fast path, no pow()).
            // sor_active && !sor_auto_v && omega_fixed_v>0 → fixed omega.
            // sor_active && sor_auto_v → per-margin adaptive omega[k].
            double eff_omega;
            if (!sor_active) {
                eff_omega = 1.0;
            } else if (!sor_auto_v && omega_fixed_v > 0.0) {
                eff_omega = omega_fixed_v;
            } else {
                eff_omega = sor_omega[k];
            }
            for (int j = 0; j < nj; j++) {
                if (!(S_lin[j] >= empty_threshold_abs) ||
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
                    // SOR under-relaxation in linear space.
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
            // X_cur scatter-multiply: contiguous write, gather read — SIMD-safe.
            // Guards converted to blend to avoid branch in vectorized loop.
#pragma omp simd
            for (int c = 0; c < ct.M_cell; c++) {
                const int j = gk[c];
                if (j >= 0 && j < nj && X_init[c] > 0.0)
                    X_cur[c] *= rescale_lin[j];
            }
            return false;
        };

        // ────────────────────────────────────────────────────────────────────
        // SRAA-m: apply_single_margin_log hoisted to homotopy scope so
        // f_eval_lf can dispatch it when use_linear=false. Body unchanged.
        // ────────────────────────────────────────────────────────────────────
        auto apply_single_margin_log = [&](int k) -> bool {
            // effective omega for this margin (log-space path).
            double eff_omega_log;
            if (!sor_active) {
                eff_omega_log = 1.0;
            } else if (!sor_auto_v && omega_fixed_v > 0.0) {
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
                    // pdzx: O(K) inner loop replaced by cell_lf shadow + 1 subtraction.
                    // cell_lf[c] = Σ_m lf[cat_offset[m] + g_per_cell[m][c]] (g>=0 only).
                    // We need Σ_{m!=k} that sum, so subtract margin-k's contribution
                    // iff g_k(c) >= 0 (mirrors cell_lf accumulation guard).
                    double s = log_X_init[c] + log_W[c] + cell_lf[c];
                    int gk = ct.g_per_cell[k][c];
                    if (gk >= 0) s -= lf[cat_offset[k] + gk];
                    lv[r] = s;
                    if (s > lv_max) lv_max = s;
                }
                if (!std::isfinite(lv_max)) { record_empty(k, j); continue; }
                double sum = 0.0;
                for (size_t r = 0; r < lv.size(); r++) {
                    if (std::isfinite(lv[r])) sum += std::exp(lv[r] - lv_max);
                }
                double log_S_kj = lv_max + std::log(sum);
                if (!std::isfinite(log_S_kj) || log_S_kj < log_empty_threshold) {
                    record_empty(k, j);
                    continue;
                }
                record_nonempty(k, j);
                double log_target = std::log(st.targets[k][j] * ct.W_input);
                // in log-space, SOR applies as a fractional step:
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
                    // J2: Jacobi gates the scattered cell_lf writes; rebuild fires after sweep.
                    // lf[off+j] write above is UNGATED — log factors must propagate for rebuild.
                    double delta = lf_new - lf_old;
                    if (!st.jacobi_log && std::fabs(delta) > 1e-12) {
                        for (int c : cells_by_margin_cat[cat_offset[k] + j])
                            cell_lf[c] += delta;
                    }
                }
            }
            return false;  // log path does not trip overflow mid-sweep
        };

        // ────────────────────────────────────────────────────────────────────
        // SRAA-m fixed-point map for the linear path of this homotopy level.
        // Hoisted out of the inner for-loop so the SRAA while-loop can invoke
        // it. Captures lf, f_lin, cell_lf, X_cur, ct, X_init, st.K, cat_offset,
        // apply_single_margin_linear, S_lin by [&].
        //
        // Receives a flat lf iterate (size cat_offset[st.K]); rebuilds derived
        // state via unpack_lf; runs one round-robin BCD pass over K margins;
        // packs updated lf back into `flat`; returns errRp on success or +inf
        // on overflow. On overflow, lf is mid-sweep (partial); SRAA safeguard
        // rejects via err_AA=inf>err_plain and reverts via swap with F_cur.
        // ────────────────────────────────────────────────────────────────────
        auto f_eval_lf = [&](std::vector<double>& flat) -> double {
            unpack_lf(flat, lf, f_lin, cell_lf, X_cur, ct, X_init, log_X_init,
                      st.K, cat_offset, cell_lf_hwm);

            bool overflow = false;
            if (use_linear) {
                for (int k = 0; k < st.K && !overflow; k++) {
                    if (apply_single_margin_linear(k)) overflow = true;
                }
            } else {
                // Log-path setup: seed X_tilde = X_cur, rebuild log_W from W.
                if (X_tilde.size() != static_cast<size_t>(ct.M_cell)) {
                    X_tilde.assign(ct.M_cell, 0.0);
                }
                std::copy(X_cur.begin(), X_cur.end(), X_tilde.begin());
                lbw::bulk_log(W.data(), log_W.data(), ct.M_cell);

                for (int k = 0; k < st.K && !overflow; k++) {
                    if (apply_single_margin_log(k)) overflow = true;  // always false for log path
                }
                // J2: Jacobi rebuild after K-margin sweep in SRAA f_eval_lf.
                if (st.jacobi_log) {
                    std::fill(cell_lf.begin(), cell_lf.end(), 0.0);
                    for (int k = 0; k < st.K; k++) {
                        const int off = cat_offset[k];
                        const int* gk = ct.g_per_cell[k].data();
                        for (int c = 0; c < ct.M_cell; c++) {
                            int g = gk[c];
                            if (g >= 0) cell_lf[c] += lf[off + g];
                        }
                    }
                }
                // Refresh X_tilde from updated cell_lf.
                for (int c = 0; c < ct.M_cell; c++) {
                    X_tilde[c] = (X_init[c] > 0.0) ? X_init[c] * std::exp(cell_lf[c]) : 0.0;
                }
            }
            // Always pack — on overflow, lf is partially updated; packing
            // preserves the SRAA invariant that `flat` reflects the current
            // iterate after the call.
            pack_lf(lf, flat);
            if (overflow) return std::numeric_limits<double>::infinity();

            const std::vector<double>& X_eval = use_linear ? X_cur : X_tilde;
            double W_total = 0.0;
            for (int c = 0; c < ct.M_cell; c++) W_total += X_eval[c];
            if (!(W_total > 0.0)) return std::numeric_limits<double>::infinity();
            double errRp = 0.0;
            for (int k = 0; k < st.K; k++) {
                const int nj = st.cat_counts[k];
                std::fill(S_lin.begin(), S_lin.begin() + nj, 0.0);
                const int* gk = ct.g_per_cell[k].data();
                for (int c = 0; c < ct.M_cell; c++) {
                    int j = gk[c];
                    if (j >= 0 && j < nj) S_lin[j] += X_eval[c];
                }
                for (int j = 0; j < nj; j++) {
                    double e = std::fabs(S_lin[j] / W_total - st.targets[k][j]);
                    if (e > errRp) errRp = e;
                }
            }
            return errRp;
        };

        // ────────────────────────────────────────────────────────────────────
        // SRAA-m outer-loop branch (B3): replaces the inner for-loop on the
        // linear path when accelerate=TRUE. The for-loop below is gated on
        // !sraa_active_lvl.
        // ────────────────────────────────────────────────────────────────────
        // SRAA-m is path-agnostic: f_eval_lf dispatches on use_linear internally.
        const bool sraa_active_lvl = st.accelerate;
        // Greedy downgrade — hoist to homotopy scope so it fires for both
        // linear and log SRAA paths (was inside for-loop, unreachable when
        // sraa_active_lvl=true after LL3 dropped the && use_linear gate).
        if (sraa_active_lvl && st.scheduler.mode == SchedulerMode::GREEDY) {
            if (st.verbose >= 1)
                st.log("[oris] greedy scheduler disabled under SRAA-m; using round_robin");
        }
        lbw::SRAAState oris_sraa;
        std::vector<double> lf_flat;
        std::vector<double> lf_best;
        const std::vector<double> dummy_L;
        const std::vector<double> dummy_U;
        int sraa_outer_stall_count = 0;
        double sraa_best_errRp     = std::numeric_limits<double>::infinity();
        double nat_metric_prev_sraa = std::numeric_limits<double>::infinity();
        int    nat_iter_prev_sraa   = -1;
        // Shared scratch — SRAA and flat-BCD are mutually exclusive; one buffer serves both.
        std::vector<double> w_ratio_scratch(ct.M_cell);
        // ════════════════════ SRAA-m accelerated path ════════════════════
        if (sraa_active_lvl) {
            oris_sraa.init(total_cats, lbw::kSRAAm);
            lf_flat.assign(total_cats, 0.0);
            lf_best.assign(total_cats, 0.0);
            pack_lf(lf, lf_flat);
            // Seed F_cur ONCE before the loop; sraa_step's Step 1 evaluates F
            // at F_cur. Subsequent steps carry F_cur forward via swap — do
            // NOT reset F_cur inside the loop.
            oris_sraa.F_cur = lf_flat;

            int  f_evals_used = 0;
            int  iter_sraa    = 0;   // SRAA-local iteration counter for SOR burnin
            bool converged    = false;
            while (f_evals_used < budget_lvl && !converged) {
                // Re-seed F_cur = current iterate before each step so that
                // step 1 of sraa_step evaluates f_eval(current_lf), not
                // f_eval(stale_F_cur). Without this reseed, after the first
                // swap F_cur holds old_lf and step 1 recomputes F(old_lf)
                // → R_k=0 → phantom convergence in 1 iteration.
                // (Mirrors raking.cpp:362: rk_sraa.F_cur = X.)
                oris_sraa.F_cur = lf_flat;
                auto r = lbw::sraa_step(f_eval_lf, lf_flat, dummy_L, dummy_U,
                                        oris_sraa, /*apply_clamp=*/false);
                f_evals_used += r.f_evals;
                iter_sraa    += r.f_evals;
                res.base.iterations = total_iters + f_evals_used;
                res.base.max_error  = r.err_rp;

                // B-narrow SOR coexistence: enable SOR adaptation on plain
                // (non-AA-accepted) SRAA steps, disable on AA-accepted steps
                // (their trajectory is non-monotone from extrapolation; SOR
                // would fight AA). Uses global errRp from r.err_rp as the
                // monotonicity proxy — coarser than per-margin errRp_k but
                // preserves the dampening effect needed to stabilize max_err.
                sor_auto_v = sor_base && !r.aa_accepted;
                if (sor_auto_v && sor_active && iter_sraa >= sor_burnin_v) {
                    const double curr_errRp = r.err_rp;
                    bool decreasing = (curr_errRp < sor_prev_errRp[0]);
                    bool sign_flip  = !decreasing && sor_prev_decreasing[0];
                    if (sign_flip) {
                        // Oscillation: damp regardless of mode.
                        for (int k = 0; k < st.K; k++) {
                            sor_omega[k] = std::max(omega_min_v,
                                sor_omega[k] * kSorOscillationDamp);
                            if (sor_omega[k] < sor_min_omega)
                                sor_min_omega = sor_omega[k];
                        }
                        sor_n_damped++;
                    } else if (decreasing) {
                        if (omega_mode_v == 2) {
                            // Spectral mode: Lehmann 2022 optimal omega from residual ratio.
                            // SRAA uses global errRp (index [0]) — single theta2 for all k.
                            // Ceiling = kSorSpectralCeiling (1.99); omega_max is fixed-mode only.
                            double theta2  = estimate_theta2(sor_prev_errRp[0], curr_errRp);
                            double omega_s = omega_from_theta2(theta2, kSorSpectralCeiling);
                            for (int k = 0; k < st.K; k++)
                                sor_omega[k] = omega_s;
                        } else if (omega_mode_v == 1) {
                            // Fixed mode: jump to omega_max on monotone convergence.
                            for (int k = 0; k < st.K; k++)
                                sor_omega[k] = omega_max_v;
                        } else {
                            // Heuristic mode (0): cautious multiplicative grow.
                            for (int k = 0; k < st.K; k++) {
                                sor_omega[k] = std::min(omega_max_v,
                                    sor_omega[k] * kSorRecoveryGrowth);
                            }
                        }
                    }
                    sor_prev_decreasing[0] = decreasing;
                    sor_prev_errRp[0]      = curr_errRp;
                }

                // Outer stall guard — revert lf to lf_best and clear AA history
                // when errRp degrades for kSRAAOuterStallWindow steps.
                // Uses errRp (not best.best_metric) since best-iterate now tracks
                // marginal_kl on the kErrCheckInterval gate, not errRp per-step.
                sraa_best_errRp = std::min(sraa_best_errRp, r.err_rp);
                if (r.err_rp > sraa_best_errRp * (1.0 + lbw::kSRAAOuterSlack)) {
                    if (++sraa_outer_stall_count >= lbw::kSRAAOuterStallWindow) {
                        lf_flat = lf_best;
                        unpack_lf(lf_flat, lf, f_lin, cell_lf, X_cur, ct, X_init,
                                  log_X_init, st.K, cat_offset, cell_lf_hwm);
                        oris_sraa.clear();
                        oris_sraa.F_cur = lf_flat;
                        sraa_outer_stall_count = 0;
                    }
                } else {
                    sraa_outer_stall_count = 0;
                }

                // Convergence check — mirrors the for-loop's kErrCheckInterval
                // gating. f_eval_lf returns errRp only; full metrics (marginal_kl,
                // kl, chi2) are computed here at check-interval boundaries so the
                // configured convergence metric (e.g. MARGINAL_KL for oris) is
                // available. Between check intervals the convergence flag is false.
                if (f_evals_used == 1 ||
                    f_evals_used % kErrCheckInterval == 0 ||
                    f_evals_used >= budget_lvl) {
                    double W_total = 0.0;
                    for (int c = 0; c < ct.M_cell; c++) W_total += X_cur[c];
                    if (W_total > 0.0) {
                        auto cm = lbw::compute_cell_metrics(st, ct, X_cur, W_total,
                                                            S_lin);
                        // MARGINAL_KL is not stored in CellMetrics.marginal_kl
                        // (returns 0.0 → spurious "trivially converged" on first
                        // check). CellMetrics.kl IS Σ_k T_kj*log(T_kj/S_kj) =
                        // the marginal calibration KL. Remap for the struct call.
                        lbw::CalibConvergenceCfg sraa_cfg = st.convergence_cfg;
                        if (sraa_cfg.metric == lbw::CalibMetric::MARGINAL_KL)
                            sraa_cfg.metric = lbw::CalibMetric::KL;
                        converged = lbw::check_convergence(sraa_cfg, cm,
                                                           prev_metric_for_rule,
                                                           st.tol_abs);
                        // jd2f: intermediate homotopy levels allow warm-jump when
                        // errRp is loose, mirroring the non-SRAA branch
                        // (converged = errRp < tol_lvl at L2038). Without it the
                        // SRAA path waited for full convergence at every level,
                        // slowing homotopy progression. Performance only — the
                        // final level (lvl == N_levels-1) still uses full tol.
                        if (!converged && lvl < N_levels - 1)
                            converged = (cm.errRp < tol_lvl);
                        // Natural-metric best-iterate update (inside W_total > 0 block).
                        // lf_best captured at check intervals only — correct: marginal_kl
                        // unavailable between checks. Fewer compute_weight_kl calls than errRp path.
                        const double nat_metric = lbw::select_metric(sraa_cfg.metric, cm);
                        if (res.base.metric_first_check == std::numeric_limits<double>::infinity())
                            res.base.metric_first_check = nat_metric;
                        if (std::isfinite(nat_metric_prev_sraa)) {
                            res.base.metric_prev_check = nat_metric_prev_sraa;
                            res.base.prev_check_iter   = nat_iter_prev_sraa;
                        }
                        nat_metric_prev_sraa = nat_metric;
                        nat_iter_prev_sraa   = res.base.iterations;
                        if (std::isfinite(nat_metric) && nat_metric < best.best_metric) {
                            auto& w_ratio = w_ratio_scratch;
                            for (int c = 0; c < ct.M_cell; c++)
                                w_ratio[c] = (X_init[c] > 0.0) ? X_cur[c] / X_init[c] : 0.0;
                            best.update(nat_metric,
                                        lbw::compute_weight_kl(X_cur, X_init, ct.M_cell, st.n,
                                                               kl_ratio_buf.data(), kl_weight_buf.data()),
                                        res.base.iterations, w_ratio);
                            lf_best          = lf_flat;
                            lf_best_snap     = lf_flat;
                            lf_best_snap_set = true;
                        }
                    }
                }

                if (st.verbose >= 1) {
                    char msg[256];
                    std::snprintf(msg, sizeof(msg),
                                  "ORIS[sraa] f_evals=%d errRp=%.3e aa=%d",
                                  f_evals_used, r.err_rp,
                                  (int)r.aa_accepted);
                    st.log(msg);
                }
            }
            res.aa_accepted_count = oris_sraa.aa_accepted_count;
            total_iters += f_evals_used;
            if (converged) {
                level_converged = true;
                // Set convergence status — mirrors the mark_converged() call
                // in the non-SRAA for-loop. Without this, status stays at
                // RK_ERR_NOCONV and is later mis-classified as RK_ERR_BUDGET.
                if (lvl == N_levels - 1) {
                    lbw::mark_converged(res, st.convergence_cfg, res.base.iterations);
                }
            }
            // Clamp cell masses to [L_cell, U_cell] before syncing to X.
            // SRAA lf extrapolation is unconstrained and can overshoot cell
            // capacity bounds; without clamping, per-obs weights explode
            // (e.g. wmax=75 on stepstone with max_weight=5).
            //
            // vtjf: mass-preserving clamp. A plain clamp sheds the over-cap excess,
            // so Σ_c X < n and the next homotopy level's warm start begins below
            // the SRAA fixed point. Instead pick a scalar r so that
            //   Σ_c clamp(r·src[c], L, U) = n.
            // f(r) is monotone non-decreasing, so bisect. If n is outside
            // [Σ L_cell, Σ U_cell] (mass-saturated/infeasible) r saturates and X
            // pins at the reachable bound (best effort) — the honest outcome.
            {
                const std::vector<double>& src_cells = use_linear ? X_cur
                    : (!X_tilde.empty() ? X_tilde : X_cur);
                const double target = static_cast<double>(st.n);
                auto clamped_sum = [&](double r) {
                    double s = 0.0;
                    for (int c = 0; c < ct.M_cell; c++)
                        s += std::clamp(r * src_cells[c], L_cell[c], U_cell[c]);
                    return s;
                };
                double r = 1.0;
                if (clamped_sum(0.0) >= target) {
                    r = 0.0;                       // even all-at-L overshoots n
                } else {
                    double r_lo = 0.0, r_hi = 1.0;
                    int grow = 0;
                    while (clamped_sum(r_hi) < target && grow < 200) { r_hi *= 2.0; ++grow; }
                    if (clamped_sum(r_hi) < target) {
                        r = r_hi;                  // mass-saturated: pin near U
                    } else {
                        for (int it = 0; it < 100; it++) {
                            r = 0.5 * (r_lo + r_hi);
                            const double s = clamped_sum(r);
                            if (std::abs(s - target) <= 1e-12 * target) break;
                            if (s < target) r_lo = r; else r_hi = r;
                        }
                    }
                }
                for (int c = 0; c < ct.M_cell; c++) {
                    X[c] = std::clamp(r * src_cells[c], L_cell[c], U_cell[c]);
                }
            }
        }

        double nat_metric_prev_nonavec = std::numeric_limits<double>::infinity();
        int    nat_iter_prev_nonavec   = -1;
        // T20 g9jm: hoist scratch vectors out of flat-BCD iter loop — eliminates K-sized alloc per iter.
        std::vector<double> per_margin_err(st.K);
        std::vector<double> per_k_errRp_cache(st.K, 0.0);
        bool per_k_errRp_valid = false;
        // ════════════════════ Non-accelerated flat BCD path ════════════════════
        if (!sraa_active_lvl) {
        for (int iter_in_lvl = 1; iter_in_lvl <= budget_lvl; iter_in_lvl++) {
        const int iter = total_iters + iter_in_lvl;
        per_k_errRp_valid = false;  // T20: reset stale-cache flag at iter start
        res.base.iterations = iter;
        alpha = compute_alpha();
        if (alpha < res.min_alpha_seen) res.min_alpha_seen = alpha;

        // J2: Jacobi freeze. Log path only; linear path uses live X_cur leave-one-out.
        if (st.jacobi_log && !use_linear) {
            std::copy(cell_lf.begin(), cell_lf.end(), cell_lf_frozen.begin());
        }

        // Margin sweep: branched by path.
        bool overflow_trip = false;

        // Do NOT read st.max_weight here — homotopy levels pass
        // current_max_weight indirectly via the shared U_cell already built.

        // Per-margin residual errRp_k = max_j |S_lin[j]/W_total - targets[k][j]|.
        // Uses the current X_cur (linear path) or rebuilt S via cells_by_margin_cat (log).
        // B12: return infinity when W_total<=0 so greedy scheduler treats this margin
        // as highest priority rather than silently signalling false perfect convergence.
        auto compute_margin_errRp_linear = [&](int k) -> double {
            const int nj = st.cat_counts[k];
            std::fill(S_lin.begin(), S_lin.begin() + nj, 0.0);
            const int* gk = ct.g_per_cell[k].data();
            double W_total = 0.0;
            for (int c = 0; c < ct.M_cell; c++) W_total += X_cur[c];
            if (W_total <= 0.0) return std::numeric_limits<double>::infinity();
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

        // B12: same sentinel contract as compute_margin_errRp_linear above.
        auto compute_margin_errRp_log = [&](int k) -> double {
            if (X_tilde.empty()) return std::numeric_limits<double>::infinity();
            double W_total = 0.0;
            for (int c = 0; c < ct.M_cell; c++) W_total += X_tilde[c];
            if (W_total <= 0.0) return std::numeric_limits<double>::infinity();
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

        // f_eval_lf is hoisted to the homotopy-level scope above (see B3 SRAA
        // wiring). Inside this for-loop it remains in scope but is unused on
        // the non-SRAA path.

        bool use_greedy = (st.scheduler.mode == SchedulerMode::GREEDY);
        if (st.accelerate && use_greedy) {
            use_greedy = false;
            res.sraa_demoted = true;
            if (st.verbose >= 1)
                st.log("[oris] greedy scheduler disabled under SRAA-m; using round_robin");
        }

        if (use_linear) {
            // WU-2 prefactored linear-space sweep (spec rev 5 §5).
            // Sequential (raking-style) access to X_cur for cache efficiency.
            // Each margin k: 2 sequential passes over M_cell cells = O(K*M_cell) total.
            // Per bucket (k,j): S_kj = sum X_cur[c]/f_kj accumulated in pass 1;
            //   X_cur[c] rescaled by f_new/f_old in pass 2. Identical to bucket loop.

            if (use_greedy) {
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
                // j <= cat_counts[k]: include NA bucket (cat_offset has +1 per margin).
                // Without NA shift, cells NA for a margin would have cell_lf decremented
                // by full shift but lf[k][NA] unchanged — invariant violated.
                for (int k = 0; k < st.K; k++) {
                    for (int j = 0; j <= st.cat_counts[k]; j++)
                        lf[cat_offset[k] + j] += lf_correction;
                    if (st.verbose >= 2) {
                        char lf_msg[128];
                        std::snprintf(lf_msg, sizeof(lf_msg),
                            "T1B lf_correction[%d]=%.6f iter=%d", k, lf_correction, iter);
                        st.log(lf_msg);
                    }
                }
                lbw::bulk_scaled_exp(1.0, lf.data(), f_lin.data(), total_cats);
                // B8: per-cell K_active scaling. cell_lf[c] sums lf entries only over
                // active (g_per_cell[k][c] >= 0) margins. After lf[k][·] += lf_correction
                // for all k (including NA bucket), cell_lf[c] decreased by exactly
                // K_active(c) * |lf_correction| = (K_active/K) * shift, NOT shift.
                // Likewise X_cur[c] = X_init[c] * Π_{active k} f_lin[k][g_k(c)] scales
                // by exp(K_active * lf_correction), not exp(-shift).
                for (int c = 0; c < ct.M_cell; c++) {
                    int K_active = 0;
                    for (int k = 0; k < st.K; k++)
                        if (ct.g_per_cell[k][c] >= 0) K_active++;
                    const double cell_corr = static_cast<double>(K_active) * lf_correction;
                    cell_lf[c] += cell_corr;
                    X_cur[c]   *= std::exp(cell_corr);
                }
                // Reset to lowest so hwm repopulates honestly from next positive delta,
                // avoiding spurious re-trigger when true max is exactly at threshold.
                cell_lf_hwm = std::numeric_limits<double>::lowest();
                if (st.use_admm_capacity) std::fill(u.begin(), u.end(), 0.0);
                if (alm_active) std::fill(lambda_cell.begin(), lambda_cell.end(), 0.0);
                if (st.verbose >= 2) {
                    char msg[128];
                    std::snprintf(msg, sizeof(msg), "ORIS T1.B renorm shift=%.2e", shift);
                    st.log(msg);
                }
            }

            if (overflow_trip && !linear_fallback_used) {
                // One-shot fallback: reset all solver state, switch to log-space,
                // restart outer loop from iter 0. State-clean list per spec rev 5 §5.
                linear_fallback_used = true;
                use_linear = false;
                // SRAA history from the linear path is stale after fallback.
                if (sraa_active_lvl) {
                    res.aa_accepted_count = oris_sraa.aa_accepted_count;
                    oris_sraa.clear();
                }
                std::fill(lf.begin(), lf.end(), 0.0);
                std::fill(cell_lf.begin(), cell_lf.end(), 0.0);
                cell_lf_hwm = std::numeric_limits<double>::lowest();
                if (st.use_admm_capacity) std::fill(u.begin(), u.end(), 0.0);
                if (alm_active) std::fill(lambda_cell.begin(), lambda_cell.end(), 0.0);
                std::fill(f_lin.begin(), f_lin.end(), 1.0);
                std::fill(X_cur.begin(), X_cur.end(), 0.0);
                std::fill(W.begin(), W.end(), 1.0);
                std::fill(log_W.begin(), log_W.end(), 0.0);
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
                // reset X_prev after fallback — X semantics changed (log-path).
                for (int c = 0; c < ct.M_cell; c++) X_prev[c] = X[c];
                // Reset best-iterate: pre-fallback snapshot from degenerate linear-space
                best.reset();
                sraa_best_errRp     = std::numeric_limits<double>::infinity();
                nat_metric_prev_sraa = std::numeric_limits<double>::infinity();
                nat_iter_prev_sraa   = -1;
                if (st.verbose >= 1) {
                    st.log("ORIS: linear-space overflow trip; fallback to log-space.");
                }
                continue;  // skip the post-sweep X_tilde / capacity / errRp blocks this round
            }
        } else {
            // Log-space sweep via helper. Greedy branch needs X_tilde for residual
            // probe; rebuild it cheaply (residual-only; no capacity side-effect).
            if (use_greedy) {
                // Rebuild X_tilde from current lf for residual scoring.
                if (X_tilde.empty()) X_tilde.assign(ct.M_cell, 0.0);
                // S3 Site A: vectorized X_tilde rebuild via bulk_exp_clipped.
#pragma omp simd
                for (int c = 0; c < ct.M_cell; c++)
                    s_buf[c] = (X_init[c] <= 0.0) ? -kLogClip : log_X_init[c] + cell_lf[c];
                lbw::bulk_exp_clipped(s_buf.data(), X_tilde.data(), ct.M_cell, kLogClip);
                for (int c = 0; c < ct.M_cell; c++)
                    if (X_init[c] <= 0.0) X_tilde[c] = 0.0;
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
                    // S3 Site B: vectorized X_tilde refresh via bulk_exp_clipped.
#pragma omp simd
                    for (int c = 0; c < ct.M_cell; c++)
                        s_buf[c] = (X_init[c] <= 0.0) ? -kLogClip : log_X_init[c] + cell_lf[c];
                    lbw::bulk_exp_clipped(s_buf.data(), X_tilde.data(), ct.M_cell, kLogClip);
                    for (int c = 0; c < ct.M_cell; c++)
                        if (X_init[c] <= 0.0) X_tilde[c] = 0.0;
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
            // J2: Jacobi rebuild. O(K * M_cell) sequential; cache-friendly.
            // Rebuilds cell_lf from current lf[] and refreshes cell_lf_hwm.
            if (st.jacobi_log) {
                std::fill(cell_lf.begin(), cell_lf.end(), 0.0);
                for (int k = 0; k < st.K; k++) {
                    const int off = cat_offset[k];
                    const int* gk = ct.g_per_cell[k].data();
                    for (int c = 0; c < ct.M_cell; c++) {
                        int g = gk[c];
                        if (g >= 0) cell_lf[c] += lf[off + g];
                    }
                }
                // Refresh hwm from live cell_lf (one O(M_cell) pass).
                cell_lf_hwm = std::numeric_limits<double>::lowest();
                for (int c = 0; c < ct.M_cell; c++) {
                    double val = cell_lf[c] + log_X_init[c];
                    if (std::isfinite(val) && val > cell_lf_hwm) cell_lf_hwm = val;
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
                if (alm_active && X_tilde_c > 0.0) {
                    // ALM: linearized Newton step, rho = mu*X_tilde balances KL Hessian.
                    // X = X_tilde * (1 - lambda + mu*z) / (1 + rho); lambda += mu*(X-z).
                    //
                    // Derivation (ORIS uses the UN-normalized KL / I-divergence
                    // D = X*log(X/X_tilde) - X + X_tilde, so dD/dX = log(X/X_tilde),
                    // NOT log(X/X_tilde)+1). Stationarity of
                    //   D + lambda*X + (mu/2)(X-z)^2  is  log(X/X_tilde)+lambda+mu(X-z)=0.
                    // Linearize log(X/X_tilde) ~= X/X_tilde - 1:
                    //   (X/X_tilde - 1) + lambda + mu(X-z) = 0
                    //   X(1 + mu*X_tilde) = X_tilde*(1 - lambda + mu*z)
                    //   X = X_tilde*(1 - lambda + mu*z)/(1 + rho).
                    // This is exactly the line below. leafblower-7emq proposed an
                    // extra "-rho" in the numerator, but that derives from the
                    // NORMALIZED KL (dD/dX = log+1), which is not ORIS's generator.
                    // Verified: both forms share the fixed point X=clamp(X_tilde,L,U)
                    // and convergence rate is a wash (2-cell scipy-style sim). Do NOT
                    // add -rho here.
                    const double z   = std::clamp(X_tilde_c, L_cell[c], U_cell[c]);
                    const double rho = st.alm.capacity_mu * X_tilde_c;
                    double X_alm = X_tilde_c * (1.0 - lambda_cell[c] + st.alm.capacity_mu * z) / (1.0 + rho);
                    if (!std::isfinite(X_alm) || X_alm <= 0.0) X_alm = z;  // NaN guard
                    X[c] = X_alm; W[c] = X_alm / X_tilde_c; X_cur[c] = X_alm;
                    const double lambda_cap = 10.0 * capacity_mu_base * st.max_weight;
                    lambda_cell[c] += st.alm.capacity_mu * (X_alm - z);
                    lambda_cell[c]  = std::clamp(lambda_cell[c], -lambda_cap, lambda_cap);
                } else {
                    // Hard clamp (oris default or X_tilde_c <= 0).
                    double xc = std::clamp(X_tilde_c, L_cell[c], U_cell[c]);
                    X[c] = xc; W[c] = xc / X_tilde_c; X_cur[c] = xc;
                }
                res.n_xcur_writes_per_iter_last++;
                if (X[c] != X_tilde_c) n_cap++;
            }
            res.n_cap_active = n_cap;
            if (overflow_detected) {
                // Full state reset on mid-loop break; partial writes to W/X/X_cur undone.
                std::fill(X_cur.begin(),   X_cur.end(),   0.0);
                std::fill(W.begin(),       W.end(),       1.0);
                std::fill(log_W.begin(),   log_W.end(),   0.0);
                std::fill(X.begin(),       X.end(),       0.0);
                std::fill(X_tilde.begin(), X_tilde.end(), 0.0);
                std::fill(lf.begin(),      lf.end(),      0.0);
                std::fill(cell_lf.begin(), cell_lf.end(), 0.0);
                cell_lf_hwm = std::numeric_limits<double>::lowest();
                std::fill(f_lin.begin(),   f_lin.end(),   1.0);
                std::fill(infeas_streak.begin(), infeas_streak.end(), 0);
                if (alm_active) std::fill(lambda_cell.begin(), lambda_cell.end(), 0.0);
                res.n_xcur_writes_per_iter_last = 0;
                use_linear = false;
                // SRAA history from linear path is stale after path flip.
                if (sraa_active_lvl && !lf_flat.empty()) {
                    res.aa_accepted_count = oris_sraa.aa_accepted_count;
                    oris_sraa.clear();
                    pack_lf(lf, lf_flat);
                    oris_sraa.F_cur = lf_flat;
                }
                linear_fallback_used = true;
                // reset X_prev after fallback — X semantics changed (log-path).
                for (int c = 0; c < ct.M_cell; c++) X_prev[c] = X[c];
                if (st.verbose >= 1) st.log("ORIS: linear-space overflow trip; fallback to log-space.");
                continue;
            }
        } else {
            // Log-path: X_tilde + capacity + X_cur unchanged from current implementation.
            bool overflow_detected = false;
            double max_log_X_tilde = -std::numeric_limits<double>::infinity();
            if (X_tilde.empty()) X_tilde.assign(ct.M_cell, 0.0);
            // S3 Site C: Pass 1 — SIMD, fill s_buf, track overflow/max.
            // continue replaced by ternary blend; reductions declared for vectorizer.
#pragma omp simd reduction(max:max_log_X_tilde) reduction(||:overflow_detected)
            for (int c = 0; c < ct.M_cell; c++) {
                // T2.A: single-stream exp via cell_lf (was K=20 DRAM streams)
                const bool active = (X_init[c] > 0.0);
                const double s = active ? (log_X_init[c] + cell_lf[c]) : -kLogClip;
                if (active && s > max_log_X_tilde) max_log_X_tilde = s;
                if (active && s > kLogClip && U_cell[c] >= 1e299) overflow_detected = true;
                s_buf[c] = s;
            }
            // Pass 2: vectorized exp with clipping.
            lbw::bulk_exp_clipped(s_buf.data(), X_tilde.data(), ct.M_cell, kLogClip);
            for (int c = 0; c < ct.M_cell; c++)
                if (X_init[c] <= 0.0) X_tilde[c] = 0.0;
            if (overflow_detected) {
                res.base.status = RK_ERR_NOCONV;
                res.base.max_error = std::numeric_limits<double>::infinity();
                if (st.verbose >= 2) {
                    char msg[256];
                    std::snprintf(msg, sizeof(msg),
                                  "ORIS: log-factor overflow (max_log_X_tilde=%.1f > 700) "
                                  "indicates ill-conditioning; try looser max_weight or "
                                  "tighter tol_abs, or method=raking.", max_log_X_tilde);
                    st.log(msg);
                }
                break;
            }
            // Capacity block: X[c] = clamp(X_tilde[c], L_c, U_c); W[c] updated for next iter.
            int n_cap = 0;
            for (int c = 0; c < ct.M_cell; c++) {
                if (alm_active && X_tilde[c] > 0.0) {
                    const double z   = std::clamp(X_tilde[c], L_cell[c], U_cell[c]);
                    const double rho = st.alm.capacity_mu * X_tilde[c];
                    double X_alm = X_tilde[c] * (1.0 - lambda_cell[c] + st.alm.capacity_mu * z) / (1.0 + rho);
                    if (!std::isfinite(X_alm) || X_alm <= 0.0) X_alm = z;
                    X[c] = X_alm;
                    W[c] = X_alm / X_tilde[c];
                    const double lambda_cap = 10.0 * capacity_mu_base * st.max_weight;
                    lambda_cell[c] += st.alm.capacity_mu * (X_alm - z);
                    lambda_cell[c]  = std::clamp(lambda_cell[c], -lambda_cap, lambda_cap);
                } else {
                    double xc = std::clamp(X_tilde[c], L_cell[c], U_cell[c]);
                    X[c] = xc;
                    if (X_tilde[c] > 0.0) {
                        W[c] = xc / X_tilde[c];
                    } else {
                        W[c] = 1.0;
                    }
                }
                if (W[c] != 1.0) n_cap++;
            }
            res.n_cap_active = n_cap;
            // S2: precompute log_W for next iteration's apply_single_margin_log calls.
            lbw::bulk_log(W.data(), log_W.data(), ct.M_cell);
        }

        // ALM adaptive mu growth: if max bound violation persists above tolerance,
        // grow mu geometrically (capped). Triggered after kAlmPersistenceThreshold
        // consecutive iterations with violation > tol_primal.
        if (alm_active) {
            constexpr int    kAlmPersistenceThreshold = 5;
            constexpr double kAlmGrowthFactor         = 2.0;
            constexpr double kAlmMaxScale             = 1000.0;
            const double mean_n_per_cell = static_cast<double>(st.n) / std::max(1, ct.M_cell);
            const double tol_primal      = 0.01 * st.max_weight * mean_n_per_cell;

            double max_violation = 0.0;
            for (int c = 0; c < ct.M_cell; c++) {
                const double v = std::max(X[c] - U_cell[c], L_cell[c] - X[c]);
                if (std::isfinite(v)) max_violation = std::max(max_violation, v);
            }

            if (max_violation > tol_primal) {
                if (++alm_violation_streak >= kAlmPersistenceThreshold &&
                    capacity_mu_adaptive < capacity_mu_base * kAlmMaxScale) {
                    capacity_mu_adaptive *= kAlmGrowthFactor;
                    st.alm.capacity_mu = eta_i_current * capacity_mu_adaptive;
                    // SRAA history stale after capacity_mu change — fixed-point shifted.
                    if (sraa_active_lvl && !lf_flat.empty()) {
                        res.aa_accepted_count = oris_sraa.aa_accepted_count;
                        oris_sraa.clear();
                        pack_lf(lf, lf_flat);
                        oris_sraa.F_cur = lf_flat;
                    }
                    res.alm_n_growth_events++;
                    alm_violation_streak = 0;
                    if (st.verbose >= 2) {
                        char msg[128];
                        std::snprintf(msg, sizeof(msg), "[oris_soft] mu growth: %.4e", capacity_mu_adaptive);
                        st.log(msg);
                    }
                }
            } else {
                alm_violation_streak = 0;
            }
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
                    double errRp_k = 0.0;                          // 773f.7
                    for (int j = 0; j < nj; j++) {
                        double e = std::fabs(S_lin[j] / W_total - st.targets[k][j]);
                        if (e > errRp) errRp = e;
                        if (e > errRp_k) errRp_k = e;             // 773f.7: capture per-margin
                    }
                    per_k_errRp_cache[k] = errRp_k;               // 773f.7
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
            per_k_errRp_valid = (use_linear && W_total > 0.0);    // 773f.7
            res.marginal_kl_at_iter = marg_kl;
            res.base.max_error = errRp;

            // w_ratio: alias to function-scope scratch (hoisted above sraa_active_lvl branch).
            // Populated lazily inside each gate only when a new best is found.
            auto& w_ratio = w_ratio_scratch;

            // WU-E / g4oj: BLOCK 1 — MAX_ERR best-iterate (errRp always valid here,
            // outside need_extra_metrics gate). Tracks min errRp when MAX_ERR is active.
            if (st.convergence_cfg.metric == lbw::CalibMetric::MAX_ERR) {
                if (errRp < best.best_metric) {
                    for (int c = 0; c < ct.M_cell; c++)
                        w_ratio[c] = (X_init[c] > 0.0) ? X[c] / X_init[c] : 0.0;
                    best.update(errRp,
                                lbw::compute_weight_kl(X, X_init, ct.M_cell, st.n, kl_ratio_buf.data(), kl_weight_buf.data()),
                                iter, w_ratio);
                }
            }
            // BLOCK 1b — MARGINAL_KL best-iterate (marg_kl always valid alongside errRp).
            // Tracks min marginal KL when MARGINAL_KL is active.
            if (st.convergence_cfg.metric == lbw::CalibMetric::MARGINAL_KL) {
                if (iter == 1 && res.base.metric_first_check == std::numeric_limits<double>::infinity())
                    res.base.metric_first_check = res.marginal_kl_at_iter;
                if (std::isfinite(nat_metric_prev_nonavec)) {
                    res.base.metric_prev_check = nat_metric_prev_nonavec;
                    res.base.prev_check_iter   = nat_iter_prev_nonavec;
                }
                nat_metric_prev_nonavec = res.marginal_kl_at_iter;
                nat_iter_prev_nonavec   = iter;
                if (res.marginal_kl_at_iter < best.best_metric) {
                    for (int c = 0; c < ct.M_cell; c++)
                        w_ratio[c] = (X_init[c] > 0.0) ? X[c] / X_init[c] : 0.0;
                    best.update(res.marginal_kl_at_iter,
                                lbw::compute_weight_kl(X, X_init, ct.M_cell, st.n, kl_ratio_buf.data(), kl_weight_buf.data()),
                                iter, w_ratio);
                }
            }

            // per-margin omega adaptation (both linear and log paths; auto mode only;
            // suppressed during burnin to let the infeas-streak damping settle first).
            if (sor_active && sor_auto_v && iter >= sor_burnin_v) {
                if (W_total > 0.0) {
                    for (int k = 0; k < st.K; k++) {
                        const int nj_k = st.cat_counts[k];
                        double errRp_k;
                        if (per_k_errRp_valid) {
                            // 773f.7: reuse per-margin errRp captured during convergence sweep.
                            errRp_k = per_k_errRp_cache[k];
                        } else {
                            // Log path: re-accumulate (unchanged behavior).
                            std::fill(S_lin.begin(), S_lin.begin() + nj_k, 0.0);
                            const int* gk_s = ct.g_per_cell[k].data();
                            for (int c = 0; c < ct.M_cell; c++) {
                                int j = gk_s[c];
                                if (j >= 0 && j < nj_k) S_lin[j] += X[c];
                            }
                            errRp_k = 0.0;
                            for (int j = 0; j < nj_k; j++) {
                                double e = std::fabs(S_lin[j] / W_total - st.targets[k][j]);
                                if (e > errRp_k) errRp_k = e;
                            }
                        }
                        bool decreasing = (errRp_k < sor_prev_errRp[k]);
                        bool sign_flip  = !decreasing && sor_prev_decreasing[k];
                        if (sign_flip) {
                            // Oscillation detected: damp regardless of mode.
                            sor_omega[k] = std::max(omega_min_v, sor_omega[k] * kSorOscillationDamp);
                            sor_n_damped++;
                        } else if (decreasing) {
                            if (omega_mode_v == 2) {
                                // Spectral mode: Lehmann 2022 per-margin optimal omega.
                                // Ceiling = kSorSpectralCeiling (1.99); omega_max is fixed-mode only.
                                double theta2_k = estimate_theta2(sor_prev_errRp[k], errRp_k);
                                sor_omega[k] = omega_from_theta2(theta2_k, kSorSpectralCeiling);
                            } else if (omega_mode_v == 1) {
                                // Fixed mode: jump to omega_max on monotone convergence.
                                sor_omega[k] = omega_max_v;
                            } else {
                                // Heuristic mode (0): cautious multiplicative grow.
                                sor_omega[k] = std::min(omega_max_v, sor_omega[k] * kSorRecoveryGrowth);
                            }
                        }
                        if (sor_omega[k] < sor_min_omega) sor_min_omega = sor_omega[k];
                        sor_prev_decreasing[k] = decreasing;
                        sor_prev_errRp[k]      = errRp_k;
                    }
                }
            }

            // Compute pct_change (max relative shift in cell mass since last check).
            double pct_change = 0.0;
            if (W_total > 0.0) {
                for (int c = 0; c < ct.M_cell; c++) {
                    double rel = std::fabs(X[c] - X_prev[c]) / std::max(X_prev[c], 1e-12);
                    if (rel > pct_change) pct_change = rel;
                }
            }

            // l1_weight = Σ|ΔX| / W_input (normalized absolute shift in cell mass).
            double l1_weight = 0.0;
            for (int c = 0; c < ct.M_cell; c++)
                l1_weight += std::fabs(X[c] - X_prev[c]);
            if (ct.W_input > 0.0) l1_weight /= ct.W_input;

            // extra metrics gate: skip O(K*M_cell) passes when the active
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
                 metric == lbw::CalibMetric::L1_WEIGHT   ||
                 iter_in_lvl == budget_lvl);
            // Note: MARGINAL_KL is in need_extra_metrics to ensure grake_norm, kl, etc.
            // are populated in the result struct at convergence (marg_kl itself is already
            // computed in the errRp loop above at zero extra cost, but the full metrics
            // pass is needed to populate the other result fields for diagnostics).
            double grake_norm = 0.0;

            double mean_err_sum = 0.0;
            double kl_max       = 0.0;
            double chi2_total   = 0.0;
            if (need_extra_metrics && W_total > 0.0) {
                // 773f.4: fuse extra-metrics sweep into single compute_cell_metrics call.
                // S_lin is sized to max_cat — satisfies bucket >= max(cat_counts[k]) requirement.
                // X is the post-capacity cell mass vector (valid for both linear and log paths).
                // Field mapping (verified against manual loop above):
                //   cm.mean_err  = mean_sum / K  ↔ mean_err_sum / K  (same formula)
                //   cm.kl        = max_k kl_k     ↔ kl_max            (same)
                //   cm.chi2      = sum chi2        ↔ chi2_total        (same)
                //   cm.grake_norm= max nm           ↔ grake_norm        (same)
                // errRp already set from first pass above; cm.errRp is identical — discard.
                // marg_kl is NOT in CellMetrics — retained from the first pass unchanged.
                const lbw::CellMetrics cm = lbw::compute_cell_metrics(st, ct, X, W_total, S_lin);
                mean_err_sum = cm.mean_err * static_cast<double>(st.K);
                kl_max       = cm.kl;
                chi2_total   = cm.chi2;
                grake_norm   = cm.grake_norm;
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
                    if (iter == 1 && res.base.metric_first_check == std::numeric_limits<double>::infinity())
                        res.base.metric_first_check = curr_best;
                    if (std::isfinite(nat_metric_prev_nonavec)) {
                        res.base.metric_prev_check = nat_metric_prev_nonavec;
                        res.base.prev_check_iter   = nat_iter_prev_nonavec;
                    }
                    nat_metric_prev_nonavec = curr_best;
                    nat_iter_prev_nonavec   = iter;
                    if (std::isfinite(curr_best) && curr_best < best.best_metric) {
                        for (int c = 0; c < ct.M_cell; c++)
                            w_ratio[c] = (X_init[c] > 0.0) ? X[c] / X_init[c] : 0.0;
                        best.update(curr_best,
                                    lbw::compute_weight_kl(X, X_init, ct.M_cell, st.n, kl_ratio_buf.data(), kl_weight_buf.data()),
                                    iter, w_ratio);
                    }
                }
            }
            double mean_err = (st.K > 0) ? (mean_err_sum / static_cast<double>(st.K)) : 0.0;

            // Store metrics in result struct (unconditional — intermediate checks
            // store 0 for gated metrics; final iter always populates all fields).
            res.base.l1_weight_change = l1_weight;    // real Σ|ΔX|/W_input (replaces pct_change stub)
            res.base.grake_norm       = grake_norm;   // max_kj normalized margin residual
            res.base.mean_error       = mean_err;
            res.base.kl               = kl_max;
            res.base.chi2             = chi2_total;

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
                              "ORIS iter %d: errRp=%.3e", iter, errRp);
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
                    static constexpr double kLn10 = 2.302585092994046;  // std::log(10.0)
                    std::snprintf(msg, sizeof(msg),
                                  "  margin=%d: log10(f) range [%.2f, %.2f]",
                                  k + 1,
                                  lf_min / kLn10,
                                  lf_max / kLn10);
                    st.log(msg);
                }
            }

            // metric + rule dispatch.
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
                if (converged && converged_abs) absolute_tol_fired = true;

                // Intermediate homotopy levels: allow warm-jump when errRp is loose.
                if (!converged && lvl < N_levels - 1) {
                    converged = (errRp < tol_lvl);
                }

                if (converged) {
                    // Early exit from this homotopy level. If this is the final level,
                    // set terminal status; else warm-jump to the next (tighter) level.
                    level_converged = true;
                    if (lvl == N_levels - 1) {
                        res.base.status = structural_infeas_pairs.empty() ? RK_OK : RK_ERR_INFEAS;
                        // za9r: pin the convergence-firing iter here rather than
                        // letting oris_finalize read res.base.iterations at exit
                        // (robust if any post-convergence work advances the counter).
                        res.base.convergence_iter = res.base.iterations;
                    }
                    res.final_alpha = alpha;
                    break;
                }
            }
        }
        res.final_alpha = alpha;
    }  // end for (iter_in_lvl)
        }  // end if (!sraa_active_lvl)
        // Log-space overflow: oris inner loop sets res.base.status = RK_ERR_NOCONV
        // and res.base.max_error = +inf via `break`. Must not proceed to next level.
        if (!std::isfinite(res.base.max_error)) {
            homotopy_break = true;
        }
        total_iters = res.base.iterations;  // global iter counter; picks up partial level
        if (level_converged && lvl == N_levels - 1) {
            homotopy_break = true;
        }
    }  // end homotopy level loop

    // ════════════════════ Post-loop: ALM projection + obs expansion + bounds ════════════════════
    oris_finalize(st, res, ct, X, X_init, L_cell, U_cell,
                   alm_active, capacity_mu_adaptive, lambda_cell,
                   best, absolute_tol_fired, structural_infeas_pairs,
                   sor_min_omega, sor_n_damped, probe_samples);
    return res;
}

} // namespace lbw
