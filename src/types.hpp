#pragma once
#include <cstdint>
#include <limits>
#include <vector>
#include "leafblower.h"   /* rk_bounds_mode_t */
#ifndef LBW_NO_R
#  include <R_ext/Print.h>
#else
#  include <cstdio>
#endif

namespace lbw {

// ── Overlay config structs ──────────────────────────────────────────────────
struct HomotopyConfigLbw {
    int    n_levels        = 1;     // 1 = disabled (single level = current behaviour)
    double start_factor    = 1.0;   // starting max_weight multiplier
    double end_factor      = 1.0;   // ending max_weight multiplier
    double budget_split_p  = 0.5;   // Chizat-inspired budget split
    // enabled is derived: (n_levels > 1) — no separate field to desync
};

enum class SchedulerMode   : int { ROUND_ROBIN = 0, GREEDY = 1 };
enum class EtaScheduleMode : int { FIXED = 0, TANG_DYNAMIC = 1 };

struct SchedulerConfigLbw {
    SchedulerMode mode                     = SchedulerMode::ROUND_ROBIN;
    double        residual_recheck_fraction = 0.1;  // internal; not exposed to ABI or R
};

struct EtaScheduleConfigLbw {
    EtaScheduleMode mode          = EtaScheduleMode::FIXED;
    double          eta_start     = 1.0;
    double          eta_end       = 1.0;
    double          schedule_power = 0.5;
};
// ── End overlay config structs ─────────────────────────────────────────────

enum class CalibMetric : int {
    MAX_ERR    = 0,
    MEAN_ERR   = 1,
    KL         = 2,
    CHI2       = 3,
    GRAKE_NORM   = 4,
    L1_WEIGHT    = 5,
    MARGINAL_KL  = 6    // Σ_k Σ_j t_kj log(t_kj / achieved_kj) — calibration quality
};

enum class CalibRule : int {
    THRESHOLD   = 0,
    IMPROVEMENT = 1,
    PLATEAU     = 2
};
enum class CalibStopWhen : int { ANY = 0, ALL = 1 };

struct CalibConvergenceCfg {
    double        pct_tol      = 0.001;
    double        absolute_tol = 0.0;
    CalibMetric   metric       = CalibMetric::MAX_ERR;
    CalibRule     rule         = CalibRule::IMPROVEMENT;
    CalibStopWhen stop_when    = CalibStopWhen::ANY;
};

struct CalibSorCfg {
    bool   enabled       = true;
    bool   auto_adapt    = true;
    double omega_init    = 1.0;
    double omega_min     = 0.3;
    double omega_max     = 1.5;  // recovery ceiling; 1.5 is in the proven (0,2) window (Thibault 2021)
    double omega_fixed   = -1.0;  // sentinel: use auto
    int    burnin        = 20;  // iterations before SOR adaptation starts; ORIS only (raking ignores)
    // omega_mode_id: 0=heuristic (0.7 damp/1.05 grow), 1=fixed (omega_max),
    //                2=iterate-change (free-coord ||dX_free||^2 estimator, default; e18t.9 SHIP)
    int    omega_mode_id = 2;
    bool   corun_aa      = false;  // e65t.1: allow SOR on AA steps; default false = exclusion
                                   // DEAD (e65t.1 NO-GO): co-run never enabled; default false. See bd leafblower-e65t.1.
};

struct ALMConfig {
    double lambda      = 0.0;  // dual variable for sum(w)=n; only read when mu > 0
    double mu          = 0.0;  // penalty coefficient; 0.0 = ALM inactive
    double capacity_mu = 0.0;  // oris_soft ALM penalty (capacity box constraint); 0.0 = inactive
};

// ── CalibResult: shared base fields present in every solver result struct ────
// Solver structs embed CalibResult and add solver-specific extras (message,
// M_cell, n_factorizations, ALM diagnostics, ORIS internals, etc.).
// Invariant: no ORIS-private diagnostics here.
struct CalibResult {
    int    status              = RK_ERR_NOCONV;
    int    iterations          = 0;
    double max_error           = 1.0;   // 1.0 is the dominant existing default
    double mean_error          = 0.0;
    double kl                  = 0.0;
    double chi2                = 0.0;
    double l1_weight_change    = 0.0;
    double grake_norm          = 0.0;
    int    convergence_metric  = 0;
    int    convergence_rule    = 1;
    double convergence_tol     = 0.001;
    int    convergence_iter    = -1;
    double best_error          = std::numeric_limits<double>::infinity();
    int    best_iter           = 0;
    double metric_first_check  = std::numeric_limits<double>::infinity();
    double metric_prev_check   = std::numeric_limits<double>::infinity();
    int    prev_check_iter     = -1;
    std::vector<double> best_weights;
    double convergence_solver_objective = std::numeric_limits<double>::infinity();
    int    convergence_minimized_metric = 0;
    // stall_kind: 0=no stall, 1=wchange (weight-change plateau, SRAA path),
    //             2=kl (metric/KL plateau, plain-IPF path).
    // Set by the solver at the RK_ERR_STALL emission site; read by r_bridge and harvest.
    // Replaces the accelerate_bool heuristic in harvest.R (leafblower-8eod).
    int    stall_kind                   = 0;
};
// ── End CalibResult ───────────────────────────────────────────────────────────

struct CalibState {
    int n;
    int K;
    double* weights;                   // caller-owned, n elements
    const int32_t** group_ids;        // caller-owned, K pointers to n int32_t elements
    const int* cat_counts;            // caller-owned, K elements
    const double** targets;           // caller-owned, K pointers to cat_counts[k] doubles
    double min_weight;
    double max_weight;
    double tol_abs;
    int inner_max_iter;
    int outer_max_iter;
    int verbose;
    bool oris_auto_selected   = false;  // true iff AUTO routing selected ORIS; used for verbose prefix
    bool use_admm_capacity    = false;  // oris_soft: ADMM P1.1; default false = hard clamp
    ALMConfig alm;  // ALM penalty state: lambda, mu, capacity_mu
    rk_bounds_mode_t bounds_mode = RK_BOUNDS_CELL;  /* P3.1: per-obs vs cell-aggregate bounds */
    // ── Overlay config ──
    HomotopyConfigLbw    homotopy;
    SchedulerConfigLbw   scheduler;
    EtaScheduleConfigLbw eta_schedule;
    CalibConvergenceCfg  convergence_cfg;
    CalibSorCfg          sor_cfg;
    bool                 accelerate = false;  // SQUAREM outer loop for raking
    bool                 jacobi_log = false;  // log path: freeze cell_lf at iter start (Jacobi semantics)
    double               newton_tsvd_ratio = 1e-8;  // newton_kl: TSVD truncation ratio (Epic-H WH-e); <=0 falls back to 1e-8
    double               ridge_lambda      = 0.0;   // Tikhonov ridge on dual λ: newton_kl H_pre[k,k]+=, greg N[j,j]+=; 0=off
    double               gk_omega          = 1.0;   // e65t.2: greenkhorn over-relaxation exponent
    // ── End overlay config ──
    void (*log_fn)(const char* msg, void* ctx);
    void* log_ctx;

    void log(const char* msg) const {
        if (verbose <= 0) return;
        if (log_fn) {
            log_fn(msg, log_ctx);
        } else {
#ifndef LBW_NO_R
            REprintf("%s\n", msg);
#else
            fprintf(stderr, "%s\n", msg);
#endif
        }
    }
};

} // namespace lbw
