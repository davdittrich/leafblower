#pragma once
#include <cstdint>
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
    bool   enabled         = false; // master toggle; set true when n_levels > 1
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
    double omega_fixed   = -1.0;  // sentinel: use auto
    int    burnin        = 20;
};

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
    int lbfgs_m;
    int verbose;
    bool ieppa_auto_selected  = false;  // true iff AUTO routing selected iEPPA; used for verbose prefix
    bool use_admm_capacity    = false;  // ieppa_soft: ADMM P1.1; default false = hard clamp
    double alm_lambda = 0.0;  // dual variable for sum(w)=n; only read when alm_mu > 0
    double alm_mu     = 0.0;  // penalty coefficient; 0.0 = ALM inactive
    rk_bounds_mode_t bounds_mode = RK_BOUNDS_CELL;  /* P3.1: per-obs vs cell-aggregate bounds */
    // ── Overlay config ──
    HomotopyConfigLbw    homotopy;
    SchedulerConfigLbw   scheduler;
    EtaScheduleConfigLbw eta_schedule;
    CalibConvergenceCfg  convergence_cfg;
    CalibSorCfg          sor_cfg;
    bool                 accelerate = false;  // SQUAREM outer loop for raking
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
