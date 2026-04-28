#ifndef LEAFBLOWER_H
#define LEAFBLOWER_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>  /* size_t */
#include <stdint.h>  /* int32_t, int64_t */

/* Bounds enforcement mode (appended field in rk_params_t; default RK_BOUNDS_CELL) */
typedef enum {
    RK_BOUNDS_CELL = 0,  /* cell-aggregate bounds (current default; may violate per-obs) */
    RK_BOUNDS_UNIT = 1   /* per-observation strict bounds via intra-cell water-filling */
} rk_bounds_mode_t;

/* ── Overlay enums ── */
typedef enum { RK_SCHED_ROUND_ROBIN = 0, RK_SCHED_GREEDY = 1 } rk_scheduler_t;
typedef enum { RK_ETA_FIXED = 0, RK_ETA_TANG_DYNAMIC = 1 } rk_eta_mode_t;

/* Homotopy config (default = single level = identity) */
typedef struct {
    int    n_levels;
    double start_factor;
    double end_factor;
    double budget_split_p;
    int    enabled;
} rk_homotopy_cfg_t;
/* ── End overlay enums ── */

/* ── Return codes ── */
#define RK_OK           0  /* Converged: improvement criterion satisfied */
#define RK_ERR_NOCONV   1  /* Legacy alias — no longer emitted by new solvers */
#define RK_ERR_INFEAS   2  /* Infeasible: empty cell with positive target */
#define RK_ERR_BADARG   3  /* Invalid argument */
#define RK_ERR_BUDGET   4  /* Budget exhausted while loss still decreasing; increase max_iterations */
#define RK_ERR_STALL    5  /* Loss function plateau — at constrained optimum; weights are valid */

/* ── Algorithm selector ── */
typedef enum {
    RK_ALG_AUTO   = 0,
    RK_ALG_IEPPA  = 1,
    RK_ALG_LBFGSB = 2,
    RK_ALG_RAKING    = 3,
    RK_ALG_SINKHORN  = 4,
    RK_ALG_CHEBYSHEV = 5,
    RK_ALG_GREG      = 6,
    RK_ALG_GRAKE      = 7,
    RK_ALG_IEPPA_SOFT = 8    /* ieppa + ADMM soft capacity enforcement */
} rk_algorithm_t;

/* ── Calibration parameters ── */
typedef struct {
    double          min_weight;      /* default 0.0 */
    double          max_weight;      /* default 5.0 */
    int             inner_max_iter;  /* inner BCD cap per outer iter, default 500 */
    int             outer_max_iter;  /* outer EPP / L-BFGS max iters, default 50.
                                       * Note: R bridge sets outer_max_iter = inner_max_iter;
                                       * independent control requires direct C API use. */
    double          tol_abs;         /* convergence tolerance, default 1e-6.
                                       * iEPPA: max absolute primal error max_k max_j |S_kj/W - tau_kj|.
                                       * L-BFGS-B: normalized dual gradient norm maxAbs(grad)/W. */
    rk_algorithm_t  algorithm;       /* default RK_ALG_AUTO */
    int             verbose;         /* 0=silent, 1=progress, 2=debug */
    double          epsilon;         /* deprecated: no longer read by any solver; kept for ABI compat */
    int             lbfgs_m;         /* L-BFGS history size, default 10 */
    void            (*log_fn)(const char* msg, void* ctx);
    void*           log_ctx;
    rk_bounds_mode_t bounds_mode;  /* default RK_BOUNDS_CELL */
    /* ── Overlay knobs ── */
    rk_homotopy_cfg_t homotopy;         /* default: n_levels=1, enabled=0 */
    rk_scheduler_t    scheduler;        /* default RK_SCHED_ROUND_ROBIN */
    rk_eta_mode_t     eta_mode;         /* default RK_ETA_FIXED */
    double            eta_start;
    double            eta_end;
    double            eta_schedule_power;
    /* ── End overlay knobs ── */
    /* ── Convergence config ── */
    double pct_tol;          /* threshold for IMPROVEMENT/PLATEAU rules (default 0.001) */
    double absolute_tol;     /* threshold for THRESHOLD rule + stop_when secondary (default 0.0) */
    int    metric;           /* CalibMetric: 0=MAX_ERR 1=MEAN_ERR 2=KL 3=CHI2 4=GRAKE_NORM 5=L1_WEIGHT */
    int    rule;             /* CalibRule: 0=THRESHOLD 1=IMPROVEMENT 2=PLATEAU */
    int    stop_when;        /* 0=ANY 1=ALL */
    /* ── SOR config (iEPPA only; ignored by raking/lbfgsb) ── */
    int    sor_enabled;
    int    sor_auto;
    double sor_omega_init;
    double sor_omega_min;
    double sor_omega_fixed;  /* -1.0 = use auto */
    int    sor_burnin;
    /* ── End convergence/SOR config ── */
} rk_params_t;

/* ── Result ── */
typedef struct {
    int             status;          /* RK_OK / RK_ERR_* */
    int             iterations;      /* outer iterations completed */
    double          max_error;       /* max calibration error at last iterate */
    rk_algorithm_t  algorithm_used;  /* actual algorithm run (never RK_ALG_AUTO) */
    char            message[256];    /* null-terminated status message */
    int             n_xcur_writes_per_iter_linear;  /* P1.1 diagnostic */
    double          min_alpha_seen;                 /* P2.1: min alpha over all sweeps; 1.0 if never damped */
    double          final_alpha;                    /* P2.1: alpha at solver exit */
    int             n_bounds_violated;  /* cell-mode diagnostic: count of w_i outside bounds (no action) */
    int             n_bounds_clamped;   /* unit-mode action: count of w_i clamped after water-fill exhausted */
    /* ── Overlay diagnostics ── */
    int             homotopy_levels_used;   /* 0 = homotopy disabled */
    double          homotopy_final_factor;  /* max_weight multiplier at last level */
    int             greedy_sweeps_taken;    /* greedy scheduler sweeps per last inner pass */
    double          eta_final;             /* alm_mu multiplier at exit; 0.0 = N/A */
    /* ── End overlay diagnostics ── */
    /* ── Extended quality metrics ── */
    double mean_error;       /* L1-over-margins */
    double kl;               /* max KL across margins */
    double chi2;             /* total chi-square */
    double l1_weight_change;   /* L1 normalized weight change Σ|Δw|/n (renamed from pct_change) */
    double grake_norm;         /* survey::grake normalized residual max_k|misfit|/(1+|pop|) */
    int    convergence_metric; /* CalibMetric at exit (0=MAX_ERR...) */
    int    convergence_rule;   /* CalibRule at exit (0=THRESHOLD...) */
    double convergence_tol;    /* threshold that fired */
    int    convergence_iter;   /* iteration at convergence (-1 if max_iter) */
    double best_error;         /* errRp at best iterate */
    int    best_iter;
    double sor_min_omega;    /* iEPPA only; non-iEPPA = 1.0 */
    int    sor_n_damped;     /* iEPPA only; non-iEPPA = 0 */
    double convergence_solver_objective;  /* solver's mathematical objective at best_iter */
    int    convergence_minimized_metric; /* CalibMetric: which metric was minimized */
    /* ── End extended quality metrics ── */
} rk_result_t;

/* Fill *p with safe defaults */
void rk_params_init(rk_params_t* p);

/* Zero-initialize *r (memset); call before passing to rk_calibrate */
void rk_result_init(rk_result_t* r);

/*
 * Calibrate survey weights in-place.
 *   n          : number of observations
 *   K          : number of calibration margins
 *   weights    : [n] on input = design weights; on output = calibrated weights
 *   group_ids  : [K] array of pointers; group_ids[k][i] in {-1, 0..cat_counts[k]-1}
 *   cat_counts : [K] number of categories per margin
 *   targets    : [K] array of pointers; targets[k][j] = target proportion for cat j
 *   params     : calibration parameters (NULL = use defaults)
 *   result     : output result struct (NULL = ignore)
 * Returns RK_OK, RK_ERR_NOCONV, RK_ERR_INFEAS, or RK_ERR_BADARG.
 *
 * C++17 note: the C++ compilation unit (c_api.cpp) applies [[nodiscard]] to
 * this function via a wrapper; the C header stays C99-clean.
 */
int rk_calibrate(
    int n, int K,
    double* weights,
    const int32_t** group_ids,
    const int* cat_counts,
    const double** targets,
    const rk_params_t* params,
    rk_result_t* result
);

#ifdef __cplusplus
}
#endif

/* ABI tripwires. If EXPECTED_RK_PARAMS_BYTES fails, a new field was added —
 * update this value after auditing ABI consumers. */
#ifdef __cplusplus
static_assert(RK_ALG_AUTO == 0, "memset(0) default must equal RK_ALG_AUTO");
/* rk_result_t tripwire. Linux x86_64, verified 2026-04-24: 448 bytes. */
#define EXPECTED_RK_RESULT_BYTES 448
static_assert(sizeof(rk_result_t) == EXPECTED_RK_RESULT_BYTES,
    "rk_result_t size changed; update EXPECTED_RK_RESULT_BYTES and ABI consumers");
/* Compute sizeof(rk_params_t) on the target platform at implementation time
 * and hard-code it here. Record the value in a comment. Example:
 *   Linux x86_64 GCC 13, verified 2026-04-24: 72 bytes.
 * After measuring, replace the placeholder below with the actual value. */
/* ABI layout (2026-04-24): added overlay fields after bounds_mode.
 *   rk_homotopy_cfg_t (n_levels int + 4B pad + 3 doubles + enabled int + 4B pad = 40B)
 *   rk_scheduler_t (int, 4B) + rk_eta_mode_t (int, 4B)
 *   eta_start (double, 8B) + eta_end (double, 8B) + eta_schedule_power (double, 8B)
 * WU-A (2026-04-25): convergence redesign — criterion→metric+rule in rk_params_t:
 *   pct_tol (double, 8B) + absolute_tol (double, 8B)
 *   metric (int, 4B) + rule (int, 4B) + stop_when (int, 4B)
 *   sor_enabled (int, 4B) + sor_auto (int, 4B)
 *   sor_omega_init (double, 8B) + sor_omega_min (double, 8B) + sor_omega_fixed (double, 8B)
 *   sor_burnin (int, 4B) + 4B pad
 * Total: 224B. Verified 2026-04-25 Linux x86_64. */
#define EXPECTED_RK_PARAMS_BYTES 224
static_assert(sizeof(rk_params_t) == EXPECTED_RK_PARAMS_BYTES,
              "rk_params_t size changed; check ABI consumers");
#endif

#endif /* LEAFBLOWER_H */
