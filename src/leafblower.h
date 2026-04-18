#ifndef LEAFBLOWER_H
#define LEAFBLOWER_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>  /* size_t */
#include <stdint.h>  /* int32_t, int64_t */

/* ── Return codes ── */
#define RK_OK         0  /* Success */
#define RK_ERR_NOCONV 1  /* Did not converge within outer_max_iter */
#define RK_ERR_INFEAS 2  /* Infeasible: empty cell with positive target */
#define RK_ERR_BADARG 3  /* Invalid argument */

/* ── Algorithm selector ── */
typedef enum {
    RK_ALG_AUTO   = 0,
    RK_ALG_IEPPA  = 1,
    RK_ALG_LBFGSB = 2
} rk_algorithm_t;

/* ── Calibration parameters ── */
typedef struct {
    double          min_weight;      /* default 0.0 */
    double          max_weight;      /* default 5.0 */
    int             inner_max_iter;  /* inner BCD cap per outer iter, default 500 */
    int             outer_max_iter;  /* outer EPP / L-BFGS max iters, default 50 */
    double          tol_abs;         /* convergence tolerance, default 1e-6 */
    rk_algorithm_t  algorithm;       /* default RK_ALG_AUTO */
    int             verbose;         /* 0=silent, 1=progress, 2=debug */
    double          epsilon;         /* iEPPA entropic parameter, default 0.05 */
    int             lbfgs_m;         /* L-BFGS history size, default 10 */
    void            (*log_fn)(const char* msg, void* ctx);
    void*           log_ctx;
} rk_params_t;

/* ── Result ── */
typedef struct {
    int             status;          /* RK_OK / RK_ERR_* */
    int             iterations;      /* outer iterations completed */
    double          max_error;       /* max calibration error at last iterate */
    rk_algorithm_t  algorithm_used;  /* actual algorithm run (never RK_ALG_AUTO) */
    char            message[256];    /* null-terminated status message */
} rk_result_t;

/* Fill *p with safe defaults */
void rk_params_init(rk_params_t* p);

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

#endif /* LEAFBLOWER_H */
