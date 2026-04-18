#include "leafblower.h"
#include "types.hpp"
#include <cstring>
#include <cstdio>
#include <cmath>
#include <climits>
#include <cstdint>

// C++17 [[nodiscard]] on rk_calibrate — silently ignored on C++14
#if __cplusplus >= 201703L
  #define LBW_NODISCARD [[nodiscard]]
#else
  #define LBW_NODISCARD
#endif

extern "C" {

void rk_params_init(rk_params_t* p) {
    if (!p) return;
    p->min_weight    = 0.0;
    p->max_weight    = 5.0;
    p->inner_max_iter = 500;
    p->outer_max_iter = 50;
    p->tol_abs       = 1e-6;
    p->algorithm     = RK_ALG_AUTO;
    p->verbose       = 0;
    p->epsilon       = 0.05;
    p->lbfgs_m       = 10;
    p->log_fn        = nullptr;
    p->log_ctx       = nullptr;
}

static int validate_inputs(int n, int K,
                            const double* weights,
                            const int32_t** group_ids,
                            const int* cat_counts,
                            const double** targets,
                            const rk_params_t* p,
                            rk_result_t* result) {
    auto err = [&](const char* msg) -> int {
        if (result) {
            result->status = RK_ERR_BADARG;
            snprintf(result->message, 256, "%s", msg);
        }
        return RK_ERR_BADARG;
    };

    if (!weights)   return err("weights pointer is NULL");
    if (!group_ids) return err("group_ids pointer is NULL");
    if (!cat_counts)return err("cat_counts pointer is NULL");
    if (!targets)   return err("targets pointer is NULL");
    if (n <= 0)     return err("n must be > 0");
    if (K <= 0)     return err("K must be > 0");

    if (p->min_weight >= p->max_weight)
        return err("min_weight must be strictly less than max_weight");

    // Logit singularity guard — each checked independently
    bool use_logit = (p->min_weight > 0.0) && std::isfinite(p->max_weight);
    if (use_logit) {
        if (p->min_weight == 1.0)
            return err("logit link undefined: min_weight=1 makes denominator (1-L)=0");
        if (p->max_weight == 1.0)
            return err("logit link undefined: max_weight=1 makes denominator (U-1)=0");
    }

    // cat_counts checks + overflow guard
    size_t total_cats = 0;
    for (int k = 0; k < K; k++) {
        if (cat_counts[k] <= 0)
            return err("cat_counts[k] must be > 0 for all k");
        if (cat_counts[k] > n)
            return err("cat_counts[k] > n: more categories than observations");
        total_cats += (size_t)cat_counts[k];
    }
    if ((size_t)n * total_cats > SIZE_MAX / 2)
        return err("problem too large for platform size_t");

    // Initial weight checks
    double total_w = 0.0;
    for (int i = 0; i < n; i++) {
        if (!std::isfinite(weights[i]))
            return err("NaN or Inf in initial weights[]");
        total_w += weights[i];
    }
    if (total_w < 1e-15)
        return err("total weight is zero or negative");

    // targets checks
    for (int k = 0; k < K; k++) {
        if (!targets[k]) return err("targets[k] is NULL");
        double sum = 0.0;
        for (int j = 0; j < cat_counts[k]; j++) {
            if (!std::isfinite(targets[k][j]))
                return err("NaN or Inf in targets[]");
            if (targets[k][j] < 0.0)
                return err("targets[k][j] < 0");
            sum += targets[k][j];
        }
        if (std::fabs(sum - 1.0) > 1e-8)
            return err("targets[k] does not sum to 1 (within 1e-8)");
    }

    // group_ids range validation — full O(n*K) pass before any weight modification
    for (int k = 0; k < K; k++) {
        if (!group_ids[k]) return err("group_ids[k] is NULL");
        for (int i = 0; i < n; i++) {
            int g = group_ids[k][i];
            if (g < -1)
                return err("group_ids[k][i] < -1: only -1 (NA) is valid");
            if (g >= cat_counts[k])
                return err("group_ids[k][i] >= cat_counts[k]");
        }
    }

    return RK_OK;
}

LBW_NODISCARD int rk_calibrate(int n, int K,
                                double* weights,
                                const int32_t** group_ids,
                                const int* cat_counts,
                                const double** targets,
                                const rk_params_t* params,
                                rk_result_t* result) {
    rk_params_t defaults;
    rk_params_init(&defaults);
    const rk_params_t* p = params ? params : &defaults;

    if (result) {
        result->status = RK_OK;
        result->iterations = 0;
        result->max_error = 0.0;
        result->algorithm_used = RK_ALG_AUTO;
        result->message[0] = '\0';
    }

    int rc = validate_inputs(n, K, weights, group_ids, cat_counts, targets, p, result);
    if (rc != RK_OK) return rc;

    // Stub: algorithm dispatch not yet implemented
    if (result) {
        result->status = RK_ERR_NOCONV;
        snprintf(result->message, 256, "rk_calibrate not fully implemented");
    }
    return RK_ERR_NOCONV;
}

} // extern "C"
