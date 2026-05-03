#pragma once
#include "leafblower.h"
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <vector>

namespace lbw {

// Validate inputs for rk_calibrate. Returns RK_OK or RK_ERR_BADARG.
// On RK_ERR_BADARG, sets result->status and result->message if result != nullptr.
// This is the single source of truth for input validation — called from both
// c_api.cpp (rk_calibrate C ABI path) and r_bridge.cpp (direct C++ path).
inline int validate_calibrate_inputs(int n, int K,
    const double* weights, const int32_t** group_ids,
    const int* cat_counts, const double** targets,
    const rk_params_t* p, rk_result_t* result, rk_algorithm_t alg)
{
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
    if (K > 64)     return err("K exceeds maximum (64); too many margin columns");

    if (p->min_weight >= p->max_weight)
        return err("min_weight must be strictly less than max_weight");

    // Logit singularity guard — only applies to L-BFGS-B.
    // iEPPA uses multiplicative IPF scaling; no logit link, no singularity.
    // logit_scale = (U-L)/((U-1)*(1-L)); denominator → 0 when L or U near 1,
    // producing logit_scale ~1/eps → overflow/cancellation in L-BFGS-B.
    if (alg == RK_ALG_LBFGSB) {
        const double kSingularityEps = 1e-6;
        if (std::fabs(p->min_weight - 1.0) < kSingularityEps)
            return err("logit link undefined: min_weight near 1 makes denominator (1-L)~0");
        if (std::fabs(p->max_weight - 1.0) < kSingularityEps)
            return err("logit link undefined: max_weight near 1 makes denominator (U-1)~0");
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
    if (total_cats > 0 && (size_t)n > SIZE_MAX / 2 / (size_t)total_cats)
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
        if (std::fabs(sum - 1.0) > 1e-6)
            return err("targets[k] does not sum to 1 (within 1e-6)");
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

    // B13: NA-only category with positive target → structural INFEAS.
    // O(n×K) two-pass bitmap: replaces O(K×ΣC×n) triple-nested scan.
    // Pass 1: mark which (k,j) combinations appear in data — O(n×K)
    std::vector<std::vector<bool>> seen(K);
    for (int k = 0; k < K; k++) seen[k].assign(cat_counts[k], false);
    for (int i = 0; i < n; i++)
        for (int k = 0; k < K; k++) {
            int g = group_ids[k][i];
            if (g >= 0 && g < cat_counts[k]) seen[k][g] = true;
        }
    // Pass 2: check targeted but not-in-data categories — O(K×max_cats)
    for (int k = 0; k < K; k++) {
        for (int j = 0; j < cat_counts[k]; j++) {
            if (!seen[k][j] && targets[k][j] > 1e-12) {
                char msg[256];
                std::snprintf(msg, sizeof(msg),
                    "margin %d, category %d: target=%.6g but 0 observations assigned "
                    "(all NA or missing); structurally infeasible problem",
                    k, j, targets[k][j]);
                if (result) {
                    result->status = RK_ERR_INFEAS;
                    std::snprintf(result->message, 256, "%s", msg);
                }
                return RK_ERR_INFEAS;
            }
        }
    }

    if (!std::isfinite(p->tol_abs) || p->tol_abs <= 0.0)
        return err("tol_abs must be finite and positive");

    /* Enum range guards: prevent silent UB from out-of-range static_cast */
    if (p->metric < 0 || p->metric > 6)
        return err("metric out of range [0,6]: 0=MAX_ERR 1=MEAN_ERR 2=KL 3=CHI2 4=GRAKE_NORM 5=L1_WEIGHT 6=MARGINAL_KL");
    if (p->rule < 0 || p->rule > 2)
        return err("rule out of range [0,2]: 0=THRESHOLD 1=IMPROVEMENT 2=PLATEAU");
    if (p->stop_when < 0 || p->stop_when > 1)
        return err("stop_when out of range [0,1]: 0=ANY 1=ALL");

    return RK_OK;
}

} // namespace lbw
