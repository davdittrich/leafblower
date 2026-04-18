#include "leafblower.h"
#include "types.hpp"
#include "lbfgsb_solver.hpp"
#include "ieppa.hpp"
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
    bool use_logit = std::isfinite(p->max_weight);
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

static rk_algorithm_t select_algorithm(int n, int K,
                                        const int* cat_counts,
                                        const rk_params_t* p,
                                        int64_t& complexity_out) {
    complexity_out = INT64_C(0);
    for (int k = 0; k < K; k++) complexity_out += (int64_t)n * cat_counts[k];
    if (p->algorithm != RK_ALG_AUTO) return p->algorithm;
    if (complexity_out > 500000L || std::isfinite(p->max_weight) || p->min_weight > 0.0)
        return RK_ALG_IEPPA;
    return RK_ALG_LBFGSB;
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

    int64_t complexity = INT64_C(0);
    rk_algorithm_t alg = select_algorithm(n, K, cat_counts, p, complexity);

    // Build CalibState
    lbw::CalibState st;
    st.n = n; st.K = K;
    st.weights = weights;
    st.group_ids = group_ids;
    st.cat_counts = cat_counts;
    st.targets = targets;
    st.min_weight    = p->min_weight;
    st.max_weight    = p->max_weight;
    st.tol_abs       = p->tol_abs;
    st.inner_max_iter = p->inner_max_iter;
    st.outer_max_iter = p->outer_max_iter;
    st.epsilon       = p->epsilon;
    st.lbfgs_m       = p->lbfgs_m;
    st.verbose       = p->verbose;
    st.log_fn        = p->log_fn;
    st.log_ctx       = p->log_ctx;
    st.total_cats    = 0;
    for (int k = 0; k < K; k++) st.total_cats += cat_counts[k];

    // Verbose routing report (complexity already computed by select_algorithm)
    if (p->verbose >= 1) {
        char msg[256];
        if (alg == RK_ALG_IEPPA)
            snprintf(msg, 256, "Auto-selected iEPPA: complexity=%lld, max_weight=%.2f, min_weight=%.2f",
                     (long long)complexity, p->max_weight, p->min_weight);
        else
            snprintf(msg, 256, "Auto-selected L-BFGS-B: complexity=%lld <= 500000, max_weight=%.2f, min_weight=%.2f",
                     (long long)complexity, p->max_weight, p->min_weight);
        st.log(msg);
    }

    int status;
    int iterations;
    double max_error;
    rk_algorithm_t used;

    if (alg == RK_ALG_LBFGSB) {
        auto res = lbw::lbfgsb_solve(st);
        status = res.status;
        iterations = res.iterations;
        max_error = res.max_error;
        used = RK_ALG_LBFGSB;
    } else {
        // iEPPA Sinkhorn-BCD solver
        auto res = lbw::ieppa_solve(st);
        status = res.status;
        iterations = res.iterations;
        max_error = res.max_error;
        used = RK_ALG_IEPPA;
    }

    if (result) {
        result->status = status;
        result->iterations = iterations;
        result->max_error = max_error;
        result->algorithm_used = used;
        if (result->message[0] == '\0') {
            snprintf(result->message, 256,
                     "%s: %d iters, max_error=%.2e",
                     used == RK_ALG_LBFGSB ? "L-BFGS-B" : "iEPPA",
                     iterations, max_error);
        }
    }
    return status;
}

} // extern "C"
