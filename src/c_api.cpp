#include "leafblower.h"
#include "types.hpp"
#include "lbfgsb_solver.hpp"
#include "ieppa.hpp"   // faithful (new)
#include "raking.hpp"  // renamed hybrid
#include <cstring>
#include <cstdio>
#include <cmath>
#include <climits>
#include <cstdint>
#include <limits>

// C++17 [[nodiscard]] on rk_calibrate — silently ignored on C++14
#if __cplusplus >= 201703L
  #define LBW_NODISCARD [[nodiscard]]
#else
  #define LBW_NODISCARD
#endif

extern "C" {

void rk_params_init(rk_params_t* p) {
    if (!p) return;
    memset(p, 0, sizeof(*p));
    p->min_weight    = 0.0;
    p->max_weight    = 5.0;
    p->inner_max_iter = 500;
    p->outer_max_iter = 50;
    p->tol_abs       = 1e-6;
    p->algorithm     = RK_ALG_AUTO;
    p->verbose       = 0;
    p->epsilon       = 0.0;   /* deprecated: not read by any solver; see leafblower.h */
    p->lbfgs_m       = 10;
    p->log_fn        = nullptr;
    p->log_ctx       = nullptr;
    p->bounds_mode   = RK_BOUNDS_CELL;  /* explicit for clarity; memset=0 already gives this */
    p->homotopy.n_levels       = 1;
    p->homotopy.start_factor   = 1.0;
    p->homotopy.end_factor     = 1.0;
    p->homotopy.budget_split_p = 0.5;
    p->homotopy.enabled        = 0;
    p->scheduler               = RK_SCHED_ROUND_ROBIN;
    p->eta_mode                = RK_ETA_FIXED;
    p->eta_start               = 1.0;
    p->eta_end                 = 1.0;
    p->eta_schedule_power      = 0.5;
    p->pct_tol                 = 0.001;
    p->absolute_tol            = 0.0;
    p->criterion               = 0;  /* PCT */
    p->stop_when               = 0;  /* ANY */
    p->sor_enabled             = 1;
    p->sor_auto                = 1;
    p->sor_omega_init          = 1.0;
    p->sor_omega_min           = 0.3;
    p->sor_omega_fixed         = -1.0;
    p->sor_burnin              = 20;
}

void rk_result_init(rk_result_t* r) {
    if (!r) return;
    memset(r, 0, sizeof(*r));
    r->best_error    = std::numeric_limits<double>::infinity();  /* Inf sentinel; R sees Inf not finite 1e308 */
    r->sor_min_omega = 1.0;    /* non-iEPPA default */
}

static int validate_inputs(int n, int K,
                            const double* weights,
                            const int32_t** group_ids,
                            const int* cat_counts,
                            const double** targets,
                            const rk_params_t* p,
                            rk_result_t* result,
                            rk_algorithm_t alg) {
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

    if (!std::isfinite(p->tol_abs) || p->tol_abs <= 0.0)
        return err("tol_abs must be finite and positive");

    /* Enum range guards: prevent silent UB from out-of-range static_cast */
    if (p->criterion < 0 || p->criterion > 4)
        return err("criterion out of range [0,4]: 0=PCT 1=ABS 2=KL 3=MEAN_ABS 4=CHI2");
    if (p->stop_when < 0 || p->stop_when > 1)
        return err("stop_when out of range [0,1]: 0=ANY 1=ALL");

    return RK_OK;
}


// Return contract: on RK_OK or RK_ERR_NOCONV, the output weight vector
// satisfies Σ weights[i] = n (where n = number of observations). Solvers
// enforce this invariant internally at exit (see src/ieppa.cpp,
// src/raking.cpp, src/lbfgsb_solver.cpp post-exit normalize blocks).
// Third-party callers should NOT apply their own sum/mean normalization —
// doing so silently invalidates bounds_mode="unit" strict-bounds guarantees,
// because ieppa's water-fill clamps would be re-pushed above max_weight.
// On RK_ERR_INFEAS / RK_ERR_BADARG the output weights are undefined.
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
        rk_result_init(result);
    }

    // Resolve algorithm before validation so the singularity guard knows which
    // link function will be used. Guard against null cat_counts or invalid K/n
    // — validate_inputs will reject those cases with a proper error message.
    rk_algorithm_t alg;
    bool auto_selected = false;
    if (cat_counts && K > 0 && n > 0) {
        switch (p->algorithm) {
            case RK_ALG_LBFGSB: alg = RK_ALG_LBFGSB; break;
            case RK_ALG_RAKING: alg = RK_ALG_RAKING; break;
            case RK_ALG_IEPPA:  alg = RK_ALG_IEPPA;  break;
            case RK_ALG_AUTO:
            default:
                alg = RK_ALG_IEPPA;  // AUTO → faithful iEPPA. Benchmark-driven refinement TBD.
                auto_selected = true;
                break;
        }
    } else {
        alg = p->algorithm;
    }

    int rc = validate_inputs(n, K, weights, group_ids, cat_counts, targets, p, result, alg);
    if (rc != RK_OK) return rc;

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
    st.lbfgs_m       = p->lbfgs_m;
    st.verbose       = p->verbose;
    st.ieppa_auto_selected = auto_selected;  // read by ieppa_solve for verbose prefix
    st.bounds_mode   = p->bounds_mode;
    st.log_fn        = p->log_fn;
    st.log_ctx       = p->log_ctx;
    // Thread overlay config
    st.homotopy.n_levels        = p->homotopy.n_levels;
    st.homotopy.start_factor    = p->homotopy.start_factor;
    st.homotopy.end_factor      = p->homotopy.end_factor;
    st.homotopy.budget_split_p  = p->homotopy.budget_split_p;
    st.homotopy.enabled         = (p->homotopy.enabled != 0) || (p->homotopy.n_levels > 1);
    st.scheduler.mode           = (p->scheduler == RK_SCHED_GREEDY)
                                    ? lbw::SchedulerMode::GREEDY
                                    : lbw::SchedulerMode::ROUND_ROBIN;
    st.eta_schedule.mode        = (p->eta_mode == RK_ETA_TANG_DYNAMIC)
                                    ? lbw::EtaScheduleMode::TANG_DYNAMIC
                                    : lbw::EtaScheduleMode::FIXED;
    st.eta_schedule.eta_start     = p->eta_start;
    st.eta_schedule.eta_end       = p->eta_end;
    st.eta_schedule.schedule_power = p->eta_schedule_power;
    st.convergence_cfg.pct_tol      = p->pct_tol;
    st.convergence_cfg.absolute_tol = p->absolute_tol;
    st.convergence_cfg.criterion    = static_cast<lbw::CalibCriterion>(p->criterion);
    st.convergence_cfg.stop_when    = static_cast<lbw::CalibStopWhen>(p->stop_when);
    st.sor_cfg.enabled              = (p->sor_enabled != 0);
    st.sor_cfg.auto_adapt           = (p->sor_auto != 0);
    st.sor_cfg.omega_init           = p->sor_omega_init;
    st.sor_cfg.omega_min            = p->sor_omega_min;
    st.sor_cfg.omega_fixed          = p->sor_omega_fixed;
    st.sor_cfg.burnin               = p->sor_burnin;
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
        if (result) {
            result->mean_error  = res.mean_error;
            result->kl          = res.kl;
            result->chi2        = res.chi2;
            result->pct_change  = res.pct_change;
            result->best_error  = res.best_error;
            result->best_iter   = res.best_iter;
            /* sor_min_omega, sor_n_damped remain at rk_result_init defaults (1.0, 0) */
        }
    } else if (alg == RK_ALG_RAKING) {
        // Classical raking: IPF + Dykstra box + Dykstra hyperplane (renamed from iEPPA)
        auto res = lbw::raking_solve(st);
        status = res.status;
        iterations = res.iterations;
        max_error = res.max_error;
        used = RK_ALG_RAKING;
        if (result) {
            result->mean_error  = res.mean_error;
            result->kl          = res.kl;
            result->chi2        = res.chi2;
            result->pct_change  = res.pct_change;
            result->best_error  = res.best_error;
            result->best_iter   = res.best_iter;
            /* sor_min_omega, sor_n_damped remain at rk_result_init defaults (1.0, 0) */
        }
    } else {
        // Default / IEPPA: paper-faithful algBCD at C=0 (new src/ieppa.cpp)
        auto res = lbw::ieppa_solve(st);
        status = res.status;
        iterations = res.iterations;
        max_error = res.max_error;
        used = RK_ALG_IEPPA;
        if (result) {
            result->n_xcur_writes_per_iter_linear = res.n_xcur_writes_per_iter_linear;
            result->min_alpha_seen  = res.min_alpha_seen;
            result->final_alpha     = res.final_alpha;
            result->n_bounds_violated = res.n_bounds_violated;
            result->n_bounds_clamped  = res.n_bounds_clamped;
            result->homotopy_levels_used  = res.homotopy_levels_used;
            result->homotopy_final_factor = res.homotopy_final_factor;
            result->greedy_sweeps_taken   = res.greedy_sweeps_taken;
            result->eta_final             = res.eta_final;
            result->mean_error    = res.mean_error;
            result->kl            = res.kl;
            result->chi2          = res.chi2;
            result->pct_change    = res.pct_change;
            result->best_error    = res.best_error;
            result->best_iter     = res.best_iter;
            result->sor_min_omega = res.sor_min_omega;
            result->sor_n_damped  = res.sor_n_damped;
        }
    }

    if (result) {
        result->status = status;
        result->iterations = iterations;
        result->max_error = max_error;
        result->algorithm_used = used;
        if (result->message[0] == '\0') {
            const char* name = (used == RK_ALG_LBFGSB) ? "L-BFGS-B"
                             : (used == RK_ALG_RAKING) ? "raking"
                             : "iEPPA";
            snprintf(result->message, 256, "%s: %d iters, max_error=%.2e",
                     name, iterations, max_error);
        }
    }
    return status;
}

} // extern "C"
