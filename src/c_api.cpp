#include "leafblower.h"
#include "validation.hpp"
#include "types.hpp"
#include "lbfgsb_solver.hpp"
#include "ieppa.hpp"      // faithful (new)
#include "raking.hpp"     // renamed hybrid
#include "sinkhorn.hpp"   // KL Bregman Dykstra
#include "greg.hpp"       // Newton QP chi2 (GREG, Deville-Sarnal 1992)
#include "chebyshev.hpp"  // Chebyshev/GRAKE LP-based solvers
#include "cell_table.hpp" // estimate_M_cell for AUTO routing
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

template <typename R>
static void pack_solver_result(rk_result_t* dst, const R& src, rk_algorithm_t alg) noexcept {
    if (!dst) return;
    dst->status                       = src.status;
    dst->iterations                   = src.iterations;
    dst->max_error                    = src.max_error;
    dst->convergence_metric           = src.convergence_metric;
    dst->convergence_rule             = src.convergence_rule;
    dst->convergence_tol              = src.convergence_tol;
    dst->convergence_iter             = src.convergence_iter;
    dst->convergence_solver_objective        = src.convergence_solver_objective;
    dst->convergence_minimized_metric = src.convergence_minimized_metric;
    dst->best_error                   = src.best_error;
    dst->best_iter                    = src.best_iter;
    dst->mean_error                   = src.mean_error;
    dst->kl                           = src.kl;
    dst->chi2                         = src.chi2;
    dst->grake_norm                   = src.grake_norm;
    dst->l1_weight_change             = src.l1_weight_change;
    dst->algorithm_used               = alg;
    std::strncpy(dst->message, src.message, sizeof(dst->message) - 1);
}

static void pack_lbfgsb_result(rk_result_t* dst, const lbw::LBFGSResult& src) noexcept {
    if (!dst) return;
    dst->mean_error                   = src.mean_error;
    dst->kl                           = src.kl;
    dst->chi2                         = src.chi2;
    dst->l1_weight_change             = src.l1_weight_change;
    dst->grake_norm                   = src.grake_norm;
    dst->convergence_metric           = src.convergence_metric;
    dst->convergence_rule             = src.convergence_rule;
    dst->convergence_tol              = src.convergence_tol;
    dst->convergence_iter             = src.convergence_iter;
    dst->convergence_solver_objective        = src.convergence_solver_objective;
    dst->convergence_minimized_metric = src.convergence_minimized_metric;
    dst->best_error                   = src.best_error;
    dst->best_iter                    = src.best_iter;
    /* sor_min_omega, sor_n_damped remain at rk_result_init defaults (1.0, 0) */
}

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
    p->metric                  = 0;  /* MAX_ERR */
    p->rule                    = 1;  /* IMPROVEMENT */
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
    r->best_error         = std::numeric_limits<double>::infinity();  /* Inf sentinel; R sees Inf not finite 1e308 */
    r->convergence_solver_objective = std::numeric_limits<double>::infinity();  /* Inf sentinel, consistent with best_error */
    r->sor_min_omega      = 1.0;    /* non-iEPPA default */
    r->convergence_rule                 = 1;      /* IMPROVEMENT */
    r->convergence_tol                  = 0.001;
    r->convergence_iter                 = -1;     /* -1 = did not converge */
    r->convergence_minimized_metric     = 0;
}

static int validate_inputs(int n, int K,
                            const double* weights,
                            const int32_t** group_ids,
                            const int* cat_counts,
                            const double** targets,
                            const rk_params_t* p,
                            rk_result_t* result,
                            rk_algorithm_t alg) {
    return lbw::validate_calibrate_inputs(
        n, K, weights, group_ids, cat_counts, targets, p, result, alg);
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
            case RK_ALG_LBFGSB:   alg = RK_ALG_LBFGSB; break;
            case RK_ALG_RAKING:   alg = RK_ALG_RAKING; break;
            case RK_ALG_IEPPA:    alg = RK_ALG_IEPPA;  break;
            case RK_ALG_SINKHORN: alg = RK_ALG_SINKHORN; break;
            case RK_ALG_GREG:    alg = RK_ALG_GREG; break;
            case RK_ALG_CHEBYSHEV: alg = RK_ALG_CHEBYSHEV; break;
            case RK_ALG_GRAKE:    alg = RK_ALG_GRAKE; break;
            case RK_ALG_AUTO:
            default: {
                // Route to raking when cell table is nearly incompressible (M_cell/n > 0.9).
                // At high ratios iEPPA has no compression benefit; raking is equivalent and simpler.
                int M_cell_est = lbw::estimate_M_cell(n, K, group_ids, cat_counts);
                // Exact integer comparison: M_cell_est / n > 0.9  ↔  M_cell_est * 10 > n * 9
                alg = (static_cast<int64_t>(M_cell_est) * 10 > static_cast<int64_t>(n) * 9)
                      ? RK_ALG_RAKING : RK_ALG_IEPPA;
                auto_selected = true;
                break;
            }
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
    st.convergence_cfg.metric       = static_cast<lbw::CalibMetric>(p->metric);
    st.convergence_cfg.rule         = static_cast<lbw::CalibRule>(p->rule);
    st.convergence_cfg.stop_when    = static_cast<lbw::CalibStopWhen>(p->stop_when);
    st.sor_cfg.enabled              = (p->sor_enabled != 0);
    st.sor_cfg.auto_adapt           = (p->sor_auto != 0);
    st.sor_cfg.omega_init           = p->sor_omega_init;
    st.sor_cfg.omega_min            = p->sor_omega_min;
    st.sor_cfg.omega_fixed          = p->sor_omega_fixed;
    st.sor_cfg.burnin               = p->sor_burnin;
    // Only the auto-fallback path needs this; skip O(n) copy for explicit method calls.
    const std::vector<double> weights_backup = (p->algorithm == RK_ALG_AUTO)
        ? std::vector<double>(weights, weights + n)
        : std::vector<double>();
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
        pack_lbfgsb_result(result, res);
    } else if (alg == RK_ALG_RAKING) {
        // Classical raking: IPF + Dykstra box + Dykstra hyperplane (renamed from iEPPA)
        auto res = lbw::raking_solve(st);
        status = res.status;
        iterations = res.iterations;
        max_error = res.max_error;
        used = RK_ALG_RAKING;
        if (result) {
            result->mean_error          = res.mean_error;
            result->kl                  = res.kl;
            result->chi2                = res.chi2;
            result->l1_weight_change    = res.l1_weight_change;
            result->grake_norm          = res.grake_norm;
            result->convergence_metric  = res.convergence_metric;
            result->convergence_rule    = res.convergence_rule;
            result->convergence_tol     = res.convergence_tol;
            result->convergence_iter                = res.convergence_iter;
            result->convergence_solver_objective           = res.convergence_solver_objective;
            result->convergence_minimized_metric    = res.convergence_minimized_metric;
            result->best_error          = res.best_error;
            result->best_iter           = res.best_iter;
            /* sor_min_omega, sor_n_damped remain at rk_result_init defaults (1.0, 0) */
        }
    } else if (alg == RK_ALG_SINKHORN) {
        auto sres = lbw::sinkhorn_solve(st);
        pack_solver_result(result, sres, alg);
        return sres.status;
    } else if (alg == RK_ALG_GREG) {
        auto gres = lbw::greg_solve(st);
        pack_solver_result(result, gres, alg);
        return gres.status;
    } else {
        if (alg == RK_ALG_CHEBYSHEV) {
            auto r = lbw::chebyshev_ipm(st, lbw::LpVariant::CHEBYSHEV);
            pack_solver_result(result, r, alg);
            return r.status;
        } else if (alg == RK_ALG_GRAKE) {
            auto r = lbw::chebyshev_ipm(st, lbw::LpVariant::GRAKE);
            pack_solver_result(result, r, alg);
            return r.status;
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
            result->mean_error          = res.mean_error;
            result->kl                  = res.kl;
            result->chi2                = res.chi2;
            result->l1_weight_change    = res.l1_weight_change;
            result->grake_norm          = res.grake_norm;
            result->convergence_metric  = res.convergence_metric;
            result->convergence_rule    = res.convergence_rule;
            result->convergence_tol     = res.convergence_tol;
            result->convergence_iter                = res.convergence_iter;
            result->convergence_solver_objective           = res.convergence_solver_objective;
            result->convergence_minimized_metric    = res.convergence_minimized_metric;
            result->best_error          = res.best_error;
            result->best_iter           = res.best_iter;
            result->sor_min_omega       = res.sor_min_omega;
            result->sor_n_damped        = res.sor_n_damped;
        }
        }
    }

    // Auto-fallback: if primary solver NOCONVs, retry with L-BFGS-B
    if (status == RK_ERR_NOCONV && p->algorithm == RK_ALG_AUTO) {
        if (p->verbose >= 1)
            st.log("auto: primary solver NOCONV; retrying with L-BFGS-B");
        // Restore original weights (only mutated field in CalibState)
        std::copy(weights_backup.begin(), weights_backup.end(), weights);
        auto fb = lbw::lbfgsb_solve(st);
        status     = fb.status;
        iterations = fb.iterations;
        max_error  = fb.max_error;
        used       = RK_ALG_LBFGSB;
        pack_lbfgsb_result(result, fb);
        /* sor_min_omega, sor_n_damped remain at rk_result_init defaults (1.0, 0) */
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
