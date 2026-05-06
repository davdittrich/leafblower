#include "leafblower.h"
#include "validation.hpp"
#include "types.hpp"
#include "ieppa.hpp"      // faithful (new)
#include "raking.hpp"     // renamed hybrid
#include "sinkhorn.hpp"   // KL Bregman Dykstra
#include "greg.hpp"       // Newton QP chi2 (GREG, Deville-Sarnal 1992)
#include "chebyshev.hpp"  // Chebyshev LP-based solver
#include "cell_table.hpp" // estimate_M_cell for AUTO routing
#include "greenkhorn.hpp"
#include "logit_calib.hpp"
#include "newton_calib.hpp"
#include <cstring>
#include <cstdio>
#include <cmath>
#include <climits>
#include <cstdint>
#include <limits>
#include <algorithm>

#define LBW_NODISCARD [[nodiscard]]

template <typename R>
static void pack_solver_result(rk_result_t* dst, const R& src, rk_algorithm_t alg) noexcept {
    if (!dst) return;
    dst->status                       = src.base.status;
    dst->iterations                   = src.base.iterations;
    dst->max_error                    = src.base.max_error;
    dst->convergence_metric           = src.base.convergence_metric;
    dst->convergence_rule             = src.base.convergence_rule;
    dst->convergence_tol              = src.base.convergence_tol;
    dst->convergence_iter             = src.base.convergence_iter;
    dst->convergence_solver_objective        = src.base.convergence_solver_objective;
    dst->convergence_minimized_metric = src.base.convergence_minimized_metric;
    dst->best_error                   = src.base.best_error;
    dst->best_iter                    = src.base.best_iter;
    dst->mean_error                   = src.base.mean_error;
    dst->kl                           = src.base.kl;
    dst->chi2                         = src.base.chi2;
    dst->grake_norm                   = src.base.grake_norm;
    dst->l1_weight_change             = src.base.l1_weight_change;
    dst->algorithm_used               = alg;
    dst->metric_first_check           = src.base.metric_first_check;
    dst->metric_prev_check            = src.base.metric_prev_check;
    dst->prev_check_iter              = src.base.prev_check_iter;
    std::strncpy(dst->message, src.message, sizeof(dst->message) - 1);
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
    p->_reserved_lbfgs_m = 0;
    p->log_fn        = nullptr;
    p->log_ctx       = nullptr;
    p->bounds_mode   = RK_BOUNDS_CELL;  /* explicit for clarity; memset=0 already gives this */
    p->homotopy.n_levels       = 1;
    p->homotopy.start_factor   = 1.0;
    p->homotopy.end_factor     = 1.0;
    p->homotopy.budget_split_p = 0.5;
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
    p->newton_tsvd_ratio       = 1e-8;  /* Epic-H WH-e: newton_kl TSVD truncation default */
    p->ridge_lambda            = 0.0;   /* Tikhonov ridge on dual λ; 0.0 = off */
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
    r->metric_first_check  = std::numeric_limits<double>::infinity();  /* Inf sentinel: not yet captured */
    r->metric_prev_check   = std::numeric_limits<double>::infinity();  /* Inf sentinel: not yet captured */
    r->prev_check_iter     = -1;                                       /* -1: not yet captured */
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
// src/raking.cpp post-exit normalize blocks).
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
    bool wh_g_severe_skew_accelerate = false;  // Epic-H WH-g: AUTO target-skew gate
    if (cat_counts && K > 0 && n > 0) {
        switch (p->algorithm) {
            case RK_ALG_RAKING:   alg = RK_ALG_RAKING; break;
            case RK_ALG_IEPPA:      alg = RK_ALG_IEPPA;      break;
            case RK_ALG_IEPPA_SOFT: alg = RK_ALG_IEPPA_SOFT; break;
            case RK_ALG_SINKHORN:   alg = RK_ALG_SINKHORN;   break;
            case RK_ALG_GREG:    alg = RK_ALG_GREG; break;
            case RK_ALG_CHEBYSHEV:   alg = RK_ALG_CHEBYSHEV;   break;
            case RK_ALG_GREENKHORN:  alg = RK_ALG_GREENKHORN;  break;
            case RK_ALG_LOGIT:       alg = RK_ALG_LOGIT;       break;
            case RK_ALG_NEWTON_KL:   alg = RK_ALG_NEWTON_KL;   break;
            case RK_ALG_AUTO:
            default: {
                // Route based on cell table compression ratio, dimension, and target skew.
                // K<5 OR M_cell/n<0.9   → RK_ALG_IEPPA / RK_ALG_RAKING (unchanged)
                // K≥5, M_cell/n≥0.9, target_skew ≤ 5 → RK_ALG_NEWTON_KL (moderate skew)
                // K≥5, M_cell/n≥0.9, target_skew > 5 → RK_ALG_IEPPA + accelerate
                //                                       (Epic-H WH-g: severe skew)
                int M_cell_est = lbw::estimate_M_cell(n, K, group_ids, cat_counts);
                // Exact integer comparison: M_cell_est / n >= 0.9  ↔  M_cell_est * 10 >= n * 9
                if (static_cast<int64_t>(M_cell_est) * 10 >= static_cast<int64_t>(n) * 9) {
                    if (K >= 5) {
                        // Compute target_skew = max_T / max(min_T, 1e-12).
                        double max_target = 0.0, min_target = 1.0;
                        for (int k = 0; k < K; ++k) {
                            for (int j = 0; j < cat_counts[k]; ++j) {
                                const double t = targets[k][j];
                                if (t > max_target) max_target = t;
                                if (t > 1e-12 && t < min_target) min_target = t;
                            }
                        }
                        const double target_skew = max_target / std::max(min_target, 1e-12);
                        const bool severe_skew = (target_skew > 5.0);
                        if (severe_skew) {
                            // Epic-H WH-g: kk1204 K=20 evidence — Newton-KL stalls at gap≈6.24e-2
                            // on severe-skew dual landscape; iEPPA+SRAA converges instead.
                            alg = RK_ALG_IEPPA;
                            wh_g_severe_skew_accelerate = true;
                        } else {
                            alg = RK_ALG_NEWTON_KL;
                        }
                    } else {
                        alg = RK_ALG_RAKING;
                    }
                } else {
                    alg = RK_ALG_IEPPA;
                }
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
    st.newton_tsvd_ratio            = p->newton_tsvd_ratio;  /* Epic-H WH-e */
    if (wh_g_severe_skew_accelerate) {
        // Epic-H WH-g: severe-skew K≥5 AUTO routes to ieppa with SRAA enabled.
        st.accelerate = true;
    }
    if (p->accelerate)                  /* PY-2: explicit opt-in to SRAA */
        st.accelerate = true;
    if (p->alm_penalty > 0.0)           /* PY-2: ALM penalty coefficient */
        st.alm.mu = p->alm_penalty;
    st.ridge_lambda = p->ridge_lambda;  /* Tikhonov ridge on dual λ */
    // Only the auto-fallback path needs this; skip O(n) copy for explicit method calls.
    const std::vector<double> weights_backup = (p->algorithm == RK_ALG_AUTO)
        ? std::vector<double>(weights, weights + n)
        : std::vector<double>();
    int status;
    int iterations;
    double max_error;
    rk_algorithm_t used;

    if (alg == RK_ALG_RAKING) {
        // Classical raking: IPF + Dykstra box + Dykstra hyperplane (renamed from iEPPA)
        auto res = lbw::raking_solve(st);
        status = res.base.status;
        iterations = res.base.iterations;
        max_error = res.base.max_error;
        used = RK_ALG_RAKING;
        if (result) {
            result->mean_error          = res.base.mean_error;
            result->kl                  = res.base.kl;
            result->chi2                = res.base.chi2;
            result->l1_weight_change    = res.base.l1_weight_change;
            result->grake_norm          = res.base.grake_norm;
            result->convergence_metric  = res.base.convergence_metric;
            result->convergence_rule    = res.base.convergence_rule;
            result->convergence_tol     = res.base.convergence_tol;
            result->convergence_iter                = res.base.convergence_iter;
            result->convergence_solver_objective           = res.base.convergence_solver_objective;
            result->convergence_minimized_metric    = res.base.convergence_minimized_metric;
            result->best_error          = res.base.best_error;
            result->best_iter           = res.base.best_iter;
            result->metric_first_check  = res.base.metric_first_check;
            result->metric_prev_check   = res.base.metric_prev_check;
            result->prev_check_iter     = res.base.prev_check_iter;
            /* sor_min_omega, sor_n_damped remain at rk_result_init defaults (1.0, 0) */
        }
    } else if (alg == RK_ALG_SINKHORN) {
        auto sres = lbw::sinkhorn_solve(st);
        pack_solver_result(result, sres, alg);
        return sres.base.status;
    } else if (alg == RK_ALG_GREG) {
        auto gres = lbw::greg_solve(st);
        pack_solver_result(result, gres, alg);
        return gres.base.status;
    } else if (alg == RK_ALG_GREENKHORN) {
        /* Direct C API callers bypass R-layer validation.
           Caller must ensure min_weight < max_weight. */
        auto res = lbw::greenkhorn_solve(st);
        // Solver stores calibrated weights only in best_weights, not in st.weights.
        // Copy to caller buffer on all exits (RK_OK, NOCONV, BUDGET); mirrors r_bridge.cpp:806-810.
        if (!res.base.best_weights.empty() &&
            static_cast<int>(res.base.best_weights.size()) == n)
            std::copy(res.base.best_weights.begin(),
                      res.base.best_weights.end(), weights);
        pack_solver_result(result, res, alg);
        used = RK_ALG_GREENKHORN;
        status = res.base.status;
        iterations = res.base.iterations;
        max_error = res.base.max_error;
    } else if (alg == RK_ALG_LOGIT) {
        /* Direct C API callers bypass R-layer validation.
           Caller must ensure max_weight is finite and > min_weight. */
        auto res = lbw::logit_calibrate(st);
        // Solver stores calibrated weights only in best_weights, not in st.weights.
        // Copy to caller buffer on all exits (RK_OK, NOCONV, BUDGET); mirrors r_bridge.cpp:806-810.
        if (!res.base.best_weights.empty() &&
            static_cast<int>(res.base.best_weights.size()) == n)
            std::copy(res.base.best_weights.begin(),
                      res.base.best_weights.end(), weights);
        pack_solver_result(result, res, alg);
        used = RK_ALG_LOGIT;
        status = res.base.status;
        iterations = res.base.iterations;
        max_error = res.base.max_error;
    } else if (alg == RK_ALG_NEWTON_KL) {
        auto nkr = lbw::newton_calibrate(st);
        if (result) {
            result->status     = nkr.base.status;
            result->iterations = nkr.base.iterations;
            result->max_error  = nkr.base.max_error;
            result->convergence_metric  = nkr.base.convergence_metric;
            result->convergence_rule    = nkr.base.convergence_rule;
            result->convergence_tol     = nkr.base.convergence_tol;
            result->convergence_iter    = nkr.base.convergence_iter;
            result->best_error          = nkr.base.best_error;
            result->best_iter           = nkr.base.best_iter;
            result->algorithm_used      = RK_ALG_NEWTON_KL;
            result->metric_first_check  = nkr.base.metric_first_check;
            result->metric_prev_check   = nkr.base.metric_prev_check;
            result->prev_check_iter     = nkr.base.prev_check_iter;
            std::strncpy(result->message, nkr.message, sizeof(result->message) - 1);
        }
        for (int i = 0; i < n; i++) weights[i] = st.weights[i];
        return nkr.base.status;
    } else {
        if (alg == RK_ALG_CHEBYSHEV) {
            // ieppa warm-start; mirrors r_bridge.cpp:628-657.
            std::vector<double> w_warm;
            double delta_warm = -1.0;
            {   // scoped: weights_copy must not outlive st_warm (dangling ptr)
                std::vector<double> weights_copy(weights, weights + n);
                lbw::CalibState st_warm = st;
                st_warm.weights = weights_copy.data();
                st_warm.inner_max_iter = std::max(50, std::min(100, st.inner_max_iter / 10));
                auto ieppa_res = lbw::ieppa_solve(st_warm);
                if (!ieppa_res.base.best_weights.empty() &&
                    static_cast<int>(ieppa_res.base.best_weights.size()) == n &&
                    std::isfinite(ieppa_res.base.max_error)) {
                    w_warm     = std::move(ieppa_res.base.best_weights);
                    delta_warm = ieppa_res.base.max_error * 1.5;
                }
            }
            auto r = lbw::chebyshev_ipm(st, w_warm, delta_warm);
            pack_solver_result(result, r, alg);
            return r.base.status;
        } else if (alg == RK_ALG_IEPPA_SOFT) {
            st.use_admm_capacity = true;
            /* capacity_penalty for ieppa_soft: direct C API callers bypass R-layer validation.
               Contract: p.capacity_penalty <= 0.0 selects auto (M_cell/n from estimate_M_cell);
               positive value is used directly. Callers must validate range externally. */
            if (p->capacity_penalty > 0.0) {
                st.alm.capacity_mu = p->capacity_penalty;
            } else {
                int M_cell_est = lbw::estimate_M_cell(n, K, group_ids, cat_counts);
                st.alm.capacity_mu = (n > 0) ? static_cast<double>(M_cell_est) / n : 1.0;
            }
            auto res = lbw::ieppa_solve(st);
            status = res.base.status;
            iterations = res.base.iterations;
            max_error = res.base.max_error;
            used = RK_ALG_IEPPA_SOFT;
            if (result) {
                result->n_xcur_writes_per_iter_last = res.n_xcur_writes_per_iter_last;
                result->min_alpha_seen  = res.min_alpha_seen;
                result->final_alpha     = res.final_alpha;
                result->n_bounds_violated = res.n_bounds_violated;
                result->n_bounds_clamped  = res.n_bounds_clamped;
                result->homotopy_levels_used  = res.homotopy_levels_used;
                result->homotopy_final_factor = res.homotopy_final_factor;
                result->greedy_sweeps_taken   = res.greedy_sweeps_taken;
                result->eta_final             = res.eta_final;
                result->mean_error          = res.base.mean_error;
                result->kl                  = res.base.kl;
                result->chi2                = res.base.chi2;
                result->l1_weight_change    = res.base.l1_weight_change;
                result->grake_norm          = res.base.grake_norm;
                result->convergence_metric  = res.base.convergence_metric;
                result->convergence_rule    = res.base.convergence_rule;
                result->convergence_tol     = res.base.convergence_tol;
                result->convergence_iter                = res.base.convergence_iter;
                result->convergence_solver_objective    = res.base.convergence_solver_objective;
                result->convergence_minimized_metric    = res.base.convergence_minimized_metric;
                result->best_error          = res.base.best_error;
                result->best_iter           = res.base.best_iter;
                result->metric_first_check  = res.base.metric_first_check;
                result->metric_prev_check   = res.base.metric_prev_check;
                result->prev_check_iter     = res.base.prev_check_iter;
                result->sor_min_omega       = res.sor_min_omega;
                result->sor_n_damped        = res.sor_n_damped;
                result->alm_capacity_mu_final = res.alm_capacity_mu_final;
                result->alm_n_growth_events   = res.alm_n_growth_events;
                result->alm_max_dual_norm     = res.alm_max_dual_norm;
                result->alm_sum_drift         = res.alm_sum_drift;
            }
        } else {
            // Default / IEPPA: paper-faithful algBCD at C=0 (new src/ieppa.cpp)
            auto res = lbw::ieppa_solve(st);
            status = res.base.status;
            iterations = res.base.iterations;
            max_error = res.base.max_error;
            used = RK_ALG_IEPPA;
        if (result) {
            result->n_xcur_writes_per_iter_last = res.n_xcur_writes_per_iter_last;
            result->min_alpha_seen  = res.min_alpha_seen;
            result->final_alpha     = res.final_alpha;
            result->n_bounds_violated = res.n_bounds_violated;
            result->n_bounds_clamped  = res.n_bounds_clamped;
            result->homotopy_levels_used  = res.homotopy_levels_used;
            result->homotopy_final_factor = res.homotopy_final_factor;
            result->greedy_sweeps_taken   = res.greedy_sweeps_taken;
            result->eta_final             = res.eta_final;
            result->mean_error          = res.base.mean_error;
            result->kl                  = res.base.kl;
            result->chi2                = res.base.chi2;
            result->l1_weight_change    = res.base.l1_weight_change;
            result->grake_norm          = res.base.grake_norm;
            result->convergence_metric  = res.base.convergence_metric;
            result->convergence_rule    = res.base.convergence_rule;
            result->convergence_tol     = res.base.convergence_tol;
            result->convergence_iter                = res.base.convergence_iter;
            result->convergence_solver_objective           = res.base.convergence_solver_objective;
            result->convergence_minimized_metric    = res.base.convergence_minimized_metric;
            result->best_error          = res.base.best_error;
            result->best_iter           = res.base.best_iter;
            result->metric_first_check  = res.base.metric_first_check;
            result->metric_prev_check   = res.base.metric_prev_check;
            result->prev_check_iter     = res.base.prev_check_iter;
            result->sor_min_omega       = res.sor_min_omega;
            result->sor_n_damped        = res.sor_n_damped;
            result->alm_capacity_mu_final = res.alm_capacity_mu_final;
            result->alm_n_growth_events   = res.alm_n_growth_events;
            result->alm_max_dual_norm     = res.alm_max_dual_norm;
            result->alm_sum_drift         = res.alm_sum_drift;
        }
        }
    }

    // Auto-fallback: if primary solver NOCONVs, retry with newton_kl
    if (status == RK_ERR_NOCONV && p->algorithm == RK_ALG_AUTO) {
        if (p->verbose >= 1)
            st.log("auto: primary NOCONV; retrying with newton_kl");
        std::copy(weights_backup.begin(), weights_backup.end(), weights);
        auto fb = lbw::newton_calibrate(st);
        status     = fb.base.status;
        iterations = fb.base.iterations;
        max_error  = fb.base.max_error;
        used       = RK_ALG_NEWTON_KL;
        if (result) {
            result->status            = fb.base.status;
            result->iterations        = fb.base.iterations;
            result->max_error         = fb.base.max_error;
            result->algorithm_used    = RK_ALG_NEWTON_KL;
            result->convergence_metric = fb.base.convergence_metric;
            result->convergence_rule   = fb.base.convergence_rule;
            result->convergence_tol    = fb.base.convergence_tol;
            result->convergence_iter   = fb.base.convergence_iter;
            result->best_error         = fb.base.best_error;
            result->best_iter          = fb.base.best_iter;
            result->metric_first_check = fb.base.metric_first_check;
            result->metric_prev_check  = fb.base.metric_prev_check;
            result->prev_check_iter    = fb.base.prev_check_iter;
            result->mean_error         = fb.base.mean_error;
            result->kl                 = fb.base.kl;
            result->chi2               = fb.base.chi2;
            result->l1_weight_change   = fb.base.l1_weight_change;
            std::strncpy(result->message, fb.message, sizeof(result->message)-1);
        }
        for (int i = 0; i < n; i++) weights[i] = st.weights[i];
    }

    if (result) {
        result->status = status;
        result->iterations = iterations;
        result->max_error = max_error;
        result->algorithm_used = used;
        if (result->message[0] == '\0') {
            const char* name = (used == RK_ALG_RAKING)      ? "raking"
                             : (used == RK_ALG_IEPPA_SOFT)  ? "iEPPA-soft"
                             : (used == RK_ALG_GREENKHORN)  ? "greenkhorn"
                             : (used == RK_ALG_LOGIT)        ? "logit"
                             : (used == RK_ALG_NEWTON_KL)   ? "newton_kl"
                             : "iEPPA";
            snprintf(result->message, 256, "%s: %d iters, max_error=%.2e",
                     name, iterations, max_error);
        }
    }
    return status;
}

} // extern "C"
