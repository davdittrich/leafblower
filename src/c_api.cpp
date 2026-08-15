#include "leafblower.h"
#include "validation.hpp"
#include "types.hpp"
#include "oris.hpp"       // ORIS (Over-Relaxed Iterative Scaling)
#include "raking.hpp"     // renamed hybrid
#include "sinkhorn.hpp"   // KL Bregman Dykstra
#include "greg.hpp"       // Newton QP chi2 (GREG, Deville-Sarnal 1992)
#include "chebyshev.hpp"  // Chebyshev LP-based solver
#include "cell_table.hpp" // estimate_M_cell for AUTO routing
#include "calib_dispatch.hpp" // SC1: shared dispatch table (leafblower-rywn)
#include "greenkhorn.hpp"
#include "logit_calib.hpp"
#include "newton_calib.hpp"
#include "design_effect.hpp"
#include <cstring>
#include <cstdio>
#include <cmath>
#include <climits>
#include <cstdint>
#include <limits>
#include <algorithm>

static_assert(RK_ALG_ORIS == 1 && RK_ALG_ORIS_SOFT == 8,
              "ORIS enum values frozen");

#define LBW_NODISCARD [[nodiscard]]

// Single source of truth: algorithm index → display name.
// Indices match rk_algorithm_t enum values in leafblower.h.
// Gaps (2, 7) use "(reserved)" so out-of-enum values are visible.
static constexpr const char* kAlgNames[] = {
    "auto",         // 0  RK_ALG_AUTO
    "ORIS",         // 1  RK_ALG_ORIS
    "(reserved)",   // 2  removed LBFGSB
    "raking",       // 3  RK_ALG_RAKING
    "sinkhorn",     // 4  RK_ALG_SINKHORN
    "chebyshev",    // 5  RK_ALG_CHEBYSHEV
    "greg",         // 6  RK_ALG_GREG
    "(gap)",        // 7  unused gap between GREG=6 and ORIS_SOFT=8
    "ORIS-soft",    // 8  RK_ALG_ORIS_SOFT
    "greenkhorn",   // 9  RK_ALG_GREENKHORN
    "logit",        // 10 RK_ALG_LOGIT
    "newton_kl"     // 11 RK_ALG_NEWTON_KL
};
static_assert(RK_ALG_NEWTON_KL == 11,
    "kAlgNames: update if enum max changes");
static_assert(sizeof(kAlgNames)/sizeof(kAlgNames[0]) == 12,
    "kAlgNames must cover all 12 enum slots 0..11");

// CR-D11 (j7x8.11): void_t detection for a top-level n_bounds_violated field, so
// the shared pack template surfaces the bound-violation diagnostic (parity with
// r_bridge's has_n_bounds) instead of leaving Python/C callers at a stale 0 for
// raking/greg/sinkhorn under the no-clamp cell contract.
template <class T, class = void> struct has_n_bounds_c : std::false_type {};
template <class T>
struct has_n_bounds_c<T, std::void_t<decltype(std::declval<T&>().n_bounds_violated)>> : std::true_type {};

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
    dst->message[sizeof(dst->message) - 1] = '\0';  // CXX.4: strncpy may not NUL-terminate
    // CR-D11 (j7x8.11): surface the bound-violation diagnostic for results that
    // carry it (raking/greg/sinkhorn now do); oris keeps its manual set below.
    if constexpr (has_n_bounds_c<R>::value) {
        dst->n_bounds_violated = src.n_bounds_violated;
        dst->n_bounds_clamped  = src.n_bounds_clamped;
    }
}

// xc1s.1: shared ORIS result-field pack. The ORIS and ORIS_SOFT dispatch
// branches set an identical 38-field block (incl. the alm_* diagnostics).
// ORISResult has no `message`, so this cannot use the pack_solver_result
// template (which unconditionally copies src.message); the common tail
// synthesizes the message when it is left empty.
static void pack_oris_result_c(rk_result_t* dst, const lbw::ORISResult& res) noexcept {
    if (!dst) return;
    dst->n_xcur_writes_per_iter_last = res.n_xcur_writes_per_iter_last;
    dst->min_alpha_seen  = res.min_alpha_seen;
    dst->final_alpha     = res.final_alpha;
    dst->n_bounds_violated = res.n_bounds_violated;
    dst->n_bounds_clamped  = res.n_bounds_clamped;
    dst->homotopy_levels_used  = res.homotopy_levels_used;
    dst->homotopy_final_factor = res.homotopy_final_factor;
    dst->greedy_sweeps_taken   = res.greedy_sweeps_taken;
    dst->eta_final             = res.eta_final;
    dst->mean_error          = res.base.mean_error;
    dst->kl                  = res.base.kl;
    dst->chi2                = res.base.chi2;
    dst->l1_weight_change    = res.base.l1_weight_change;
    dst->grake_norm          = res.base.grake_norm;
    dst->convergence_metric  = res.base.convergence_metric;
    dst->convergence_rule    = res.base.convergence_rule;
    dst->convergence_tol     = res.base.convergence_tol;
    dst->convergence_iter                = res.base.convergence_iter;
    dst->convergence_solver_objective    = res.base.convergence_solver_objective;
    dst->convergence_minimized_metric    = res.base.convergence_minimized_metric;
    dst->best_error          = res.base.best_error;
    dst->best_iter           = res.base.best_iter;
    dst->metric_first_check  = res.base.metric_first_check;
    dst->metric_prev_check   = res.base.metric_prev_check;
    dst->prev_check_iter     = res.base.prev_check_iter;
    dst->sor_min_omega       = res.sor_min_omega;
    dst->sor_n_damped        = res.sor_n_damped;
    dst->sor_omega_mean      = res.sor_omega_mean;
    dst->sor_any_latched     = res.sor_any_latched;
    dst->sor_n_pinned_fb     = res.sor_n_pinned_fb;
    dst->sor_n_warmup_fb     = res.sor_n_warmup_fb;
    dst->sor_n_conv_fb       = res.sor_n_conv_fb;
    dst->sor_n_resid_grew    = res.sor_n_resid_grew;
    dst->sor_n_monotone_cd   = res.sor_n_monotone_cd;
    dst->alm_capacity_mu_final = res.alm_capacity_mu_final;
    dst->alm_n_growth_events   = res.alm_n_growth_events;
    dst->alm_max_dual_norm     = res.alm_max_dual_norm;
    dst->alm_sum_drift         = res.alm_sum_drift;
}

// CR-D6 (j7x8.6): single newton_kl result pack used by BOTH the explicit
// RK_ALG_NEWTON_KL branch and the AUTO-fallback branch, so newton_kl produces
// identical fields regardless of how it is reached (also closes xc1s.1.1's
// explicit-vs-fallback base-field divergence). Base fields go through the shared
// pack; newton's TSVD/LM diagnostics have no rk_result_t counterpart (see inline
// note below); and the ORIS-only fields are reset to documented non-ORIS defaults
// (leafblower.h: 1.0/0) so
// a fallback after an ORIS primary does not surface stale sor_*/alm_*/homotopy_*.
// The 1.0/0 defaults also match r_bridge's explicit-newton output → R/Python parity.
template <typename R>
static void pack_newton_result_c(rk_result_t* dst, const R& nkr, rk_algorithm_t alg) noexcept {
    if (!dst) return;
    pack_solver_result(dst, nkr, alg);   // base fields + message + n_bounds (trait)
    // (n_projected_dims / lm_mu_final are NewtonCalibResult-only; rk_result_t does
    //  not carry them, so nothing newton-specific to copy into the C result here.)
    // ORIS-only diagnostics → documented non-ORIS defaults
    dst->n_xcur_writes_per_iter_last = 0;
    dst->min_alpha_seen        = 1.0;
    dst->final_alpha           = 1.0;
    dst->homotopy_levels_used  = 0;
    dst->homotopy_final_factor = 1.0;
    dst->greedy_sweeps_taken   = 0;
    dst->eta_final             = 0.0;
    dst->sor_min_omega     = 1.0;
    dst->sor_n_damped      = 0;
    dst->sor_omega_mean    = 1.0;
    dst->sor_any_latched   = 0;
    dst->sor_n_pinned_fb   = 0;
    dst->sor_n_warmup_fb   = 0;
    dst->sor_n_conv_fb     = 0;
    dst->sor_n_resid_grew  = 0;
    dst->sor_n_monotone_cd = 0;
    dst->alm_capacity_mu_final = 0.0;
    dst->alm_n_growth_events   = 0;
    dst->alm_max_dual_norm     = 0.0;
    dst->alm_sum_drift         = 0.0;
}

// SC1 (leafblower-rywn): narrow the shared, unconstrained lbw::DispatchResult
// (calib_dispatch.hpp) into the ABI-frozen rk_result_t. The 4 superset-only
// fields (n_projected_dims, lm_mu_final, sraa_demoted, plus the obs-level
// best_weights vector) have no rk_result_t counterpart and are deliberately
// NOT copied here (matches Python's existing, already-accepted behavior).
static void pack_dispatch_result_c(rk_result_t* dst, const lbw::DispatchResult& dres) noexcept {
    if (!dst) return;
    dst->status                       = dres.status;
    dst->iterations                   = dres.iterations;
    dst->max_error                    = dres.max_error;
    dst->algorithm_used               = dres.alg_used;
    dst->mean_error                   = dres.mean_error;
    dst->kl                           = dres.kl;
    dst->chi2                         = dres.chi2;
    dst->l1_weight_change             = dres.l1_weight_change;
    dst->grake_norm                   = dres.grake_norm;
    dst->convergence_metric           = dres.convergence_metric;
    dst->convergence_rule             = dres.convergence_rule;
    dst->convergence_tol              = dres.convergence_tol;
    dst->convergence_iter             = dres.convergence_iter;
    dst->convergence_solver_objective = dres.convergence_solver_objective;
    dst->convergence_minimized_metric = dres.convergence_minimized_metric;
    dst->best_error                   = dres.best_error;
    dst->best_iter                    = dres.best_iter;
    dst->metric_first_check           = dres.metric_first_check;
    dst->metric_prev_check            = dres.metric_prev_check;
    dst->prev_check_iter              = dres.prev_check_iter;
    dst->n_bounds_violated            = dres.n_bounds_violated;
    dst->n_bounds_clamped             = dres.n_bounds_clamped;
    std::strncpy(dst->message, dres.solver_message, sizeof(dst->message) - 1);
    dst->message[sizeof(dst->message) - 1] = '\0';
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
    p->sor_omega_max           = 1.5;
    p->sor_burnin              = 20;
    p->sor_omega_mode_id       = 2;   /* default: iterate-change ||dX_free||^2 (e18t.9 SHIP) */
    p->_pad_aa_33kz            = 0;
    p->newton_tsvd_ratio       = 1e-8;  /* Epic-H WH-e: newton_kl TSVD truncation default */
    p->ridge_lambda            = 0.0;   /* Tikhonov ridge on dual λ; 0.0 = off */
}

void rk_result_init(rk_result_t* r) {
    if (!r) return;
    memset(r, 0, sizeof(*r));
    r->best_error         = std::numeric_limits<double>::infinity();  /* Inf sentinel; R sees Inf not finite 1e308 */
    r->convergence_solver_objective = std::numeric_limits<double>::infinity();  /* Inf sentinel, consistent with best_error */
    r->sor_min_omega      = 1.0;    /* non-ORIS default */
    r->sor_omega_mean     = 1.0;    /* non-ORIS default */
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
// enforce this invariant internally at exit (see src/oris.cpp,
// src/raking.cpp post-exit normalize blocks).
// Third-party callers should NOT apply their own sum/mean normalization —
// doing so silently invalidates bounds_mode="unit" strict-bounds guarantees,
// because oris's water-fill clamps would be re-pushed above max_weight.
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
    // link function will be used. Guard against null cat_counts/group_ids/targets
    // or invalid K/n — AUTO routing dereferences group_ids (estimate_M_cell) and
    // targets[k][j] before validate_inputs runs (CR-D1/j7x8.1: direct C-ABI callers
    // with algorithm=AUTO + NULL pointers segfaulted here; pybind pre-validates but
    // other FFI callers do not). On a failed guard we fall through to the else and
    // let validate_inputs reject the NULL/invalid inputs with RK_ERR_BADARG.
    rk_algorithm_t alg;
    bool auto_selected = false;
    bool wh_g_severe_skew_accelerate = false;  // Epic-H WH-g: AUTO target-skew gate
    if (cat_counts && group_ids && targets && K > 0 && n > 0) {
        switch (p->algorithm) {
            case RK_ALG_RAKING:   alg = RK_ALG_RAKING; break;
            case RK_ALG_ORIS:      alg = RK_ALG_ORIS;      break;
            case RK_ALG_ORIS_SOFT: alg = RK_ALG_ORIS_SOFT; break;
            case RK_ALG_SINKHORN:   alg = RK_ALG_SINKHORN;   break;
            case RK_ALG_GREG:    alg = RK_ALG_GREG; break;
            case RK_ALG_CHEBYSHEV:   alg = RK_ALG_CHEBYSHEV;   break;
            case RK_ALG_GREENKHORN:  alg = RK_ALG_GREENKHORN;  break;
            case RK_ALG_LOGIT:       alg = RK_ALG_LOGIT;       break;
            case RK_ALG_NEWTON_KL:   alg = RK_ALG_NEWTON_KL;   break;
            case RK_ALG_AUTO:
            default: {
                // Route based on cell table compression ratio, dimension, and target skew.
                // K<5 OR M_cell/n<0.9   → RK_ALG_ORIS / RK_ALG_RAKING (unchanged)
                // K≥5, M_cell/n≥0.9, target_skew ≤ 5 → RK_ALG_NEWTON_KL (moderate skew)
                // K≥5, M_cell/n≥0.9, target_skew > 5 → RK_ALG_ORIS + accelerate
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
                            // on severe-skew dual landscape; ORIS+SRAA converges instead.
                            alg = RK_ALG_ORIS;
                            wh_g_severe_skew_accelerate = true;
                        } else {
                            alg = RK_ALG_NEWTON_KL;
                        }
                    } else {
                        alg = RK_ALG_RAKING;
                    }
                } else {
                    alg = RK_ALG_ORIS;
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

    // j7x8.21: reject non-finite / non-positive homotopy factors at the c_api boundary
    // (mirrors the r_bridge .Call check). A NaN/Inf factor flows into the geometric
    // schedule k_start*pow(k_end/k_start,frac) → NaN current_max_weight; a NaN
    // start_factor also silently disables the CR-D20 tightest-level guard.
    if (!std::isfinite(p->homotopy.start_factor) || p->homotopy.start_factor <= 0.0 ||
        !std::isfinite(p->homotopy.end_factor)   || p->homotopy.end_factor   <= 0.0) {
        if (result)
            std::snprintf(result->message, sizeof(result->message),
                "homotopy start/end factor must be finite and > 0");
        return RK_ERR_BADARG;
    }

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
    st.oris_auto_selected = auto_selected;  // read by oris_solve for verbose prefix
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
    st.sor_cfg.omega_max            = p->sor_omega_max;
    st.sor_cfg.burnin               = p->sor_burnin;
    st.sor_cfg.omega_mode_id        = p->sor_omega_mode_id;
    st.newton_tsvd_ratio            = p->newton_tsvd_ratio;  /* Epic-H WH-e */
    if (wh_g_severe_skew_accelerate) {
        // Epic-H WH-g: severe-skew K≥5 AUTO routes to oris with SRAA enabled.
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
        // SC1 (leafblower-rywn): routed through the shared dispatch table
        // instead of calling the raking solver + a manual field copy
        // directly. sraa_demoted is a superset-only DispatchResult field
        // with no rk_result_t counterpart (pack_dispatch_result_c) — Python's
        // result dict is unaffected, exactly as before this migration.
        lbw::DispatchResult dres_rk;
        lbw::dispatch_solver(alg, st, dres_rk);
        pack_dispatch_result_c(result, dres_rk);
        status = dres_rk.status;
        iterations = dres_rk.iterations;
        max_error = dres_rk.max_error;
        used = RK_ALG_RAKING;
    } else if (alg == RK_ALG_SINKHORN) {
        // SC1 (leafblower-rywn): routed through the shared dispatch table
        // instead of calling lbw::sinkhorn_solve + pack_solver_result directly.
        lbw::DispatchResult dres;
        lbw::dispatch_solver(alg, st, dres);
        pack_dispatch_result_c(result, dres);
        return dres.status;
    } else if (alg == RK_ALG_GREG) {
        // SC1 (leafblower-rywn): routed through the shared dispatch table
        // instead of calling lbw::greg_solve + pack_solver_result directly.
        lbw::DispatchResult dres;
        lbw::dispatch_solver(alg, st, dres);
        pack_dispatch_result_c(result, dres);
        return dres.status;
    } else if (alg == RK_ALG_GREENKHORN) {
        /* Direct C API callers bypass R-layer validation.
           Caller must ensure min_weight < max_weight. */
        // SC1 (leafblower-rywn): routed through the shared dispatch table
        // instead of calling lbw::greenkhorn_solve + pack_solver_result directly.
        lbw::DispatchResult dres_gk;
        lbw::dispatch_solver(alg, st, dres_gk);
        // Solver stores calibrated weights only in best_weights, not in st.weights.
        // Copy to caller buffer on all exits (RK_OK, NOCONV, BUDGET); mirrors
        // r_bridge.cpp's centralized greenkhorn/logit weights copy-back after
        // the solver dispatch chain.
        if (!dres_gk.best_weights.empty() &&
            static_cast<int>(dres_gk.best_weights.size()) == n)
            std::copy(dres_gk.best_weights.begin(),
                      dres_gk.best_weights.end(), weights);
        pack_dispatch_result_c(result, dres_gk);
        used = RK_ALG_GREENKHORN;
        status = dres_gk.status;
        iterations = dres_gk.iterations;
        max_error = dres_gk.max_error;
    } else if (alg == RK_ALG_LOGIT) {
        /* Direct C API callers bypass R-layer validation.
           Caller must ensure max_weight is finite and > min_weight. */
        // SC1 (leafblower-rywn): routed through the shared dispatch table
        // instead of calling lbw::logit_calibrate + pack_solver_result directly.
        lbw::DispatchResult dres_lg;
        lbw::dispatch_solver(alg, st, dres_lg);
        // Solver stores calibrated weights only in best_weights, not in st.weights.
        // Copy to caller buffer on all exits (RK_OK, NOCONV, BUDGET); mirrors
        // r_bridge.cpp's centralized greenkhorn/logit weights copy-back after
        // the solver dispatch chain.
        if (!dres_lg.best_weights.empty() &&
            static_cast<int>(dres_lg.best_weights.size()) == n)
            std::copy(dres_lg.best_weights.begin(),
                      dres_lg.best_weights.end(), weights);
        pack_dispatch_result_c(result, dres_lg);
        used = RK_ALG_LOGIT;
        status = dres_lg.status;
        iterations = dres_lg.iterations;
        max_error = dres_lg.max_error;
    } else if (alg == RK_ALG_NEWTON_KL) {
        auto nkr = lbw::newton_calibrate(st);
        // CR-D6 (j7x8.6): route through the shared newton pack — previously this
        // branch hand-copied a partial field set (omitting mean_error/kl/chi2/
        // l1_weight_change/grake_norm/convergence_solver_objective/
        // convergence_minimized_metric), diverging from the AUTO-fallback branch.
        pack_newton_result_c(result, nkr, RK_ALG_NEWTON_KL);
        for (int i = 0; i < n; i++) weights[i] = st.weights[i];
        return nkr.base.status;
    } else {
        if (alg == RK_ALG_CHEBYSHEV) {
            // kxna.23: reject inner_max_iter < 1 BEFORE the ~5-iter ORIS warm-start
            // below (the shared-table solver arm guards this too, but only AFTER
            // the wasted warm-start solve). Mirror the solver's exact message.
            if (st.inner_max_iter < 1) {
                lbw::ChebyshevResult r;
                r.base.status = RK_ERR_BADARG;
                std::snprintf(r.message, sizeof(r.message),
                    "chebyshev: inner_max_iter (%d) must >= 1", st.inner_max_iter);
                pack_solver_result(result, r, alg);
                return r.base.status;
            }
            // SC1 (leafblower-rywn): routed through the shared dispatch table
            // instead of running the oris warm-start + solver call +
            // pack_solver_result inline (mirrors r_bridge.cpp's "chebyshev"
            // branch). best_weights is deliberately NOT copied into the
            // caller's `weights` buffer here: the solver already mutates
            // st.weights (aliased to `weights`) in place on every path that
            // reaches it, matching sinkhorn/greg/raking's in-place contract;
            // rk_result_t has no best_weights field (pack_dispatch_result_c).
            lbw::DispatchResult dres_cheb;
            lbw::dispatch_solver(alg, st, dres_cheb);
            pack_dispatch_result_c(result, dres_cheb);
            return dres_cheb.status;
        } else if (alg == RK_ALG_ORIS_SOFT) {
            st.use_admm_capacity = true;
            /* capacity_penalty for oris_soft: direct C API callers bypass R-layer validation.
               Contract: p.capacity_penalty <= 0.0 selects auto (M_cell/n from estimate_M_cell);
               positive value is used directly. Callers must validate range externally. */
            if (p->capacity_penalty > 0.0) {
                st.alm.capacity_mu = p->capacity_penalty;
            } else {
                int M_cell_est = lbw::estimate_M_cell(n, K, group_ids, cat_counts);
                st.alm.capacity_mu = (n > 0) ? static_cast<double>(M_cell_est) / n : 1.0;
            }
            auto res = lbw::oris_solve(st);
            status = res.base.status;
            iterations = res.base.iterations;
            max_error = res.base.max_error;
            used = RK_ALG_ORIS_SOFT;
            if (result) pack_oris_result_c(result, res);
        } else {
            // Default / ORIS: paper-faithful algBCD at C=0 (src/oris.cpp)
            auto res = lbw::oris_solve(st);
            status = res.base.status;
            iterations = res.base.iterations;
            max_error = res.base.max_error;
            used = RK_ALG_ORIS;
        if (result) pack_oris_result_c(result, res);
        }
    }

    // Auto-fallback: if primary solver NOCONV or BUDGET, retry with newton_kl.
    // Parity with r_bridge.cpp:558 — both NOCONV and BUDGET indicate "still
    // improving, ran out of iterations" and should fall back. STALL stays put
    // (at constrained optimum, fallback won't help).
    if ((status == RK_ERR_NOCONV || status == RK_ERR_BUDGET)
        && p->algorithm == RK_ALG_AUTO) {
        if (p->verbose >= 1)
            st.log("auto: primary NOCONV; retrying with newton_kl");
        std::copy(weights_backup.begin(), weights_backup.end(), weights);
        auto fb = lbw::newton_calibrate(st);
        status     = fb.base.status;
        iterations = fb.base.iterations;
        max_error  = fb.base.max_error;
        used       = RK_ALG_NEWTON_KL;
        // CR-D6 (j7x8.6): same shared pack as the explicit branch → identical fields,
        // plus it resets the ORIS-only diagnostics the abandoned primary populated via
        // pack_oris_result_c (grake_norm/convergence_solver_objective/
        // convergence_minimized_metric were also previously dropped here).
        pack_newton_result_c(result, fb, RK_ALG_NEWTON_KL);
        for (int i = 0; i < n; i++) weights[i] = st.weights[i];
    }

    if (result) {
        result->status = status;
        result->iterations = iterations;
        result->max_error = max_error;
        result->algorithm_used = used;
        if (result->message[0] == '\0') {
            const int idx = static_cast<int>(used);
            const char* name = (idx >= 0 && idx < static_cast<int>(
                sizeof(kAlgNames)/sizeof(kAlgNames[0])))
                ? kAlgNames[idx] : "unknown";
            snprintf(result->message, sizeof(result->message),  // CXX.4: was literal 256
                     "%s: %d iters, max_error=%.2e",
                     name, iterations, max_error);
        }
    }
    return status;
}

int rk_design_effect(
    const double* weights,
    const double* outcome,
    const int*    data_codes,
    const int*    cat_counts,
    int n, int K,
    rk_design_effect_result_t* out
) {
    if (out == nullptr || weights == nullptr) return RK_ERR_BADARG;
    return lbw::design_effect_compute(weights, outcome, data_codes, cat_counts, n, K, out);
}

} // extern "C"
