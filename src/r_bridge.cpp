#include "leafblower.h"
#include "validation.hpp"
#ifdef __GLIBC__
#include <malloc.h>
#endif
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <cstring>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>
#include "logit.hpp"
#include "cell_table.hpp"
#include "types.hpp"
#include "oris.hpp"
#include "raking.hpp"
#include "sinkhorn.hpp"
#include "greg.hpp"
#include "chebyshev.hpp"
#include "greenkhorn.hpp"
#include "logit_calib.hpp"
#include "newton_calib.hpp"

static inline double scalar_real(SEXP x, const char* name) {
    if (TYPEOF(x) != REALSXP || LENGTH(x) != 1)
        Rf_error("leafblower: '%s' must be a length-1 numeric (got type %s length %d)",
                 name, Rf_type2char(TYPEOF(x)), (int)LENGTH(x));
    return REAL(x)[0];
}
static inline int scalar_int(SEXP x, const char* name) {
    if (TYPEOF(x) != INTSXP || LENGTH(x) != 1)
        Rf_error("leafblower: '%s' must be a length-1 integer (got type %s length %d)",
                 name, Rf_type2char(TYPEOF(x)), (int)LENGTH(x));
    return INTEGER(x)[0];
}


namespace {
const std::unordered_map<std::string_view, rk_algorithm_t> kAlgMap = {
    {"oris",       RK_ALG_ORIS},
    {"oris_soft",  RK_ALG_ORIS_SOFT},
    {"raking",     RK_ALG_RAKING},
    {"greg",       RK_ALG_GREG},
    {"chebyshev",  RK_ALG_CHEBYSHEV},
    {"sinkhorn",   RK_ALG_SINKHORN},
    {"auto",       RK_ALG_AUTO},
    {"greenkhorn", RK_ALG_GREENKHORN},
    {"logit",      RK_ALG_LOGIT},
    {"newton_kl",  RK_ALG_NEWTON_KL},
};
} // anonymous namespace

extern "C" {
SEXP C_logit_F_at_zero(SEXP, SEXP);
SEXP C_logit_range_check(SEXP, SEXP, SEXP);
SEXP C_logit_Hprime_check(SEXP, SEXP, SEXP);
SEXP C_rk_calibrate(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_leafblower_cell_table_probe(SEXP, SEXP);
}

// R log trampoline: forwards CalibState.log() calls to Rprintf
static void r_log_trampoline(const char* msg, void* /*ctx*/) {
    Rprintf("[leafblower] %s\n", msg);
}

// Initialize CalibState from rk_params_t plus pre-built data buffers.
// Replaces the field-by-field copy block; encapsulates the rk_params_t -> CalibState
// projection contract so it lives in one named place.
static void init_calib_state(lbw::CalibState& st,
                             int n, int K,
                             const rk_params_t& p,
                             std::vector<double>& weights,
                             std::vector<const int32_t*>& group_ids,
                             std::vector<int>& cat_counts,
                             std::vector<const double*>& targets) {
    st.n = n; st.K = K;
    st.weights = weights.data();
    st.group_ids = const_cast<const int32_t**>(group_ids.data());
    st.cat_counts = cat_counts.data();
    st.targets = const_cast<const double**>(targets.data());
    st.min_weight    = p.min_weight;
    st.max_weight    = p.max_weight;
    st.tol_abs       = p.tol_abs;
    st.inner_max_iter = p.inner_max_iter;
    st.outer_max_iter = p.outer_max_iter;
    st.verbose       = p.verbose;
    st.bounds_mode   = p.bounds_mode;
    st.log_fn        = p.log_fn;
    st.log_ctx       = p.log_ctx;
    st.homotopy.n_levels        = p.homotopy.n_levels;
    st.homotopy.start_factor    = p.homotopy.start_factor;
    st.homotopy.end_factor      = p.homotopy.end_factor;
    st.homotopy.budget_split_p  = p.homotopy.budget_split_p;
    st.scheduler.mode           = (p.scheduler == RK_SCHED_GREEDY)
                                    ? lbw::SchedulerMode::GREEDY
                                    : lbw::SchedulerMode::ROUND_ROBIN;
    st.eta_schedule.mode        = (p.eta_mode == RK_ETA_TANG_DYNAMIC)
                                    ? lbw::EtaScheduleMode::TANG_DYNAMIC
                                    : lbw::EtaScheduleMode::FIXED;
    st.eta_schedule.eta_start     = p.eta_start;
    st.eta_schedule.eta_end       = p.eta_end;
    st.eta_schedule.schedule_power = p.eta_schedule_power;
    st.convergence_cfg.pct_tol      = p.pct_tol;
    st.convergence_cfg.absolute_tol = p.absolute_tol;
    st.convergence_cfg.metric       = static_cast<lbw::CalibMetric>(p.metric);
    st.convergence_cfg.rule         = static_cast<lbw::CalibRule>(p.rule);
    st.convergence_cfg.stop_when    = static_cast<lbw::CalibStopWhen>(p.stop_when);
    st.sor_cfg.enabled              = (p.sor_enabled != 0);
    st.sor_cfg.auto_adapt           = (p.sor_auto != 0);
    st.sor_cfg.omega_init           = p.sor_omega_init;
    st.sor_cfg.omega_min            = p.sor_omega_min;
    st.sor_cfg.omega_fixed          = p.sor_omega_fixed;
    st.sor_cfg.burnin               = p.sor_burnin;
}

extern "C" {

SEXP C_rk_design_effect(SEXP, SEXP, SEXP, SEXP, SEXP);  // defined below

void R_init_leafblower(DllInfo* dll) {
    static const R_CallMethodDef call_methods[] = {
        {"C_logit_F_at_zero",    (DL_FUNC)&C_logit_F_at_zero,    2},
        {"C_logit_range_check",  (DL_FUNC)&C_logit_range_check,  3},
        {"C_logit_Hprime_check", (DL_FUNC)&C_logit_Hprime_check, 3},
        {"C_rk_calibrate",       (DL_FUNC)&C_rk_calibrate,       37},
        {"C_leafblower_cell_table_probe", (DL_FUNC)&C_leafblower_cell_table_probe, 2},
        {"C_rk_design_effect",           (DL_FUNC)&C_rk_design_effect,           5},
        {NULL, NULL, 0}
    };
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

// Test bridge: return F(0) for given L, U
SEXP C_logit_F_at_zero(SEXP Lsxp, SEXP Usxp) {
    double L = scalar_real(Lsxp, "L");
    double U = scalar_real(Usxp, "U");
    lbw::LinkFn fn(L, U);
    return Rf_ScalarReal(fn.F(0.0));
}

// Test bridge: return F(u) for a vector of u values
SEXP C_logit_range_check(SEXP Lsxp, SEXP Usxp, SEXP usxp) {
    double L = scalar_real(Lsxp, "L");
    double U = scalar_real(Usxp, "U");
    int n = LENGTH(usxp);
    SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
    lbw::LinkFn fn(L, U);
    const double* u = REAL(usxp);
    double* res = REAL(out);
    for (int i = 0; i < n; i++) res[i] = fn.F(u[i]);
    UNPROTECT(1);
    return out;
}

// Test bridge: return |H'(u0) - F(u0)| via numerical diff (step 1e-7)
SEXP C_logit_Hprime_check(SEXP Lsxp, SEXP Usxp, SEXP u0sxp) {
    double L = scalar_real(Lsxp, "L");
    double U = scalar_real(Usxp, "U");
    double u0 = scalar_real(u0sxp, "u0");
    double h = 1e-7;
    lbw::LinkFn fn(L, U);
    double Hprime_numerical = (fn.H(u0 + h) - fn.H(u0 - h)) / (2.0 * h);
    double diff = std::fabs(Hprime_numerical - fn.F(u0));
    return Rf_ScalarReal(diff);
}

// Main calibration bridge.
// group_ids_sexp: VECSXP of K pre-encoded INTSXP (0-indexed, -1=OOV/NA)
// cat_counts_sexp: INTSXP of length K (category counts)
// targets_sexp: VECSXP of K REALSXP (proportions)
// n_obs_sexp: INTSXP scalar (number of observations)
// Returns: list(weights=double[n], result=list(status, iterations, max_error, algorithm_used, message))
SEXP C_rk_calibrate(SEXP group_ids_sexp, SEXP cat_counts_sexp,
                    SEXP targets_sexp,   SEXP n_obs_sexp,
                    SEXP min_weight_sexp, SEXP max_weight_sexp,
                    SEXP method_sexp, SEXP verbose_sexp,
                    SEXP inner_max_iter_sexp, SEXP start_weights_sexp,
                    SEXP capacity_penalty_sexp, SEXP alm_penalty_sexp,
                    SEXP tol_abs_sexp, SEXP bounds_mode_sexp,
                    SEXP homotopy_levels_sexp, SEXP homotopy_start_factor_sexp,
                    SEXP homotopy_end_factor_sexp, SEXP homotopy_budget_p_sexp,
                    SEXP scheduler_sexp, SEXP eta_schedule_sexp,
                    SEXP eta_start_sexp, SEXP eta_end_sexp,
                    SEXP eta_schedule_power_sexp,
                    /* Convergence config (WU-A) */
                    SEXP pct_tol_sexp, SEXP absolute_tol_sexp,
                    SEXP metric_sexp, SEXP rule_sexp, SEXP stop_when_sexp,
                    /* SOR config (WU-A) */
                    SEXP sor_enabled_sexp, SEXP sor_auto_sexp,
                    SEXP sor_omega_init_sexp, SEXP sor_omega_min_sexp,
                    SEXP sor_omega_fixed_sexp, SEXP sor_burnin_sexp,
                    /* SQUAREM */
                    SEXP accelerate_sexp,
                    /* Epic-H WH-e: newton_kl TSVD truncation ratio (default 1e-8 from R layer). */
                    SEXP newton_tsvd_ratio_sexp,
                    /* Tikhonov ridge on dual λ; 0.0 = off. */
                    SEXP ridge_lambda_sexp) {
    int K = LENGTH(group_ids_sexp);
    int n = scalar_int(n_obs_sexp, "n_obs");

    std::string pre_error;

    std::vector<std::vector<int32_t>> gids_storage(K);
    std::vector<const int32_t*> group_ids(K);
    std::vector<int> cat_counts(K);
    std::vector<std::vector<double>> tgt_storage(K);
    std::vector<const double*> targets(K);

    {
        if (TYPEOF(cat_counts_sexp) != INTSXP)
            pre_error = "cat_counts must be an integer vector";
        // CXX.1: bound-check container lengths vs K (=LENGTH(group_ids_sexp))
        // BEFORE any indexed read, so a malformed direct .Call cannot OOB-read.
        else if (LENGTH(cat_counts_sexp) != K)
            pre_error = "cat_counts length != number of margins (K)";
        else if (TYPEOF(targets_sexp) != VECSXP || LENGTH(targets_sexp) < K)
            pre_error = "targets must be a list of length >= number of margins (K)";
        const int* cc = pre_error.empty() ? INTEGER(cat_counts_sexp) : nullptr;
        for (int k = 0; k < K && pre_error.empty(); k++) {
            cat_counts[k] = cc[k];
            SEXP gid_vec = VECTOR_ELT(group_ids_sexp, k);
            if (TYPEOF(gid_vec) != INTSXP) {
                pre_error = "group_ids[[" + std::to_string(k + 1) + "]] must be an integer vector";
                break;
            }
            if (LENGTH(gid_vec) != n) {
                pre_error = "group_ids[[" + std::to_string(k + 1) + "]] length != n";
                break;
            }
            const int* gp = INTEGER(gid_vec);
            gids_storage[k].assign(gp, gp + n);
            group_ids[k] = gids_storage[k].data();
            SEXP tgt_vec = VECTOR_ELT(targets_sexp, k);
            if (TYPEOF(tgt_vec) != REALSXP) {
                pre_error = "targets[[" + std::to_string(k + 1) + "]] must be a numeric vector";
                break;
            }
            // CXX.1: guard the assign range — tp+cat_counts[k] must stay in bounds.
            if (LENGTH(tgt_vec) != cat_counts[k]) {
                pre_error = "targets[[" + std::to_string(k + 1) + "]] length != cat_counts[" + std::to_string(k + 1) + "]";
                break;
            }
            const double* tp = REAL(tgt_vec);
            tgt_storage[k].assign(tp, tp + cat_counts[k]);
            targets[k] = tgt_storage[k].data();
        }
    }

    // Build weights (start_weights already normalized to mean=1 by R layer)
    std::vector<double> weights(n);
    if (pre_error.empty()) {
    if (Rf_isNull(start_weights_sexp)) {
        for (int i = 0; i < n; i++) weights[i] = 1.0;
    } else {
        if (LENGTH(start_weights_sexp) != n) {
            char buf[128];
            std::snprintf(buf, sizeof(buf),
                "leafblower: start_weights length %d != n=%d",
                (int)LENGTH(start_weights_sexp), n);
            pre_error = buf;
        } else {
            const double* sw = REAL(start_weights_sexp);
            for (int i = 0; i < n; i++) weights[i] = sw[i];
        }
    }
    } // if (pre_error.empty())

    // Set calibration params
    rk_params_t p;
    rk_params_init(&p);
    p.min_weight     = scalar_real(min_weight_sexp, "min_weight");
    p.max_weight     = scalar_real(max_weight_sexp, "max_weight");
    p.verbose        = scalar_int(verbose_sexp, "verbose");
    p.inner_max_iter = scalar_int(inner_max_iter_sexp, "max_iter");
    // newton_kl: 2nd-order outer loop defaults to 50 (see newton_calib.cpp:90).
    // Other solvers: outer = inner (outer_max_iter is ignored for non-newton_kl paths).
    if (LENGTH(method_sexp) == 1 && TYPEOF(method_sexp) == STRSXP &&
        strcmp(CHAR(STRING_ELT(method_sexp, 0)), "newton_kl") == 0) {
        p.outer_max_iter = 0;  // triggers C-side 50-default in newton_calib.cpp:90
    } else {
        p.outer_max_iter = scalar_int(inner_max_iter_sexp, "max_iter");
    }
    p.tol_abs        = scalar_real(tol_abs_sexp, "tol_abs");
    p.bounds_mode    = (rk_bounds_mode_t) scalar_int(bounds_mode_sexp, "bounds_mode");
    p.log_fn         = (p.verbose > 0) ? r_log_trampoline : nullptr;
    /* Overlay knobs */
    p.homotopy.n_levels        = scalar_int(homotopy_levels_sexp, "homotopy_levels");
    p.homotopy.start_factor    = scalar_real(homotopy_start_factor_sexp, "homotopy_start_factor");
    p.homotopy.end_factor      = scalar_real(homotopy_end_factor_sexp, "homotopy_end_factor");
    p.homotopy.budget_split_p  = scalar_real(homotopy_budget_p_sexp, "homotopy_budget_p");
    {
        if (pre_error.empty() && (TYPEOF(scheduler_sexp) != STRSXP || LENGTH(scheduler_sexp) != 1))
            pre_error = "scheduler must be a length-1 character string";
        const char* sched_str = pre_error.empty() ? CHAR(STRING_ELT(scheduler_sexp, 0)) : "";
        p.scheduler = (strcmp(sched_str, "greedy") == 0) ? RK_SCHED_GREEDY : RK_SCHED_ROUND_ROBIN;
    }
    {
        if (pre_error.empty() && (TYPEOF(eta_schedule_sexp) != STRSXP || LENGTH(eta_schedule_sexp) != 1))
            pre_error = "eta_schedule must be a length-1 character string";
        const char* eta_str = pre_error.empty() ? CHAR(STRING_ELT(eta_schedule_sexp, 0)) : "";
        p.eta_mode = (strcmp(eta_str, "tang_dynamic") == 0) ? RK_ETA_TANG_DYNAMIC : RK_ETA_FIXED;
    }
    p.eta_start           = scalar_real(eta_start_sexp, "eta_start");
    p.eta_end             = scalar_real(eta_end_sexp, "eta_end");
    p.eta_schedule_power  = scalar_real(eta_schedule_power_sexp, "eta_schedule_power");
    /* Convergence config (WU-A) */
    p.pct_tol             = scalar_real(pct_tol_sexp, "pct_tol");
    p.absolute_tol        = scalar_real(absolute_tol_sexp, "absolute_tol");
    p.metric              = scalar_int(metric_sexp, "metric");
    p.rule                = scalar_int(rule_sexp, "rule");
    p.stop_when           = scalar_int(stop_when_sexp, "stop_when");
    /* SOR config (WU-A) */
    p.sor_enabled         = scalar_int(sor_enabled_sexp, "sor_enabled");
    p.sor_auto            = scalar_int(sor_auto_sexp, "sor_auto");
    p.sor_omega_init      = scalar_real(sor_omega_init_sexp, "sor_omega_init");
    p.sor_omega_min       = scalar_real(sor_omega_min_sexp, "sor_omega_min");
    p.sor_omega_fixed     = scalar_real(sor_omega_fixed_sexp, "sor_omega_fixed");
    p.sor_burnin          = scalar_int(sor_burnin_sexp, "sor_burnin");

    if (pre_error.empty() && LENGTH(method_sexp) != 1)
        pre_error = "method must be a length-1 character string";
    const char* method_str = pre_error.empty() ? CHAR(STRING_ELT(method_sexp, 0)) : "";
    const auto alg_it = pre_error.empty() ? kAlgMap.find(method_str) : kAlgMap.end();
    if (pre_error.empty())
        p.algorithm = (alg_it != kAlgMap.end()) ? alg_it->second : RK_ALG_ORIS;

    // Full input validation — shared with c_api.cpp path via validation.hpp.
    if (pre_error.empty()) {
        rk_result_t validation_result;
        rk_result_init(&validation_result);
        rk_algorithm_t alg_for_validation = (alg_it != kAlgMap.end()) ? alg_it->second : RK_ALG_RAKING;
        int vrc = lbw::validate_calibrate_inputs(
            n, K,
            weights.data(),
            const_cast<const int32_t**>(group_ids.data()),
            cat_counts.data(),
            const_cast<const double**>(targets.data()),
            &p, &validation_result, alg_for_validation);
        if (vrc != RK_OK) {
            char buf[300];
            std::snprintf(buf, sizeof(buf),
                "leafblower: invalid arguments \342\200\224 %s", validation_result.message);
            pre_error = buf;
        }
    }

    // WU-E: call C++ solvers directly (bypasses flat C ABI) to access best_weights vector.
    // Build CalibState mirroring c_api.cpp:rk_calibrate() setup.
    lbw::CalibState st;
    init_calib_state(st, n, K, p, weights, group_ids, cat_counts, targets);
    // Epic-H WH-e: newton_kl TSVD truncation ratio. <=0 falls back to internal default 1e-8.
    st.newton_tsvd_ratio = scalar_real(newton_tsvd_ratio_sexp, "newton_tsvd_ratio");
    // Tikhonov ridge on dual λ; 0.0 = off.
    st.ridge_lambda = scalar_real(ridge_lambda_sexp, "ridge_lambda");
    st.oris_auto_selected          = false;  // R bridge always resolves method explicitly
    st.alm.lambda = 0.0;
    st.alm.mu     = 0.0;
    st.accelerate = (scalar_int(accelerate_sexp, "accelerate") != 0);

    // Resolve capacity_mu for oris_soft: harmonized to estimate_M_cell path matching c_api.cpp.
    // leafblower-yh0l.4: T3 benchmark identified Py as winner on fulldata; differentiator was
    // R using exact M_cell/n via build_cell_table (~0.073 at K=9), Py using estimate_M_cell
    // which caps at n for K>8 (capacity_mu=1.0). ~14x scaling difference drove iter-0 divergence.
    // capacity_penalty <= 0.0 (or NULL) selects auto via estimate_M_cell (matches c_api.cpp:380);
    // positive value is used directly. ALM block is gated by st.use_admm_capacity,
    // so st.alm.capacity_mu is harmless for non-oris_soft callers.
    {
        const double cp_val = Rf_isNull(capacity_penalty_sexp)
            ? -1.0
            : (LENGTH(capacity_penalty_sexp) == 1 ? REAL(capacity_penalty_sexp)[0] : -1.0);
        if (cp_val > 0.0) {
            st.alm.capacity_mu = cp_val;
        } else if (pre_error.empty()) {
            // CXX.1: estimate_M_cell dereferences group_ids/cat_counts; skip when a
            // length-validation error is already pending (the throw at the deferred
            // check below converts pre_error into a graceful R error). Filling these
            // with malformed input would OOB-read.
            int M_cell_est = lbw::estimate_M_cell(n, K,
                const_cast<const int32_t**>(group_ids.data()),
                cat_counts.data());
            st.alm.capacity_mu = (n > 0) ? static_cast<double>(M_cell_est) / n : 1.0;
        }
    }

    // alm_penalty: objective ALM penalty coefficient (st.alm.mu). Positive value activates.
    {
        const double alm_penalty_val = Rf_isNull(alm_penalty_sexp)
            ? -1.0
            : (LENGTH(alm_penalty_sexp) == 1 ? REAL(alm_penalty_sexp)[0] : -1.0);
        st.alm.mu = (alm_penalty_val > 0.0) ? alm_penalty_val : 0.0;
    }

    // Scalar fields mirrored from rk_result_t for compatibility with downstream assembly.
    int    res_status     = RK_ERR_NOCONV;
    int    res_iterations = 0;
    double res_max_error  = 1.0;
    int    res_alg_used   = (int)RK_ALG_ORIS;
    char   res_message[256] = "";
    int    res_n_xcur_writes = 0;
    double res_min_alpha  = 1.0;
    double res_final_alpha = 1.0;
    int    res_n_bounds_violated = 0;
    int    res_n_bounds_clamped  = 0;
    int    res_homotopy_levels_used  = 0;
    double res_homotopy_final_factor = 1.0;
    int    res_greedy_sweeps_taken   = 0;
    double res_eta_final             = 0.0;
    double res_mean_error        = 0.0;
    double res_kl                = 0.0;
    double res_chi2              = 0.0;
    double res_l1_weight_change  = 0.0;
    double res_grake_norm        = 0.0;
    int    res_conv_metric       = 0;
    int    res_conv_rule         = 1;
    double res_conv_tol          = 0.001;
    int    res_conv_iter         = -1;
    double res_best_error        = std::numeric_limits<double>::infinity();
    int    res_best_iter   = 0;
    double res_sor_min_omega = 1.0;
    int    res_sor_n_damped  = 0;
    double res_conv_objective          = 0.0;
    int    res_conv_minimized_metric   = 0;
    /* ALM diagnostics (populated only in oris_soft dispatch; zero elsewhere) */
    double res_alm_capacity_mu_final   = 0.0;
    int    res_alm_n_growth_events     = 0;
    double res_alm_max_dual_norm       = 0.0;
    double res_alm_sum_drift           = 0.0;
    /* Acceleration (SRAA) diagnostic (oris/oris_soft only; zero elsewhere) */
    int    res_aa_accepted_count       = 0;
    /* SRAA scheduler-demotion flag (oris/raking only; FALSE elsewhere).
       TRUE iff accelerate=TRUE AND greedy scheduler was demoted to round_robin. */
    int    res_sraa_demoted            = 0;
    /* Newton-KL TSVD diagnostic (Epic-Dβ WL-1; non-zero only for newton_kl) */
    int    res_n_projected_dims        = 0;
    /* Newton-KL Levenberg-Marquardt diagnostic (newton_kl only; zero elsewhere) */
    double res_lm_mu_final             = 0.0;
    double res_metric_first_check      = std::numeric_limits<double>::infinity();
    double res_metric_prev_check       = std::numeric_limits<double>::infinity();
    int    res_prev_check_iter         = -1;
    /* stall_kind: 0=no stall, 1=wchange (SRAA path), 2=kl (plain-IPF path).
       Populated by the solver at the RK_ERR_STALL emission site (leafblower-8eod). */
    int    res_stall_kind              = 0;
    std::vector<double> res_best_weights;  // obs-level, length n

    // DRY helper: pack the 8 convergence-diagnostic fields shared by all solvers.
    auto pack_solver_result = [&](const auto& res) {
        res_l1_weight_change         = res.base.l1_weight_change;
        res_grake_norm               = res.base.grake_norm;
        res_conv_metric              = res.base.convergence_metric;
        res_conv_rule                = res.base.convergence_rule;
        res_conv_tol                 = res.base.convergence_tol;
        res_conv_iter                = res.base.convergence_iter;
        res_conv_objective           = res.base.convergence_solver_objective;
        res_conv_minimized_metric    = res.base.convergence_minimized_metric;
        res_mean_error               = res.base.mean_error;
        res_kl                       = res.base.kl;
        res_chi2                     = res.base.chi2;
        res_best_error               = res.base.best_error;
        res_best_iter                = res.base.best_iter;
        res_metric_first_check       = res.base.metric_first_check;
        res_metric_prev_check        = res.base.metric_prev_check;
        res_prev_check_iter          = res.base.prev_check_iter;
        res_stall_kind               = res.base.stall_kind;
    };

    std::string solver_error;
    {  // Rf_error longjmp skips C++ dtors; destroy all RAII objects before calling it (R-exts §5.5)
    try {
    // A1b: pre-dispatch validation errors are deferred to here so all std::vectors
    // above are still in scope; throw converts them into solver_error and they are
    // destroyed at '}' below before Rf_error fires (R-exts §5.5).
    if (!pre_error.empty()) throw std::runtime_error(pre_error);
#ifdef __GLIBC__
    malloc_trim(0);  // compact heap before solver — reduces fragmentation-induced LLC cache misses
#endif
    if (strcmp(method_str, "raking") == 0) {
        auto res = lbw::raking_solve(st);
        res_status     = res.base.status;
        res_iterations = res.base.iterations;
        res_max_error  = res.base.max_error;
        res_alg_used   = (int)RK_ALG_RAKING;
        pack_solver_result(res);
        res_sraa_demoted = res.sraa_demoted ? 1 : 0;
        res_best_weights = std::move(res.base.best_weights);
    } else if (strcmp(method_str, "auto") == 0) {
        // AUTO routing (Epic-H WH-g):
        //   K<5 OR M_cell/n ≤ 0.9 → raking / ORIS (compression-based, unchanged)
        //   K≥5, M_cell/n > 0.9, target_skew ≤ 5 → newton_kl  (moderate skew)
        //   K≥5, M_cell/n > 0.9, target_skew  > 5 → oris+SRAA (severe skew)
        //   target_skew = max(targets) / max(min(targets), 1e-12)
        int M_cell_est = lbw::estimate_M_cell(n, K,
            const_cast<const int32_t**>(group_ids.data()),
            cat_counts.data());
        // Exact integer comparison: M_cell_est / n >= 0.9  ↔  M_cell_est * 10 >= n * 9
        // PAR.1: use >= to match c_api.cpp:190; > diverged at the exact 0.9 boundary.
        const bool zero_compression = (static_cast<int64_t>(M_cell_est) * 10 >= static_cast<int64_t>(n) * 9);
        // Save for auto-fallback: only st.weights is mutated by solvers in-place
        const std::vector<double> weights_backup(weights);
        if (zero_compression && K >= 5) {
            // Compute target_skew on the AUTO dispatch path.
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
                st.oris_auto_selected = true;
                st.accelerate = true;
                auto res = lbw::oris_solve(st);
                res_status     = res.base.status;
                res_iterations = res.base.iterations;
                res_max_error  = res.base.max_error;
                res_alg_used   = (int)RK_ALG_ORIS;
                res_n_xcur_writes         = res.n_xcur_writes_per_iter_last;
                res_min_alpha             = res.min_alpha_seen;
                res_final_alpha           = res.final_alpha;
                res_n_bounds_violated     = res.n_bounds_violated;
                res_n_bounds_clamped      = res.n_bounds_clamped;
                res_homotopy_levels_used  = res.homotopy_levels_used;
                res_homotopy_final_factor = res.homotopy_final_factor;
                res_greedy_sweeps_taken   = res.greedy_sweeps_taken;
                res_eta_final             = res.eta_final;
                pack_solver_result(res);
                res_sor_min_omega    = res.sor_min_omega;
                res_sor_n_damped     = res.sor_n_damped;
                res_aa_accepted_count     = res.aa_accepted_count;
                res_sraa_demoted          = res.sraa_demoted ? 1 : 0;
                res_best_weights = std::move(res.base.best_weights);
            } else {
                auto res = lbw::newton_calibrate(st);
                pack_solver_result(res);
                res_status     = res.base.status;
                res_iterations = res.base.iterations;
                res_max_error  = res.base.max_error;
                res_alg_used   = (int)RK_ALG_NEWTON_KL;
                res_n_projected_dims = res.n_projected_dims;
                res_lm_mu_final      = res.lm_mu_final;
                res_n_bounds_violated = res.n_bounds_violated;
                if (!res.base.best_weights.empty())
                    res_best_weights = std::move(res.base.best_weights);
                else
                    res_best_weights.assign(st.n, 0.0);
            }
        } else if (zero_compression) {
            // K<5, zero-compression: raking
            auto res = lbw::raking_solve(st);
            res_status     = res.base.status;
            res_iterations = res.base.iterations;
            res_max_error  = res.base.max_error;
            res_alg_used   = (int)RK_ALG_RAKING;
            pack_solver_result(res);
            res_sraa_demoted = res.sraa_demoted ? 1 : 0;
            res_best_weights = std::move(res.base.best_weights);
        } else {
            // Compressed regime: ORIS (any K)
            st.oris_auto_selected = true;
            auto res = lbw::oris_solve(st);
            res_status     = res.base.status;
            res_iterations = res.base.iterations;
            res_max_error  = res.base.max_error;
            res_alg_used   = (int)RK_ALG_ORIS;
            res_n_xcur_writes         = res.n_xcur_writes_per_iter_last;
            res_min_alpha             = res.min_alpha_seen;
            res_final_alpha           = res.final_alpha;
            res_n_bounds_violated     = res.n_bounds_violated;
            res_n_bounds_clamped      = res.n_bounds_clamped;
            res_homotopy_levels_used  = res.homotopy_levels_used;
            res_homotopy_final_factor = res.homotopy_final_factor;
            res_greedy_sweeps_taken   = res.greedy_sweeps_taken;
            res_eta_final             = res.eta_final;
            pack_solver_result(res);
            res_sor_min_omega    = res.sor_min_omega;
            res_sor_n_damped     = res.sor_n_damped;
            res_aa_accepted_count     = res.aa_accepted_count;
            res_sraa_demoted          = res.sraa_demoted ? 1 : 0;
            res_best_weights = std::move(res.base.best_weights);
        }
        // Auto-fallback: if primary solver NOCONVs or exhausts budget (still
        // improving), retry with newton_kl.  STALL(5) is excluded: the solver
        // is at the constrained optimum and fallback cannot improve it.
        if (res_status == RK_ERR_NOCONV || res_status == RK_ERR_BUDGET) {
            if (st.verbose >= 1)
                st.log("auto: primary solver NOCONV/BUDGET; retrying with newton_kl");
            // Restore original weights (only mutated field in CalibState)
            std::copy(weights_backup.begin(), weights_backup.end(), weights.begin());
            st.oris_auto_selected = false;
            auto fb = lbw::newton_calibrate(st);
            res_status     = fb.base.status;
            res_iterations = fb.base.iterations;
            res_max_error  = fb.base.max_error;
            res_alg_used   = (int)RK_ALG_NEWTON_KL;
            res_n_projected_dims = fb.n_projected_dims;
            res_lm_mu_final      = fb.lm_mu_final;
            pack_solver_result(fb);
            if (!fb.base.best_weights.empty())
                res_best_weights = std::move(fb.base.best_weights);
            else
                res_best_weights.assign(st.n, 0.0);
        }
    } else if (strcmp(method_str, "sinkhorn") == 0) {
        auto res = lbw::sinkhorn_solve(st);
        pack_solver_result(res);
        res_status     = res.base.status;
        res_iterations = res.base.iterations;
        res_max_error  = res.base.max_error;
        res_alg_used   = static_cast<int>(RK_ALG_SINKHORN);
        if (!res.base.best_weights.empty())
            res_best_weights = std::move(res.base.best_weights);
        else
            res_best_weights.assign(st.n, 0.0);
    } else if (strcmp(method_str, "greg") == 0) {
        auto res = lbw::greg_solve(st);
        pack_solver_result(res);
        res_status     = res.base.status;
        res_iterations = res.base.iterations;
        res_max_error  = res.base.max_error;
        res_alg_used   = static_cast<int>(RK_ALG_GREG);
        if (!res.base.best_weights.empty())
            res_best_weights = std::move(res.base.best_weights);
        else
            res_best_weights.assign(st.n, 0.0);
    } else if (strcmp(method_str, "greenkhorn") == 0) {
        // validate: min_weight < max_weight
        if (st.min_weight >= st.max_weight && st.max_weight > 0.0)
            throw std::runtime_error("greenkhorn requires min_weight < max_weight");
        auto res = lbw::greenkhorn_solve(st);
        pack_solver_result(res);
        res_status     = res.base.status;
        res_iterations = res.base.iterations;
        res_max_error  = res.base.max_error;
        res_alg_used   = (int)RK_ALG_GREENKHORN;
        if (!res.base.best_weights.empty())
            res_best_weights = std::move(res.base.best_weights);
        else
            res_best_weights.assign(st.n, 0.0);
    } else if (strcmp(method_str, "logit") == 0) {
        auto res = lbw::logit_calibrate(st);
        pack_solver_result(res);
        res_status     = res.base.status;
        res_iterations = res.base.iterations;
        res_max_error  = res.base.max_error;
        res_alg_used   = (int)RK_ALG_LOGIT;
        if (!res.base.best_weights.empty())
            res_best_weights = std::move(res.base.best_weights);
        else
            res_best_weights.assign(st.n, 0.0);
    } else if (strcmp(method_str, "newton_kl") == 0) {
        auto res = lbw::newton_calibrate(st);
        pack_solver_result(res);
        res_status     = res.base.status;
        res_iterations = res.base.iterations;
        res_max_error  = res.base.max_error;
        res_alg_used   = (int)RK_ALG_NEWTON_KL;
        res_n_projected_dims  = res.n_projected_dims;
        res_lm_mu_final       = res.lm_mu_final;
        res_n_bounds_violated = res.n_bounds_violated;  // surface violation count (leafblower-73d7)
        if (!res.base.best_weights.empty())
            res_best_weights = std::move(res.base.best_weights);
        else
            res_best_weights.assign(st.n, 0.0);  // sentinel zeros: violation guard left best_weights empty
    } else {
        // Dispatch for chebyshev, oris_soft, and default oris.

        // Run oris warm-start BEFORE chebyshev dispatch.
        std::vector<double> w_warm_obs;
        double delta_warm = -1.0;
        if (strcmp(method_str, "chebyshev") == 0) {
            // SAFETY: weights_copy protects st.weights from oris_solve mutation.
            std::vector<double> weights_copy(st.weights, st.weights + st.n);
            lbw::CalibState st_warm = st;
            st_warm.weights = weights_copy.data();
            st_warm.inner_max_iter = std::max(5, std::min(100, st.inner_max_iter / 10));
            auto oris_res = lbw::oris_solve(st_warm);
            // st_warm must NOT escape this block (dangling pointer after weights_copy destroyed).
            if (!oris_res.base.best_weights.empty() &&
                static_cast<int>(oris_res.base.best_weights.size()) == st.n &&
                std::isfinite(oris_res.base.max_error)) {
                w_warm_obs = std::move(oris_res.base.best_weights);
                delta_warm = oris_res.base.max_error * 1.5;
            }
        }

        auto dispatch_cheb = [&](int alg_code) {
            const std::vector<double>& warm_ref = w_warm_obs;
            const double d_warm = delta_warm;
            auto res = lbw::chebyshev_ipm(st, warm_ref, d_warm);
            pack_solver_result(res);
            res_status     = res.base.status;
            res_iterations = res.base.iterations;
            res_max_error  = res.base.max_error;
            res_alg_used   = alg_code;
            if (!res.base.best_weights.empty())
                res_best_weights = std::move(res.base.best_weights);
            else
                res_best_weights.assign(st.n, 0.0);
        };
        if (strcmp(method_str, "chebyshev") == 0) {
            dispatch_cheb(static_cast<int>(RK_ALG_CHEBYSHEV));
        } else if (strcmp(method_str, "oris_soft") == 0) {
            st.oris_auto_selected = false;
            st.use_admm_capacity   = true;
            auto res = lbw::oris_solve(st);
            res_status     = res.base.status;
            res_iterations = res.base.iterations;
            res_max_error  = res.base.max_error;
            res_alg_used   = (int)RK_ALG_ORIS_SOFT;
            res_n_xcur_writes         = res.n_xcur_writes_per_iter_last;
            res_min_alpha             = res.min_alpha_seen;
            res_final_alpha           = res.final_alpha;
            res_n_bounds_violated     = res.n_bounds_violated;
            res_n_bounds_clamped      = res.n_bounds_clamped;
            res_homotopy_levels_used  = res.homotopy_levels_used;
            res_homotopy_final_factor = res.homotopy_final_factor;
            res_greedy_sweeps_taken   = res.greedy_sweeps_taken;
            res_eta_final             = res.eta_final;
            pack_solver_result(res);
            res_sor_min_omega    = res.sor_min_omega;
            res_sor_n_damped     = res.sor_n_damped;
            res_alm_capacity_mu_final = res.alm_capacity_mu_final;
            res_alm_n_growth_events   = res.alm_n_growth_events;
            res_alm_max_dual_norm     = res.alm_max_dual_norm;
            res_alm_sum_drift         = res.alm_sum_drift;
            res_aa_accepted_count     = res.aa_accepted_count;
            res_sraa_demoted          = res.sraa_demoted ? 1 : 0;
            res_best_weights = std::move(res.base.best_weights);
        } else {
            // Default / oris
            st.oris_auto_selected = (strcmp(method_str, "oris") != 0);
            auto res = lbw::oris_solve(st);
            res_status     = res.base.status;
            res_iterations = res.base.iterations;
            res_max_error  = res.base.max_error;
            res_alg_used   = (int)RK_ALG_ORIS;
            res_n_xcur_writes         = res.n_xcur_writes_per_iter_last;
            res_min_alpha             = res.min_alpha_seen;
            res_final_alpha           = res.final_alpha;
            res_n_bounds_violated     = res.n_bounds_violated;
            res_n_bounds_clamped      = res.n_bounds_clamped;
            res_homotopy_levels_used  = res.homotopy_levels_used;
            res_homotopy_final_factor = res.homotopy_final_factor;
            res_greedy_sweeps_taken   = res.greedy_sweeps_taken;
            res_eta_final             = res.eta_final;
            pack_solver_result(res);
            res_sor_min_omega    = res.sor_min_omega;
            res_sor_n_damped     = res.sor_n_damped;
            res_aa_accepted_count     = res.aa_accepted_count;
            res_sraa_demoted          = res.sraa_demoted ? 1 : 0;
            res_best_weights = std::move(res.base.best_weights);
        }
    }
    } catch (const std::exception& e) {
        solver_error = e.what();
    } catch (...) {
        solver_error = "unknown exception";
    }
    }
    if (!solver_error.empty()) {
        // pre_error messages are already fully-qualified; solver errors get the
        // "internal solver error" prefix so callers can distinguish the two.
        if (!pre_error.empty())
            Rf_error("%s", solver_error.c_str());
        else
            Rf_error("leafblower: internal solver error \xe2\x80\x94 %s",
                     solver_error.c_str());
    }

    // Single source of truth for rk_algorithm_t → R-visible name.
    // Indices match enum values in leafblower.h. Update both together.
    static const char* kAlgNames[] = {
        "",           // 0 = RK_ALG_AUTO
        "oris",       // 1 = RK_ALG_ORIS
        "",           // 2 = (removed lbfgsb slot)
        "raking",     // 3 = RK_ALG_RAKING
        "sinkhorn",   // 4 = RK_ALG_SINKHORN
        "chebyshev",  // 5 = RK_ALG_CHEBYSHEV
        "greg",       // 6 = RK_ALG_GREG
        "",           // 7 = deprecated GRAKE
        "oris_soft",  // 8 = RK_ALG_ORIS_SOFT
        "greenkhorn", // 9 = RK_ALG_GREENKHORN
        "logit",      // 10 = RK_ALG_LOGIT
        "newton_kl",  // 11 = RK_ALG_NEWTON_KL
    };
    static const int kAlgNamesLen = 12;
    static_assert(RK_ALG_NEWTON_KL == 11, "kAlgNames table needs update on enum change");
    const char* alg_name_cstr = (res_alg_used >= 0 && res_alg_used < kAlgNamesLen)
        ? kAlgNames[res_alg_used]
        : "unknown";
    std::snprintf(res_message, 256, "%s: %d iters, max_error=%.2e",
                  alg_name_cstr, res_iterations, res_max_error);

    // greenkhorn and logit do not modify st.weights in-place; copy calibrated
    // weights into the weights vector so raw$weights in harvest.R is correct.
    // (raking/oris already write to st.weights in-place — don't copy.)
    if (!res_best_weights.empty() && (int)res_best_weights.size() == n &&
        (res_alg_used == static_cast<int>(RK_ALG_GREENKHORN) ||
         res_alg_used == static_cast<int>(RK_ALG_LOGIT))) {
        std::copy(res_best_weights.begin(), res_best_weights.end(), weights.begin());
    }

    // Build return list: list(weights=numeric[n], result=list(42 fields))
    // PROTECT out first so wts/res_list/res_names sit above it on the stack;
    // each is UNPROTECTed immediately after adoption, in LIFO order.
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));  // owns wts + res_list below

    SEXP wts = PROTECT(Rf_allocVector(REALSXP, n));
    memcpy(REAL(wts), weights.data(), (size_t)n * sizeof(double));

    constexpr int N_RESULT_FIELDS = 42;
    SEXP res_list  = PROTECT(Rf_allocVector(VECSXP,  N_RESULT_FIELDS));  // 14 prior + 8 scalars + best_weights + 7 convergence fields + 4 ALM diagnostics + 1 SRAA diagnostic + 1 Newton-KL TSVD diagnostic + 1 Newton-KL LM diagnostic + 1 metric_first_check + 1 metric_prev_check + 1 prev_check_iter + 1 sraa_demoted + 1 convergence_stall_kind
    SEXP res_names = PROTECT(Rf_allocVector(STRSXP,  N_RESULT_FIELDS));
    SET_STRING_ELT(res_names, 0, Rf_mkChar("status"));
    SET_STRING_ELT(res_names, 1, Rf_mkChar("iterations"));
    SET_STRING_ELT(res_names, 2, Rf_mkChar("max_error"));
    SET_STRING_ELT(res_names, 3, Rf_mkChar("algorithm_used"));
    SET_STRING_ELT(res_names, 4, Rf_mkChar("message"));
    SET_STRING_ELT(res_names, 5, Rf_mkChar("n_xcur_writes_per_iter_last"));
    SET_STRING_ELT(res_names, 6, Rf_mkChar("min_alpha_seen"));
    SET_STRING_ELT(res_names, 7, Rf_mkChar("final_alpha"));
    SET_STRING_ELT(res_names, 8, Rf_mkChar("n_bounds_violated"));
    SET_STRING_ELT(res_names, 9, Rf_mkChar("n_bounds_clamped"));
    SET_VECTOR_ELT(res_list, 0, Rf_ScalarInteger(res_status));
    SET_VECTOR_ELT(res_list, 1, Rf_ScalarInteger(res_iterations));
    SET_VECTOR_ELT(res_list, 2, Rf_ScalarReal(res_max_error));
    SET_VECTOR_ELT(res_list, 3, Rf_mkString(alg_name_cstr));
    SET_VECTOR_ELT(res_list, 4, Rf_mkString(res_message));
    SET_VECTOR_ELT(res_list, 5, Rf_ScalarInteger(res_n_xcur_writes));
    SET_VECTOR_ELT(res_list, 6, Rf_ScalarReal(res_min_alpha));
    SET_VECTOR_ELT(res_list, 7, Rf_ScalarReal(res_final_alpha));
    SET_VECTOR_ELT(res_list, 8, Rf_ScalarInteger(res_n_bounds_violated));
    SET_VECTOR_ELT(res_list, 9, Rf_ScalarInteger(res_n_bounds_clamped));
    SET_STRING_ELT(res_names, 10, Rf_mkChar("homotopy_levels_used"));
    SET_STRING_ELT(res_names, 11, Rf_mkChar("homotopy_final_factor"));
    SET_STRING_ELT(res_names, 12, Rf_mkChar("greedy_sweeps_taken"));
    SET_STRING_ELT(res_names, 13, Rf_mkChar("eta_final"));
    SET_VECTOR_ELT(res_list, 10, Rf_ScalarInteger(res_homotopy_levels_used));
    SET_VECTOR_ELT(res_list, 11, Rf_ScalarReal(res_homotopy_final_factor));
    SET_VECTOR_ELT(res_list, 12, Rf_ScalarInteger(res_greedy_sweeps_taken));
    SET_VECTOR_ELT(res_list, 13, Rf_ScalarReal(res_eta_final));
    /* Quality metric scalars (indices 14-21) */
    SET_STRING_ELT(res_names, 14, Rf_mkChar("mean_error"));
    SET_STRING_ELT(res_names, 15, Rf_mkChar("kl"));
    SET_STRING_ELT(res_names, 16, Rf_mkChar("chi2"));
    SET_STRING_ELT(res_names, 17, Rf_mkChar("l1_weight_change"));
    SET_STRING_ELT(res_names, 18, Rf_mkChar("best_error"));
    SET_STRING_ELT(res_names, 19, Rf_mkChar("best_iter"));
    SET_STRING_ELT(res_names, 20, Rf_mkChar("sor_min_omega"));
    SET_STRING_ELT(res_names, 21, Rf_mkChar("sor_n_damped"));
    SET_VECTOR_ELT(res_list, 14, Rf_ScalarReal(res_mean_error));
    SET_VECTOR_ELT(res_list, 15, Rf_ScalarReal(res_kl));
    SET_VECTOR_ELT(res_list, 16, Rf_ScalarReal(res_chi2));
    SET_VECTOR_ELT(res_list, 17, Rf_ScalarReal(res_l1_weight_change));
    SET_VECTOR_ELT(res_list, 18, Rf_ScalarReal(res_best_error));
    SET_VECTOR_ELT(res_list, 19, Rf_ScalarInteger(res_best_iter));
    SET_VECTOR_ELT(res_list, 20, Rf_ScalarReal(res_sor_min_omega));
    SET_VECTOR_ELT(res_list, 21, Rf_ScalarInteger(res_sor_n_damped));
    /* Element 22: best_weights REALSXP (WU-E) */
    SET_STRING_ELT(res_names, 22, Rf_mkChar("best_weights"));
    {
        int bw_n = (int)res_best_weights.size();
        SEXP bw_sxp = PROTECT(Rf_allocVector(REALSXP, bw_n));
        double* bw  = REAL(bw_sxp);
        for (int i = 0; i < bw_n; i++) bw[i] = res_best_weights[i];
        SET_VECTOR_ELT(res_list, 22, bw_sxp);
        UNPROTECT(1);  // bw_sxp adopted by res_list
    }
    /* Elements 23-27: convergence diagnostics (WU-A) */
    SET_STRING_ELT(res_names, 23, Rf_mkChar("grake_norm"));
    SET_STRING_ELT(res_names, 24, Rf_mkChar("convergence_metric"));
    SET_STRING_ELT(res_names, 25, Rf_mkChar("convergence_rule"));
    SET_STRING_ELT(res_names, 26, Rf_mkChar("convergence_tol"));
    SET_STRING_ELT(res_names, 27, Rf_mkChar("convergence_iter"));
    SET_VECTOR_ELT(res_list, 23, Rf_ScalarReal(res_grake_norm));
    SET_VECTOR_ELT(res_list, 24, Rf_ScalarInteger(res_conv_metric));
    SET_VECTOR_ELT(res_list, 25, Rf_ScalarInteger(res_conv_rule));
    SET_VECTOR_ELT(res_list, 26, Rf_ScalarReal(res_conv_tol));
    SET_VECTOR_ELT(res_list, 27, Rf_ScalarInteger(res_conv_iter));
    /* Elements 28-29: convergence_objective and convergence_minimized_metric (Task 1) */
    SET_STRING_ELT(res_names, 28, Rf_mkChar("solver_objective"));
    SET_STRING_ELT(res_names, 29, Rf_mkChar("convergence_minimized_metric"));
    SET_VECTOR_ELT(res_list,  28, Rf_ScalarReal(res_conv_objective));
    SET_VECTOR_ELT(res_list,  29, Rf_ScalarInteger(res_conv_minimized_metric));
    /* Elements 30-33: ALM diagnostics (non-zero only for oris_soft) */
    SET_STRING_ELT(res_names, 30, Rf_mkChar("alm_capacity_mu_final"));
    SET_STRING_ELT(res_names, 31, Rf_mkChar("alm_n_growth_events"));
    SET_STRING_ELT(res_names, 32, Rf_mkChar("alm_max_dual_norm"));
    SET_STRING_ELT(res_names, 33, Rf_mkChar("alm_sum_drift"));
    SET_VECTOR_ELT(res_list,  30, Rf_ScalarReal(res_alm_capacity_mu_final));
    SET_VECTOR_ELT(res_list,  31, Rf_ScalarInteger(res_alm_n_growth_events));
    SET_VECTOR_ELT(res_list,  32, Rf_ScalarReal(res_alm_max_dual_norm));
    SET_VECTOR_ELT(res_list,  33, Rf_ScalarReal(res_alm_sum_drift));
    /* Element 34: SRAA acceleration diagnostic (oris/oris_soft only; zero elsewhere) */
    SET_STRING_ELT(res_names, 34, Rf_mkChar("aa_accepted_count"));
    SET_VECTOR_ELT(res_list,  34, Rf_ScalarInteger(res_aa_accepted_count));
    /* Element 35: Newton-KL TSVD diagnostic (newton_kl only; zero elsewhere) */
    SET_STRING_ELT(res_names, 35, Rf_mkChar("n_projected_dims"));
    SET_VECTOR_ELT(res_list,  35, Rf_ScalarInteger(res_n_projected_dims));
    /* Element 36: Newton-KL Levenberg-Marquardt diagnostic (newton_kl only; zero elsewhere) */
    SET_STRING_ELT(res_names, 36, Rf_mkChar("lm_mu_final"));
    SET_VECTOR_ELT(res_list,  36, Rf_ScalarReal(res_lm_mu_final));
    /* Element 37: first-check metric value (oris only; Inf elsewhere) */
    SET_STRING_ELT(res_names, 37, Rf_mkChar("metric_first_check"));
    SET_VECTOR_ELT(res_list,  37, Rf_ScalarReal(res_metric_first_check));
    SET_STRING_ELT(res_names, 38, Rf_mkChar("metric_prev_check"));
    SET_VECTOR_ELT(res_list,  38, Rf_ScalarReal(res_metric_prev_check));
    SET_STRING_ELT(res_names, 39, Rf_mkChar("prev_check_iter"));
    SET_VECTOR_ELT(res_list,  39, Rf_ScalarInteger(res_prev_check_iter));
    /* Element 40: SRAA scheduler-demotion flag (oris/raking only; FALSE elsewhere) */
    SET_STRING_ELT(res_names, 40, Rf_mkChar("sraa_demoted"));
    SET_VECTOR_ELT(res_list,  40, Rf_ScalarLogical(res_sraa_demoted));
    /* Element 41: solver-emitted stall kind (leafblower-8eod).
       0=no stall, 1=wchange (SRAA path), 2=kl (plain-IPF path).
       Set at RK_ERR_STALL emission site; replaces accelerate_bool heuristic in harvest.R. */
    SET_STRING_ELT(res_names, 41, Rf_mkChar("convergence_stall_kind"));
    SET_VECTOR_ELT(res_list,  41, Rf_ScalarInteger(res_stall_kind));
    Rf_setAttrib(res_list, R_NamesSymbol, res_names);
    UNPROTECT(1);  // res_names adopted by res_list

    SET_VECTOR_ELT(out, 1, res_list);
    UNPROTECT(1);  // res_list adopted by out
    SET_VECTOR_ELT(out, 0, wts);
    UNPROTECT(1);  // wts adopted by out
    {
        SEXP out_names = PROTECT(Rf_allocVector(STRSXP,  2));
        SET_STRING_ELT(out_names, 0, Rf_mkChar("weights"));
        SET_STRING_ELT(out_names, 1, Rf_mkChar("result"));
        Rf_setAttrib(out, R_NamesSymbol, out_names);
        UNPROTECT(1);  // out_names adopted by out
    }
    UNPROTECT(1);  // out: return transfers ownership to caller
    return out;
}

// test-only: exposes CellTable internals for unit tests
extern "C" SEXP C_leafblower_cell_table_probe(SEXP r_group_ids_list, SEXP r_n) {
    int n = scalar_int(r_n, "r_n");
    int K = Rf_length(r_group_ids_list);
    if (K > lbw::K_MAX) {
        Rf_error("K (%d) exceeds K_MAX (%d)", K, lbw::K_MAX);
    }
    // Extract pointers
    std::vector<const int32_t*> gid_ptrs(K);
    std::vector<int> cat_counts(K);
    for (int k = 0; k < K; k++) {
        SEXP v = VECTOR_ELT(r_group_ids_list, k);
        gid_ptrs[k] = (const int32_t*) INTEGER(v);
        // Derive cat_counts from max value + 1 (excluding -1)
        int max_g = -1;
        for (int i = 0; i < n; i++) {
            int g = gid_ptrs[k][i];
            if (g > max_g) max_g = g;
        }
        cat_counts[k] = (max_g < 0) ? 1 : (max_g + 1);
    }
    std::vector<double> uniform_weights(n, 1.0);
    lbw::CellTable ct;
    int rc = lbw::build_cell_table(n, K, gid_ptrs.data(), cat_counts.data(),
                                    uniform_weights.data(), ct);
    if (rc != 0) Rf_error("build_cell_table failed (rc=%d)", rc);

    // Build return list: list(M_cell, cell_of, n_per_cell)
    SEXP ret = PROTECT(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ret, 0, Rf_ScalarInteger(ct.M_cell));
    SEXP cell_of_sexp = PROTECT(Rf_allocVector(INTSXP, n));
    std::memcpy(INTEGER(cell_of_sexp), ct.cell_of.data(), (size_t)n * sizeof(int));
    SET_VECTOR_ELT(ret, 1, cell_of_sexp);
    SEXP npc = PROTECT(Rf_allocVector(INTSXP, ct.M_cell));
    std::memcpy(INTEGER(npc), ct.n_per_cell.data(), (size_t)ct.M_cell * sizeof(int));
    SET_VECTOR_ELT(ret, 2, npc);
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, Rf_mkChar("M_cell"));
    SET_STRING_ELT(names, 1, Rf_mkChar("cell_of"));
    SET_STRING_ELT(names, 2, Rf_mkChar("n_per_cell"));
    Rf_setAttrib(ret, R_NamesSymbol, names);
    UNPROTECT(4);  // ret + cell_of_sexp + npc + names
    return ret;
}

SEXP C_rk_design_effect(SEXP weights_sexp, SEXP outcome_sexp,
                         SEXP data_codes_sexp, SEXP cat_counts_sexp,
                         SEXP K_sexp) {
    if (TYPEOF(weights_sexp) != REALSXP)
        Rf_error("design_effect: 'weights' must be REALSXP (got TYPEOF=%d)", TYPEOF(weights_sexp));
    if (TYPEOF(K_sexp) != INTSXP || Rf_length(K_sexp) != 1)
        Rf_error("design_effect: 'K' must be INTSXP scalar");

    const int n = Rf_length(weights_sexp);
    const int K = INTEGER(K_sexp)[0];
    const double* weights = REAL(weights_sexp);

    const double* outcome = nullptr;
    if (TYPEOF(outcome_sexp) == REALSXP) {
        if (Rf_length(outcome_sexp) != n)
            Rf_error("design_effect: 'outcome' length (%d) must equal length(weights) (%d)",
                     Rf_length(outcome_sexp), n);
        outcome = REAL(outcome_sexp);
    } else if (TYPEOF(outcome_sexp) != NILSXP) {
        Rf_error("design_effect: 'outcome' must be REALSXP or NULL (got TYPEOF=%d)",
                 TYPEOF(outcome_sexp));
    }

    const int* data_codes = nullptr;
    const int* cat_counts = nullptr;
    if (K > 0 && outcome != nullptr) {
        if (TYPEOF(data_codes_sexp) != INTSXP)
            Rf_error("design_effect: 'data_codes' must be INTSXP (got TYPEOF=%d)",
                     TYPEOF(data_codes_sexp));
        if (Rf_length(data_codes_sexp) != n * K)
            Rf_error("design_effect: length(data_codes) must be n*K = %d (got %d)",
                     n * K, Rf_length(data_codes_sexp));
        if (TYPEOF(cat_counts_sexp) != INTSXP || Rf_length(cat_counts_sexp) != K)
            Rf_error("design_effect: 'cat_counts' must be INTSXP of length K=%d", K);
        data_codes = INTEGER(data_codes_sexp);
        cat_counts = INTEGER(cat_counts_sexp);
    }

    rk_design_effect_result_t out;
    const int status = rk_design_effect(weights, outcome, data_codes, cat_counts, n, K, &out);
    if (status != RK_OK)
        Rf_error("design_effect: %s", out.message);

    SEXP res       = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP res_names = PROTECT(Rf_allocVector(STRSXP, 4));
    SET_VECTOR_ELT(res, 0, Rf_ScalarReal(out.deff_K));
    SET_STRING_ELT(res_names, 0, Rf_mkChar("deff_K"));
    SET_VECTOR_ELT(res, 1, Rf_ScalarReal(out.deff_H));
    SET_STRING_ELT(res_names, 1, Rf_mkChar("deff_H"));
    SET_VECTOR_ELT(res, 2, Rf_ScalarInteger(out.rank_def));
    SET_STRING_ELT(res_names, 2, Rf_mkChar("rank_def"));
    SET_VECTOR_ELT(res, 3, Rf_mkString(out.message));
    SET_STRING_ELT(res_names, 3, Rf_mkChar("message"));
    Rf_setAttrib(res, R_NamesSymbol, res_names);
    UNPROTECT(2);  // res_names adopted by res; res returned (protected by caller)
    return res;
}

} // extern "C"
