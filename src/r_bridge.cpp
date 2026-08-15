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
#include "calib_dispatch.hpp" // SC1: shared dispatch table (leafblower-rywn)
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

// CR-D9 (j7x8.9): deferred-throw variants for the C_rk_calibrate hot path. The
// direct-Rf_error versions above longjmp past live C++ destructors (group_ids,
// tgt_storage, weights, res_best_weights, pre_error/solver_error) -- a leak on
// every bad-scalar arg. These set pre_error (first-error-wins) and return a
// harmless 0 instead; the single deferred error exit at the end of
// C_rk_calibrate releases all heap-backed locals before its Rf_error, so this
// path is leak-free. The throwing 2-arg forms are kept only for the tiny
// test-bridge helpers with no live RAII.
static inline double scalar_real(SEXP x, const char* name, std::string& pre_error) {
    if (TYPEOF(x) != REALSXP || LENGTH(x) != 1) {
        if (pre_error.empty()) {
            char buf[128];
            std::snprintf(buf, sizeof(buf),
                "leafblower: '%s' must be a length-1 numeric (got type %s length %d)",
                name, Rf_type2char(TYPEOF(x)), (int)LENGTH(x));
            pre_error = buf;
        }
        return 0.0;
    }
    return REAL(x)[0];
}
static inline int scalar_int(SEXP x, const char* name, std::string& pre_error) {
    if (TYPEOF(x) != INTSXP || LENGTH(x) != 1) {
        if (pre_error.empty()) {
            char buf[128];
            std::snprintf(buf, sizeof(buf),
                "leafblower: '%s' must be a length-1 integer (got type %s length %d)",
                name, Rf_type2char(TYPEOF(x)), (int)LENGTH(x));
            pre_error = buf;
        }
        return 0;
    }
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
                    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
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
    st.sor_cfg.omega_max            = p.sor_omega_max;
    st.sor_cfg.omega_fixed          = p.sor_omega_fixed;
    st.sor_cfg.burnin               = p.sor_burnin;
    st.sor_cfg.omega_mode_id        = p.sor_omega_mode_id;
}

extern "C" {

SEXP C_rk_design_effect(SEXP, SEXP, SEXP, SEXP, SEXP);  // defined below

void R_init_leafblower(DllInfo* dll) {
    static const R_CallMethodDef call_methods[] = {
        {"C_logit_F_at_zero",    (DL_FUNC)&C_logit_F_at_zero,    2},
        {"C_logit_range_check",  (DL_FUNC)&C_logit_range_check,  3},
        {"C_logit_Hprime_check", (DL_FUNC)&C_logit_Hprime_check, 3},
        {"C_rk_calibrate",       (DL_FUNC)&C_rk_calibrate,       39},
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
                    SEXP sor_omega_max_sexp,
                    SEXP sor_omega_fixed_sexp, SEXP sor_burnin_sexp,
                    SEXP sor_omega_mode_id_sexp,
                    /* SQUAREM */
                    SEXP accelerate_sexp,
                    /* Epic-H WH-e: newton_kl TSVD truncation ratio (default 1e-8 from R layer). */
                    SEXP newton_tsvd_ratio_sexp,
                    /* Tikhonov ridge on dual λ; 0.0 = off. */
                    SEXP ridge_lambda_sexp) {
    std::string pre_error;
    if (TYPEOF(group_ids_sexp) != VECSXP)
        pre_error = "group_ids must be a list";
    int K = pre_error.empty() ? LENGTH(group_ids_sexp) : 0;
    int n = scalar_int(n_obs_sexp, "n_obs", pre_error);
    // CR-D7 (j7x8.7): range-check n BEFORE any n-sized allocation. A direct .Call
    // with n_obs=-1 or NA_integer_ (=INT_MIN) reaches the std::vector<double>
    // weights(n) below — OUTSIDE the solver try block — so a bad_alloc/length_error
    // there would be uncaught → std::terminate → R session abort. Deferring via
    // pre_error (not Rf_error) keeps the CR-D9 RAII-clean-unwind contract (every
    // pre-dispatch error unwinds the live std::string/vectors before one throw); the
    // guarded alloc below then stays empty so nothing throws before the deferred throw.
    if (pre_error.empty() && n <= 0)
        pre_error = "n_obs must be a positive integer";
    // CR-D9a (j7x8.17): with K=0 (group_ids=list()) there are no per-margin
    // LENGTH(gid)!=n cross-checks to constrain n, so a huge positive n_obs would
    // still reach the pre-try weights(n) alloc -> bad_alloc -> std::terminate. A
    // zero-margin calibration is degenerate anyway; reject it via the same deferred
    // path so the guarded alloc below stays empty.
    if (pre_error.empty() && K <= 0)
        pre_error = "at least one margin is required (group_ids is empty)";

    // CR-H6 (xc1s.6): alias R's INTSXP data directly instead of a K×n deep copy.
    // Audit (all 8 solvers + cell_table/validation/newton) confirms group_ids is
    // read-only (int g = group_ids[k][i]; zero writes, zero const_cast to mutable
    // int32). The R INTSXP args stay PROTECT'd for the whole .Call and solvers run
    // synchronously, so the aliased pointers outlive their use.
    // is_same (not just sizeof) makes the reinterpret_cast a provable no-op —
    // closes any strict-aliasing question. Holds on every R platform (glibc/musl/
    // macOS/MinGW all typedef int32_t as int).
    static_assert(std::is_same<int, std::int32_t>::value,
                  "R INTSXP (C int) must be int32_t to alias as const int32_t*");
    std::vector<const int32_t*> group_ids(K);
    std::vector<int> cat_counts(K);
    std::vector<std::vector<double>> tgt_storage(K);
    std::vector<const double*> targets(K);

    if (pre_error.empty()) {  // skip container validation if group_ids already failed (first-error-wins)
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
            const int* gp = INTEGER_RO(gid_vec);
            group_ids[k] = reinterpret_cast<const int32_t*>(gp);  // alias, no copy
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

    // Build weights (start_weights already normalized to mean=1 by R layer).
    // CR-D7: size 0 when a pre_error is pending (e.g. n<=0) so this pre-try alloc
    // never receives a negative/huge n — the deferred throw at :~623 handles it.
    std::vector<double> weights(pre_error.empty() ? n : 0);
    if (pre_error.empty()) {
    if (Rf_isNull(start_weights_sexp)) {
        for (int i = 0; i < n; i++) weights[i] = 1.0;
    } else {
        if (TYPEOF(start_weights_sexp) != REALSXP || LENGTH(start_weights_sexp) != n) {
            char buf[128];
            std::snprintf(buf, sizeof(buf),
                "leafblower: start_weights must be numeric of length n=%d (got type %s length %d)",
                n, Rf_type2char(TYPEOF(start_weights_sexp)), (int)LENGTH(start_weights_sexp));
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
    p.min_weight     = scalar_real(min_weight_sexp, "min_weight", pre_error);
    p.max_weight     = scalar_real(max_weight_sexp, "max_weight", pre_error);
    p.verbose        = scalar_int(verbose_sexp, "verbose", pre_error);
    p.inner_max_iter = scalar_int(inner_max_iter_sexp, "max_iter", pre_error);
    // newton_kl: 2nd-order outer loop defaults to 50 (see newton_calib.cpp:90).
    // Other solvers: outer = inner (outer_max_iter is ignored for non-newton_kl paths).
    if (LENGTH(method_sexp) == 1 && TYPEOF(method_sexp) == STRSXP &&
        strcmp(CHAR(STRING_ELT(method_sexp, 0)), "newton_kl") == 0) {
        p.outer_max_iter = 0;  // triggers C-side 50-default in newton_calib.cpp:90
    } else {
        p.outer_max_iter = scalar_int(inner_max_iter_sexp, "max_iter", pre_error);
    }
    p.tol_abs        = scalar_real(tol_abs_sexp, "tol_abs", pre_error);
    p.bounds_mode    = (rk_bounds_mode_t) scalar_int(bounds_mode_sexp, "bounds_mode", pre_error);
    p.log_fn         = (p.verbose > 0) ? r_log_trampoline : nullptr;
    /* Overlay knobs */
    p.homotopy.n_levels        = scalar_int(homotopy_levels_sexp, "homotopy_levels", pre_error);
    p.homotopy.start_factor    = scalar_real(homotopy_start_factor_sexp, "homotopy_start_factor", pre_error);
    p.homotopy.end_factor      = scalar_real(homotopy_end_factor_sexp, "homotopy_end_factor", pre_error);
    // j7x8.21: scalar_real validates only TYPEOF/LENGTH. A non-finite / non-positive
    // homotopy factor flows into the geometric schedule k_start*pow(k_end/k_start,frac)
    // → NaN current_max_weight (oris.cpp), and a NaN start_factor silently disables the
    // CR-D20 tightest-level guard (min(NaN,end)=NaN). Reject via the deferred-throw
    // pre_error idiom (leak-free per CR-D9).
    if (pre_error.empty() &&
        (!std::isfinite(p.homotopy.start_factor) || p.homotopy.start_factor <= 0.0 ||
         !std::isfinite(p.homotopy.end_factor)   || p.homotopy.end_factor   <= 0.0))
        pre_error = "leafblower: 'homotopy_start_factor' and 'homotopy_end_factor' "
                    "must be finite and > 0";
    p.homotopy.budget_split_p  = scalar_real(homotopy_budget_p_sexp, "homotopy_budget_p", pre_error);
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
    p.eta_start           = scalar_real(eta_start_sexp, "eta_start", pre_error);
    p.eta_end             = scalar_real(eta_end_sexp, "eta_end", pre_error);
    p.eta_schedule_power  = scalar_real(eta_schedule_power_sexp, "eta_schedule_power", pre_error);
    /* Convergence config (WU-A) */
    p.pct_tol             = scalar_real(pct_tol_sexp, "pct_tol", pre_error);
    p.absolute_tol        = scalar_real(absolute_tol_sexp, "absolute_tol", pre_error);
    p.metric              = scalar_int(metric_sexp, "metric", pre_error);
    p.rule                = scalar_int(rule_sexp, "rule", pre_error);
    p.stop_when           = scalar_int(stop_when_sexp, "stop_when", pre_error);
    /* SOR config (WU-A) */
    p.sor_enabled         = scalar_int(sor_enabled_sexp, "sor_enabled", pre_error);
    p.sor_auto            = scalar_int(sor_auto_sexp, "sor_auto", pre_error);
    p.sor_omega_init      = scalar_real(sor_omega_init_sexp, "sor_omega_init", pre_error);
    p.sor_omega_min       = scalar_real(sor_omega_min_sexp, "sor_omega_min", pre_error);
    p.sor_omega_max       = scalar_real(sor_omega_max_sexp, "sor_omega_max", pre_error);
    p.sor_omega_fixed     = scalar_real(sor_omega_fixed_sexp, "sor_omega_fixed", pre_error);
    p.sor_burnin          = scalar_int(sor_burnin_sexp, "sor_burnin", pre_error);
    p.sor_omega_mode_id   = scalar_int(sor_omega_mode_id_sexp, "sor_omega_mode_id", pre_error);

    if (pre_error.empty() && (TYPEOF(method_sexp) != STRSXP || LENGTH(method_sexp) != 1))
        pre_error = "method must be a length-1 character string";
    const char* method_str = pre_error.empty() ? CHAR(STRING_ELT(method_sexp, 0)) : "";
    const auto alg_it = pre_error.empty() ? kAlgMap.find(method_str) : kAlgMap.end();
    // CR-D10 (j7x8.10): reject unrecognized method instead of silently running
    // ORIS. Direct .Call bypasses R's match.arg(), so an unknown string (typo
    // "grg") previously validated as RAKING but executed ORIS — a silent
    // contract mismatch. Name the offending string in a graceful pre_error.
    if (pre_error.empty()) {
        if (alg_it == kAlgMap.end())
            pre_error = std::string("method: unrecognized algorithm \"") + method_str + "\"";
        else
            p.algorithm = alg_it->second;
    }

    // Full input validation — shared with c_api.cpp path via validation.hpp.
    if (pre_error.empty()) {
        rk_result_t validation_result;
        rk_result_init(&validation_result);
        // After CR-D10 (j7x8.10), pre_error.empty() here implies the method was
        // recognized (a miss sets pre_error above), so p.algorithm is resolved.
        rk_algorithm_t alg_for_validation = p.algorithm;
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
    st.newton_tsvd_ratio = scalar_real(newton_tsvd_ratio_sexp, "newton_tsvd_ratio", pre_error);
    // Tikhonov ridge on dual λ; 0.0 = off.
    st.ridge_lambda = scalar_real(ridge_lambda_sexp, "ridge_lambda", pre_error);
    st.oris_auto_selected          = false;  // R bridge always resolves method explicitly
    st.alm.lambda = 0.0;
    st.alm.mu     = 0.0;
    st.accelerate = (scalar_int(accelerate_sexp, "accelerate", pre_error) != 0);

    // Resolve capacity_mu for oris_soft: harmonized to estimate_M_cell path matching c_api.cpp.
    // leafblower-yh0l.4: T3 benchmark identified Py as winner on fulldata; differentiator was
    // R using exact M_cell/n via build_cell_table (~0.073 at K=9), Py using estimate_M_cell
    // which caps at n for K>8 (capacity_mu=1.0). ~14x scaling difference drove iter-0 divergence.
    // capacity_penalty <= 0.0 (or NULL) selects auto via estimate_M_cell (matches c_api.cpp:380);
    // positive value is used directly. ALM block is gated by st.use_admm_capacity,
    // so st.alm.capacity_mu is harmless for non-oris_soft callers.
    // eb79.15/SC1 (plan 07): lazy cache for lbw::resolve_m_cell_est (which
    // itself wraps estimate_M_cell, O(n·K)). This site (A) and route_auto's
    // AUTO-routing site below (B) share this ONE cache; on the method="auto"
    // + capacity_penalty≤0 path both blocks run, but the underlying estimate
    // computes once. Sentinel -1 = uncomputed; whichever site runs first
    // stores its result, the other reuses it. A's CXX.1 length-validity skip
    // is NOT hoisted into resolve_m_cell_est, so B still computes fresh when
    // A's branch was skipped.
    int m_cell_est_cache = -1;
    {
        const double cp_val = Rf_isNull(capacity_penalty_sexp)
            ? -1.0
            : (TYPEOF(capacity_penalty_sexp) == REALSXP && LENGTH(capacity_penalty_sexp) == 1 ? REAL(capacity_penalty_sexp)[0] : -1.0);
        if (cp_val > 0.0) {
            st.alm.capacity_mu = cp_val;
        } else if (pre_error.empty()) {
            // CXX.1: estimate_M_cell dereferences group_ids/cat_counts; skip when a
            // length-validation error is already pending (the throw at the deferred
            // check below converts pre_error into a graceful R error). Filling these
            // with malformed input would OOB-read.
            const int m_cell_est = lbw::resolve_m_cell_est(n, K,
                const_cast<const int32_t**>(group_ids.data()),
                cat_counts.data(), m_cell_est_cache);
            st.alm.capacity_mu = (n > 0) ? static_cast<double>(m_cell_est) / n : 1.0;
        }
    }

    // alm_penalty: objective ALM penalty coefficient (st.alm.mu). Positive value activates.
    {
        const double alm_penalty_val = Rf_isNull(alm_penalty_sexp)
            ? -1.0
            : (TYPEOF(alm_penalty_sexp) == REALSXP && LENGTH(alm_penalty_sexp) == 1 ? REAL(alm_penalty_sexp)[0] : -1.0);
        st.alm.mu = (alm_penalty_val > 0.0) ? alm_penalty_val : 0.0;
    }

    // Scalar fields mirrored from rk_result_t for compatibility with downstream assembly.
    int    res_status     = RK_ERR_NOCONV;
    int    res_iterations = 0;
    double res_max_error  = 1.0;
    int    res_alg_used   = (int)RK_ALG_ORIS;
    char   res_message[256] = "";
    char   res_solver_message[256] = "";  // eb79.25: captured solver res.message (error statuses)
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
    double res_sor_min_omega  = 1.0;
    int    res_sor_n_damped   = 0;
    double res_sor_omega_mean = 1.0;
    int    res_sor_any_latched  = 0;
    int    res_sor_n_pinned_fb  = 0;
    int    res_sor_n_warmup_fb  = 0;
    int    res_sor_n_conv_fb    = 0;
    int    res_sor_n_resid_grew = 0;
    int    res_sor_n_monotone_cd = 0;
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
    // SC1 (leafblower-rywn): scratch struct for the shared dispatch table
    // (calib_dispatch.hpp). Solvers not yet migrated (D-01) leave this
    // default-constructed and unused; migrated branches populate it via
    // lbw::dispatch_solver then copy into the res_* locals above.
    // Function-scope (not branch-scope) so its heap-backed best_weights
    // member is covered by both Rf_error swap-release blocks below, per the
    // RAII convention (RESEARCH.md Pitfall 4).
    lbw::DispatchResult dres;

    // SC1 (plan 07): pack_solver_result/pack_oris_result (the DRY lambdas
    // that used to feed the AUTO branch's severe-skew/moderate-skew/
    // compressed sub-branches) are retired — every AUTO outcome now flows
    // through lbw::dispatch_solver + a direct dres.* -> res_* copy, same
    // shape as every other already-migrated branch below. has_message/
    // has_n_bounds (top of file) become unused with them.

    std::string solver_error;
    {  // Rf_error longjmp skips C++ dtors; destroy all RAII objects before calling it (R-exts §5.5)
    try {
    // A1b: pre-dispatch validation errors are deferred to here so all std::vectors
    // above are still in scope; throw converts them into solver_error and they are
    // destroyed at '}' below before Rf_error fires (R-exts §5.5). CR-D9 (j7x8.9):
    // this now includes bad-scalar-arg errors — scalar_real/scalar_int(...,
    // pre_error) set pre_error instead of calling Rf_error directly, so EVERY
    // pre-dispatch error path unwinds RAII cleanly (no leak).
    if (!pre_error.empty()) throw std::runtime_error(pre_error);
#ifdef __GLIBC__
    malloc_trim(0);  // compact heap before solver — reduces fragmentation-induced LLC cache misses
#endif
    // SC1 (leafblower-rywn, plan 07 Task 3): p.algorithm was already resolved
    // and validated against kAlgMap above (CR-D10, near method_str's
    // resolution) — C_rk_calibrate now dispatches purely off that resolved
    // enum, the same shape c_api.cpp's rk_calibrate() already used, instead
    // of re-deriving it via a strcmp(method_str, ...) chain. AUTO gets
    // lbw::route_auto() plus up to two lbw::dispatch_solver() calls (primary
    // + fallback); every explicit method gets exactly one.
    const rk_algorithm_t alg = p.algorithm;
    if (alg == RK_ALG_GREENKHORN && st.min_weight >= st.max_weight && st.max_weight > 0.0)
        throw std::runtime_error("greenkhorn requires min_weight < max_weight");

    if (alg == RK_ALG_AUTO) {
        // AUTO routing (Epic-H WH-g): the routing DECISION now lives in the
        // single lbw::route_auto() (calib_dispatch.hpp, SC1 plan 07), shared
        // with c_api.cpp's algorithm-resolution switch. This branch calls
        // it, dispatches once through the shared table, and (on NOCONV/
        // BUDGET) dispatches a SECOND time as the newton_kl fallback — the
        // fallback is now a second dispatch_solver call, not nested routing.
        //   K<5 OR M_cell/n ≤ 0.9 → raking / ORIS (compression-based, unchanged)
        //   K≥5, M_cell/n > 0.9, target_skew ≤ 5 → newton_kl  (moderate skew)
        //   K≥5, M_cell/n > 0.9, target_skew  > 5 → oris+SRAA (severe skew)
        // Save for auto-fallback: only st.weights is mutated by solvers in-place.
        const std::vector<double> weights_backup(weights);
        // CR-D5 (j7x8.5): the severe-skew branch forces st.accelerate=true for ORIS;
        // capture the user's original value to restore before a newton_kl fallback.
        const bool accel_backup = st.accelerate;
        // eb79.15: reuse the site-A m_cell_est_cache when it already ran; route_auto
        // computes fresh (and stores back into the cache) otherwise.
        lbw::AutoRouteResult route = lbw::route_auto(
            st.n, st.K, st.group_ids, st.cat_counts, st.targets, m_cell_est_cache);
        // types.hpp: "true iff AUTO routing selected ORIS" — precise per-algorithm
        // semantics (R's pre-existing, fixture-pinned behavior; adopted by c_api.cpp
        // too in this plan, which previously set this unconditionally for any AUTO
        // pick — see leafblower-* filed alongside this migration).
        st.oris_auto_selected = (route.algorithm == RK_ALG_ORIS);
        if (route.force_accelerate) st.accelerate = true;  // severe-skew ORIS+SRAA
        lbw::dispatch_solver(route.algorithm, st, dres);
        res_status             = dres.status;
        res_iterations         = dres.iterations;
        res_max_error          = dres.max_error;
        res_alg_used           = static_cast<int>(dres.alg_used);
        res_mean_error         = dres.mean_error;
        res_kl                 = dres.kl;
        res_chi2               = dres.chi2;
        res_l1_weight_change   = dres.l1_weight_change;
        res_grake_norm         = dres.grake_norm;
        res_conv_metric        = dres.convergence_metric;
        res_conv_rule          = dres.convergence_rule;
        res_conv_tol           = dres.convergence_tol;
        res_conv_iter          = dres.convergence_iter;
        res_conv_objective     = dres.convergence_solver_objective;
        res_conv_minimized_metric = dres.convergence_minimized_metric;
        res_best_error         = dres.best_error;
        res_best_iter          = dres.best_iter;
        res_metric_first_check = dres.metric_first_check;
        res_metric_prev_check  = dres.metric_prev_check;
        res_prev_check_iter    = dres.prev_check_iter;
        res_stall_kind         = dres.stall_kind;
        res_n_bounds_violated  = dres.n_bounds_violated;
        res_n_bounds_clamped   = dres.n_bounds_clamped;
        std::snprintf(res_solver_message, sizeof(res_solver_message), "%s", dres.solver_message);
        res_n_xcur_writes         = dres.n_xcur_writes_per_iter_last;
        res_min_alpha              = dres.min_alpha_seen;
        res_final_alpha            = dres.final_alpha;
        res_homotopy_levels_used   = dres.homotopy_levels_used;
        res_homotopy_final_factor  = dres.homotopy_final_factor;
        res_greedy_sweeps_taken    = dres.greedy_sweeps_taken;
        res_eta_final              = dres.eta_final;
        res_sor_min_omega      = dres.sor_min_omega;
        res_sor_n_damped       = dres.sor_n_damped;
        res_sor_omega_mean     = dres.sor_omega_mean;
        res_sor_any_latched    = dres.sor_any_latched;
        res_sor_n_pinned_fb    = dres.sor_n_pinned_fb;
        res_sor_n_warmup_fb    = dres.sor_n_warmup_fb;
        res_sor_n_conv_fb      = dres.sor_n_conv_fb;
        res_sor_n_resid_grew   = dres.sor_n_resid_grew;
        res_sor_n_monotone_cd  = dres.sor_n_monotone_cd;
        res_aa_accepted_count  = dres.aa_accepted_count;
        res_sraa_demoted       = dres.sraa_demoted ? 1 : 0;
        res_n_projected_dims   = dres.n_projected_dims;
        res_lm_mu_final        = dres.lm_mu_final;
        res_best_weights       = std::move(dres.best_weights);

        // Auto-fallback: if primary solver NOCONVs or exhausts budget (still
        // improving), retry with newton_kl.  STALL(5) is excluded: the solver
        // is at the constrained optimum and fallback cannot improve it.
        if (res_status == RK_ERR_NOCONV || res_status == RK_ERR_BUDGET) {
            if (st.verbose >= 1)
                st.log("auto: primary solver NOCONV/BUDGET; retrying with newton_kl");
            // Restore original weights (only mutated field in CalibState)
            std::copy(weights_backup.begin(), weights_backup.end(), weights.begin());
            st.oris_auto_selected = false;
            st.accelerate = accel_backup;  // CR-D5: undo severe-skew ORIS override
            lbw::dispatch_solver(RK_ALG_NEWTON_KL, st, dres);
            res_status             = dres.status;
            res_iterations         = dres.iterations;
            res_max_error          = dres.max_error;
            res_alg_used           = static_cast<int>(dres.alg_used);
            res_mean_error         = dres.mean_error;
            res_kl                 = dres.kl;
            res_chi2               = dres.chi2;
            res_l1_weight_change   = dres.l1_weight_change;
            res_grake_norm         = dres.grake_norm;
            res_conv_metric        = dres.convergence_metric;
            res_conv_rule          = dres.convergence_rule;
            res_conv_tol           = dres.convergence_tol;
            res_conv_iter          = dres.convergence_iter;
            res_conv_objective     = dres.convergence_solver_objective;
            res_conv_minimized_metric = dres.convergence_minimized_metric;
            res_best_error         = dres.best_error;
            res_best_iter          = dres.best_iter;
            res_metric_first_check = dres.metric_first_check;
            res_metric_prev_check  = dres.metric_prev_check;
            res_prev_check_iter    = dres.prev_check_iter;
            res_stall_kind         = dres.stall_kind;
            res_n_bounds_violated  = dres.n_bounds_violated;
            res_n_bounds_clamped   = dres.n_bounds_clamped;
            std::snprintf(res_solver_message, sizeof(res_solver_message), "%s", dres.solver_message);
            res_n_projected_dims   = dres.n_projected_dims;
            res_lm_mu_final        = dres.lm_mu_final;
            // CR-D5 (j7x8.5): when the abandoned primary was ORIS, dres still
            // carries its ORIS-only diagnostics (dispatch_solver's NEWTON_KL
            // arm does not touch them); newton owns none of them. Reset to
            // their documented non-ORIS defaults (leafblower.h) so the
            // exported result reflects the winning solver, not stale ORIS
            // state. (n_bounds_violated/clamped ARE reset above — newton's
            // own real values, already copied from the fresh dres.)
            res_n_xcur_writes         = 0;
            res_min_alpha             = 1.0;
            res_final_alpha           = 1.0;
            res_homotopy_levels_used  = 0;
            res_homotopy_final_factor = 1.0;
            res_greedy_sweeps_taken   = 0;
            res_eta_final             = 0.0;
            res_sor_min_omega     = 1.0;
            res_sor_n_damped      = 0;
            res_sor_omega_mean    = 1.0;
            res_sor_any_latched   = 0;
            res_sor_n_pinned_fb   = 0;
            res_sor_n_warmup_fb   = 0;
            res_sor_n_conv_fb     = 0;
            res_sor_n_resid_grew  = 0;
            res_sor_n_monotone_cd = 0;
            res_aa_accepted_count = 0;
            res_sraa_demoted      = 0;
            res_best_weights      = std::move(dres.best_weights);
        }
    } else {
        // Any explicit method (never AUTO-selected): types.hpp's
        // oris_auto_selected contract ("true iff AUTO routing selected
        // ORIS") is always false here. use_admm_capacity is set inside
        // dispatch_solver's RK_ALG_ORIS_SOFT arm (calib_dispatch.hpp, plan
        // 07 Task 3) — no per-method assignment needed on this side.
        st.oris_auto_selected = false;
        lbw::dispatch_solver(alg, st, dres);
        res_status             = dres.status;
        res_iterations         = dres.iterations;
        res_max_error          = dres.max_error;
        res_alg_used           = static_cast<int>(dres.alg_used);
        res_mean_error         = dres.mean_error;
        res_kl                 = dres.kl;
        res_chi2               = dres.chi2;
        res_l1_weight_change   = dres.l1_weight_change;
        res_grake_norm         = dres.grake_norm;
        res_conv_metric        = dres.convergence_metric;
        res_conv_rule          = dres.convergence_rule;
        res_conv_tol           = dres.convergence_tol;
        res_conv_iter          = dres.convergence_iter;
        res_conv_objective     = dres.convergence_solver_objective;
        res_conv_minimized_metric = dres.convergence_minimized_metric;
        res_best_error         = dres.best_error;
        res_best_iter          = dres.best_iter;
        res_metric_first_check = dres.metric_first_check;
        res_metric_prev_check  = dres.metric_prev_check;
        res_prev_check_iter    = dres.prev_check_iter;
        res_stall_kind         = dres.stall_kind;
        res_n_bounds_violated  = dres.n_bounds_violated;
        res_n_bounds_clamped   = dres.n_bounds_clamped;
        // eb79.25/CR-D11: dres.solver_message is already correctly populated
        // (or cleared — RakingResult/ORISResult carry no message field, and
        // dispatch_solver's RAKING/ORIS/ORIS_SOFT arms clear it there) per
        // algorithm by dispatch_solver itself; one snprintf covers every arm.
        std::snprintf(res_solver_message, sizeof(res_solver_message), "%s", dres.solver_message);
        // ORIS/ORIS_SOFT-only diagnostics (n_xcur_writes .. sor_n_monotone_cd):
        // dispatch_solver leaves them at DispatchResult's default-constructed
        // values for every OTHER algorithm — identical to the res_* locals'
        // own top-of-function defaults (documented on the struct itself) —
        // so copying them unconditionally is a no-op there and correct here.
        res_n_xcur_writes         = dres.n_xcur_writes_per_iter_last;
        res_min_alpha              = dres.min_alpha_seen;
        res_final_alpha            = dres.final_alpha;
        res_homotopy_levels_used   = dres.homotopy_levels_used;
        res_homotopy_final_factor  = dres.homotopy_final_factor;
        res_greedy_sweeps_taken    = dres.greedy_sweeps_taken;
        res_eta_final              = dres.eta_final;
        res_sor_min_omega      = dres.sor_min_omega;
        res_sor_n_damped       = dres.sor_n_damped;
        res_sor_omega_mean     = dres.sor_omega_mean;
        res_sor_any_latched    = dres.sor_any_latched;
        res_sor_n_pinned_fb    = dres.sor_n_pinned_fb;
        res_sor_n_warmup_fb    = dres.sor_n_warmup_fb;
        res_sor_n_conv_fb      = dres.sor_n_conv_fb;
        res_sor_n_resid_grew   = dres.sor_n_resid_grew;
        res_sor_n_monotone_cd  = dres.sor_n_monotone_cd;
        res_aa_accepted_count  = dres.aa_accepted_count;
        res_sraa_demoted       = dres.sraa_demoted ? 1 : 0;
        // newton_kl-only diagnostics: same no-op-elsewhere argument.
        res_n_projected_dims   = dres.n_projected_dims;
        res_lm_mu_final        = dres.lm_mu_final;
        // oris_soft-only ALM diagnostics: same no-op-elsewhere argument.
        res_alm_capacity_mu_final = dres.alm_capacity_mu_final;
        res_alm_n_growth_events   = dres.alm_n_growth_events;
        res_alm_max_dual_norm     = dres.alm_max_dual_norm;
        res_alm_sum_drift         = dres.alm_sum_drift;
        // best_weights: dispatch_solver's CHEBYSHEV/NEWTON_KL arms already
        // apply their own empty-check + zero-fill sentinel internally;
        // sinkhorn.cpp/raking.cpp guarantee a length-st.n vector via their
        // own internal fallback (assign(st.n, 0.0) on the violation path),
        // and greg/greenkhorn/logit_calib/oris never leave it empty on any
        // path (verified by reading each solver's best_weights assignment
        // sites) — a plain move reproduces every algorithm's pre-migration
        // behavior exactly.
        res_best_weights       = std::move(dres.best_weights);
    }
    } catch (const std::exception& e) {
        solver_error = e.what();
    } catch (...) {
        solver_error = "unknown exception";
    }
    }
    if (!solver_error.empty()) {
        // CR-D9 (j7x8.9): make this the single leak-free error exit. Rf_error
        // longjmps past C++ dtors (R-exts §5.5), so copy the message to a stack
        // buffer and explicitly release EVERY heap-backed function-scope local
        // before it fires. CalibState/rk_params_t own no heap (verified), so the
        // vectors + strings below are the complete set. Only runs on the error
        // path (solver_error non-empty) — the success path keeps them intact for
        // building the result. pre_error messages are already fully-qualified;
        // solver errors get the "internal solver error" prefix.
        char msg[512];
        if (!pre_error.empty())
            std::snprintf(msg, sizeof(msg), "%s", solver_error.c_str());
        else
            std::snprintf(msg, sizeof(msg),
                          "leafblower: internal solver error \xe2\x80\x94 %s",
                          solver_error.c_str());
        std::string().swap(pre_error);
        std::string().swap(solver_error);
        std::vector<const int32_t*>().swap(group_ids);
        std::vector<int>().swap(cat_counts);
        std::vector<std::vector<double>>().swap(tgt_storage);
        std::vector<const double*>().swap(targets);
        std::vector<double>().swap(weights);
        std::vector<double>().swap(res_best_weights);
        std::vector<double>().swap(dres.best_weights);  // SC1: DispatchResult's heap member
        Rf_error("%s", msg);
    }

    // SC1 (plan 07): the single lbw::kAlgNames table of record
    // (calib_dispatch.hpp) replaces this file's own copy.
    const char* alg_name_cstr = (res_alg_used >= 0 && res_alg_used < lbw::kAlgNamesLen)
        ? lbw::kAlgNames[res_alg_used]
        : "unknown";
    // eb79.25: for ERROR statuses, surface the solver's own message (e.g. logit's
    // structural-INFEAS margin name, a specific BADARG reason) when it set one — instead
    // of the generic 'alg: N iters, max_error' summary that discarded it. Empty-message
    // error paths (sinkhorn/greg INFEAS/BADARG) and non-error statuses keep the summary.
    if ((res_status == RK_ERR_INFEAS || res_status == RK_ERR_BADARG)
        && res_solver_message[0] != '\0') {
        std::snprintf(res_message, 256, "%s", res_solver_message);
    } else {
        std::snprintf(res_message, 256, "%s: %d iters, max_error=%.2e",
                      alg_name_cstr, res_iterations, res_max_error);
    }

    // greenkhorn and logit do not modify st.weights in-place; copy calibrated
    // weights into the weights vector so raw$weights in harvest.R is correct.
    // (raking/oris already write to st.weights in-place — don't copy.)
    if (!res_best_weights.empty() && (int)res_best_weights.size() == n &&
        (res_alg_used == static_cast<int>(RK_ALG_GREENKHORN) ||
         res_alg_used == static_cast<int>(RK_ALG_LOGIT))) {
        std::copy(res_best_weights.begin(), res_best_weights.end(), weights.begin());
    }

    // CR-D9b (j7x8.18): the result-list build below calls Rf_allocVector/Rf_mkChar,
    // any of which can Rf_error (OOM) and longjmp past every function-scope C++ local
    // -> leak. Only `weights` (->wts) and `res_best_weights` (->best_weights field)
    // still feed the result; the input-side locals are dead. Swap-release them now
    // (mirrors the error-site release above), shrinking the OOM-leak surface to just
    // the two still-live result vectors (weights, res_best_weights).
    std::string().swap(pre_error);
    std::string().swap(solver_error);
    std::vector<const int32_t*>().swap(group_ids);
    std::vector<int>().swap(cat_counts);
    std::vector<std::vector<double>>().swap(tgt_storage);
    std::vector<const double*>().swap(targets);
    std::vector<double>().swap(dres.best_weights);  // SC1: DispatchResult is dead past this point

    // Build return list: list(weights=numeric[n], result=list(42 fields))
    // PROTECT out first so wts/res_list/res_names sit above it on the stack;
    // each is UNPROTECTed immediately after adoption, in LIFO order.
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));  // owns wts + res_list below

    SEXP wts = PROTECT(Rf_allocVector(REALSXP, n));
    memcpy(REAL(wts), weights.data(), (size_t)n * sizeof(double));

    constexpr int N_RESULT_FIELDS = 49;
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
    /* Elements 22-28: ORIS SOR observability fields (e18t.4) */
    SET_STRING_ELT(res_names, 22, Rf_mkChar("sor_omega_mean"));
    SET_STRING_ELT(res_names, 23, Rf_mkChar("sor_any_latched"));
    SET_STRING_ELT(res_names, 24, Rf_mkChar("sor_n_pinned_fb"));
    SET_STRING_ELT(res_names, 25, Rf_mkChar("sor_n_warmup_fb"));
    SET_STRING_ELT(res_names, 26, Rf_mkChar("sor_n_conv_fb"));
    SET_STRING_ELT(res_names, 27, Rf_mkChar("sor_n_resid_grew"));
    SET_STRING_ELT(res_names, 28, Rf_mkChar("sor_n_monotone_cd"));
    SET_VECTOR_ELT(res_list,  22, Rf_ScalarReal(res_sor_omega_mean));
    SET_VECTOR_ELT(res_list,  23, Rf_ScalarInteger(res_sor_any_latched));
    SET_VECTOR_ELT(res_list,  24, Rf_ScalarInteger(res_sor_n_pinned_fb));
    SET_VECTOR_ELT(res_list,  25, Rf_ScalarInteger(res_sor_n_warmup_fb));
    SET_VECTOR_ELT(res_list,  26, Rf_ScalarInteger(res_sor_n_conv_fb));
    SET_VECTOR_ELT(res_list,  27, Rf_ScalarInteger(res_sor_n_resid_grew));
    SET_VECTOR_ELT(res_list,  28, Rf_ScalarInteger(res_sor_n_monotone_cd));
    /* Element 29: best_weights REALSXP (WU-E) */
    SET_STRING_ELT(res_names, 29, Rf_mkChar("best_weights"));
    {
        int bw_n = (int)res_best_weights.size();
        SEXP bw_sxp = PROTECT(Rf_allocVector(REALSXP, bw_n));
        double* bw  = REAL(bw_sxp);
        // eb79.14: memcpy the contiguous REALSXP from the contiguous std::vector<double>
        // (matches the nearby weights memcpy); 0-length is a no-op.
        if (bw_n > 0) std::memcpy(bw, res_best_weights.data(), (size_t)bw_n * sizeof(double));
        SET_VECTOR_ELT(res_list, 29, bw_sxp);
        UNPROTECT(1);  // bw_sxp adopted by res_list
    }
    /* Elements 30-34: convergence diagnostics (WU-A) */
    SET_STRING_ELT(res_names, 30, Rf_mkChar("grake_norm"));
    SET_STRING_ELT(res_names, 31, Rf_mkChar("convergence_metric"));
    SET_STRING_ELT(res_names, 32, Rf_mkChar("convergence_rule"));
    SET_STRING_ELT(res_names, 33, Rf_mkChar("convergence_tol"));
    SET_STRING_ELT(res_names, 34, Rf_mkChar("convergence_iter"));
    SET_VECTOR_ELT(res_list, 30, Rf_ScalarReal(res_grake_norm));
    SET_VECTOR_ELT(res_list, 31, Rf_ScalarInteger(res_conv_metric));
    SET_VECTOR_ELT(res_list, 32, Rf_ScalarInteger(res_conv_rule));
    SET_VECTOR_ELT(res_list, 33, Rf_ScalarReal(res_conv_tol));
    SET_VECTOR_ELT(res_list, 34, Rf_ScalarInteger(res_conv_iter));
    /* Elements 35-36: convergence_objective and convergence_minimized_metric (Task 1) */
    SET_STRING_ELT(res_names, 35, Rf_mkChar("solver_objective"));
    SET_STRING_ELT(res_names, 36, Rf_mkChar("convergence_minimized_metric"));
    SET_VECTOR_ELT(res_list,  35, Rf_ScalarReal(res_conv_objective));
    SET_VECTOR_ELT(res_list,  36, Rf_ScalarInteger(res_conv_minimized_metric));
    /* Elements 37-40: ALM diagnostics (non-zero only for oris_soft) */
    SET_STRING_ELT(res_names, 37, Rf_mkChar("alm_capacity_mu_final"));
    SET_STRING_ELT(res_names, 38, Rf_mkChar("alm_n_growth_events"));
    SET_STRING_ELT(res_names, 39, Rf_mkChar("alm_max_dual_norm"));
    SET_STRING_ELT(res_names, 40, Rf_mkChar("alm_sum_drift"));
    SET_VECTOR_ELT(res_list,  37, Rf_ScalarReal(res_alm_capacity_mu_final));
    SET_VECTOR_ELT(res_list,  38, Rf_ScalarInteger(res_alm_n_growth_events));
    SET_VECTOR_ELT(res_list,  39, Rf_ScalarReal(res_alm_max_dual_norm));
    SET_VECTOR_ELT(res_list,  40, Rf_ScalarReal(res_alm_sum_drift));
    /* Element 41: SRAA acceleration diagnostic (oris/oris_soft only; zero elsewhere) */
    SET_STRING_ELT(res_names, 41, Rf_mkChar("aa_accepted_count"));
    SET_VECTOR_ELT(res_list,  41, Rf_ScalarInteger(res_aa_accepted_count));
    /* Element 42: Newton-KL TSVD diagnostic (newton_kl only; zero elsewhere) */
    SET_STRING_ELT(res_names, 42, Rf_mkChar("n_projected_dims"));
    SET_VECTOR_ELT(res_list,  42, Rf_ScalarInteger(res_n_projected_dims));
    /* Element 43: Newton-KL Levenberg-Marquardt diagnostic (newton_kl only; zero elsewhere) */
    SET_STRING_ELT(res_names, 43, Rf_mkChar("lm_mu_final"));
    SET_VECTOR_ELT(res_list,  43, Rf_ScalarReal(res_lm_mu_final));
    /* Element 44: first-check metric value (oris only; Inf elsewhere) */
    SET_STRING_ELT(res_names, 44, Rf_mkChar("metric_first_check"));
    SET_VECTOR_ELT(res_list,  44, Rf_ScalarReal(res_metric_first_check));
    SET_STRING_ELT(res_names, 45, Rf_mkChar("metric_prev_check"));
    SET_VECTOR_ELT(res_list,  45, Rf_ScalarReal(res_metric_prev_check));
    SET_STRING_ELT(res_names, 46, Rf_mkChar("prev_check_iter"));
    SET_VECTOR_ELT(res_list,  46, Rf_ScalarInteger(res_prev_check_iter));
    /* Element 47: SRAA scheduler-demotion flag (oris/raking only; FALSE elsewhere) */
    SET_STRING_ELT(res_names, 47, Rf_mkChar("sraa_demoted"));
    SET_VECTOR_ELT(res_list,  47, Rf_ScalarLogical(res_sraa_demoted));
    /* Element 48: solver-emitted stall kind (leafblower-8eod).
       0=no stall, 1=wchange (SRAA path), 2=kl (plain-IPF path).
       Set at RK_ERR_STALL emission site; replaces accelerate_bool heuristic in harvest.R. */
    SET_STRING_ELT(res_names, 48, Rf_mkChar("convergence_stall_kind"));
    SET_VECTOR_ELT(res_list,  48, Rf_ScalarInteger(res_stall_kind));
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
    if (TYPEOF(r_group_ids_list) != VECSXP)
        Rf_error("group_ids must be a list");
    int n = scalar_int(r_n, "r_n");
    int K = Rf_length(r_group_ids_list);
    // CR-D8 (j7x8.8): validate ALL inputs in a leak-free pre-pass (nothing
    // allocated yet) before any dereference. This registered entry point is
    // .Call-able directly, bypassing the R wrapper's coercion — a short or
    // non-INTSXP vector previously OOB-read (the loop reads INTEGER(v)[i] for
    // i<n) or hit an accessor error, and a negative/huge n fed a std::vector
    // allocation. Mirrors the main entry's CXX.1 discipline.
    if (n < 0)
        Rf_error("n must be non-negative, got %d", n);
    if (K == 0 && n > 0)
        Rf_error("group_ids is empty but n (%d) > 0", n);
    if (K > lbw::K_MAX)
        Rf_error("K (%d) exceeds K_MAX (%d)", K, lbw::K_MAX);
    for (int k = 0; k < K; k++) {
        SEXP v = VECTOR_ELT(r_group_ids_list, k);
        if (TYPEOF(v) != INTSXP)
            Rf_error("group_ids[[%d]] must be an integer vector", k + 1);
        if (Rf_length(v) < n)
            Rf_error("group_ids[[%d]] length (%d) < n (%d)", k + 1,
                     (int) Rf_length(v), n);
    }
    // Extract pointers (all inputs validated above).
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
        // j7x8.24: n*K is a long-vector length. Rf_length returns int and TRUNCATES a
        // vector longer than INT_MAX (exactly the overflow case), and int n*K overflows.
        // Use XLENGTH (R_xlen_t) on both sides; (R_xlen_t)n * K (n,K <= INT_MAX) fits the
        // signed 64-bit R_xlen_t without overflow.
        const R_xlen_t expected_len = static_cast<R_xlen_t>(n) * static_cast<R_xlen_t>(K);
        if (XLENGTH(data_codes_sexp) != expected_len)
            Rf_error("design_effect: length(data_codes) must be n*K = %lld (got %lld)",
                     static_cast<long long>(expected_len),
                     static_cast<long long>(XLENGTH(data_codes_sexp)));
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
