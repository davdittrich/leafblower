#include "leafblower.h"
#include "validation.hpp"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <cstring>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>
#include "logit.hpp"
#include "cell_table.hpp"
#include "types.hpp"
#include "ieppa.hpp"
#include "raking.hpp"
#include "sinkhorn.hpp"
#include "greg.hpp"
#include "chebyshev.hpp"
#include "lbfgsb_solver.hpp"
#include "greenkhorn.hpp"
#include "logit_calib.hpp"
#include "newton_calib.hpp"

namespace {
const std::unordered_map<std::string_view, rk_algorithm_t> kAlgMap = {
    {"ieppa",      RK_ALG_IEPPA},
    {"ieppa_soft", RK_ALG_IEPPA_SOFT},
    {"lbfgsb",     RK_ALG_LBFGSB},
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
                    SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_leafblower_cell_table_probe(SEXP, SEXP);
}

// R log trampoline: forwards CalibState.log() calls to Rprintf
static void r_log_trampoline(const char* msg, void* /*ctx*/) {
    Rprintf("[leafblower] %s\n", msg);
}

extern "C" {

void R_init_leafblower(DllInfo* dll) {
    static const R_CallMethodDef call_methods[] = {
        {"C_logit_F_at_zero",    (DL_FUNC)&C_logit_F_at_zero,    2},
        {"C_logit_range_check",  (DL_FUNC)&C_logit_range_check,  3},
        {"C_logit_Hprime_check", (DL_FUNC)&C_logit_Hprime_check, 3},
        {"C_rk_calibrate",       (DL_FUNC)&C_rk_calibrate,       35},
        {"C_leafblower_cell_table_probe", (DL_FUNC)&C_leafblower_cell_table_probe, 2},
        {NULL, NULL, 0}
    };
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

// Test bridge: return F(0) for given L, U
SEXP C_logit_F_at_zero(SEXP Lsxp, SEXP Usxp) {
    double L = REAL(Lsxp)[0];
    double U = REAL(Usxp)[0];
    lbw::LinkFn fn(L, U);
    return Rf_ScalarReal(fn.F(0.0));
}

// Test bridge: return F(u) for a vector of u values
SEXP C_logit_range_check(SEXP Lsxp, SEXP Usxp, SEXP usxp) {
    double L = REAL(Lsxp)[0];
    double U = REAL(Usxp)[0];
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
    double L = REAL(Lsxp)[0];
    double U = REAL(Usxp)[0];
    double u0 = REAL(u0sxp)[0];
    double h = 1e-7;
    lbw::LinkFn fn(L, U);
    double Hprime_numerical = (fn.H(u0 + h) - fn.H(u0 - h)) / (2.0 * h);
    double diff = std::fabs(Hprime_numerical - fn.F(u0));
    return Rf_ScalarReal(diff);
}

// Main calibration bridge.
// data_sexp: data.frame with factor/character columns
// target_sexp: named list of named numeric vectors (variable -> proportions)
// Returns: list(weights=double[n], result=list(status, iterations, max_error, algorithm_used, message))
SEXP C_rk_calibrate(SEXP data_sexp, SEXP target_sexp,
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
                    /* Jacobi log-path sweep flag (passed by R harvest(); currently unused in
                     * bridge but consumed by R registry to keep positional alignment with newer
                     * args appended below). */
                    SEXP jacobi_sweep_sexp,
                    /* Epic-H WH-e: newton_kl TSVD truncation ratio (default 1e-8 from R layer). */
                    SEXP newton_tsvd_ratio_sexp) {
    (void)jacobi_sweep_sexp;  // Silenced: already-existing pass-through, no consumer yet.
    int n = Rf_nrows(VECTOR_ELT(data_sexp, 0));
    SEXP target_names = Rf_getAttrib(target_sexp, R_NamesSymbol);
    int K = LENGTH(target_sexp);

    // Build group_ids and targets from data.frame factor/character columns
    std::vector<std::vector<int32_t>> gids_storage(K);
    std::vector<const int32_t*> group_ids(K);
    std::vector<int> cat_counts(K);
    std::vector<std::vector<double>> tgt_storage(K);
    std::vector<const double*> targets(K);

    SEXP col_names = Rf_getAttrib(data_sexp, R_NamesSymbol);

    // Build column name → index map once (O(C)) for O(1) per-margin lookup.
    std::unordered_map<std::string, int> col_map;
    col_map.reserve(LENGTH(col_names));
    for (int c = 0; c < LENGTH(col_names); c++)
        col_map[CHAR(STRING_ELT(col_names, c))] = c;

    for (int k = 0; k < K; k++) {
        const char* varname = CHAR(STRING_ELT(target_names, k));
        // Find column in data frame
        auto col_it = col_map.find(varname);
        if (col_it == col_map.end())
            Rf_error("Variable '%s' not found in data", varname);
        int col_idx = col_it->second;

        SEXP col     = VECTOR_ELT(data_sexp, col_idx);
        SEXP tgt_vec = VECTOR_ELT(target_sexp, k);
        SEXP tgt_names = Rf_getAttrib(tgt_vec, R_NamesSymbol);
        int ncat = LENGTH(tgt_vec);
        cat_counts[k] = ncat;

        tgt_storage[k].resize(ncat);
        for (int j = 0; j < ncat; j++) tgt_storage[k][j] = REAL(tgt_vec)[j];
        targets[k] = tgt_storage[k].data();

        // O(1) level-name -> index map
        std::unordered_map<std::string, int> level_to_idx;
        level_to_idx.reserve(ncat);
        for (int j = 0; j < ncat; j++)
            level_to_idx[CHAR(STRING_ELT(tgt_names, j))] = j;

        // Encode factor/character -> 0-indexed int (NA -> -1, OOV -> -1)
        gids_storage[k].resize(n);
        if (Rf_isFactor(col)) {
            SEXP flevels = Rf_getAttrib(col, R_LevelsSymbol);
            const int* codes = INTEGER(col);
            for (int i = 0; i < n; i++) {
                if (codes[i] == NA_INTEGER) { gids_storage[k][i] = -1; continue; }
                const char* lv = CHAR(STRING_ELT(flevels, codes[i] - 1));
                auto it = level_to_idx.find(lv);
                gids_storage[k][i] = (it != level_to_idx.end()) ? it->second : -1;
            }
        } else if (TYPEOF(col) == STRSXP) {
            for (int i = 0; i < n; i++) {
                if (STRING_ELT(col, i) == NA_STRING) { gids_storage[k][i] = -1; continue; }
                const char* lv = CHAR(STRING_ELT(col, i));
                auto it = level_to_idx.find(lv);
                gids_storage[k][i] = (it != level_to_idx.end()) ? it->second : -1;
            }
        } else {
            Rf_error("Column '%s' must be a factor or character vector", varname);
        }
        group_ids[k] = gids_storage[k].data();
    }

    // Build weights (start_weights already normalized to mean=1 by R layer)
    std::vector<double> weights(n);
    if (Rf_isNull(start_weights_sexp)) {
        for (int i = 0; i < n; i++) weights[i] = 1.0;
    } else {
        if (LENGTH(start_weights_sexp) != n)
            Rf_error("leafblower: start_weights length %d != n=%d",
                     (int)LENGTH(start_weights_sexp), n);
        const double* sw = REAL(start_weights_sexp);
        for (int i = 0; i < n; i++) weights[i] = sw[i];
    }

    // Set calibration params
    rk_params_t p;
    rk_params_init(&p);
    p.min_weight     = REAL(min_weight_sexp)[0];
    p.max_weight     = REAL(max_weight_sexp)[0];
    p.verbose        = INTEGER(verbose_sexp)[0];
    p.inner_max_iter = INTEGER(inner_max_iter_sexp)[0];
    // outer_max_iter = max_iterations: user controls both iEPPA inner BCD and L-BFGS-B
    // outer step budget via the same parameter. With the O(n) Wolfe inner loop, 500
    // outer steps costs ~500 O(K*n) grad evals — acceptable; typical convergence is <50.
    p.outer_max_iter = INTEGER(inner_max_iter_sexp)[0];
    p.tol_abs        = REAL(tol_abs_sexp)[0];
    p.bounds_mode    = (rk_bounds_mode_t) INTEGER(bounds_mode_sexp)[0];
    p.log_fn         = (p.verbose > 0) ? r_log_trampoline : nullptr;
    /* Overlay knobs */
    p.homotopy.n_levels        = INTEGER(homotopy_levels_sexp)[0];
    p.homotopy.start_factor    = REAL(homotopy_start_factor_sexp)[0];
    p.homotopy.end_factor      = REAL(homotopy_end_factor_sexp)[0];
    p.homotopy.budget_split_p  = REAL(homotopy_budget_p_sexp)[0];
    {
        const char* sched_str = CHAR(STRING_ELT(scheduler_sexp, 0));
        p.scheduler = (strcmp(sched_str, "greedy") == 0) ? RK_SCHED_GREEDY : RK_SCHED_ROUND_ROBIN;
    }
    {
        const char* eta_str = CHAR(STRING_ELT(eta_schedule_sexp, 0));
        p.eta_mode = (strcmp(eta_str, "tang_dynamic") == 0) ? RK_ETA_TANG_DYNAMIC : RK_ETA_FIXED;
    }
    p.eta_start           = REAL(eta_start_sexp)[0];
    p.eta_end             = REAL(eta_end_sexp)[0];
    p.eta_schedule_power  = REAL(eta_schedule_power_sexp)[0];
    /* Convergence config (WU-A) */
    p.pct_tol             = REAL(pct_tol_sexp)[0];
    p.absolute_tol        = REAL(absolute_tol_sexp)[0];
    p.metric              = INTEGER(metric_sexp)[0];
    p.rule                = INTEGER(rule_sexp)[0];
    p.stop_when           = INTEGER(stop_when_sexp)[0];
    /* SOR config (WU-A) */
    p.sor_enabled         = INTEGER(sor_enabled_sexp)[0];
    p.sor_auto            = INTEGER(sor_auto_sexp)[0];
    p.sor_omega_init      = REAL(sor_omega_init_sexp)[0];
    p.sor_omega_min       = REAL(sor_omega_min_sexp)[0];
    p.sor_omega_fixed     = REAL(sor_omega_fixed_sexp)[0];
    p.sor_burnin          = INTEGER(sor_burnin_sexp)[0];

    if (LENGTH(method_sexp) != 1)
        Rf_error("method must be a length-1 character string");
    const char* method_str = CHAR(STRING_ELT(method_sexp, 0));
    const auto alg_it = kAlgMap.find(method_str);
    p.algorithm = (alg_it != kAlgMap.end()) ? alg_it->second : RK_ALG_IEPPA;

    // Full input validation — shared with c_api.cpp path via validation.hpp.
    {
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
        if (vrc != RK_OK)
            Rf_error("leafblower: invalid arguments \342\200\224 %s",
                     validation_result.message);
    }

    // WU-E: call C++ solvers directly (bypasses flat C ABI) to access best_weights vector.
    // Build CalibState mirroring c_api.cpp:rk_calibrate() setup.
    lbw::CalibState st;
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
    st.lbfgs_m       = p.lbfgs_m;
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
    // Epic-H WH-e: newton_kl TSVD truncation ratio. <=0 falls back to internal default 1e-8.
    st.newton_tsvd_ratio = REAL(newton_tsvd_ratio_sexp)[0];
    st.ieppa_auto_selected          = false;  // R bridge always resolves method explicitly
    st.alm.lambda = 0.0;
    st.alm.mu     = 0.0;
    st.accelerate = (INTEGER(accelerate_sexp)[0] != 0);

    // Resolve capacity_mu for ieppa_soft: build cell table to obtain auto value.
    // Done here (not inside solver) so r_bridge.cpp controls the resolution contract.
    // capacity_penalty <= 0.0 (or NULL) selects auto (M_cell/n from cell_table);
    // positive value is used directly. ALM block is gated by st.use_admm_capacity,
    // so st.alm.capacity_mu is harmless for non-ieppa_soft callers.
    {
        lbw::CellTable ct_tmp;
        int ct_rc = lbw::build_cell_table(n, K,
            const_cast<const int32_t**>(group_ids.data()),
            cat_counts.data(),
            weights.data(),
            ct_tmp);
        if (ct_rc != 0) {
            // build_cell_table failure is non-fatal here — capacity_mu falls back to 1.0.
            // Infeasibility (empty cell with positive target) is caught by validate_inputs above.
            st.alm.capacity_mu = 1.0;
        } else {
            const double cp_val = Rf_isNull(capacity_penalty_sexp)
                ? -1.0
                : (LENGTH(capacity_penalty_sexp) == 1 ? REAL(capacity_penalty_sexp)[0] : -1.0);
            st.alm.capacity_mu = (cp_val <= 0.0) ? ct_tmp.capacity_mu_auto : cp_val;
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
    int    res_alg_used   = (int)RK_ALG_IEPPA;
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
    /* ALM diagnostics (populated only in ieppa_soft dispatch; zero elsewhere) */
    double res_alm_capacity_mu_final   = 0.0;
    int    res_alm_n_growth_events     = 0;
    double res_alm_max_dual_norm       = 0.0;
    double res_alm_sum_drift           = 0.0;
    /* Acceleration (SRAA) diagnostic (ieppa/ieppa_soft only; zero elsewhere) */
    int    res_aa_accepted_count       = 0;
    /* Newton-KL TSVD diagnostic (Epic-Dβ WL-1; non-zero only for newton_kl) */
    int    res_n_projected_dims        = 0;
    /* Newton-KL Levenberg-Marquardt diagnostic (newton_kl only; zero elsewhere) */
    double res_lm_mu_final             = 0.0;
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
    };

    try {
    if (strcmp(method_str, "lbfgsb") == 0) {
        auto res = lbw::lbfgsb_solve(st);
        res_status     = res.base.status;
        res_iterations = res.base.iterations;
        res_max_error  = res.base.max_error;
        res_alg_used   = (int)RK_ALG_LBFGSB;
        pack_solver_result(res);
        res_best_weights = std::move(res.base.best_weights);
    } else if (strcmp(method_str, "raking") == 0) {
        auto res = lbw::raking_solve(st);
        res_status     = res.base.status;
        res_iterations = res.base.iterations;
        res_max_error  = res.base.max_error;
        res_alg_used   = (int)RK_ALG_RAKING;
        pack_solver_result(res);
        res_best_weights = std::move(res.base.best_weights);
    } else if (strcmp(method_str, "auto") == 0) {
        // AUTO routing (Epic-H WH-g):
        //   K<5 OR M_cell/n ≤ 0.9 → raking / iEPPA (compression-based, unchanged)
        //   K≥5, M_cell/n > 0.9, target_skew ≤ 5 → newton_kl  (moderate skew)
        //   K≥5, M_cell/n > 0.9, target_skew  > 5 → ieppa+SRAA (severe skew)
        //   target_skew = max(targets) / max(min(targets), 1e-12)
        int M_cell_est = lbw::estimate_M_cell(n, K,
            const_cast<const int32_t**>(group_ids.data()),
            cat_counts.data());
        // Exact integer comparison: M_cell_est / n > 0.9  ↔  M_cell_est * 10 > n * 9
        const bool zero_compression = (static_cast<int64_t>(M_cell_est) * 10 > static_cast<int64_t>(n) * 9);
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
                // on severe-skew dual landscape; iEPPA+SRAA converges instead.
                st.ieppa_auto_selected = true;
                st.accelerate = true;
                auto res = lbw::ieppa_solve(st);
                res_status     = res.base.status;
                res_iterations = res.base.iterations;
                res_max_error  = res.base.max_error;
                res_alg_used   = (int)RK_ALG_IEPPA;
                res_n_xcur_writes         = res.n_xcur_writes_per_iter_linear;
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
            res_best_weights = std::move(res.base.best_weights);
        } else {
            // Compressed regime: iEPPA (any K)
            st.ieppa_auto_selected = true;
            auto res = lbw::ieppa_solve(st);
            res_status     = res.base.status;
            res_iterations = res.base.iterations;
            res_max_error  = res.base.max_error;
            res_alg_used   = (int)RK_ALG_IEPPA;
            res_n_xcur_writes         = res.n_xcur_writes_per_iter_linear;
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
            res_best_weights = std::move(res.base.best_weights);
        }
        // Auto-fallback: if primary solver NOCONVs or exhausts budget (still
        // improving), retry with L-BFGS-B.  STALL(5) is excluded: the solver
        // is at the constrained optimum and fallback cannot improve it.
        if (res_status == RK_ERR_NOCONV || res_status == RK_ERR_BUDGET) {
            if (st.verbose >= 1)
                st.log("auto: primary solver NOCONV/BUDGET; retrying with L-BFGS-B");
            // Restore original weights (only mutated field in CalibState)
            std::copy(weights_backup.begin(), weights_backup.end(), weights.begin());
            st.ieppa_auto_selected = false;
            auto fb = lbw::lbfgsb_solve(st);
            res_status     = fb.base.status;
            res_iterations = fb.base.iterations;
            res_max_error  = fb.base.max_error;
            res_alg_used   = (int)RK_ALG_LBFGSB;
            pack_solver_result(fb);
            res_best_weights = std::move(fb.base.best_weights);
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
            Rf_error("leafblower: greenkhorn requires min_weight < max_weight");
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
        res_n_projected_dims = res.n_projected_dims;
        res_lm_mu_final      = res.lm_mu_final;
        if (!res.base.best_weights.empty())
            res_best_weights = std::move(res.base.best_weights);
        else
            res_best_weights.assign(st.n, 0.0);
    } else {
        // Dispatch for chebyshev, ieppa_soft, and default ieppa.

        // Run ieppa warm-start BEFORE chebyshev dispatch.
        std::vector<double> w_warm_obs;
        double delta_warm = -1.0;
        if (strcmp(method_str, "chebyshev") == 0) {
            // SAFETY: weights_copy protects st.weights from ieppa_solve mutation.
            std::vector<double> weights_copy(st.weights, st.weights + st.n);
            lbw::CalibState st_warm = st;
            st_warm.weights = weights_copy.data();
            st_warm.inner_max_iter = std::max(5, std::min(100, st.inner_max_iter / 10));
            auto ieppa_res = lbw::ieppa_solve(st_warm);
            // st_warm must NOT escape this block (dangling pointer after weights_copy destroyed).
            if (!ieppa_res.base.best_weights.empty() &&
                static_cast<int>(ieppa_res.base.best_weights.size()) == st.n &&
                std::isfinite(ieppa_res.base.max_error)) {
                w_warm_obs = std::move(ieppa_res.base.best_weights);
                delta_warm = ieppa_res.base.max_error * 1.5;
            }
        }

        auto dispatch_cheb = [&](lbw::LpVariant variant, int alg_code) {
            const std::vector<double>& warm_ref = w_warm_obs;
            const double d_warm = delta_warm;
            auto res = lbw::chebyshev_ipm(st, variant, warm_ref, d_warm);
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
            dispatch_cheb(lbw::LpVariant::CHEBYSHEV, static_cast<int>(RK_ALG_CHEBYSHEV));
        } else if (strcmp(method_str, "ieppa_soft") == 0) {
            st.ieppa_auto_selected = false;
            st.use_admm_capacity   = true;
            auto res = lbw::ieppa_solve(st);
            res_status     = res.base.status;
            res_iterations = res.base.iterations;
            res_max_error  = res.base.max_error;
            res_alg_used   = (int)RK_ALG_IEPPA_SOFT;
            res_n_xcur_writes         = res.n_xcur_writes_per_iter_linear;
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
            res_best_weights = std::move(res.base.best_weights);
        } else {
            // Default / ieppa
            st.ieppa_auto_selected = (strcmp(method_str, "ieppa") != 0);
            auto res = lbw::ieppa_solve(st);
            res_status     = res.base.status;
            res_iterations = res.base.iterations;
            res_max_error  = res.base.max_error;
            res_alg_used   = (int)RK_ALG_IEPPA;
            res_n_xcur_writes         = res.n_xcur_writes_per_iter_linear;
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
            res_best_weights = std::move(res.base.best_weights);
        }
    }
    } catch (const std::exception& e) {
        // N=0 PROTECTs before this try block — no UNPROTECT needed.
        // Rf_error longjmps; must not be inside a try that owns RAII objects.
        Rf_error("leafblower: internal solver error — %s", e.what());
    } catch (...) {
        // Catch non-std::exception throws (raw types, third-party exceptions).
        Rf_error("leafblower: internal solver error — unknown exception type");
    }

    const char* alg_name = (res_alg_used == (int)RK_ALG_LBFGSB)                       ? "L-BFGS-B"
                         : (res_alg_used == (int)RK_ALG_RAKING)                       ? "raking"
                         : (res_alg_used == static_cast<int>(RK_ALG_SINKHORN))        ? "sinkhorn"
                         : (res_alg_used == static_cast<int>(RK_ALG_GREG))            ? "greg"
                         : (res_alg_used == static_cast<int>(RK_ALG_CHEBYSHEV))       ? "chebyshev"
                         : (res_alg_used == static_cast<int>(RK_ALG_GREENKHORN))      ? "greenkhorn"
                         : (res_alg_used == static_cast<int>(RK_ALG_LOGIT))           ? "logit"
                         : (res_alg_used == static_cast<int>(RK_ALG_NEWTON_KL))      ? "newton_kl"
                         : "iEPPA";
    std::snprintf(res_message, 256, "%s: %d iters, max_error=%.2e",
                  alg_name, res_iterations, res_max_error);

    // greenkhorn and logit do not modify st.weights in-place; copy calibrated
    // weights into the weights vector so raw$weights in harvest.R is correct.
    // (raking/ieppa/lbfgsb already write to st.weights in-place — don't copy.)
    if (!res_best_weights.empty() && (int)res_best_weights.size() == n &&
        (res_alg_used == static_cast<int>(RK_ALG_GREENKHORN) ||
         res_alg_used == static_cast<int>(RK_ALG_LOGIT))) {
        std::copy(res_best_weights.begin(), res_best_weights.end(), weights.begin());
    }

    // Build return list: list(weights=numeric[n], result=list(23 fields))
    SEXP wts = PROTECT(Rf_allocVector(REALSXP, n));
    memcpy(REAL(wts), weights.data(), (size_t)n * sizeof(double));

    SEXP res_list  = PROTECT(Rf_allocVector(VECSXP,  37));  // 14 prior + 8 scalars + best_weights + 7 convergence fields + 4 ALM diagnostics + 1 SRAA diagnostic + 1 Newton-KL TSVD diagnostic + 1 Newton-KL LM diagnostic
    SEXP res_names = PROTECT(Rf_allocVector(STRSXP,  37));
    SET_STRING_ELT(res_names, 0, Rf_mkChar("status"));
    SET_STRING_ELT(res_names, 1, Rf_mkChar("iterations"));
    SET_STRING_ELT(res_names, 2, Rf_mkChar("max_error"));
    SET_STRING_ELT(res_names, 3, Rf_mkChar("algorithm_used"));
    SET_STRING_ELT(res_names, 4, Rf_mkChar("message"));
    SET_STRING_ELT(res_names, 5, Rf_mkChar("n_xcur_writes_per_iter_linear"));
    SET_STRING_ELT(res_names, 6, Rf_mkChar("min_alpha_seen"));
    SET_STRING_ELT(res_names, 7, Rf_mkChar("final_alpha"));
    SET_STRING_ELT(res_names, 8, Rf_mkChar("n_bounds_violated"));
    SET_STRING_ELT(res_names, 9, Rf_mkChar("n_bounds_clamped"));
    SET_VECTOR_ELT(res_list, 0, Rf_ScalarInteger(res_status));
    SET_VECTOR_ELT(res_list, 1, Rf_ScalarInteger(res_iterations));
    SET_VECTOR_ELT(res_list, 2, Rf_ScalarReal(res_max_error));
    SET_VECTOR_ELT(res_list, 3, Rf_ScalarInteger(res_alg_used));
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
    /* Elements 30-33: ALM diagnostics (non-zero only for ieppa_soft) */
    SET_STRING_ELT(res_names, 30, Rf_mkChar("alm_capacity_mu_final"));
    SET_STRING_ELT(res_names, 31, Rf_mkChar("alm_n_growth_events"));
    SET_STRING_ELT(res_names, 32, Rf_mkChar("alm_max_dual_norm"));
    SET_STRING_ELT(res_names, 33, Rf_mkChar("alm_sum_drift"));
    SET_VECTOR_ELT(res_list,  30, Rf_ScalarReal(res_alm_capacity_mu_final));
    SET_VECTOR_ELT(res_list,  31, Rf_ScalarInteger(res_alm_n_growth_events));
    SET_VECTOR_ELT(res_list,  32, Rf_ScalarReal(res_alm_max_dual_norm));
    SET_VECTOR_ELT(res_list,  33, Rf_ScalarReal(res_alm_sum_drift));
    /* Element 34: SRAA acceleration diagnostic (ieppa/ieppa_soft only; zero elsewhere) */
    SET_STRING_ELT(res_names, 34, Rf_mkChar("aa_accepted_count"));
    SET_VECTOR_ELT(res_list,  34, Rf_ScalarInteger(res_aa_accepted_count));
    /* Element 35: Newton-KL TSVD diagnostic (newton_kl only; zero elsewhere) */
    SET_STRING_ELT(res_names, 35, Rf_mkChar("n_projected_dims"));
    SET_VECTOR_ELT(res_list,  35, Rf_ScalarInteger(res_n_projected_dims));
    /* Element 36: Newton-KL Levenberg-Marquardt diagnostic (newton_kl only; zero elsewhere) */
    SET_STRING_ELT(res_names, 36, Rf_mkChar("lm_mu_final"));
    SET_VECTOR_ELT(res_list,  36, Rf_ScalarReal(res_lm_mu_final));
    Rf_setAttrib(res_list, R_NamesSymbol, res_names);

    SEXP out       = PROTECT(Rf_allocVector(VECSXP,  2));
    SEXP out_names = PROTECT(Rf_allocVector(STRSXP,  2));
    SET_STRING_ELT(out_names, 0, Rf_mkChar("weights"));
    SET_STRING_ELT(out_names, 1, Rf_mkChar("result"));
    SET_VECTOR_ELT(out, 0, wts);
    SET_VECTOR_ELT(out, 1, res_list);
    Rf_setAttrib(out, R_NamesSymbol, out_names);
    UNPROTECT(5);
    return out;
}

// test-only: exposes CellTable internals for unit tests
extern "C" SEXP C_leafblower_cell_table_probe(SEXP r_group_ids_list, SEXP r_n) {
    int n = INTEGER(r_n)[0];
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

} // extern "C"
