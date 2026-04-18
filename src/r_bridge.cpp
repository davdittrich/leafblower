#include "leafblower.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>
#include "logit.hpp"

extern "C" {
SEXP C_logit_F_at_zero(SEXP, SEXP);
SEXP C_logit_range_check(SEXP, SEXP, SEXP);
SEXP C_logit_Hprime_check(SEXP, SEXP, SEXP);
SEXP C_rk_calibrate(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
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
        {"C_rk_calibrate",       (DL_FUNC)&C_rk_calibrate,       8},
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
                    SEXP inner_max_iter_sexp, SEXP start_weights_sexp) {
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

    for (int k = 0; k < K; k++) {
        const char* varname = CHAR(STRING_ELT(target_names, k));
        // Find column in data frame
        int col_idx = -1;
        for (int c = 0; c < LENGTH(col_names); c++) {
            if (strcmp(CHAR(STRING_ELT(col_names, c)), varname) == 0) {
                col_idx = c; break;
            }
        }
        if (col_idx < 0)
            Rf_error("Variable '%s' not found in data", varname);

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
    p.log_fn         = (p.verbose > 0) ? r_log_trampoline : nullptr;

    if (LENGTH(method_sexp) != 1)
        Rf_error("method must be a length-1 character string");
    const char* method_str = CHAR(STRING_ELT(method_sexp, 0));
    if      (strcmp(method_str, "ieppa")  == 0) p.algorithm = RK_ALG_IEPPA;
    else if (strcmp(method_str, "lbfgsb") == 0) p.algorithm = RK_ALG_LBFGSB;
    else                                          p.algorithm = RK_ALG_AUTO;

    rk_result_t result;
    rk_calibrate(n, K, weights.data(),
                 group_ids.data(),
                 cat_counts.data(),
                 (const double**)targets.data(),
                 &p, &result);

    // Build return list: list(weights=numeric[n], result=list(5 fields))
    SEXP wts = PROTECT(Rf_allocVector(REALSXP, n));
    memcpy(REAL(wts), weights.data(), (size_t)n * sizeof(double));

    SEXP res_list  = PROTECT(Rf_allocVector(VECSXP,  5));
    SEXP res_names = PROTECT(Rf_allocVector(STRSXP,  5));
    SET_STRING_ELT(res_names, 0, Rf_mkChar("status"));
    SET_STRING_ELT(res_names, 1, Rf_mkChar("iterations"));
    SET_STRING_ELT(res_names, 2, Rf_mkChar("max_error"));
    SET_STRING_ELT(res_names, 3, Rf_mkChar("algorithm_used"));
    SET_STRING_ELT(res_names, 4, Rf_mkChar("message"));
    SET_VECTOR_ELT(res_list, 0, Rf_ScalarInteger(result.status));
    SET_VECTOR_ELT(res_list, 1, Rf_ScalarInteger(result.iterations));
    SET_VECTOR_ELT(res_list, 2, Rf_ScalarReal(result.max_error));
    SET_VECTOR_ELT(res_list, 3, Rf_ScalarInteger((int)result.algorithm_used));
    SET_VECTOR_ELT(res_list, 4, Rf_mkString(result.message));
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

} // extern "C"
