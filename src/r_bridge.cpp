#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "logit.hpp"

extern "C" {

// Forward declarations for all .Call entries defined in this file.
// (C_rk_calibrate is added in Task 7; its declaration is added there too.)
SEXP C_logit_F_at_zero(SEXP, SEXP);
SEXP C_logit_range_check(SEXP, SEXP, SEXP);
SEXP C_logit_Hprime_check(SEXP, SEXP, SEXP);

// R_init_leafblower: called automatically by R when the shared library is loaded.
// Registers all .Call entry points so they can be found by name.
// NOTE: updated in Task 7 to include C_rk_calibrate.
void R_init_leafblower(DllInfo* dll) {
    static const R_CallMethodDef call_methods[] = {
        {"C_logit_F_at_zero",    (DL_FUNC)&C_logit_F_at_zero,    2},
        {"C_logit_range_check",  (DL_FUNC)&C_logit_range_check,  3},
        {"C_logit_Hprime_check", (DL_FUNC)&C_logit_Hprime_check, 3},
        {NULL, NULL, 0}
    };
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);  // prevent lookup of unregistered symbols
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

} // extern "C"
