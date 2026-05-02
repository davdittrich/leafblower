// Epic-J WU-1 stub
#define STRICT_R_HEADERS 1
#define R_NO_REMAP 1

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include "cp_calib.hpp"
#include "ipm_calib.hpp"

extern "C" SEXP cp_solve_R(
    SEXP A_csr,
    SEXP b,
    SEXP d,
    SEXP lo,
    SEXP hi,
    SEXP max_iterations,
    SEXP capture_trace,
    SEXP seed
) {
  int n = Rf_length(d);

  // Build output list
  SEXP result;
  PROTECT(result = Rf_allocVector(VECSXP, 6));

  // weights: stub zeros
  SEXP weights;
  PROTECT(weights = Rf_allocVector(REALSXP, n));
  for (int i = 0; i < n; i++) {
    REAL(weights)[i] = 1.0;
  }
  SET_VECTOR_ELT(result, 0, weights);

  // status_code: 99 (stub)
  SEXP status_code;
  PROTECT(status_code = Rf_allocVector(INTSXP, 1));
  INTEGER(status_code)[0] = 99;
  SET_VECTOR_ELT(result, 1, status_code);

  // status_msg: "WU-1 stub"
  SEXP status_msg;
  PROTECT(status_msg = Rf_allocVector(STRSXP, 1));
  SET_STRING_ELT(status_msg, 0, Rf_mkChar("WU-1 stub"));
  SET_VECTOR_ELT(result, 2, status_msg);

  // iterations: 0
  SEXP iterations;
  PROTECT(iterations = Rf_allocVector(INTSXP, 1));
  INTEGER(iterations)[0] = 0;
  SET_VECTOR_ELT(result, 3, iterations);

  // wall_time_ms: 0.0
  SEXP wall_time;
  PROTECT(wall_time = Rf_allocVector(REALSXP, 1));
  REAL(wall_time)[0] = 0.0;
  SET_VECTOR_ELT(result, 4, wall_time);

  // trace: empty matrix with 0 rows, 6 columns
  SEXP trace;
  PROTECT(trace = Rf_allocMatrix(REALSXP, 0, 6));
  SET_VECTOR_ELT(result, 5, trace);

  // Set names
  SEXP names;
  PROTECT(names = Rf_allocVector(STRSXP, 6));
  SET_STRING_ELT(names, 0, Rf_mkChar("weights"));
  SET_STRING_ELT(names, 1, Rf_mkChar("status_code"));
  SET_STRING_ELT(names, 2, Rf_mkChar("status_msg"));
  SET_STRING_ELT(names, 3, Rf_mkChar("iterations"));
  SET_STRING_ELT(names, 4, Rf_mkChar("wall_time_ms"));
  SET_STRING_ELT(names, 5, Rf_mkChar("trace"));
  Rf_setAttrib(result, R_NamesSymbol, names);

  UNPROTECT(8);
  return result;
}

extern "C" SEXP ipm_solve_R(
    SEXP A_csr,
    SEXP b,
    SEXP d,
    SEXP lo,
    SEXP hi,
    SEXP max_iterations,
    SEXP capture_trace,
    SEXP seed
) {
  int n = Rf_length(d);

  // Build output list
  SEXP result;
  PROTECT(result = Rf_allocVector(VECSXP, 6));

  // weights: stub ones
  SEXP weights;
  PROTECT(weights = Rf_allocVector(REALSXP, n));
  for (int i = 0; i < n; i++) {
    REAL(weights)[i] = 1.0;
  }
  SET_VECTOR_ELT(result, 0, weights);

  // status_code: 99 (stub)
  SEXP status_code;
  PROTECT(status_code = Rf_allocVector(INTSXP, 1));
  INTEGER(status_code)[0] = 99;
  SET_VECTOR_ELT(result, 1, status_code);

  // status_msg: "WU-1 stub"
  SEXP status_msg;
  PROTECT(status_msg = Rf_allocVector(STRSXP, 1));
  SET_STRING_ELT(status_msg, 0, Rf_mkChar("WU-1 stub"));
  SET_VECTOR_ELT(result, 2, status_msg);

  // iterations: 0
  SEXP iterations;
  PROTECT(iterations = Rf_allocVector(INTSXP, 1));
  INTEGER(iterations)[0] = 0;
  SET_VECTOR_ELT(result, 3, iterations);

  // wall_time_ms: 0.0
  SEXP wall_time;
  PROTECT(wall_time = Rf_allocVector(REALSXP, 1));
  REAL(wall_time)[0] = 0.0;
  SET_VECTOR_ELT(result, 4, wall_time);

  // trace: empty matrix with 0 rows, 6 columns
  SEXP trace;
  PROTECT(trace = Rf_allocMatrix(REALSXP, 0, 6));
  SET_VECTOR_ELT(result, 5, trace);

  // Set names
  SEXP names;
  PROTECT(names = Rf_allocVector(STRSXP, 6));
  SET_STRING_ELT(names, 0, Rf_mkChar("weights"));
  SET_STRING_ELT(names, 1, Rf_mkChar("status_code"));
  SET_STRING_ELT(names, 2, Rf_mkChar("status_msg"));
  SET_STRING_ELT(names, 3, Rf_mkChar("iterations"));
  SET_STRING_ELT(names, 4, Rf_mkChar("wall_time_ms"));
  SET_STRING_ELT(names, 5, Rf_mkChar("trace"));
  Rf_setAttrib(result, R_NamesSymbol, names);

  UNPROTECT(8);
  return result;
}
