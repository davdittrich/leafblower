// Epic-J WU-1 stub + WU-2 body
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
    if (!Rf_isVectorList(A_csr) || Rf_length(A_csr) < 5) {
        Rf_error("A_csr must be a list of length 5 (p, j, n_row, n_col, x)");
    }
    
    SEXP A_p = VECTOR_ELT(A_csr, 0);
    SEXP A_j = VECTOR_ELT(A_csr, 1);
    SEXP A_n_row = VECTOR_ELT(A_csr, 2);
    SEXP A_n_col = VECTOR_ELT(A_csr, 3);
    SEXP A_x = VECTOR_ELT(A_csr, 4);

    /* Type validation BEFORE any pointer dereference (memory safety). */
    if (!Rf_isInteger(A_n_row) || Rf_length(A_n_row) < 1)
        Rf_error("A_csr$n_row must be INTEGER scalar");
    if (!Rf_isInteger(A_n_col) || Rf_length(A_n_col) < 1)
        Rf_error("A_csr$n_col must be INTEGER scalar");
    if (!Rf_isInteger(max_iterations) || Rf_length(max_iterations) < 1)
        Rf_error("max_iterations must be INTEGER scalar");
    if (!Rf_isLogical(capture_trace) || Rf_length(capture_trace) < 1)
        Rf_error("capture_trace must be LOGICAL scalar");
    if (!Rf_isInteger(A_p) || !Rf_isInteger(A_j))
        Rf_error("A_csr$p and A_csr$j must be INTEGER");
    if (!Rf_isReal(b) || !Rf_isReal(d) || !Rf_isReal(lo) || !Rf_isReal(hi) || !Rf_isReal(A_x))
        Rf_error("b, d, lo, hi, A_csr$x must be REAL");

    int n_row = INTEGER(A_n_row)[0];
    int n_col = INTEGER(A_n_col)[0];

    if (Rf_length(d) != n_row) Rf_error("length(d) != n_row");
    if (Rf_length(b) != n_col) Rf_error("length(b) != n_col");
    if (Rf_length(lo) != n_row) Rf_error("length(lo) != n_row");
    if (Rf_length(hi) != n_row) Rf_error("length(hi) != n_row");

    if (n_row > 100000000) Rf_error("research shim: n_row exceeded 1e8");
    if (n_col > 1000000) Rf_error("research shim: n_col exceeded 1e6");

    int max_iters = INTEGER(max_iterations)[0];
    if (max_iters > 100000) Rf_error("research shim: max_iterations exceeded 1e5");

    bool do_trace = LOGICAL(capture_trace)[0];

    const double* p_lo = REAL(lo);
    const double* p_hi = REAL(hi);
    const double delta = 1e-8;
    bool infeasible = false;
    for (int i = 0; i < n_row; i++) {
        if (p_lo[i] >= p_hi[i] - 2.0 * delta) {
            infeasible = true;
            break;
        }
    }

    CPResult cpres;
    if (infeasible) {
        cpres.weights.assign(n_row, 1.0);
        cpres.status_code = 3;
        cpres.status_msg = "infeasible bounds";
        cpres.iterations = 0;
        cpres.wall_time_ms = 0.0;
    } else {
        cpres = cp_calibrate(
            n_row, n_col,
            INTEGER(A_p), INTEGER(A_j), REAL(A_x),
            REAL(b), REAL(d), p_lo, p_hi,
            max_iters, do_trace
        );
    }

    int num_protects = 0;
    SEXP result;
    PROTECT(result = Rf_allocVector(VECSXP, 6));
    num_protects++;

    SEXP out_weights;
    PROTECT(out_weights = Rf_allocVector(REALSXP, n_row));
    num_protects++;
    for (int i = 0; i < n_row; i++) {
        REAL(out_weights)[i] = cpres.weights[i];
    }
    SET_VECTOR_ELT(result, 0, out_weights);

    SEXP out_status_code;
    PROTECT(out_status_code = Rf_allocVector(INTSXP, 1));
    num_protects++;
    INTEGER(out_status_code)[0] = cpres.status_code;
    SET_VECTOR_ELT(result, 1, out_status_code);

    SEXP out_status_msg;
    PROTECT(out_status_msg = Rf_allocVector(STRSXP, 1));
    num_protects++;
    SET_STRING_ELT(out_status_msg, 0, Rf_mkChar(cpres.status_msg.c_str()));
    SET_VECTOR_ELT(result, 2, out_status_msg);

    SEXP out_iters;
    PROTECT(out_iters = Rf_allocVector(INTSXP, 1));
    num_protects++;
    INTEGER(out_iters)[0] = cpres.iterations;
    SET_VECTOR_ELT(result, 3, out_iters);

    SEXP out_wall_time;
    PROTECT(out_wall_time = Rf_allocVector(REALSXP, 1));
    num_protects++;
    REAL(out_wall_time)[0] = cpres.wall_time_ms;
    SET_VECTOR_ELT(result, 4, out_wall_time);

    // Trace
    int num_cols = 6;
    int num_rows = cpres.trace_data.size() / num_cols;
    if (num_rows > 1000) {
        Rf_error("trace exceeded 1000 rows limit");
    }
    SEXP out_trace;
    PROTECT(out_trace = Rf_allocMatrix(REALSXP, num_rows, num_cols));
    num_protects++;
    
    double* trace_ptr = REAL(out_trace);
    for (int r = 0; r < num_rows; r++) {
        for (int c = 0; c < num_cols; c++) {
            trace_ptr[c * num_rows + r] = cpres.trace_data[r * num_cols + c];
        }
    }

    SEXP trace_colnames;
    PROTECT(trace_colnames = Rf_allocVector(STRSXP, num_cols));
    num_protects++;
    SET_STRING_ELT(trace_colnames, 0, Rf_mkChar("iter"));
    SET_STRING_ELT(trace_colnames, 1, Rf_mkChar("time_ms"));
    SET_STRING_ELT(trace_colnames, 2, Rf_mkChar("max_err_last"));
    SET_STRING_ELT(trace_colnames, 3, Rf_mkChar("max_err_ergodic"));
    SET_STRING_ELT(trace_colnames, 4, Rf_mkChar("primal_resid"));
    SET_STRING_ELT(trace_colnames, 5, Rf_mkChar("primal_stationarity_proxy"));

    SEXP trace_dimnames;
    PROTECT(trace_dimnames = Rf_allocVector(VECSXP, 2));
    num_protects++;
    SET_VECTOR_ELT(trace_dimnames, 0, R_NilValue);
    SET_VECTOR_ELT(trace_dimnames, 1, trace_colnames);
    Rf_setAttrib(out_trace, R_DimNamesSymbol, trace_dimnames);
    
    SET_VECTOR_ELT(result, 5, out_trace);

    SEXP names;
    PROTECT(names = Rf_allocVector(STRSXP, 6));
    num_protects++;
    SET_STRING_ELT(names, 0, Rf_mkChar("weights"));
    SET_STRING_ELT(names, 1, Rf_mkChar("status_code"));
    SET_STRING_ELT(names, 2, Rf_mkChar("status_msg"));
    SET_STRING_ELT(names, 3, Rf_mkChar("iterations"));
    SET_STRING_ELT(names, 4, Rf_mkChar("wall_time_ms"));
    SET_STRING_ELT(names, 5, Rf_mkChar("trace"));
    Rf_setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(num_protects);
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

  SEXP result;
  PROTECT(result = Rf_allocVector(VECSXP, 6));

  SEXP weights;
  PROTECT(weights = Rf_allocVector(REALSXP, n));
  for (int i = 0; i < n; i++) {
    REAL(weights)[i] = 1.0;
  }
  SET_VECTOR_ELT(result, 0, weights);

  SEXP status_code;
  PROTECT(status_code = Rf_allocVector(INTSXP, 1));
  INTEGER(status_code)[0] = 99;
  SET_VECTOR_ELT(result, 1, status_code);

  SEXP status_msg;
  PROTECT(status_msg = Rf_allocVector(STRSXP, 1));
  SET_STRING_ELT(status_msg, 0, Rf_mkChar("WU-1 stub"));
  SET_VECTOR_ELT(result, 2, status_msg);

  SEXP iterations;
  PROTECT(iterations = Rf_allocVector(INTSXP, 1));
  INTEGER(iterations)[0] = 0;
  SET_VECTOR_ELT(result, 3, iterations);

  SEXP wall_time;
  PROTECT(wall_time = Rf_allocVector(REALSXP, 1));
  REAL(wall_time)[0] = 0.0;
  SET_VECTOR_ELT(result, 4, wall_time);

  SEXP trace;
  PROTECT(trace = Rf_allocMatrix(REALSXP, 0, 6));
  SET_VECTOR_ELT(result, 5, trace);

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
