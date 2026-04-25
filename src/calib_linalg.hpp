#pragma once
#include "calib_validate.hpp"   // provides kNCatsTotalMax — DO NOT redefine
#include "cell_table.hpp"
#include <cstddef>

namespace lbw {

/**
 * Compute N = A × diag(D) × Aᵀ where A is the (n_cats_total × M_cell)
 * marginal incidence matrix.
 * N is n_cats_total × n_cats_total, stored row-major.
 * Returns RK_OK or RK_ERR_BADARG (if n_cats_total > kNCatsTotalMax).
 * NOTE: implementation in Plan D (calib_linalg.cpp).
 */
int compute_normal_equations(const CellTable& ct,
                              const double* D,
                              double* N,
                              const int* cat_offset,
                              size_t n_cats_total);

/**
 * In-place LDLT factorization of a symmetric positive semidefinite n×n matrix.
 * eps_perturb: minimum diagonal value after factorization (Gill-Murray stability).
 * Returns RK_OK or RK_ERR_BADARG if n > kNCatsTotalMax.
 * NOTE: implementation in Plan D (calib_linalg.cpp).
 */
int ldlt_factor_inplace(double* A, size_t n, double eps_perturb);

/**
 * Solve A×x = b using a previously LDLT-factored matrix.
 * Overwrites b with the solution.
 * NOTE: implementation in Plan D (calib_linalg.cpp).
 */
void ldlt_solve(const double* L, const double* d_diag, double* b, size_t n);

} // namespace lbw
