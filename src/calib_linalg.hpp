#pragma once
#include "calib_validate.hpp"   // provides kNCatsTotalMax — DO NOT redefine
#include "cell_table.hpp"
#include <cstddef>

namespace lbw {

/**
 * Compute N = A × diag(D) × Aᵀ where A is the (n_cats_total × M_cell)
 * marginal incidence matrix.
 * N is n_cats_total × n_cats_total, stored row-major.
 * K: number of margins (ct has no K field; caller must supply it).
 * Returns RK_OK or RK_ERR_BADARG (if n_cats_total > kNCatsTotalMax).
 */
int compute_normal_equations(const CellTable& ct,
                              const double* D,
                              double* N,
                              const int* cat_offset,
                              int K,
                              size_t n_cats_total);

/**
 * Compute N = A_red × diag(D) × A_red^T where A_red excludes reference margins.
 * full_to_red[m] == -1 for reference margins (skipped); >= 0 for the reduced row/col index.
 * nct_red: number of non-reference constraints.
 * N output: nct_red × nct_red, row-major. No bounds check on nct_red (caller guards via nct).
 */
int compute_normal_equations_reduced(
    const CellTable& ct,
    const double* D,
    double* N,
    const int* cat_offset,
    int K,
    size_t nct_red,
    const int* full_to_red
) noexcept;

/**
 * In-place LDLT factorization of a symmetric positive semidefinite n×n matrix.
 * eps_perturb: minimum diagonal value after factorization (Gill-Murray stability).
 * Returns RK_OK or RK_ERR_BADARG if n > kNCatsTotalMax.
 */
int cholesky_factor_inplace(double* A, size_t n, double eps_perturb);

/**
 * Solve A×x = b using a previously LDLT-factored matrix (in-place combined storage).
 * A is the factored matrix from cholesky_factor_inplace: lower triangle holds L,
 * diagonal holds d.  Overwrites b with the solution.
 */
void cholesky_solve(const double* A, size_t n, double* b);

} // namespace lbw
