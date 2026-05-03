#include "calib_linalg.hpp"
#include "leafblower.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>

#ifndef LBW_NO_R
#  include <R_ext/Error.h>
#  include <R_ext/Lapack.h>
#  include <R_ext/RS.h>   // F77_CALL
#  ifndef FCONE
#    define FCONE
#  endif
#else
// System LAPACK declarations for Python build
extern "C" {
    void dpotrf_(char* uplo, int* n, double* A, int* lda, int* info);
    void dpotrs_(char* uplo, int* n, int* nrhs, double* A, int* lda,
                 double* B, int* ldb, int* info);
}
#  define F77_CALL(x) x ## _
#  ifndef FCONE
#    define FCONE
#  endif
#endif

namespace lbw {

int compute_normal_equations(const CellTable& ct,
                              const double* D,
                              double* N,
                              const int* cat_offset,
                              int K,
                              size_t n_cats_total)
{
    if (n_cats_total > static_cast<size_t>(kNCatsTotalMax)) return RK_ERR_BADARG;
    std::fill(N, N + n_cats_total * n_cats_total, 0.0);
    for (int c = 0; c < ct.M_cell; c++) {
        if (D[c] <= 0.0) continue;
        for (int k1 = 0; k1 < K; k1++) {
            int j1 = ct.g_per_cell[k1][c];
            if (j1 < 0) continue;
            size_t row = static_cast<size_t>(cat_offset[k1]) +
                         static_cast<size_t>(j1);
            // G2: only fill lower triangle (row >= col); dpotrf(uplo='L') reads lower only.
            for (int k2 = 0; k2 < K; k2++) {
                int j2 = ct.g_per_cell[k2][c];
                if (j2 < 0) continue;
                size_t col = static_cast<size_t>(cat_offset[k2]) +
                             static_cast<size_t>(j2);
                if (row >= col) N[row * n_cats_total + col] += D[c];
            }
        }
    }
    return RK_OK;
}

int compute_normal_equations_reduced(const CellTable& ct,
                                      const double* D,
                                      double* N,
                                      const int* cat_offset,
                                      int K,
                                      size_t nct_red,
                                      const int* full_to_red) noexcept
{
    std::fill(N, N + nct_red * nct_red, 0.0);
    for (int c = 0; c < ct.M_cell; c++) {
        if (D[c] <= 0.0) continue;
        for (int k1 = 0; k1 < K; k1++) {
            int j1 = ct.g_per_cell[k1][c];
            if (j1 < 0) continue;
            int m1 = cat_offset[k1] + j1;
            int r1 = full_to_red[m1];
            if (r1 < 0) continue;  // reference margin — skip row
            for (int k2 = 0; k2 < K; k2++) {
                int j2 = ct.g_per_cell[k2][c];
                if (j2 < 0) continue;
                int m2 = cat_offset[k2] + j2;
                int r2 = full_to_red[m2];
                if (r2 < 0) continue;  // reference margin — skip col
                N[static_cast<size_t>(r1) * nct_red + static_cast<size_t>(r2)] += D[c];
            }
        }
    }
    return RK_OK;
}

// ─────────────────────────────────────────────────────────────────────────────
// G1: LAPACK Cholesky factor + solve.
//
// Replaces the prior O(n^3) scalar LDLT with BLAS-3 dpotrf/dpotrs.
//
// Storage convention:
//   The matrix A is symmetric and is written full by compute_normal_equations*.
//   The buffer is row-major in C++; LAPACK reads it as column-major. Because A
//   is symmetric, A_rowmajor == A_colmajor, so the factorization is correct
//   regardless. We use uplo='L' for both factor and solve so dpotrs reads the
//   same triangle dpotrf wrote. The factor (Cholesky L in column-major) lives
//   in the row-major upper triangle of the buffer — irrelevant to callers, who
//   only ever pass the buffer back into ldlt_solve.
//
// GMW dependency:
//   dpotrf requires positive-definite input. C1 introduced the Gill-Murray-
//   Wright column-norm bound; G1 preserves PD by perturbing the diagonal up-
//   front using a single bump derived from the same gamma/xi scale, plus the
//   caller-supplied eps_perturb floor. We do not retry dpotrf on info>0:
//   dpotrf overwrites the lower triangle in place, so a retry would need a
//   full copy. Instead we pre-bump aggressively (max(eps_perturb,
//   sqrt(eps_machine) * (gamma+beta2))) which is safe relative to numerical
//   noise yet small enough to preserve solution accuracy. Truly indefinite
//   inputs return RK_ERR_BADARG and the caller falls back to its error path.
// ─────────────────────────────────────────────────────────────────────────────
int ldlt_factor_inplace(double* A, size_t n, double eps_perturb)
{
    if (n > static_cast<size_t>(kNCatsTotalMax)) return RK_ERR_BADARG;
    if (n == 0) return RK_OK;

    // Gill-Murray-Wright bound parameter (same definition as the C1 fix):
    //   beta^2 = max(gamma, xi / sqrt(n^2 - 1), DBL_EPSILON)
    // gamma = max |A_ii|, xi = max |A_ij| for i != j.
    double gamma = 0.0, xi = 0.0;
    for (size_t i = 0; i < n; ++i) {
        gamma = std::max(gamma, std::abs(A[i * n + i]));
        for (size_t j = 0; j < i; ++j)
            xi = std::max(xi, std::abs(A[i * n + j]));
    }
    const double n2m1 = static_cast<double>(n) * static_cast<double>(n) - 1.0;
    const double beta2 = std::max({
        gamma,
        xi / std::sqrt(std::max(1.0, n2m1)),
        std::numeric_limits<double>::epsilon()
    });

    // Up-front diagonal bump. We use sqrt(machine_eps) ~ 1.49e-8 as the relative
    // floor, scaled by (gamma + beta2) which captures the matrix's dominant
    // magnitude. This is large enough to absorb roundoff on PSD inputs that
    // sit at the edge of indefiniteness, but small enough not to perturb the
    // solution beyond solver tolerances (kEpsLdlt = 1e-10 is a typical caller
    // floor; we honor that as a lower bound).
    const double sqrt_eps = std::sqrt(std::numeric_limits<double>::epsilon());
    const double bump = std::max(eps_perturb, sqrt_eps * (gamma + beta2));
    for (size_t i = 0; i < n; ++i) A[i * n + i] += bump;

    char uplo = 'L';
    int n_int = static_cast<int>(n);
    int info = 0;

    F77_CALL(dpotrf)(&uplo, &n_int, A, &n_int, &info FCONE);

    if (info != 0) return RK_ERR_BADARG;
    return RK_OK;
}

void ldlt_solve(const double* A, size_t n, double* b)
{
    if (n == 0) return;
    char uplo = 'L';
    int n_int = static_cast<int>(n);
    int nrhs = 1, info = 0;

    // dpotrs requires a non-const pointer for A; the call does not modify A.
    double* A_mut = const_cast<double*>(A);
    F77_CALL(dpotrs)(&uplo, &n_int, &nrhs, A_mut, &n_int, b, &n_int, &info FCONE);
    if (info != 0) {
#ifndef LBW_NO_R
        Rf_error("leafblower: dpotrs failed with info=%d", info);
#else
        throw std::runtime_error("leafblower: dpotrs failed");
#endif
    }
}

} // namespace lbw
