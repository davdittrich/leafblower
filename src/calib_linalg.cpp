#include "calib_linalg.hpp"
#include "leafblower.h"
#include <R_ext/Error.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

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
            for (int k2 = 0; k2 < K; k2++) {
                int j2 = ct.g_per_cell[k2][c];
                if (j2 < 0) continue;
                size_t col = static_cast<size_t>(cat_offset[k2]) +
                             static_cast<size_t>(j2);
                N[row * n_cats_total + col] += D[c];
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

int ldlt_factor_inplace(double* A, size_t n, double eps_perturb)
{
    if (n > static_cast<size_t>(kNCatsTotalMax)) return RK_ERR_BADARG;

    // Gill-Murray-Wright (1981) bound parameter:
    //   beta^2 = max(gamma, xi / sqrt(n^2 - 1), DBL_EPSILON)
    // where gamma = max |A_ii|, xi = max |A_ij| for i != j.
    // The off-diagonal scan reads the strict lower triangle (A is symmetric).
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

    for (size_t j = 0; j < n; j++) {
        // Compute diagonal d_j = A[j,j] - sum_{k<j} d_k * L[j,k]^2
        double d_j = A[j * n + j];
        for (size_t k = 0; k < j; k++) {
            double l_jk = A[j * n + k];
            d_j -= A[k * n + k] * l_jk * l_jk;
        }

        // Compute the unscaled column entries s_ij = A[i,j] - sum_{k<j} d_k L[i,k] L[j,k]
        // and track theta_j = max_{i>j} |s_ij| in the same pass (no extra read).
        // We stash s_ij into A[i*n+j] temporarily; the final L[i,j] = s_ij / d_j is
        // assigned below once d_j has been bounded.
        double theta_j = 0.0;
        for (size_t i = j + 1; i < n; i++) {
            double s = A[i * n + j];
            for (size_t k = 0; k < j; k++)
                s -= A[k * n + k] * A[i * n + k] * A[j * n + k];
            A[i * n + j] = s;
            theta_j = std::max(theta_j, std::abs(s));
        }

        // Gill-Murray-Wright modified LDLT: bound d_j by max(|d_j|, theta_j^2 / beta^2, eps_perturb).
        // The theta_j^2 / beta^2 term prevents L[i,j] = s_ij / d_j from blowing up when
        // a small d_j coincides with O(1) off-diagonal entries.
        d_j = std::max({std::abs(d_j), (theta_j * theta_j) / beta2, eps_perturb});
        A[j * n + j] = d_j;

        // Finalize column j of L: L[i,j] = s_ij / d_j (s_ij currently stored at A[i,j]).
        for (size_t i = j + 1; i < n; i++)
            A[i * n + j] /= d_j;
    }
    return RK_OK;
}

void ldlt_solve(const double* A, size_t n, double* b)
{
    // Forward substitution: L y = b  (L is unit lower triangular)
    for (size_t i = 1; i < n; i++)
        for (size_t j = 0; j < i; j++)
            b[i] -= A[i * n + j] * b[j];

    // Diagonal scaling: D z = y
    for (size_t i = 0; i < n; i++) {
        if (std::fabs(A[i * n + i]) < 1e-300)
            Rf_error("leafblower: LDLT singular at index %d", static_cast<int>(i));
        b[i] /= A[i * n + i];
    }

    // Backward substitution: L^T x = z
    // Use signed int to avoid size_t underflow when decrementing past 0
    for (int i = static_cast<int>(n) - 1; i >= 0; i--) {
        for (int j = i + 1; j < static_cast<int>(n); j++)
            b[i] -= A[static_cast<size_t>(j) * n +
                      static_cast<size_t>(i)] * b[j];
    }
}

} // namespace lbw
