#include "calib_linalg.hpp"
#include "leafblower.h"
#include <algorithm>
#include <cmath>
#include <cstring>

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

int ldlt_factor_inplace(double* A, size_t n, double eps_perturb)
{
    if (n > static_cast<size_t>(kNCatsTotalMax)) return RK_ERR_BADARG;
    for (size_t j = 0; j < n; j++) {
        // Compute diagonal d_j = A[j,j] - sum_{k<j} d_k * L[j,k]^2
        double d_j = A[j * n + j];
        for (size_t k = 0; k < j; k++) {
            double l_jk = A[j * n + k];
            d_j -= A[k * n + k] * l_jk * l_jk;
        }
        // Gill-Murray modified LDLT: clamp diagonal away from zero
        d_j = std::max(d_j, eps_perturb);
        A[j * n + j] = d_j;
        // Fill column j of L below the diagonal
        for (size_t i = j + 1; i < n; i++) {
            double s = A[i * n + j];
            for (size_t k = 0; k < j; k++)
                s -= A[k * n + k] * A[i * n + k] * A[j * n + k];
            A[i * n + j] = s / d_j;
        }
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
    for (size_t i = 0; i < n; i++)
        b[i] /= A[i * n + i];

    // Backward substitution: L^T x = z
    // Use signed int to avoid size_t underflow when decrementing past 0
    for (int i = static_cast<int>(n) - 1; i >= 0; i--) {
        for (int j = i + 1; j < static_cast<int>(n); j++)
            b[i] -= A[static_cast<size_t>(j) * n +
                      static_cast<size_t>(i)] * b[j];
    }
}

} // namespace lbw
