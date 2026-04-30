#pragma once
#include "leafblower.h"
#include "cell_table.hpp"
#include "types.hpp"

namespace lbw {

// kNCatsTotalMax is defined HERE (single definition).
// calib_linalg.hpp includes this header to share the constant.
constexpr int kNCatsTotalMax = 2048;

/**
 * Pre-entry validation for new cell-table solvers (sinkhorn, chebyshev, greg, greenkhorn).
 * Returns RK_OK, RK_ERR_BADARG, or RK_ERR_INFEAS.
 * On error, writes a human-readable message to result->message.
 *
 * Checks (in order):
 * 1. n_cats_total <= kNCatsTotalMax        → RK_ERR_BADARG  ("too many categories")
 * 2. L_c <= U_c for all cells              → RK_ERR_BADARG  with cell index
 * 3. X_init[c]==0 && L_c>0 for any cell   → RK_ERR_INFEAS  ("structural zero")
 * 4. sum(L_c) <= n <= sum(U_c)             → RK_ERR_INFEAS  ("total capacity")
 * 5. |sum(T_kj) - 1| <= 1e-6 for all k   → RK_ERR_BADARG  (normalize before calling)
 */
int calib_validate_preentry(const CellTable& ct,
                             const CalibState& st,
                             rk_result_t* result,
                             const double* X_init,
                             int n_cats_total);

} // namespace lbw
