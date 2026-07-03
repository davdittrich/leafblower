#pragma once
#include "types.hpp"
#include "cell_table.hpp"

namespace lbw {

struct GregResult {
    CalibResult base;
    // ── GREG-specific extras ──
    // GREG defaults differ from CalibResult: set in greg.cpp init block.
    //   convergence_rule=0, convergence_tol=0.0, convergence_iter=1,
    //   convergence_metric=CHI2, convergence_minimized_metric=CHI2, best_iter=1
    int    M_cell                       = 0;
    char   message[256]                 = {};
    // CR-D11 (j7x8.11): cell-mode per-obs bound-violation diagnostic (count-only)
    // / unit-mode water-fill clamp count. Surfaced via has_n_bounds in r_bridge.
    int    n_bounds_violated            = 0;
    int    n_bounds_clamped             = 0;
    // ── End extras ──
};

// Newton QP for chi2 calibration (GREG — Deville-Sarnal 1992).
// Minimizes sum_c (X[c] - X_init[c])^2 / X_init[c] subject to margin constraints + capacity.
// One Newton iteration (exact when no bounds active). Active-set for bounds: up to kMaxNewtonIters=50 iterations.
GregResult greg_solve(CalibState& st);

} // namespace lbw
