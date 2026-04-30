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
    double best_objective_seen          = std::numeric_limits<double>::infinity();
    int    M_cell                       = 0;
    char   message[256]                 = {};
    // ── End extras ──
};

// Newton QP for chi2 calibration (GREG — Deville-Sarnal 1992).
// Minimizes sum_c (X[c] - X_init[c])^2 / X_init[c] subject to margin constraints + capacity.
// One Newton iteration (exact when no bounds active). Active-set for bounds: <=10 iterations.
GregResult greg_solve(CalibState& st);

} // namespace lbw
