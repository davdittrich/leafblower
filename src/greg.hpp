#pragma once
#include "types.hpp"
#include "cell_table.hpp"
#include <vector>
#include <limits>

namespace lbw {

struct GregResult {
    int    status        = RK_ERR_NOCONV;
    int    iterations    = 0;
    double max_error     = 1.0;
    double mean_error    = 0.0;
    double kl            = 0.0;
    double chi2          = 0.0;
    double grake_norm    = 0.0;
    double l1_weight_change = 0.0;
    int    convergence_metric = static_cast<int>(CalibMetric::CHI2);
    int    convergence_rule   = 0;
    double convergence_tol    = 0.0;
    int    convergence_iter   = 1;
    double best_objective_seen          = std::numeric_limits<double>::infinity();
    double convergence_solver_objective = std::numeric_limits<double>::infinity();
    int    convergence_minimized_metric = static_cast<int>(CalibMetric::CHI2);
    double best_error   = std::numeric_limits<double>::infinity();
    int    best_iter    = 1;
    std::vector<double> best_weights;
    int    M_cell       = 0;
    char   message[256] = {};
};

// Newton QP for chi2 calibration (GREG — Deville-Sarnal 1992).
// Minimizes sum_c (X[c] - X_init[c])^2 / X_init[c] subject to margin constraints + capacity.
// One Newton iteration (exact when no bounds active). Active-set for bounds: <=10 iterations.
GregResult greg_solve(CalibState& st);

} // namespace lbw
