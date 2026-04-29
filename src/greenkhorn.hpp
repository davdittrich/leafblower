#pragma once
#include "lbw_config.h"
#include "types.hpp"
#include "leafblower.h"
#include <limits>
#include <vector>

namespace lbw {

struct GreenkornResult {
    int status;
    int iterations;
    double max_error;
    // ── Extended quality metrics (WU-A scaffold; populated in WU-B+) ──
    double mean_error          = 0.0;
    double kl                  = 0.0;
    double chi2                = 0.0;
    double l1_weight_change    = 0.0;  // WU-A: renamed from pct_change; computation in WU-B
    double grake_norm          = 0.0;  // WU-A stub; computation in WU-D
    int    convergence_metric  = 0;    // WU-A stub; CalibMetric at exit
    int    convergence_rule    = 1;    // WU-A stub; CalibRule at exit (IMPROVEMENT)
    double convergence_tol     = 0.001; // WU-A stub; threshold that fired
    int    convergence_iter    = -1;   // WU-A stub; iteration at convergence (-1=max_iter)
    double best_error          = std::numeric_limits<double>::infinity();
    int    best_iter           = 0;
    std::vector<double> best_weights;  // obs-level; length n; sum-normalized to n; empty if never checked
    double best_objective_seen          = 0.0;   // internal: weight KL at best_iter
    double convergence_solver_objective = 0.0;   // exposed: solver's mathematical objective
    int    convergence_minimized_metric = 0;     // CalibMetric: which metric was minimized
    // ── End extended quality metrics ──
    char message[256] = {0};
};

GreenkornResult greenkhorn_solve(CalibState& st);

} // namespace lbw
