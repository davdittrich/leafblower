#pragma once
#include "types.hpp"
#include "logit.hpp"
#include <limits>
#include <vector>

namespace lbw {

struct LBFGSResult {
    int    status;       // RK_OK or RK_ERR_NOCONV
    int    iterations;
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
    double convergence_objective          = 0.0;   // value of minimized metric at convergence
    int    convergence_minimized_metric   = 0;     // CalibMetric: which metric was minimized
    // ── End extended quality metrics ──
};

LBFGSResult lbfgsb_solve(CalibState& state);

} // namespace lbw
