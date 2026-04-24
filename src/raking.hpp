#pragma once
#include "types.hpp"
#include <limits>
#include <vector>
namespace lbw {
struct RakingResult {
    int status;
    int iterations;
    double max_error;
    // ── Extended quality metrics (WU-A scaffold; populated in WU-B+) ──
    double mean_error     = 0.0;
    double kl             = 0.0;
    double chi2           = 0.0;
    double pct_change     = 0.0;
    double best_error     = std::numeric_limits<double>::infinity();
    int    best_iter      = 0;
    std::vector<double> best_weights;  // obs-level; length n; sum-normalized to n; empty if never checked
    // ── End extended quality metrics ──
};
RakingResult raking_solve(CalibState& state);
} // namespace lbw
