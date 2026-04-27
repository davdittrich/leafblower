#pragma once
#include "types.hpp"
#include <limits>
#include <vector>

namespace lbw {

struct IEPPAResult {
    int status;              // RK_OK / RK_ERR_NOCONV / RK_ERR_INFEAS / RK_ERR_BADARG
    int iterations;          // outer iterations completed
    double max_error;        // final errRp
    int M_cell;              // compression info
    int n_cap_active;        // cells with W[c] != 1 at convergence
    int n_xcur_writes_per_iter_linear;  // 0 outside linear path; counter for P1.1 RED test
    double min_alpha_seen;   // min alpha over all sweeps; 1.0 if damping never engaged
    double final_alpha;      // alpha at solver exit (after last sweep)
    int n_bounds_violated;  // cell-mode diagnostic: count of w_i outside bounds (no action)
    int n_bounds_clamped;   // unit-mode action: count of w_i clamped after water-fill exhausted
    // ── Overlay diagnostics ──
    int    homotopy_levels_used  = 0;   // 0 iff homotopy disabled
    double homotopy_final_factor = 1.0; // max_weight multiplier at final level
    int    greedy_sweeps_taken   = 0;   // total greedy margin sweeps across all levels and iterations
    double eta_final             = 0.0; // alm_mu multiplier at solver exit (0 = N/A)
    // ── End overlay diagnostics ──
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
    double sor_min_omega  = 1.0;
    int    sor_n_damped   = 0;
    double best_objective_seen          = 0.0;   // internal: weight KL at best_iter
    double convergence_solver_objective = 0.0;   // exposed: solver's mathematical objective
    int    convergence_minimized_metric = 0;     // CalibMetric: which metric was minimized
    double marginal_kl_at_iter            = 0.0;   // marginal KL at current outer iteration
    // ── End extended quality metrics ──
};

// Faithful iEPPA (paper-faithful algBCD at C=0). See
// docs/superpowers/specs/2026-04-23-ieppa-faithful-design.md.
IEPPAResult ieppa_solve(CalibState& state);

} // namespace lbw
