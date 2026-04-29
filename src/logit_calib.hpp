#pragma once
#include "types.hpp"
#include "cell_table.hpp"
#include <vector>
#include <limits>

namespace lbw {

// Result struct for logit-distance calibration (Deville-Sarndal 1992).
// Mirrors GregResult fields for downstream interchangeability + adds alm_*
// telemetry slots (zeroed; logit does not run an augmented-Lagrangian loop).
struct LogitCalibResult {
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
    char   message[256] = {0};

    // ALM telemetry (logit does not use ALM — kept zero for schema parity).
    double alm_capacity_mu_final  = 0.0;
    int    alm_n_growth_events    = 0;
    double alm_max_dual_norm      = 0.0;
    double alm_sum_drift          = 0.0;
};

// Logit-distance Newton-Raphson calibration (Deville-Sarndal 1992).
// Bounds enforced by construction via logit link: w[c] = L + (U-L)*sigma(z_c).
// Reuses calib_linalg normal-equations + LDLT infrastructure (same as greg).
LogitCalibResult logit_calibrate(CalibState& st);

} // namespace lbw
