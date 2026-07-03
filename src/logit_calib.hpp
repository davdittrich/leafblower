#pragma once
#include "types.hpp"
#include "cell_table.hpp"

namespace lbw {

// Result struct for logit-distance calibration (Deville-Sarndal 1992).
// Mirrors GregResult fields for downstream interchangeability.
struct LogitCalibResult {
    CalibResult base;
    // ── Logit-specific extras ──
    // Logit defaults differ from CalibResult: set in logit_calib.cpp init block.
    //   convergence_rule=0, convergence_tol=0.0, convergence_iter=1,
    //   convergence_metric=CHI2, convergence_minimized_metric=CHI2, best_iter=1
    int    M_cell                       = 0;
    int    n_bounds_violated            = 0;   // CR-D11b (j7x8.16): cell-mode bound-active count
    int    n_bounds_clamped             = 0;   // CR-D11b (j7x8.16): unit-mode clamp count
    char   message[256]                 = {0};
    // ── End extras ──
};

// Logit-distance Newton-Raphson calibration (Deville-Sarndal 1992).
// Bounds enforced by construction via logit link: w[c] = L + (U-L)*sigma(z_c).
// Reuses calib_linalg normal-equations + LDLT infrastructure (same as greg).
LogitCalibResult logit_calibrate(CalibState& st);

} // namespace lbw
