#pragma once
#include "types.hpp"

namespace lbw {

struct NewtonCalibResult {
    CalibResult base;
    int    n_lambda   = 0;
    double dual_gap   = 0.0;
    double step_norm  = 0.0;
    double line_alpha = 1.0;
    double lm_mu_final = 0.0;
    int    n_projected_dims = 0;  // Epic-Dβ: count of TSVD-truncated null-space dims per last solve
    char   message[256] = {0};
};

NewtonCalibResult newton_calibrate(CalibState& st);

} // namespace lbw
