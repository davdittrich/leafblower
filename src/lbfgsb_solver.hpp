#pragma once
#include "types.hpp"
#include "logit.hpp"
#include <vector>

namespace lbw {

struct LBFGSResult {
    int    status;       // RK_OK or RK_ERR_NOCONV
    int    iterations;
    double max_error;
};

LBFGSResult lbfgsb_solve(CalibState& state);

} // namespace lbw
