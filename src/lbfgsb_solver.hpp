#pragma once
#include "types.hpp"
#include "logit.hpp"

namespace lbw {

struct LBFGSResult {
    CalibResult base;
    // ── L-BFGS-B-specific extras ──
    double best_objective_seen = 0.0;   // internal: weight KL at best_iter
    // ── End extras ──
};

LBFGSResult lbfgsb_solve(CalibState& state);

} // namespace lbw
