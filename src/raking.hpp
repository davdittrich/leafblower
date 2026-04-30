#pragma once
#include "types.hpp"
namespace lbw {
struct RakingResult {
    CalibResult base;
    // ── Raking-specific extras ──
    double best_objective_seen          = 0.0;   // internal: weight KL at best_iter
    // ── End extras ──
};
RakingResult raking_solve(CalibState& state);
} // namespace lbw
