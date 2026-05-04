#pragma once
#include "types.hpp"
namespace lbw {
struct RakingResult {
    CalibResult base;
    // ── Raking-specific extras ──
    double best_objective_seen          = 0.0;   // internal: weight KL at best_iter
    // SRAA scheduler-demotion flag: TRUE iff st.accelerate && scheduler==GREEDY
    // (greedy demoted to round_robin under SRAA-m). FALSE otherwise.
    bool   sraa_demoted                 = false;
    // ── End extras ──
};
RakingResult raking_solve(CalibState& state);
} // namespace lbw
