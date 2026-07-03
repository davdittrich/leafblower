#pragma once
#include <limits>
#include "types.hpp"
namespace lbw {
struct RakingResult {
    CalibResult base;
    // ── Raking-specific extras ──
    // +inf sentinel — every real KL objective is finite, so first-iter values
    // strictly improve the sentinel. 0.0 was a footgun for non-negative metrics.
    double best_objective_seen          = std::numeric_limits<double>::infinity();
    // SRAA scheduler-demotion flag: TRUE iff st.accelerate && scheduler==GREEDY
    // (greedy demoted to round_robin under SRAA-m). FALSE otherwise.
    bool   sraa_demoted                 = false;
    // CR-D11 (j7x8.11): cell-mode per-obs bound-violation diagnostic (count-only;
    // no clamp) / unit-mode water-fill clamp count. Surfaced via has_n_bounds in
    // r_bridge so a no-clamp cell-mode result reports its true violation count.
    int    n_bounds_violated            = 0;
    int    n_bounds_clamped             = 0;
    // ── End extras ──
};
RakingResult raking_solve(CalibState& state);
} // namespace lbw
