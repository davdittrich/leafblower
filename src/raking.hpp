#pragma once
#include "types.hpp"
#include "leafblower.h"
namespace lbw {
struct RakingResult {
    CalibResult base;
    // ── Raking-specific extras ──
    double best_objective_seen          = 0.0;   // internal: weight KL at best_iter
    // SRAA scheduler-demotion flag: TRUE iff st.accelerate && scheduler==GREEDY
    // (greedy demoted to round_robin under SRAA-m). FALSE otherwise.
    bool   sraa_demoted                 = false;
    // ── Hierarchical diagnostics (populated by raking_solve_hierarchical) ──
    int    hier_n_cells_total           = 0;
    int    hier_n_cells_skipped         = 0;
    int    hier_outer_iterations_used   = 0;
    double hier_outer_residual_final    = 0.0;
    int    hier_levels_used             = 0;
    // ── End extras ──
};
RakingResult raking_solve(CalibState& state);
// Two-stage hierarchical raking: dispatches on p->hierarchical_enabled.
// p == nullptr treated as hierarchical_enabled == 0 (single-stage passthrough).
RakingResult raking_solve_hierarchical(CalibState& st, const rk_params_t* p);
} // namespace lbw
