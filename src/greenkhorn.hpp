#pragma once
#include "lbw_config.h"
#include "types.hpp"
#include "leafblower.h"

namespace lbw {

struct GreenkornResult {
    CalibResult base;
    // ── Greenkhorn-specific extras ──
    double best_objective_seen = 0.0;   // internal: weight KL at best_iter
    char   message[256]        = {0};
    // ── Hierarchical diagnostics (populated by greenkhorn_solve_hierarchical) ──
    int    hier_n_cells_total           = 0;
    int    hier_n_cells_skipped         = 0;
    int    hier_outer_iterations_used   = 0;
    double hier_outer_residual_final    = 0.0;
    int    hier_levels_used             = 0;
};

GreenkornResult greenkhorn_solve(CalibState& st);

// Two-stage hierarchical Greenkhorn: dispatches on p->hierarchical_enabled.
// p == nullptr treated as hierarchical_enabled == 0 (single-stage passthrough).
// Queue isolation: greenkhorn_solve() rebuilds all state (X, W, S_flat, errRp,
// cells_per_cat) from CalibState on every call — no global or static queue state.
// Within-cell calls receive a sub-CalibState restricted to cell observations and
// fine margins, so the greedy priority queue is constructed on cell-local residuals.
GreenkornResult greenkhorn_solve_hierarchical(CalibState& st, const rk_params_t* p);

} // namespace lbw
