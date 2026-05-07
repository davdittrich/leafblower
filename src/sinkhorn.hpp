#pragma once
#include "leafblower.h"
#include "types.hpp"
#include "cell_table.hpp"

namespace lbw {

struct SinkhornResult {
    CalibResult base;
    // ── Sinkhorn-specific extras ──
    double best_objective_seen = std::numeric_limits<double>::infinity();
    int    M_cell              = 0;
    char   message[256]        = {};
    // ── End extras ──
    // convergence_solver_objective defaults to infinity in sinkhorn: set in sinkhorn.cpp init.

    // ── Hierarchical diagnostics (populated by sinkhorn_solve_hierarchical) ──
    int    hier_n_cells_total           = 0;
    int    hier_n_cells_skipped         = 0;
    int    hier_outer_iterations_used   = 0;
    double hier_outer_residual_final    = 0.0;
    int    hier_levels_used             = 0;
};

SinkhornResult sinkhorn_solve(CalibState& st);

// Two-stage hierarchical Sinkhorn: dispatches on p->hierarchical_enabled.
// p == nullptr treated as hierarchical_enabled == 0 (single-stage passthrough).
SinkhornResult sinkhorn_solve_hierarchical(CalibState& st, const rk_params_t* p);

} // namespace lbw
