#pragma once
#include "types.hpp"
#include "cell_table.hpp"

namespace lbw {

struct SinkhornResult {
    CalibResult base;
    // ── Sinkhorn-specific extras ──
    int    M_cell              = 0;
    char   message[256]        = {};
    // CR-D11 (j7x8.11): cell-mode per-obs bound-violation diagnostic (count-only)
    // / unit-mode water-fill clamp count. Surfaced via has_n_bounds in r_bridge.
    int    n_bounds_violated   = 0;
    int    n_bounds_clamped    = 0;
    // ── End extras ──
    // convergence_solver_objective defaults to infinity in sinkhorn: set in sinkhorn.cpp init.
};

SinkhornResult sinkhorn_solve(CalibState& st);

} // namespace lbw
