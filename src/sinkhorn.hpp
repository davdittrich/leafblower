#pragma once
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
};

SinkhornResult sinkhorn_solve(CalibState& st);

} // namespace lbw
