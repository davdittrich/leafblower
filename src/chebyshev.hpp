#pragma once
#include "types.hpp"
#include "cell_table.hpp"
#include <vector>

namespace lbw {

struct ChebyshevResult {
    CalibResult base;
    // ── Chebyshev-specific extras ──
    // Chebyshev defaults differ from CalibResult: set in chebyshev.cpp init block.
    //   convergence_rule=0, convergence_tol=0.0, best_iter=1,
    //   convergence_solver_objective=inf
    double best_objective_seen          = std::numeric_limits<double>::infinity();  // reserved; never written in chebyshev.cpp
    int    n_factorizations             = 0;   // Mehrotra audit counter (populated in T5)
    int    M_cell                       = 0;
    char   message[256]                 = {};
    // ── End extras ──
};

ChebyshevResult chebyshev_ipm(
    CalibState& st,
    const std::vector<double>& w_warm_obs = {},  // obs-level warm weights; empty=cold start
    double      delta_warm = -1.0                // reserved, currently unused (see chebyshev.cpp (void)delta_warm)
);

inline ChebyshevResult chebyshev_solve(CalibState& st) {
    return chebyshev_ipm(st);
}

} // namespace lbw
