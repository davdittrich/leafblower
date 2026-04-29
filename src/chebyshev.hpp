#pragma once
#include "types.hpp"
#include "cell_table.hpp"
#include <vector>
#include <limits>

namespace lbw {

enum class LpVariant { CHEBYSHEV, GRAKE };

struct ChebyshevResult {
    int    status        = RK_ERR_NOCONV;
    int    iterations    = 0;
    double max_error     = 1.0;
    double mean_error    = 0.0;
    double kl            = 0.0;
    double chi2          = 0.0;
    double grake_norm    = 0.0;
    double l1_weight_change = 0.0;
    int    convergence_metric = 0;
    int    convergence_rule   = 0;
    double convergence_tol    = 0.0;
    int    convergence_iter   = -1;
    double best_objective_seen          = std::numeric_limits<double>::infinity();
    double convergence_solver_objective = std::numeric_limits<double>::infinity();
    int    convergence_minimized_metric = 0;
    double best_error   = std::numeric_limits<double>::infinity();
    int    best_iter    = 1;
    std::vector<double> best_weights;
    int    n_factorizations = 0;  // Mehrotra audit counter (populated in T5)
    int    M_cell       = 0;
    char   message[256] = {};
};

ChebyshevResult chebyshev_ipm(
    CalibState& st,
    LpVariant   variant,
    const std::vector<double>& w_warm_obs = {},  // obs-level warm weights; empty=cold start
    double      delta_warm = -1.0                // ieppa_max_err*1.5; -1=use default 1.0
);

inline ChebyshevResult chebyshev_solve(CalibState& st) {
    return chebyshev_ipm(st, LpVariant::CHEBYSHEV);
}
inline ChebyshevResult grake_solve(CalibState& st) {
    return chebyshev_ipm(st, LpVariant::GRAKE);
}

} // namespace lbw
