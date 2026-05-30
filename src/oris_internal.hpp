#ifndef LBW_ORIS_INTERNAL_HPP
#define LBW_ORIS_INTERNAL_HPP

// Internal (non-public) declarations shared between oris.cpp and the cold
// translation units split out of it. NOT part of the public oris.hpp API.
// Created in leafblower-uu8r (cold-only monolith split): the hot oris_solve
// stays in oris.cpp; only cold, once-per-solve helpers move out.

#include <set>
#include <utility>
#include <vector>

#include "oris.hpp"            // ORISResult
#include "calib_dispatch.hpp"  // CalibState, BestIterTracker
#include "cell_table.hpp"      // CellTable

namespace lbw {

// Trajectory diagnostics (driven by LBW_TRAJECTORY_AT / LBW_TRAJECTORY_OUT env
// vars). Definitions live in oris_trajectory.cpp (uu8r.1).
std::vector<int> parse_trajectory_iters();
void write_trajectory_csv(const std::vector<std::pair<int, double>>& samples);

// Post-loop finalization — obs expansion, bounds enforcement, best-iterate
// fallback. Called once per oris_solve() after the homotopy loop exits.
// Definition in oris_finalize.cpp (uu8r.2).
void oris_finalize(
    CalibState&                               st,
    ORISResult&                               res,
    const CellTable&                          ct,
    std::vector<double>&                      X,
    const std::vector<double>&                X_init,
    const std::vector<double>&                L_cell,
    const std::vector<double>&                U_cell,
    bool                                      alm_active,
    double                                    capacity_mu_adaptive,
    const std::vector<double>&                lambda_cell,
    const BestIterTracker&                    best,
    bool                                      absolute_tol_fired,
    const std::set<std::pair<int,int>>&       structural_infeas_pairs,
    double                                    sor_min_omega,
    int                                       sor_n_damped,
    const std::vector<std::pair<int,double>>& probe_samples);

}  // namespace lbw

#endif  // LBW_ORIS_INTERNAL_HPP
