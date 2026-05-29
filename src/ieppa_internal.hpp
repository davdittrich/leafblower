#ifndef LBW_IEPPA_INTERNAL_HPP
#define LBW_IEPPA_INTERNAL_HPP

// Internal (non-public) declarations shared between ieppa.cpp and the cold
// translation units split out of it. NOT part of the public ieppa.hpp API.
// Created in leafblower-uu8r (cold-only monolith split): the hot ieppa_solve
// stays in ieppa.cpp; only cold, once-per-solve helpers move out.

#include <utility>
#include <vector>

namespace lbw {

// Trajectory diagnostics (driven by LBW_TRAJECTORY_AT / LBW_TRAJECTORY_OUT env
// vars). Definitions live in ieppa_trajectory.cpp (uu8r.1).
std::vector<int> parse_trajectory_iters();
void write_trajectory_csv(const std::vector<std::pair<int, double>>& samples);

}  // namespace lbw

#endif  // LBW_IEPPA_INTERNAL_HPP
