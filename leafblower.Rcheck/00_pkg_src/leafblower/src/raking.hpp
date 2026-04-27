#pragma once
#include "types.hpp"
namespace lbw {
struct RakingResult { int status; int iterations; double max_error; };
RakingResult raking_solve(CalibState& state);
} // namespace lbw
