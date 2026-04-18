#pragma once
#include "types.hpp"
#include "logit.hpp"
namespace lbw {
struct LBFGSResult { int status; int iterations; double max_error; };
LBFGSResult lbfgsb_solve(CalibState& state);
} // namespace lbw
