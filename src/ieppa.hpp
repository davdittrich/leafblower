#pragma once
#include "types.hpp"
namespace lbw {
struct IEPPAResult { int status; int iterations; double max_error; };
IEPPAResult ieppa_solve(CalibState& state);
} // namespace lbw
