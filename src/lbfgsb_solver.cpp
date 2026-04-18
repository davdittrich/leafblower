#include "lbfgsb_solver.hpp"
#include "leafblower.h"
namespace lbw {
LBFGSResult lbfgsb_solve(CalibState& st) {
    LBFGSResult r; r.status = 1; r.iterations = 0; r.max_error = 1.0;
    return r;
}
} // namespace lbw
