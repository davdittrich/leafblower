#include "ieppa.hpp"
#include "leafblower.h"
namespace lbw {
IEPPAResult ieppa_solve(CalibState& st) {
    IEPPAResult r; r.status = 1; r.iterations = 0; r.max_error = 1.0;
    return r;
}
} // namespace lbw
