#pragma once
#include "types.hpp"
#include <vector>

namespace lbw {

struct IEPPAResult {
    int status;              // RK_OK / RK_ERR_NOCONV / RK_ERR_INFEAS / RK_ERR_BADARG
    int iterations;          // outer iterations completed
    double max_error;        // final errRp
    int M_cell;              // compression info
    int n_cap_active;        // cells with W[c] != 1 at convergence
};

// Faithful iEPPA (paper-faithful algBCD at C=0). See
// docs/superpowers/specs/2026-04-23-ieppa-faithful-design.md.
IEPPAResult ieppa_solve(CalibState& state);

} // namespace lbw
