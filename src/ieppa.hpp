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
    int n_xcur_writes_per_iter_linear;  // 0 outside linear path; counter for P1.1 RED test
    double min_alpha_seen;   // min alpha over all sweeps; 1.0 if damping never engaged
    double final_alpha;      // alpha at solver exit (after last sweep)
    int n_bounds_violated;  // cell-mode diagnostic: count of w_i outside bounds (no action)
    int n_bounds_clamped;   // unit-mode action: count of w_i clamped after water-fill exhausted
};

// Faithful iEPPA (paper-faithful algBCD at C=0). See
// docs/superpowers/specs/2026-04-23-ieppa-faithful-design.md.
IEPPAResult ieppa_solve(CalibState& state);

} // namespace lbw
