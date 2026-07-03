#pragma once
#include <limits>
#include "lbw_config.h"
#include "types.hpp"
#include "leafblower.h"

namespace lbw {

struct GreenkornResult {
    CalibResult base;
    // ── Greenkhorn-specific extras ──
    // +inf sentinel — first-iter finite objective always strictly improves it.
    int    n_bounds_violated            = 0;   // CR-D11b (j7x8.16): cell-mode bound-active count
    int    n_bounds_clamped             = 0;   // CR-D11b (j7x8.16): unit-mode clamp count
    char   message[256]        = {0};
    // ── End extras ──
};

GreenkornResult greenkhorn_solve(CalibState& st);

} // namespace lbw
