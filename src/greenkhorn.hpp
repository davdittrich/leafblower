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
    double best_objective_seen = std::numeric_limits<double>::infinity();
    char   message[256]        = {0};
    // ── End extras ──
};

GreenkornResult greenkhorn_solve(CalibState& st);

} // namespace lbw
