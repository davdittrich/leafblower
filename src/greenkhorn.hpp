#pragma once
#include "lbw_config.h"
#include "types.hpp"
#include "leafblower.h"

namespace lbw {

struct GreenkornResult {
    CalibResult base;
    // ── Greenkhorn-specific extras ──
    double best_objective_seen = 0.0;   // internal: weight KL at best_iter
    char   message[256]        = {0};
    // ── End extras ──
};

GreenkornResult greenkhorn_solve(CalibState& st);

} // namespace lbw
