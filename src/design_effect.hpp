#pragma once

#include "leafblower.h"

namespace lbw {

// Internal C++ entry — called by C-ABI rk_design_effect (in c_api.cpp).
// Returns RK_OK / RK_ERR_BADARG; populates *out on RK_OK.
int design_effect_compute(
    const double* weights,
    const double* outcome,        // nullable
    const int*    data_codes,     // nullable
    const int*    cat_counts,     // nullable
    int n, int K,
    rk_design_effect_result_t* out
) noexcept;

}  // namespace lbw
