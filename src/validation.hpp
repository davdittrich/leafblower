#pragma once
#include "leafblower.h"
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <vector>

namespace lbw {

// Validate inputs for rk_calibrate. Returns RK_OK or RK_ERR_BADARG.
// On RK_ERR_BADARG, sets result->status and result->message if result != nullptr.
// This is the single source of truth for input validation — called from both
// c_api.cpp (rk_calibrate C ABI path) and r_bridge.cpp (direct C++ path).
int validate_calibrate_inputs(int n, int K,
    const double* weights, const int32_t** group_ids,
    const int* cat_counts, const double** targets,
    const rk_params_t* p, rk_result_t* result, rk_algorithm_t alg);

} // namespace lbw
