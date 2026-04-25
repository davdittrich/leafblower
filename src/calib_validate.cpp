#include "calib_validate.hpp"
#include <cmath>
#include <cstring>
#include <cstdio>

namespace lbw {

int calib_validate_preentry(const CellTable& ct,
                             const CalibState& st,
                             rk_result_t* result,
                             const double* X_init,
                             int n_cats_total)
{
    auto fail = [&](int code, const char* msg) -> int {
        if (result) {
            result->status = code;
            std::strncpy(result->message, msg, 255);
            result->message[255] = '\0';
        }
        return code;
    };

    // 1. n_cats_total capacity check (before any allocation)
    if (n_cats_total > kNCatsTotalMax) {
        char msg[256];
        std::snprintf(msg, sizeof(msg),
            "n_cats_total=%d exceeds limit %d; use method='ieppa' or 'raking'",
            n_cats_total, kNCatsTotalMax);
        return fail(RK_ERR_BADARG, msg);
    }

    // Compute per-cell bounds
    const double lo = st.min_weight;
    const double hi = st.max_weight;
    double sum_L = 0.0, sum_U = 0.0;

    for (int c = 0; c < ct.M_cell; c++) {
        double L_c = lo * ct.n_per_cell[c];
        double U_c = hi * ct.n_per_cell[c];

        // 2. L_c <= U_c
        if (L_c > U_c + 1e-12) {
            char msg[256];
            std::snprintf(msg, sizeof(msg),
                "cell %d: L_c=%.4g > U_c=%.4g (min_weight > max_weight)", c, L_c, U_c);
            return fail(RK_ERR_BADARG, msg);
        }

        // 3. X_init[c]==0 && L_c>0 → structural infeasibility
        if (X_init != nullptr && X_init[c] <= 0.0 && L_c > 1e-12) {
            char msg[256];
            std::snprintf(msg, sizeof(msg),
                "cell %d: X_init=0 but L_c=%.4g > 0; multiplicative solver cannot "
                "move zero cell above lower bound", c, L_c);
            return fail(RK_ERR_INFEAS, msg);
        }

        sum_L += L_c;
        sum_U += U_c;
    }

    // 4. Total capacity vs target mass n
    const double n = static_cast<double>(st.n);
    if (sum_L > n + 1e-6) {
        char msg[256];
        std::snprintf(msg, sizeof(msg),
            "sum(L_c)=%.4g > n=%.4g: lower bounds exceed total target mass", sum_L, n);
        return fail(RK_ERR_INFEAS, msg);
    }
    if (sum_U < n - 1e-6) {
        char msg[256];
        std::snprintf(msg, sizeof(msg),
            "sum(U_c)=%.4g < n=%.4g: upper bounds cannot contain total target mass", sum_U, n);
        return fail(RK_ERR_INFEAS, msg);
    }

    // 5. Target sum validation
    for (int k = 0; k < st.K; k++) {
        double s = 0.0;
        for (int j = 0; j < st.cat_counts[k]; j++) s += st.targets[k][j];
        if (std::fabs(s - 1.0) > 1e-6) {
            char msg[256];
            std::snprintf(msg, sizeof(msg),
                "margin %d targets sum to %.8f (expected 1.0±1e-6); "
                "normalize targets before calling", k, s);
            return fail(RK_ERR_BADARG, msg);
        }
    }

    return RK_OK;
}

} // namespace lbw
