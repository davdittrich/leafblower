#pragma once
#include <cstdint>
#ifndef LBW_NO_R
#  include <R_ext/Print.h>
#else
#  include <cstdio>
#endif

namespace lbw {

struct CalibState {
    int n;
    int K;
    double* weights;                   // caller-owned, n elements
    const int32_t** group_ids;        // caller-owned, K pointers to n int32_t elements
    const int* cat_counts;            // caller-owned, K elements
    const double** targets;           // caller-owned, K pointers to cat_counts[k] doubles
    double min_weight;
    double max_weight;
    double tol_abs;
    int inner_max_iter;
    int outer_max_iter;
    int lbfgs_m;
    int verbose;
    void (*log_fn)(const char* msg, void* ctx);
    void* log_ctx;

    void log(const char* msg) const {
        if (verbose <= 0) return;
        if (log_fn) {
            log_fn(msg, log_ctx);
        } else {
#ifndef LBW_NO_R
            REprintf("%s\n", msg);
#else
            fprintf(stderr, "%s\n", msg);
#endif
        }
    }
};

} // namespace lbw
