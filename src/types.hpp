#pragma once
#include <cstdint>
#include <R_ext/Print.h>

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
    double epsilon;
    int lbfgs_m;
    int verbose;
    void (*log_fn)(const char* msg, void* ctx);
    void* log_ctx;

    // Derived
    int total_cats;                    // sum of cat_counts[k]

    void log(const char* msg) const {
        if (verbose <= 0) return;
        if (log_fn) {
            log_fn(msg, log_ctx);
        } else {
            REprintf("%s\n", msg);
        }
    }
};

} // namespace lbw
