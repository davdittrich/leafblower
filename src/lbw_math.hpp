#pragma once
#include "lbw_config.h"
#include <cmath>

#if LBW_HAS_GLIBC_MVEC
#  include <immintrin.h>
extern "C" __m256d _ZGVdN4v_exp(__m256d) __attribute__((target("avx2")));
#endif

namespace lbw {

// Compute out[i] = exp(scale * u[i]) for i in [0, n).
// With glibc libmvec: 4 doubles/cycle via AVX2 _ZGVdN4v_exp.
// Otherwise: scalar std::exp.
inline void bulk_scaled_exp(double scale,
                             const double* __restrict__ u,
                             double*       __restrict__ out,
                             int n) {
#if LBW_HAS_GLIBC_MVEC
    const __m256d vs = _mm256_set1_pd(scale);
    int i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256d vu = _mm256_loadu_pd(u + i);
        __m256d res = _ZGVdN4v_exp(_mm256_mul_pd(vs, vu));
        _mm256_storeu_pd(out + i, res);
    }
    for (; i < n; ++i) out[i] = std::exp(scale * u[i]);
#else
    for (int i = 0; i < n; ++i) out[i] = std::exp(scale * u[i]);
#endif
}

} // namespace lbw
