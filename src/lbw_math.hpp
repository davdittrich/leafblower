#pragma once
#include "lbw_config.h"
#include <algorithm>
#include <cmath>

#if LBW_HAS_GLIBC_MVEC
#  include <immintrin.h>
extern "C" __m256d _ZGVdN4v_exp(__m256d) __attribute__((target("avx2")));
extern "C" __m256d _ZGVdN4v_log(__m256d) __attribute__((target("avx2")));
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

// bulk_log: out[i] = log(in[i]) for i in [0, n). in[i] must be > 0.
// Guard: non-positive values are clamped to 1e-300 floor (log(1e-300) ≈ -690,
// finite) to prevent NaN/-Inf from propagating into convergence logic.
// CR-D3 (j7x8.3): NOT __restrict__ — compute_weight_kl (calib_dispatch.hpp) calls
// this in-place as bulk_log(ratio_buf, ratio_buf, n). The map is per-element
// (out[i] = log(clamp(in[i])), same index, idempotent clamp) so aliasing is
// functionally safe, but __restrict__ made it formal UB. Removing it is free: the
// SIMD log dominates and the distinct-buffer oris callsites are unaffected.
inline void bulk_log(const double* in,
                     double*       out,
                     int n) {
#if LBW_HAS_GLIBC_MVEC
    // MVEC path: pre-scan scalar to clamp non-positives before SIMD log.
    // Avoid a separate clamped copy — write clamped values into out[], then
    // log in-place via the scalar tail. For the SIMD body, any non-positive
    // input would produce -Inf from _ZGVdN4v_log; scalar pre-scan avoids it.
    for (int i = 0; i < n; ++i) {
        out[i] = (in[i] > 0.0) ? in[i] : 1e-300;
    }
    int i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256d v = _mm256_loadu_pd(out + i);
        _mm256_storeu_pd(out + i, _ZGVdN4v_log(v));
    }
    for (; i < n; ++i) out[i] = std::log(in[i] > 0.0 ? in[i] : 1e-300);
#else
    for (int i = 0; i < n; ++i) out[i] = std::log(in[i] > 0.0 ? in[i] : 1e-300);
#endif
}

// bulk_scaled_log: out[i] = log(scale * in[i]) for i in [0, n).
// scale > 0, in[i] > 0 required.
inline void bulk_scaled_log(double scale,
                             const double* __restrict__ in,
                             double*       __restrict__ out,
                             int n) {
#if LBW_HAS_GLIBC_MVEC
    const __m256d vs = _mm256_set1_pd(scale);
    int i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256d v = _mm256_mul_pd(_mm256_loadu_pd(in + i), vs);
        _mm256_storeu_pd(out + i, _ZGVdN4v_log(v));
    }
    for (; i < n; ++i) out[i] = std::log(scale * in[i]);
#else
    for (int i = 0; i < n; ++i) out[i] = std::log(scale * in[i]);
#endif
}

// bulk_exp_clipped: out[i] = exp(clamp(in[i], -clip, clip)) for i in [0, n).
// clip must be <= 700.0 (exp(700) < DBL_MAX).
inline void bulk_exp_clipped(const double* __restrict__ in,
                              double*       __restrict__ out,
                              int n, double clip) {
#if LBW_HAS_GLIBC_MVEC
    const __m256d vclip  = _mm256_set1_pd(clip);
    const __m256d vnclip = _mm256_set1_pd(-clip);
    int i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256d v = _mm256_loadu_pd(in + i);
        v = _mm256_min_pd(_mm256_max_pd(v, vnclip), vclip);
        _mm256_storeu_pd(out + i, _ZGVdN4v_exp(v));
    }
    for (; i < n; ++i) out[i] = std::exp(std::clamp(in[i], -clip, clip));
#else
    for (int i = 0; i < n; ++i) out[i] = std::exp(std::clamp(in[i], -clip, clip));
#endif
}

} // namespace lbw
