# xc1s.16 (CR-H-MV): `-mavx2` → function multiversioning — investigation

**Date:** 2026-07-04 · **Verdict:** drafted approach NOT viable; recommend defer (won't-fix as drafted). **Priority:** P4.

## Question

Replace the global `-mavx2` build flag with per-function multiversioning
(`__attribute__((target_clones("avx2","default")))`) on the AVX2 hot kernels, to
(1) add runtime dispatch (portability: an AVX2-compiled binary run on a non-AVX2 CPU
SIGILLs today — no runtime guard) and (2) drop the global `-m` flag (future CRAN
`.check_make_vars` hygiene).

## Findings

### 1. The drafted target — `target_clones` on `bulk_scaled_exp` — is a perf regression

`bulk_scaled_exp`, `bulk_log`, `bulk_exp_clipped` are **`inline` header kernels**
(`src/lbw_math.hpp:17,43,67`) built around `_mm256_*` intrinsics + glibc libmvec
`_ZGVdN4v_exp`/`_ZGVdN4v_log`, gated on `LBW_HAS_GLIBC_MVEC`. They are called from the
**per-iteration hot loops** of the solvers:

- `src/oris.cpp:1440` `bulk_scaled_exp`; `:1486,:1507,:1654` `bulk_exp_clipped`; `:855,:1700` `bulk_log`
- `src/sinkhorn.cpp`, `src/chebyshev.cpp` — same kernels in their inner loops

A function carrying `target_clones` (or `target("avx2")`) **cannot be inlined** into a
differently-targeted caller — the compiler must emit an indirect call through the IFUNC
resolver. Converting these fully-inlined vectorized loops into an indirect call *per
invocation* on the hottest path is a direct regression. This is the same class of hazard
the project already documents for TU splitting without LTO (`CLAUDE.md`: "cross-TU calls
don't inline; move only COLD code"). These kernels are the definition of HOT.

### 2. libmvec linkage is already `target`-scoped, not flag-scoped

The vectorized-exp symbols are declared `extern "C" __m256d _ZGVdN4v_exp(__m256d)
__attribute__((target("avx2")))` (`lbw_math.hpp:8-9`). The `-mavx2` flag is needed to
*compile the `_mm256_*` intrinsics in the kernel bodies*, not to resolve the libmvec
symbol. So the flag is load-bearing for the intrinsic bodies specifically.

### 3. Severity is low under the R source-install model

`-mavx2` is delivered by a `configure` **feature-test** (`configure:73-88`) that runs on
the machine doing the compile, substituted via the opaque `@MAVX2_FLAG@` placeholder
(`src/Makevars.in:1-8`) which `tools:::.check_make_vars` does not flag. For the standard
CRAN **source** install (compile on the target), build machine == run machine, so the
feature-test correctly disables AVX2 on non-AVX2 CPUs — **no SIGILL**. The hazard is real
only for **binary redistribution** across CPUs (pre-built `.so`/Python wheel shipped to a
lesser CPU); the Python `CMakeLists.txt` gate is x86_64 + build-time feature-test, same
model. The CRAN-NOTE concern is **hypothetical** (current R greps `Makevars.in` paths[1L]
only; the placeholder is already clean).

## Viable alternatives (all heavier than the concern warrants at P4)

- **FMV at the OUTER solve function** (not the inline leaf): multiversion
  `oris_solve`/`sinkhorn_solve`/`chebyshev_ipm` so the whole hot region gets an avx2 clone
  with the leaf kernels inlined *inside* it, plus a default clone. Preserves inlining but
  duplicates large functions and requires the solver bodies to be free of other
  target-incompatible constructs — a substantial refactor + bench campaign.
- **Dedicated non-inline AVX2 TU** with a runtime-dispatched entry (`__builtin_cpu_supports`
  + a scalar fallback TU). Removes the global flag and adds dispatch, but reintroduces a
  cross-TU call for the kernels (same inlining loss as §1 unless the batch size is large
  enough to amortize — needs measurement).

## Recommendation

Defer. The drafted `target_clones`-on-the-leaf approach is a **regression** and must not
ship. The viable alternatives are architectural (outer-function FMV or a dispatched TU)
and warrant their own design + stepstone bench campaign — disproportionate to a P4
concern that is largely mitigated by the source-install feature-test. Keep the current
build-time gate. Revisit only if binary redistribution to non-AVX2 targets becomes a
supported channel, at which point outer-function FMV is the preferred SOTA path.
