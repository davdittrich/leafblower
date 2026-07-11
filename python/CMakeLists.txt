cmake_minimum_required(VERSION 3.18)
project(leafblower_python CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(pybind11 2.11 REQUIRED)
# If building manually: cmake -Dpybind11_DIR=$(python -m pybind11 --cmakedir) ..
# scikit-build-core sets this automatically via pyproject.toml.

# --- Hermetic feature detection (CR-E1) ---------------------------------------
# Reproduce ./configure's probes so the Python build generates its OWN
# lbw_config.h instead of depending on R's gitignored, machine-specific header.
# Of the three macros, only LBW_HAS_GLIBC_MVEC gates src/ code (the AVX2 exp path
# in lbw_math.hpp) — it is the parity-critical one and is probed, never hardcoded.
include(CheckCXXSourceCompiles)
include(CheckCXXSourceRuns)

# LBW_HAS_OMP_SIMD: compiler accepts -fopenmp-simd (vectorisation hints only;
# no src/ code is #if-gated on it). Mirrors configure's simd_test.cpp.
set(CMAKE_REQUIRED_FLAGS "-fopenmp-simd")
check_cxx_source_compiles(
    "int main(){double a[8]={};\n#pragma omp simd\nfor(int i=0;i<8;i++)a[i]=i*2.0;\nreturn a[0]!=0.0;}"
    LBW_SIMD_OK)
unset(CMAKE_REQUIRED_FLAGS)
if(LBW_SIMD_OK)
    set(LBW_HAS_OMP_SIMD 1)
else()
    set(LBW_HAS_OMP_SIMD 0)
endif()

# LBW_HAS_OMP: full OpenMP threading runtime. The extension links no OMP runtime
# (matching R here, whose SHLIB_OPENMP_CXXFLAGS is empty) and no src/ code is
# #if-gated on it; fixed to 0 to preserve R<->Python parity/determinism.
set(LBW_HAS_OMP 0)

# LBW_HAS_GLIBC_MVEC: glibc libmvec _ZGVdN4v_exp (AVX2, 4-wide double exp).
# THE parity-critical macro. Mirrors configure's mvec_test.cpp: compile+link+run
# with -mavx2 -lm -lmvec. Only meaningful on x86_64 (intrinsics are x86-only).
set(LBW_HAS_GLIBC_MVEC 0)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    set(CMAKE_REQUIRED_FLAGS "-mavx2")
    set(CMAKE_REQUIRED_LIBRARIES "-lm" "-lmvec")
    check_cxx_source_runs(
        "#include <immintrin.h>\nextern \"C\" __m256d _ZGVdN4v_exp(__m256d);\nint main(){__m256d x=_mm256_set1_pd(1.0);__m256d r=_ZGVdN4v_exp(x);volatile double s=((double*)&r)[0];return s>0.0?0:1;}"
        LBW_MVEC_OK)
    unset(CMAKE_REQUIRED_FLAGS)
    unset(CMAKE_REQUIRED_LIBRARIES)
    if(LBW_MVEC_OK)
        set(LBW_HAS_GLIBC_MVEC 1)
    endif()
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lbw_config.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/generated/lbw_config.h @ONLY)
message(STATUS "leafblower hermetic lbw_config.h: OMP_SIMD=${LBW_HAS_OMP_SIMD} "
               "OMP=${LBW_HAS_OMP} GLIBC_MVEC=${LBW_HAS_GLIBC_MVEC}")

# Sources: all src/*.cpp EXCEPT r_bridge.cpp (includes Rinternals.h unavailable here)
set(CORE_SOURCES
    ../src/c_api.cpp
    ../src/calib_linalg.cpp
    ../src/calib_validate.cpp
    ../src/cell_table.cpp
    ../src/chebyshev.cpp
    ../src/greenkhorn.cpp
    ../src/greg.cpp
    ../src/oris.cpp
    ../src/oris_trajectory.cpp
    ../src/trajectory.cpp  # STUDY-BRANCH-ONLY-DO-NOT-MERGE
    ../src/oris_finalize.cpp
    ../src/logit.cpp
    ../src/logit_calib.cpp
    ../src/newton_calib.cpp
    ../src/raking.cpp
    ../src/sinkhorn.cpp
    ../src/design_effect.cpp
    ../src/validation.cpp
)

pybind11_add_module(_leafblower
    leafblower/_bindings.cpp
    ${CORE_SOURCES}
)

# CR-E1: the CMake-generated hermetic header comes FIRST. It is also force-
# included below so that a stale, R-generated ../src/lbw_config.h (quote-includes
# `#include "lbw_config.h"` resolve to the source file's own directory first)
# cannot shadow it — the force-include sets the LBW_CONFIG_H guard so any later
# `#include "lbw_config.h"` finding the src/ copy is a no-op.
target_include_directories(_leafblower PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/generated ../src)
target_compile_options(_leafblower PRIVATE -include ${CMAKE_CURRENT_BINARY_DIR}/generated/lbw_config.h)
target_compile_definitions(_leafblower PRIVATE LBW_NO_R)
find_package(LAPACK REQUIRED)
target_link_libraries(_leafblower PRIVATE ${LAPACK_LIBRARIES})
# -mavx2 required: oris.cpp, sinkhorn.cpp, chebyshev.cpp use _mm256_* intrinsics
# via bulk_scaled_exp when LBW_HAS_GLIBC_MVEC=1 (default in lbw_config.h).
# Without it, gcc errors "inlining failed in call to 'always_inline'
# '_mm256_set1_pd': target specific option mismatch".
target_compile_options(_leafblower PRIVATE -O3)
# -mavx2 only on x86_64; the _mm256_* intrinsics are x86-only and #ifdef-guarded
# (LBW_HAS_GLIBC_MVEC) so ARM/non-AVX2 builds link without it.
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    target_compile_options(_leafblower PRIVATE -mavx2)
endif()

install(TARGETS _leafblower LIBRARY DESTINATION leafblower)
