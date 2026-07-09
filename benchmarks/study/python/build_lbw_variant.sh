#!/usr/bin/env bash
# WU-8b: build a leafblower Python wheel with either portable (-O2) or
# native-tuned (-O3 -march=native) compile flags, via a TEMPORARY,
# transiently-applied patch to python/CMakeLists.txt (LBW_PORTABLE option).
#
# HEAD MUST STAY PRISTINE: the patch is applied, the wheel is built and its
# compile flags are proven from the verbose build log, then the patch is
# reverted and pristine-ness of python/CMakeLists.txt is asserted before
# exit. If any step after "apply" fails, the trap below still reverts.
#
# Usage:
#   benchmarks/study/python/build_lbw_variant.sh {portable|native}
#
# After using this script, restore the venv's normal dev build (plain -O3,
# no -march=native, matching unmodified HEAD) with:
#   cd python && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
#     .venv/bin/python -m uv pip install -e . --reinstall-package leafblower
#
# NOTE ON PARITY: -march=native (the "native" variant) selects host-specific
# codegen (may exceed baseline -mavx2, e.g. FMA/AVX-512 depending on the
# build host) and is intended ONLY for local single-host speed-regression
# comparison in the WU-11 run matrix -- not for producing a distributable
# artifact. -mavx2 is unconditionally mandatory in BOTH variants (the
# _mm256_* intrinsics in oris.cpp/sinkhorn.cpp/chebyshev.cpp require it).

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
CMAKELISTS="${REPO_ROOT}/python/CMakeLists.txt"
PATCH="${REPO_ROOT}/benchmarks/study/python/portable_build.patch"
BUILD_LOG="$(mktemp -t lbw-build-XXXXXX.log)"

VARIANT="${1:-}"
if [[ "${VARIANT}" != "portable" && "${VARIANT}" != "native" ]]; then
    echo "Usage: $0 {portable|native}" >&2
    exit 1
fi

if [[ "${VARIANT}" == "portable" ]]; then
    LBW_DEFINE="LBW_PORTABLE=ON"
    O_FLAG_GREP='-O2'
else
    LBW_DEFINE="LBW_PORTABLE=OFF"
    O_FLAG_GREP='-O3 -march=native'
fi

cd "${REPO_ROOT}"

# --- Preconditions: HEAD must be pristine before we touch anything ---------
if ! git diff --quiet HEAD -- "${CMAKELISTS}"; then
    echo "FATAL: ${CMAKELISTS} is already dirty -- refusing to apply patch on top of it." >&2
    git diff -- "${CMAKELISTS}" >&2
    exit 1
fi

REVERTED=0
revert_patch() {
    if [[ "${REVERTED}" -eq 0 ]]; then
        git apply -R "${PATCH}"
        REVERTED=1
    fi
}
# Always attempt revert on exit (success or failure) so HEAD never stays dirty.
trap 'revert_patch || true' EXIT

echo "==> [${VARIANT}] applying transient patch: ${PATCH}"
git apply "${PATCH}"

echo "==> [${VARIANT}] building wheel (LBW_PORTABLE via ${LBW_DEFINE}, verbose, single-thread BLAS)"
(
    cd "${REPO_ROOT}/python"
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    uv -v pip install --python .venv/bin/python -e . --reinstall-package leafblower --no-cache \
        --config-settings="cmake.define.${LBW_DEFINE}" \
        --config-settings="build.verbose=true" \
        2>&1
) | tee "${BUILD_LOG}"

echo "==> [${VARIANT}] extracting _leafblower.cpp / .cxx compile line from build log"
COMPILE_LINE="$(grep -E '\-c .*(_bindings|oris)\.cpp' "${BUILD_LOG}" | head -1 || true)"
if [[ -z "${COMPILE_LINE}" ]]; then
    echo "FATAL: could not find a compile line in ${BUILD_LOG}" >&2
    exit 1
fi

echo "==> [${VARIANT}] compile line:"
echo "${COMPILE_LINE}"

# NOTE: CMake's implicit CMAKE_BUILD_TYPE=Release flags (-O3 -DNDEBUG) and
# pybind11_add_module()'s own default Release LTO (-flto=auto) both land
# EARLIER on the command line unconditionally -- pre-existing scikit-build-core
# / pybind11 behavior, unrelated to this patch, present in the unpatched build
# too. Our target_compile_options() -O2/-O3 is appended LAST by CMake, and gcc
# applies the LAST -O flag seen on the command line, so it is the one that
# actually governs codegen. We therefore check the LAST -O token, not mere
# substring presence.
LAST_O_FLAG="$(grep -oE -- '-O[0-3sz]' <<<"${COMPILE_LINE}" | tail -1)"
echo "==> [${VARIANT}] last (governing) -O flag: ${LAST_O_FLAG}"

if [[ "${VARIANT}" == "portable" ]]; then
    if [[ "${LAST_O_FLAG}" != "-O2" ]]; then
        echo "FATAL: governing -O flag is '${LAST_O_FLAG}', expected -O2 for portable" >&2
        exit 1
    fi
    if grep -qE -- '-march=native' <<<"${COMPILE_LINE}"; then
        echo "FATAL: -march=native unexpectedly present in portable compile line" >&2
        exit 1
    fi
else
    if [[ "${LAST_O_FLAG}" != "-O3" ]]; then
        echo "FATAL: governing -O flag is '${LAST_O_FLAG}', expected -O3 for native" >&2
        exit 1
    fi
    if ! grep -qE -- '-march=native' <<<"${COMPILE_LINE}"; then
        echo "FATAL: -march=native not found in native compile line" >&2
        exit 1
    fi
fi
if ! grep -qE -- '-mavx2' <<<"${COMPILE_LINE}"; then
    echo "FATAL: -mavx2 not found in ${VARIANT} compile line (mandatory in both variants)" >&2
    exit 1
fi

echo "==> [${VARIANT}] flags verified: build-tag=lbw-py-${VARIANT}-$(date -u +%Y%m%dT%H%M%SZ)"

echo "==> [${VARIANT}] reverting transient patch"
revert_patch

echo "==> [${VARIANT}] asserting python/CMakeLists.txt is pristine"
if ! git diff --quiet HEAD -- "${CMAKELISTS}"; then
    echo "FATAL: ${CMAKELISTS} not pristine after revert!" >&2
    git diff -- "${CMAKELISTS}" >&2
    exit 1
fi
echo "==> [${VARIANT}] PRISTINE. Build log: ${BUILD_LOG}"
