# Python build variants for the WU-11 benchmark run matrix (WU-8b)

## Why this exists

`python/CMakeLists.txt` hardcodes `-O3` unconditionally (was line 99) and
`-mavx2` on x86_64 (lines 102-104). There is no build-time toggle for a
portable (`-O2`, generic-tuning) build vs. a locally-tuned (`-march=native`)
build — the R arm gets this choice for free via `~/.R/Makevars` flags (no
source edit needed), but the Python/scikit-build-core arm has nothing
equivalent to override from outside the source tree.

WU-11 (the actual production timing run, out of scope here) wants to compare
a **portable** wheel against a **natively-tuned** wheel in the benchmark run
matrix. This requires a source-level change. Per the constraint that
`python/CMakeLists.txt` must never carry this change on `HEAD`, the change
lives as a **transient, reviewed patch** applied only for the duration of one
build, then reverted immediately.

## The patch

`portable_build.patch` (unified diff against `python/CMakeLists.txt`):

- Adds `option(LBW_PORTABLE "Portable -O2, generically-tuned build (vs -O3
  -march=native)" OFF)` near the top of the file.
- Replaces the unconditional `target_compile_options(_leafblower PRIVATE
  -O3)` with:
  - `LBW_PORTABLE=ON`  -> `-O2` (portable)
  - `LBW_PORTABLE=OFF` (default) -> `-O3 -march=native` (native)
- Leaves the `-mavx2` block (lines 102-104 pre-patch) untouched and
  unconditional in both variants: `oris.cpp`, `sinkhorn.cpp`, and
  `chebyshev.cpp` call `_mm256_*` intrinsics via `bulk_scaled_exp` whenever
  `LBW_HAS_GLIBC_MVEC=1` (the default), and gcc refuses to inline those
  intrinsics without `-mavx2` ("target specific option mismatch"). A portable
  build that drops `-mavx2` does not compile — it is not an option.
- Touches nothing else. No `src/`, `r_bridge.cpp`, `R/`, `python/leafblower/`,
  `DESCRIPTION`, or `NAMESPACE` changes.

Apply with `git apply`, revert with `git apply -R` — both directions were
verified to round-trip cleanly (`git apply --check` / `git apply -R --check`)
before use.

## Passing `LBW_PORTABLE` through scikit-build-core

scikit-build-core (>=0.8, per `python/pyproject.toml`) forwards CMake cache
defines via `--config-settings=cmake.define.<VAR>=<VALUE>` (or the `-C`
shorthand), and enables verbose build-tool output via
`--config-settings=build.verbose=true`. Verified against the installed
scikit-build-core 1.0.1 in this repo's isolated build env (visible in the
build log as `Building leafblower ... scikit_build_core-1.0.1...`):

```bash
cd python
uv pip install --python .venv/bin/python -e . --reinstall-package leafblower --no-cache \
    --config-settings="cmake.define.LBW_PORTABLE=ON" \
    --config-settings="build.verbose=true"
```

Two notes specific to this repo:
- `python -m uv` does **not** work (`uv` is not a module in the venv's
  Python — it's the standalone `/usr/bin/uv` binary). Use `uv pip install
  --python .venv/bin/python ...` instead of activating the venv, or
  `uv --python`.
- `uv -v` (verbose on `uv` itself, not just `--config-settings=build.verbose`)
  is required for uv to relay the underlying ninja/gcc compile lines to
  stdout at all — otherwise uv swallows the build backend's subprocess
  output on a successful build. `build_lbw_variant.sh` passes both.
- `--no-cache` forces a full recompile on every invocation. Without it, uv's
  build cache means re-running the same variant produces no compile lines to
  grep (cache hit, nothing to build) — defeating the flag-proof step.

## Pre-existing artifacts in the compile line (NOT from this patch)

The captured compile lines contain flags this patch did not add:

```
... -O3 -DNDEBUG ... -include .../lbw_config.h -O2 -mavx2 -flto=auto -fno-fat-lto-objects ...
```

- `-O3 -DNDEBUG` earlier in the line is CMake's implicit
  `CMAKE_CXX_FLAGS_RELEASE` (triggered by `cmake.build-type = "Release"` in
  `pyproject.toml`'s `[tool.scikit-build]`), unconditional and present in the
  **unpatched** build too (verified: baseline unpatched build's compile line
  is `... -O3 -DNDEBUG ... -O3 -mavx2 -flto=auto ...`).
- `-flto=auto -fno-fat-lto-objects` is `pybind11_add_module()`'s own default
  Release-mode LTO, also present in the unpatched baseline build. Out of
  scope for WU-8b (constraint: change only the O-level/tuning flags) and
  left untouched.

Because gcc applies the **last** `-O` flag seen on the command line, and
CMake appends `target_compile_options()` flags after the implicit
`CMAKE_CXX_FLAGS_RELEASE`, the value from `target_compile_options(_leafblower
PRIVATE -O2)` (or `-O3 -march=native`) is the one that actually governs
codegen despite the earlier `-O3` token. `build_lbw_variant.sh` verifies this
correctly: it extracts every `-O[0-3sz]` token from the compile line in
order and checks the **last** one, not mere substring presence.

## `build_lbw_variant.sh` procedure

```bash
benchmarks/study/python/build_lbw_variant.sh portable   # or: native
```

For each variant, the script:

1. Asserts `python/CMakeLists.txt` is pristine before touching anything
   (refuses to layer the patch on top of an already-dirty file).
2. `git apply`s `portable_build.patch`.
3. Builds the wheel with `LBW_PORTABLE=ON` (portable) or `LBW_PORTABLE=OFF`
   (native), single-thread BLAS (`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=
   MKL_NUM_THREADS=1`), full verbose logging.
4. Extracts the `oris.cpp` compile line from the log and asserts the
   *governing* `-O` flag matches the variant (`-O2` for portable, `-O3` for
   native), asserts `-march=native` is present only for native, and asserts
   `-mavx2` is present in **both**.
5. `git apply -R`s the patch (also happens via an `EXIT` trap if any earlier
   step fails, so a failed build never leaves `CMakeLists.txt` dirty).
6. Asserts `git diff --quiet HEAD -- python/CMakeLists.txt` and fails loudly
   otherwise.

Each invocation leaves the *installed* wheel built with that variant's flags
(useful for immediately smoke-testing it), and always leaves
`python/CMakeLists.txt` byte-identical to `HEAD`.

## Restoring the normal dev build afterwards

The script only reverts the CMakeLists patch — it does not restore the
*installed* wheel to the plain unpatched build. After using the script (or
before starting other work), rebuild the ordinary dev wheel on the pristine
(unpatched) tree:

```bash
cd python && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  uv pip install --python .venv/bin/python -e . --reinstall-package leafblower
```

This reproduces the exact pre-WU-8b behavior (`-O3`, no `-march=native`,
`-mavx2` on x86_64) since it runs against unmodified `HEAD`.

## Verified end-to-end (2026-07-09)

- Portable compile line captured, governing flag `-O2`, `-mavx2` present,
  no `-march=native`. Portable wheel imported and ran `leafblower.harvest()`
  successfully (`weights.sum() == n`).
- Native compile line captured, governing flag `-O3`, `-mavx2` present,
  `-march=native` present.
- `python/CMakeLists.txt` asserted pristine (`git diff --quiet HEAD --
  python/CMakeLists.txt`) after each variant's revert, and after the final
  default-wheel rebuild.
- `git status --short benchmarks/study/python/` shows only the 3 new files
  added by WU-8b (`portable_build.patch`, `build_lbw_variant.sh`,
  this README) — no residual diff to `python/CMakeLists.txt` or any other
  tracked file.
