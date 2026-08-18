# CRAN comments

## Submission type

This is an intermediate GitHub + r-universe release, not a CRAN web-form
submission (D-05). `R CMD check --as-cran` is nonetheless run and must pass
as the honest quality bar regardless of publication channel (D-06); the
result below is the real, observed output of that run against the cleaned
package tree, not a placeholder.

## Suggests availability correction

An earlier CI configuration (`.github/workflows/r-check.yml` before commit
`790341e`) set `_R_CHECK_FORCE_SUGGESTS_: false` and
`setup-r-dependencies`'s `dependencies: '"hard"'`, together downgrading
*every* unavailable Suggests package from a check ERROR to a silent skip --
not just the one genuinely-unresolvable entry (`autumn (>= 0.2.0)`, a
GitHub-only package at v0.1 this repository never actually calls -- only
roxygen prose and `.Rbuildignore`d `benchmarks/` reference it) the override
was added to route around. That configuration produced a green
`R CMD check` result without PracTools, survey, arrow or DiceKriging ever
being installed, so the ~25 tests those packages guard never ran in CI even
though the run reported success.

Commit `8e28df1` removed the unresolvable `autumn` Suggests entry (the
actual root cause -- it is unsatisfiable from CRAN, r-universe, or any
`Additional_repositories` feed). Commit `790341e` removed both masking
settings, restoring `R CMD check`'s default Suggests-required behavior: an
unavailable Suggests package now fails the check instead of silently
widening the skip set. The results below are from that corrected
configuration -- GitHub Actions run **32080650203** is the first CI run
under it, and its log confirms PracTools 1.7.5, DiceKriging 1.6.1,
arrow 25.0.0 and survey 4.5 installing, and the testthat run completing
with 0 skips attributable to any missing Suggests package
(`FAIL 0 | WARN 130 | SKIP 31 | PASS 1794`; all 31 skips are `On CRAN`
gates or local-only benchmark-fixture guards, none reference an
uninstalled package).

## Test environments

* R 4.6.1, x86_64-pc-linux-gnu (Arch Linux), local build, `g++ (GCC) 16.2.1`,
  compiled with `--as-cran`, full Suggests set installed (all 12 declared
  Suggests packages, including PracTools, survey, arrow and DiceKriging --
  no dependency-resolution override, `checkbashisms`/`pandoc` installed):
  **0 errors, 0 warnings, 3 NOTEs**.
* R 4.6.1, x86_64-pc-linux-gnu (Ubuntu 24.04, GitHub Actions `ubuntu-latest`),
  real CI run of `.github/workflows/r-check.yml`, run 32080650203
  (https://github.com/davdittrich/leafblower/actions/runs/32080650203),
  `gcc 13.3.0`, `--as-cran`, full Suggests set installed (the workflow's
  `setup-r-dependencies` step now resolves and installs every declared
  Suggests package -- see "Suggests availability correction" above),
  TinyTeX-built PDF manual: **0 errors, 0 warnings, 2 NOTEs**.
* R 4.6.1, x86_64-w64-mingw32 (Windows Server, GitHub Actions
  `windows-latest`, Rtools45), real CI run of
  `.github/workflows/r-check.yml`, run 32082804108
  (https://github.com/davdittrich/leafblower/actions/runs/32082804108),
  full Suggests set installed, `configure.win` generating
  `src/lbw_config.h`/`src/Makevars.win`: package compiles, links, and the
  full testthat suite passes (previously this failed to compile at all --
  `lbw_config.h: No such file or directory`, leafblower-sc9t, now closed).
  **1 warning, 2 NOTEs** -- not a clean run: `checking PDF version of
  manual ... WARNING` (the with-index LaTeX pass; the without-index retry
  succeeds and produces a valid PDF). This warning is a TinyTeX-on-Windows
  packaging gap unrelated to `configure.win` or the compile path, tracked
  separately on leafblower-fxyj -- not resolved by this change and not
  claimed as resolved here. arm64 Windows is not covered by this runner
  (GitHub-hosted `windows-latest` is x86_64 only); see r-universe's own
  build farm results below for arm64 coverage.
* r-universe's own multi-platform build farm
  (https://davdittrich.r-universe.dev/api/packages/leafblower), build
  https://github.com/r-universe/davdittrich/actions/runs/32082532844,
  commit `43894c8` (the `configure.win` fix plus both tlmgr-invocation
  fixes; the tlmgr-lookup-path hardening in the final 05-10 commit
  `6a312ce` only touches this repository's own `windows-latest` GitHub
  Actions job, not anything r-universe's build farm runs, so this result
  fully covers it). Every `_jobs` entry checked individually, not the
  aggregate `_status` field (05-06's mistake, per 05-VERIFICATION.md):
  **all 5 Windows configs are `OK`** -- `windows-devel-arm64`,
  `windows-devel-x86_64`, `windows-oldrel-x86_64`, `windows-release-arm64`,
  `windows-release-x86_64`. This is the only oracle available for arm64
  Windows (GitHub-hosted `windows-latest` is x86_64 only). All 4 Linux and
  4 macOS configs plus `source` are also `OK` on this build (the PracTools
  `deffH` NaN mismatch 05-VERIFICATION.md found on r-universe's Linux/macOS
  jobs is not reproduced here; that gap is tracked separately and is
  outside 05-10's scope, unchanged by this plan). `wasm-release` is the
  sole remaining `FAIL` -- a separate Emscripten toolchain this plan does
  not touch, filed as leafblower-soci rather than implied fixed by the
  green Windows result.
* Python 3.9-3.13 wheel matrix (GitHub Actions `ubuntu-latest` +
  `macos-14`/arm64, `.github/workflows/python-wheels.yml`): wheels build,
  pass `twine check`, and import + calibrate cleanly on all 5 versions on
  both platforms -- real CI run, not authored-but-unrun. `macos-13` (Intel)
  was dropped from the matrix: the runner never scheduled on this account
  across 25+ minutes of queued time while `ubuntu-latest` and `macos-14`
  both ran immediately, consistent with GitHub's phase-out of Intel macOS
  runners. x86_64 macOS wheel coverage is not proven by this CI matrix as a
  result (arm64 macOS is).

## R CMD check results

Real output from `R CMD build . && R CMD check --as-cran leafblower_0.1.1.tar.gz`
against the hygiene-cleaned tree (git-tracked dev artifacts removed,
additional `.Rbuildignore` patterns added for tracked and untracked
non-package files that `R CMD build` was otherwise sweeping into the
tarball), with every declared Suggests package installed and the default
missing-Suggests guard left in place (no `_R_CHECK_FORCE_SUGGESTS_`
override).

**Local (Arch Linux, full Suggests installed):**

```
Status: 3 NOTEs
```

**CI (GitHub Actions `ubuntu-latest`, run 32080650203, full Suggests
installed):**

```
Status: 2 NOTEs
```

Both environments are clean of errors and warnings under a
Suggests-complete check (D-06) -- the check that CRAN's own machinery
actually runs, not a weakened one.

* NOTE (CRAN incoming feasibility, local only): `New submission` --
  standard for any package not yet on CRAN. `_R_CHECK_CRAN_INCOMING_` is on
  locally under `--as-cran` and disabled by default in
  `r-lib/actions/check-r-package`, which is why this NOTE does not appear
  in the CI Status above. Expected and accepted on an actual CRAN
  submission.
* NOTE (compilation flags used): `-mavx2` (and, on this developer's
  machine, `-march=native` from a personal `~/.R/Makevars`) -- feature-tested
  by `configure` and only substituted into `PKG_CXXFLAGS` on hosts where it
  compiles (see "Notes on build configuration" below); load-bearing for the
  SIMD intrinsics in `oris.cpp`/`sinkhorn.cpp`/`chebyshev.cpp`. Present on
  both the local machine and CI's `ubuntu-latest`.
* NOTE (HTML version of manual): local: `package 'V8' unavailable`; CI:
  `no command 'tidy' found` / `package 'V8' unavailable`. Both messages come
  from `R CMD check`'s own check-time HTML/MathJax rendering validation
  tooling (HTML Tidy CLI, the V8 R package) -- neither is a package
  Suggests entry, neither affects package correctness, and neither is used
  to *build* the manual, only to validate its rendered HTML. The PDF manual
  check passes `OK` on both environments. `ubuntu-latest` doesn't ship HTML
  Tidy by default and the CI workflow only installs TinyTeX (for the PDF
  manual); this developer's local machine currently lacks the V8 R package.

No dependency-resolution override is set in either environment: every
declared Suggests package (`testthat`, `bench`, `lhs`, `DiceKriging`,
`ggplot2`, `rprojroot`, `survey`, `PracTools`, `arrow`, `callr`, `jsonlite`,
`withr`) is installed and its guarded tests execute for real in both
environments -- see "Suggests availability correction" above.

## Notes on build configuration

R and Python intentionally build with different optimization flags. From
`src/Makevars.in`:

> `-O3` is intentionally NOT set: R supplies the user/site `-O` level via
> `$(CXXFLAGS)` in `$(ALL_CXXFLAGS)`; a user wanting `-O3` sets it in
> `~/.R/Makevars`.

R's own build-supplied `$(LAPACK_LIBS)` and `$(BLAS_LIBS)` (`src/Makevars.in`
line 16, `PKG_LIBS = @MVEC_LIBS@ $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)`)
link against whatever LAPACK/BLAS R itself was configured against, so no
extra `SystemRequirements:` entry is needed on the R side. The Python build
has no such constraint and is a hard dependency instead
(`python/CMakeLists.txt`):

> `find_package(LAPACK REQUIRED)`

with `target_compile_options(_leafblower PRIVATE -O3)` set unconditionally.
This asymmetry is deliberate and documented in `python/CMakeLists.txt`
(citing phase-02 SC2, `leafblower-qzto`) -- CRAN's own portability check
(`tools:::.check_make_vars`) rejects `-O*` flags in `PKG_CXXFLAGS`, so R
cannot hard-set `-O3` the way the Python build does; the R/Python parity
tests treat their tolerances as the bound on how much this asymmetry may
move a result, not as something to equalize away.

`-mavx2` (the one package-set non-portable flag CRAN's check flagged above)
is required because `oris.cpp`, `sinkhorn.cpp`, and `chebyshev.cpp` use
`_mm256_*` intrinsics via `bulk_scaled_exp`; it is substituted into
`PKG_CXXFLAGS` via the `@MAVX2_FLAG@` placeholder only after `configure`
feature-tests that the host compiler accepts it (see `configure` lines
67-88), so it does not appear, and the intrinsics are not compiled, on
hosts where it would fail.
