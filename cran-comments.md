# CRAN comments

## Submission type

This is an intermediate GitHub + r-universe release, not a CRAN web-form
submission (D-05). `R CMD check --as-cran` is nonetheless run and must pass
as the honest quality bar regardless of publication channel (D-06); the
result below is the real, observed output of that run against the cleaned
package tree, not a placeholder.

## Test environments

* R 4.6.1, x86_64-pc-linux-gnu (Arch Linux), local build, `g++ (GCC) 16.2.1`,
  compiled with `--as-cran`, `checkbashisms`/`tidy`/`pandoc` installed:
  **0 errors, 0 warnings, 1 NOTE**.
* R 4.6.1, x86_64-pc-linux-gnu (Ubuntu 24.04, GitHub Actions `ubuntu-latest`),
  real CI run of `.github/workflows/r-check.yml`
  (https://github.com/davdittrich/leafblower/actions), `gcc 13.3.0`,
  `--as-cran`, TinyTeX-built PDF manual: **0 errors, 0 warnings, 2 NOTEs**.
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

Real output from `R CMD build . && R CMD check --as-cran leafblower_0.1.0.tar.gz`
against the hygiene-cleaned tree (git-tracked dev artifacts removed, 35
additional `.Rbuildignore` patterns added for tracked and untracked
non-package files that `R CMD build` was otherwise sweeping into the
tarball).

**Local (Arch Linux, `checkbashisms`/`tidy`/`pandoc` installed):**

```
Status: 1 NOTE
```

**CI (GitHub Actions `ubuntu-latest`, real run):**

```
Status: 2 NOTEs
```

Both environments are clean of errors and warnings. The local run is now
the stronger result: with `checkbashisms`/`tidy`/`pandoc` installed, the
checkbashisms WARNING and the HTML-manual NOTE this cran-comments.md
previously reported both close, leaving only the one package-controlled
NOTE below. CI still carries a second NOTE (HTML manual) because
`ubuntu-latest` doesn't ship HTML Tidy by default and the CI workflow
doesn't install it (only TinyTeX, for the PDF manual, is installed there).

* NOTE (compilation flags used): `-mavx2` -- feature-tested by `configure`
  and only substituted into `PKG_CXXFLAGS` on hosts where it compiles (see
  "Notes on build configuration" below); load-bearing for the SIMD
  intrinsics in `oris.cpp`/`sinkhorn.cpp`/`chebyshev.cpp`. Present on both
  the local machine and CI's `ubuntu-latest`. (This developer's Arch Linux
  system `Makeconf` previously added extra flags to the local run's NOTE
  text -- `-Werror=format-security`, `-Wp,-D_FORTIFY_SOURCE=3`, etc. --
  confirmed via `R_MAKEVARS_USER=/dev/null` to come from the system
  `Makeconf`, not the package's own `Makevars.in`/`configure`; only
  `-march=native` now shows locally, from this developer's personal
  `~/.R/Makevars`, same non-package-controlled category.)
* NOTE (HTML version of manual, CI only): `no command 'tidy' found` /
  `package 'V8' unavailable` -- `ubuntu-latest` doesn't ship HTML Tidy by
  default; used only for check-time HTML/MathJax validation of the manual,
  not for building it (the PDF manual check passes `OK` on both
  environments).

The CRAN-incoming-feasibility NOTE (`New submission`, `autumn` not on
CRAN) that appeared locally does not appear in the CI output above because
`r-lib/actions/check-r-package` disables `_R_CHECK_CRAN_INCOMING_` by
default; it is still expected and accepted on an actual CRAN submission
for the reason given locally: `autumn` is the upstream package
`leafblower::harvest()` is a drop-in replacement for and is not itself on
CRAN. `_R_CHECK_FORCE_SUGGESTS_=false` is set in both environments because
`autumn`, `bench`, `lhs`, `DiceKriging`, `ggplot2`, `survey`, `PracTools`,
and `arrow` are optional comparison/benchmark-only `Suggests` with no
effect on package correctness; every test referencing them is
`skip_if_not_installed()`-guarded (or lives only in the
`.Rbuildignore`d `tests/testthat/fixtures/` directory).

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
