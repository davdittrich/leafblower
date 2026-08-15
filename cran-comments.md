# CRAN comments

## Submission type

This is an intermediate GitHub + r-universe release, not a CRAN web-form
submission (D-05). `R CMD check --as-cran` is nonetheless run and must pass
as the honest quality bar regardless of publication channel (D-06); the
result below is the real, observed output of that run against the cleaned
package tree, not a placeholder.

## Test environments

* R 4.6.1, x86_64-pc-linux-gnu (Arch Linux), local build, `g++ (GCC) 16.2.1`,
  compiled with `--as-cran`.
* A CI matrix (GitHub Actions, R release/devel x Ubuntu, and a Python
  3.9-3.13 matrix on the Python side) is authored in phase 05 plans 05-03/
  05-04 and will extend coverage beyond this one local environment.

## R CMD check results

Real output from `R CMD build . && R CMD check --as-cran leafblower_0.1.0.tar.gz`
against the hygiene-cleaned tree (git-tracked dev artifacts removed, 35
additional `.Rbuildignore` patterns added for tracked and untracked
non-package files that `R CMD build` was otherwise sweeping into the
tarball):

```
Status: 1 WARNING, 3 NOTEs
```

* NOTE (CRAN incoming feasibility): `New submission`, plus `Suggests or
  Enhances not in mainstream repositories: autumn` -- expected. `autumn` is
  the upstream package `leafblower::harvest()` is a drop-in replacement for
  and is not itself on CRAN; this is the "at most the new-submission NOTE"
  target this check run is measured against.
* NOTE (compilation flags used): `-Werror=format-security`, `-Wformat`,
  `-Wp,-D_FORTIFY_SOURCE=3`, `-Wp,-D_GLIBCXX_ASSERTIONS`, `-march=x86-64`,
  `-mavx2`, `-mno-omit-leaf-frame-pointer`. Re-running with
  `R_MAKEVARS_USER=/dev/null` (bypassing this developer's personal
  `~/.R/Makevars`, which injects `-march=native -mtune=native`) confirmed
  the remaining flags come from this Arch Linux R installation's own system
  `Makeconf` -- not from this package's `Makevars.in`/`configure`, which set
  no `-march`/`-O` flags of their own (see "Notes on build configuration"
  below). `-mavx2` is the one package-controlled flag in that list and is
  load-bearing (see below); it is feature-tested by `configure` and only
  substituted in on hosts where it compiles.
* NOTE (HTML version of manual): `no command 'tidy' found` / `package 'V8'
  unavailable` -- this local machine lacks HTML Tidy and the R `V8` package,
  used only for check-time HTML/MathJax validation of the manual, not for
  building it (the PDF manual check above passed `OK`).
* WARNING (top-level files): `A complete check needs the 'checkbashisms'
  script` -- this local machine lacks `devscripts`' `checkbashisms`, used
  only to lint the package's `configure`/`cleanup` shell scripts for
  non-portable bashisms; not evidence of an actual bashism.

The HTML-manual NOTE and checkbashisms WARNING are local check-environment
tooling gaps (absent optional dependencies of `R CMD check` itself, not of
the package), not package defects; the CI matrix authored in 05-03/05-04
runs on a standard GitHub Actions image where these tools are present and is
expected to close both. `_R_CHECK_FORCE_SUGGESTS_=false` was set for this
run because `PracTools` isn't installed locally and this machine's `autumn`
install is pinned at 0.1 (`Suggests` requires `>= 0.2.0`) -- both are
`Suggests`-only, optional comparison/benchmark dependencies with no effect
on package correctness.

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
