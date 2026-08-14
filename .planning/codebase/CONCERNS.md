<!-- refreshed: 2026-08-14 -->
# Codebase Concerns

**Analysis Date:** 2026-08-14

Scope: C++17 core `src/`, R layer `R/`, Python layer `python/leafblower/`, tests
(`tests/testthat/`, `python/leafblower/test_*.py`, `tests/test_parity_weights.py`).

Notable positive: a repo-wide grep for `TODO`/`FIXME`/`HACK`/`XXX` across
`src/`, `R/`, `python/leafblower/`, `tests/` returns **zero** hits. The debt here
is not marker-comment debt — it is dual-build drift, duplicated invariant
enforcement, and a polluted repository root.

## Tech Debt

**Two independent build source lists for `src/*.cpp`:**
- Issue: R auto-globs `src/*.cpp`; the Python build does not. `python/CMakeLists.txt:60-78`
  hard-codes `CORE_SOURCES`. Adding a new `src/*.cpp` and forgetting the CMake entry
  yields a pybind11 link failure with undefined symbols — and only on the Python side,
  after R tests are already green.
- Files: `python/CMakeLists.txt:60-78`, `src/Makevars.in` (its `PKG_SOURCES` list is
  decorative — R ignores it, so it silently drifts).
- Impact: Python build breaks late; the two lists can also *silently* disagree in the
  other direction (a file in `CORE_SOURCES` but deleted from `src/`).
- Current state: the lists agree today — 17 of 18 `src/*.cpp` are listed; `r_bridge.cpp`
  is correctly excluded (`LBW_NO_R`).
- Fix approach: replace the literal list with a CMake `file(GLOB ../src/*.cpp)` minus
  `r_bridge.cpp`, or add a check-script that diffs `ls src/*.cpp` against `CORE_SOURCES`
  and runs in the test gate.

**Divergent optimization flags between the two builds:**
- Issue: `python/CMakeLists.txt:99` sets `-O3` unconditionally; the R side deliberately
  sets **no** `-O` level (CRAN's `tools:::.check_make_vars` rejects `-O` in `PKG_CXXFLAGS`),
  so R inherits the user's `~/.R/Makevars` — typically `-O2`.
- Files: `python/CMakeLists.txt:99`, `configure`, `src/Makevars.in`.
- Impact: R↔Python parity tests compare binaries built at different optimization levels.
  Parity tolerances (`1e-10`, `1e-6` for logit in `tests/test_parity_weights.py:93`) absorb
  this today, but a future FP-sensitive change can fail parity for reasons unrelated to the
  change.
- Fix approach: document the asymmetry in the parity test header, or pin the CMake build to
  match the R default.

**Duplicated unit-mode water-fill implementation:**
- Issue: `finalize_weights_buf` in `src/calib_dispatch.hpp:359-383+` implements per-cell
  water-fill, its own comment declaring it a "mirror of oris_finalize unit branch"
  (`src/oris_finalize.cpp:19`). Two copies of a numerically delicate invariant.
- Files: `src/calib_dispatch.hpp:359`, `src/oris_finalize.cpp:19`.
- Impact: a fix applied to one copy leaves `bounds_mode="unit"` wrong on the other solver
  family — exactly the class of bug the "no renormalize after water-fill" rule guards.
- Fix approach: extract one shared `water_fill_cell()` into `calib_dispatch.hpp` and call it
  from both. Cold (once-per-solve) code, so the no-LTO inlining constraint does not apply.

**Repository root polluted with tracked non-source artifacts:**
- Files (all tracked): `cell_table_92c4f45.cpp`, `cell_table_92c4f45.hpp`,
  `ieppa_92c4f45.cpp` (commit-hash-suffixed snapshot copies of `src/` files),
  `patch_raking.py`, `patch_wolfe.py` (one-off patch scripts), `test_output.log`,
  `leafblower_0.1.0.tar.gz` (a built package tarball committed into the package),
  `REVIEW_FINDINGS.md` (56 lines) and `code-review-findings.md` (117 lines), last touched
  2026-05-03.
- Impact: **CRAN/build risk** — `.Rbuildignore` excludes `.beads`, `.claude`, `docs`,
  `python`, `conftest.py`, `benchmarks` etc., but does **not** exclude any of the files
  above, so `R CMD build` ships stale `.cpp` snapshots, a nested tarball, and a log into
  the source package. The hash-suffixed `.cpp` copies also confuse grep-based navigation
  (two `cell_table` definitions, only one compiled).
- Fix approach: delete the snapshots/patch scripts/log/tarball; add `.Rbuildignore` entries
  and `.gitignore` rules for `*.tar.gz` and `test_output.log`; fold the two findings
  documents into beads tickets.

**Uncommitted build directories in the working tree:**
- `bwarn.Rcheck` (23M), `perf-maint.Rcheck` (23M), `py-parity.Rcheck` (23M),
  `py-r-parity.Rcheck` (23M), `leafblower.Rcheck` (1.1M), `node_modules`, `Cholesky`.
- Impact: ~93M of check output; `.gitignore` covers `*.Rcheck/` so they are untracked, but
  they make `git status`/`find` noisy and are trivially confused with sources.
- Fix approach: routine `make clean` step; move check output under one ignored directory.

**Reserved-but-unusable enum slots:**
- `src/leafblower.h:41-53` — slot `2` is a documented hole (`/* 2 = removed (was RK_ALG_LBFGSB) */`).
  Slot `7` is an **undocumented** hole (values jump `RK_ALG_GREG = 6` → `RK_ALG_ORIS_SOFT = 8`).
- Impact: reusing either value silently reinterprets serialized/fixture-stored algorithm codes
  as a different solver. Slot 7 has no comment warning anyone off it.
- Fix approach: add a `/* 7 = reserved */` comment mirroring the slot-2 one.

## Known Bugs

**Wrong struct comment on `homotopy_levels_used`:**
- Symptoms: `src/leafblower.h:124` documents "1 for single-pass (n_levels=1), not 0". The
  producing code `src/oris.cpp:620` sets `res.homotopy_levels_used = lvl + 1`, which is
  consistent; but the field's own default is `0` (`src/oris.hpp:25`) and two paths write a
  literal `0` (`src/c_api.cpp:156`, `src/r_bridge.cpp:768`) on the non-ORIS/failure route.
- Files: `src/leafblower.h:124`, `src/oris.hpp:25`, `src/c_api.cpp:156`, `src/r_bridge.cpp:768`.
- Trigger: any non-ORIS method, or an early-exit path — the documented "never 0" contract
  does not hold there.
- Workaround: treat `0` as "not applicable"; do not branch on `== 1` alone.
  Covered by `tests/testthat/test-homotopy-enabled-field.R`.

**`res.base.*` field access (structural footgun, currently correct):**
- Every solver result wraps `CalibResult base;` (`src/oris.hpp:8`, `src/raking.hpp:6`,
  `src/sinkhorn.hpp:8`, `src/chebyshev.hpp:9`, `src/greg.hpp:8`, `src/greenkhorn.hpp:10`,
  `src/logit_calib.hpp:10`, `src/newton_calib.hpp:7`).
- Impact: writing `res.max_error` instead of `res.base.max_error` where a same-named
  shadowing member exists compiles and silently reads the wrong field.
- Current state: the SRAA best-iterate path correctly uses `res.base.metric_first_check`
  (`src/oris.cpp:1197-1200`).

**SRAA best-iterate metric — regression re-introduced twice:**
- Correct code is at `src/oris.cpp:1196`:
  `const double nat_metric = lbw::select_metric(sraa_cfg.metric, cm);` evaluated at
  `kErrCheckInterval` boundaries (`src/oris.cpp:96`), *not* the fast `errRp` proxy.
  `src/sraa.hpp:36-42` documents the rule explicitly.
- Risk: the fast `errRp` proxy is available in the same scope and looks interchangeable;
  the bug has been re-introduced twice historically.
- Guard: `tests/testthat/test-best-iterate.R`, `test-greenkhorn-best-metric.R`,
  `test-logit-best-iterate.R`.

**Lambda `[&]` capture-at-definition in `src/raking.cpp`:**
- `last_F_metrics` and `f_eval_full_metrics` are declared at `src/raking.cpp:245-251`,
  *before* `auto F_eval = [&]` at `src/raking.cpp:261`, with an explicit comment
  ("Declared here (outer scope) so the [&] lambda capture includes it").
- Trigger: adding a new guard variable *after* line 261 and using it inside `F_eval` —
  `[&]` captures only names in scope at the definition site, so this fails to compile in the
  best case and captures a different shadowed name in the worst.

## Security Considerations

**Not a network/credential-handling codebase.** No sockets, no auth, no secrets in tree;
attack surface is untrusted numeric input into the C++ core.

**Untrusted-length input into raw pointer loops:**
- Risk: the C API takes bare `double*` + `int n` (`src/c_api.cpp`, `src/leafblower.h:56+`).
  Length mismatches are caught in the bridges, not in the core.
- Files: `src/r_bridge.cpp`, `src/calib_validate.cpp`, `src/c_api.cpp`.
- Current mitigation: `src/calib_validate.cpp` plus bridge-level checks; regression tests in
  `tests/testthat/test-bridge-length-checks.R`, `test-input-validation.R`,
  `test-cr-d7-nobs-guard.R`, `test-diagnostics-guards.R`, `test-cr-wave10-call-guards.R`.
- Recommendations: direct C-API callers (bypassing R/Python) get no length validation —
  document that `rk_*` entry points assume caller-validated buffers, or add an explicit
  precondition check in `c_api.cpp`.

**File writes from inside the Python `.so`:**
- `src/oris_trajectory.cpp:49` correctly uses `std::FILE* f = std::fopen(path, "w")`. This
  is deliberate: `std::ofstream` from inside the loaded module SIGSEGVs. A repo-wide grep
  confirms **no** `std::ofstream` in `src/`.
- Risk: the trajectory path is caller-supplied and written unvalidated. Low severity
  (local library, caller already has filesystem access), but the path is not canonicalized.
- Smoke coverage: `python/leafblower/test_trajectory_csv_smoke.py`.

## Performance Bottlenecks

**No LTO — TU splits cost inlining:**
- Problem: `-flto` is absent from both `configure` and `src/Makevars.in`. Cross-TU calls do
  not inline, so moving hot per-iteration code out of its caller's TU is a silent regression.
- Files: `src/oris.cpp` (2246 lines) is split this way on purpose — `oris_finalize.cpp` (339)
  and `oris_trajectory.cpp` hold only COLD once-per-solve code; the hot `oris_solve` stays put.
- Improvement path: any further split of `oris.cpp`/`raking.cpp`/`chebyshev.cpp` must move
  cold code only, and must be validated against the stepstone benchmark
  (`tests/testthat/fixtures/stepstone_reference_run.R`, gate in `tests/testthat/test-bench-gate.R`).

**Per-cell `std::vector<std::vector<int>>` build in unit-mode finalization:**
- `src/calib_dispatch.hpp:379-380` allocates one `std::vector<int>` per cell on every
  `finalize_weights_buf` call in `bounds_mode="unit"`.
- Cause: `M_cell` separate heap allocations for what is a CSR-shaped grouping already
  available from `CellTable`.
- Improvement path: reuse a CSR offsets/indices pair from `cell_table.hpp` instead of a
  vector-of-vectors. Cold path, so this is a memory-churn concern rather than a hot-loop one.

**Parity tests shell out to `Rscript` per case:**
- `python/leafblower/test_solver_parity.py:139-152` spawns an `Rscript -e` subprocess
  (90s timeout) per method, as does `tests/test_parity_weights.py:59`.
- Impact: the parity suite is process-spawn-bound; six-plus R interpreter startups per run.
- Improvement path: acceptable for correctness isolation; batch the methods into one
  `Rscript` invocation emitting a JSON map if runtime becomes a gate problem.

## Fragile Areas

**Weight finalization order (normalize → bounds):**
- Files: `src/calib_dispatch.hpp:359-390`, `src/oris_finalize.cpp:19+`.
- Why fragile: the sanctioned order is a single pre-bounds scale to `n`
  (`calib_dispatch.hpp:364-367`) and *then* `bounds_mode` dispatch. Renormalizing **after**
  water-fill silently breaks the `bounds_mode="unit"` clamps — the result still looks
  plausible (Σw=n holds) while per-obs bounds are violated.
- Safe modification: never add a normalization pass after line 376. The degenerate guard
  `total_w > kMinSafeTotalWeight` (1e-100, `calib_dispatch.hpp:314-318`) exists to stop
  subnormal totals from overflowing `norm` to `+inf` — do not relax it.
- Test coverage: `tests/testthat/test-clamp-contract.R`, `test-harvest-bounds-mode.R`,
  `test-oris-bounds-mode.R`, `python/leafblower/test_returned_weights_invariant.py`.

**Two solver formulas that look wrong and are not:**
- The chebyshev Mehrotra corrector's linear `y·Δs_aff` term (`src/chebyshev.cpp`) and the
  ORIS ALM Newton step `X̃(1−λ+μz)/(1+ρ)` (`src/oris.cpp`). Both are guarded by code
  comments; reviewer "fixes" to add/drop terms have come from the wrong divergence or
  derivation.
- Safe modification: verify against the actual (un-normalized-KL) generator before touching.

**`src/oris.cpp` — 2246 lines, the largest TU:**
- Contains homotopy levels, SRAA, best-iterate selection, adaptive omega, and the flat BCD
  mode-2 path in one function body (SRAA state alone spans `src/oris.cpp:936-1041`).
- Why fragile: the no-LTO constraint forbids the obvious remedy (splitting it), so the file
  can only grow. Interacting state (`sraa_best_errRp`, `next_check`, `nat_metric_prev_sraa`,
  per-level resets at `src/oris.cpp:989`, `:1007`, `:1041`) has to be reset consistently at
  every level boundary.

**`R/harvest.R` — 1269 lines, and `src/r_bridge.cpp` — 1274 lines:**
- The R↔C boundary hand-marshals a 10+ field result list by positional index
  (`src/r_bridge.cpp:1027-1031` — `SET_STRING_ELT(res_names, 10, ...)` paired with
  `SET_VECTOR_ELT(res_list, 10, ...)`).
- Why fragile: inserting a field mid-list requires renumbering every subsequent index in two
  parallel call sequences; a mismatch mislabels a result field with no compiler help.
- Test coverage: `tests/testthat/test-calibration-result-names.R`,
  `test-calib-result-consolidation.R` pin the names.

**`design_weights=` vs `weights=`:**
- `harvest()` accepts `design_weights=` for the per-observation `d_i`. There is **no**
  `weights=` argument — passing it lands in `...` and is silently ignored, producing a
  plausible unweighted result.
- Fix approach: add an explicit `...`-name check in `R/harvest.R` that errors on `weights`.

## Scaling Limits

**Single-threaded solve, BLAS threading must be pinned:**
- Current capacity: deterministic results require `OMP_NUM_THREADS`,
  `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS` all `=1` — in Python they must be set *before*
  `import numpy` (`conftest.py`).
- Limit: without the pin, R↔Python parity and stepstone benchmarks drift; with it, no
  multicore scaling is available at all.
- Scaling path: parallelism would have to be introduced inside the BCD margin loop with a
  deterministic reduction order, not by unpinning BLAS.

**AVX2 fast path is x86-only:**
- `python/CMakeLists.txt:95-104` adds `-mavx2` only when `CMAKE_SYSTEM_PROCESSOR` matches
  `x86_64|AMD64`; `_mm256_*` intrinsics in `oris.cpp`/`sinkhorn.cpp`/`chebyshev.cpp` are
  `#ifdef`-guarded by `LBW_HAS_GLIBC_MVEC`.
- Limit: ARM builds fall back to scalar `bulk_scaled_exp` with no vectorized replacement.
- Scaling path: an SVE/NEON path, or rely on the compiler's autovectorizer and measure.

## Dependencies at Risk

**LAPACK, hard `REQUIRED`:**
- `python/CMakeLists.txt:93` — `find_package(LAPACK REQUIRED)`. No fallback.
- Impact: no LAPACK, no Python wheel. The R side gets LAPACK via R itself, so the two
  builds can resolve to *different* LAPACK implementations on one machine — another source
  of parity drift not covered by any test.
- Migration plan: record the resolved LAPACK library in both builds and assert they match in
  the parity conftest, or vendor a reference path for the few routines used.

**`lbw_config.h` shadowing hazard:**
- `python/CMakeLists.txt:85-91` documents an active workaround: a stale R-generated
  `src/lbw_config.h` would win quote-include resolution, so the CMake-generated header is
  force-included first to set the guard.
- Impact: the workaround is load-bearing; removing the `-include` flag silently reintroduces
  a stale-config build. `.Rbuildignore` correctly excludes `src/lbw_config.h`.

## Missing Critical Features

**No automated check that the two build source lists agree:**
- Problem: nothing fails when `src/*.cpp` and `python/CMakeLists.txt:CORE_SOURCES` diverge —
  it surfaces as a link error at Python build time.
- Blocks: confident single-language iteration; every new solver file needs the 8-step manual
  checklist in `CLAUDE.md` followed exactly.

**No version-sync automation:**
- `DESCRIPTION` and `python/pyproject.toml` versions are bumped by hand with no check.
- Blocks: releasing R and Python artifacts that claim the same version.

## Test Coverage Gaps

**`oris_soft` (RK_ALG_ORIS_SOFT = 8) has no weight-vector parity test:**
- What's not tested: `tests/test_parity_weights.py:73-75` parametrizes over
  `greenkhorn, logit, raking, oris, sinkhorn, newton_kl`;
  `python/leafblower/test_solver_parity.py` covers `greg, newton_kl, logit, chebyshev,
  greenkhorn`. `oris_soft` appears only as a one-off fixture at
  `tests/test_parity_weights.py:193`, not in the parametrized R↔Python weight comparison.
- Files: `tests/test_parity_weights.py:73`, `python/leafblower/test_solver_parity.py:193-257`.
- Risk: an ADMM soft-capacity divergence between R and Python ships undetected.
- Priority: **High** — it is the one shipped algorithm outside the parity matrix.

**`raking` and `sinkhorn` are absent from `test_solver_parity.py`:**
- What's not tested: the convergence-rule/`max_error` parity checks (as opposed to raw
  weight-vector parity) skip both.
- Files: `python/leafblower/test_solver_parity.py`.
- Risk: default-rule resolution can diverge per method — exactly what
  `test_logit_default_rule_parity` (line 212) exists to catch, and only logit has it.
- Priority: **Medium**.

**Parity tolerance asymmetry is unexplained:**
- `tests/test_parity_weights.py:93` — `tol = 1e-6 if method == "logit" else 1e-10`. The
  three-orders-of-magnitude relaxation for logit has no comment justifying it.
- Risk: the loose tolerance may be masking a real logit divergence rather than accommodating
  a known conditioning issue.
- Priority: **Medium** — investigate and document, or tighten.

**Benchmark regression gate is opt-in:**
- Stepstone no-regression runs only under `LBW_BENCH_GATE=1`
  (`tests/testthat/test-bench-gate.R`, `.coverage-thresholds.json`).
- Risk: given the no-LTO inlining constraint, perf regressions from TU-boundary changes are
  invisible in the default gate.
- Priority: **Medium** — at minimum, require the gate on any commit touching `src/*.cpp`
  file boundaries or `CMakeLists.txt`.

**No line/branch coverage instrumentation at all:**
- Neither `covr` (R) nor `pytest-cov` (Python) is wired; the quality gate is deliberately
  behavioral (`R CMD INSTALL --preclean .` + 0 testthat FAIL + 0 pytest FAIL).
- Risk: untested branches inside the 2246-line `src/oris.cpp` are undetectable. This is a
  deliberate project choice, not an oversight — noted so it is not mistaken for one.
- Priority: **Low**.

**`.wolf/buglog.json` unavailable:**
- The file is deleted in the working tree (visible as `D` in `git status`), so no logged
  recurring-bug history could be cross-referenced for this audit. Recurring-bug evidence
  above comes from code comments (`src/sraa.hpp:36-42`, `src/raking.cpp:249`) instead.

---

*Concerns audit: 2026-08-14*
