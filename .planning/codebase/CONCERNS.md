<!-- refreshed: 2026-08-15 -->
# Codebase Concerns

**Analysis Date:** 2026-08-15

## Tech Debt

**Legacy Metrics Implementations (Benchmark Suite):**
- Issue: Legacy metric calculation functions in `benchmarks/study/common/metrics.R` and `benchmarks/study/common/test_metrics.py` hardcode `d_i=1` and silently drop divergent terms on starved categories
- Files: `benchmarks/study/common/metrics.R:71`, `benchmarks/study/common/metrics.py:93`, `benchmarks/study/common/test_metrics.R:71-93`
- Impact: Benchmark comparisons underestimate performance gaps on problems with non-uniform design weights or zero-weight categories. Epic-B (leafblower-fmot) introduced correct metric implementations but legacy stubs remain in benchmarks for reference/validation only.
- Fix approach: Keep legacy implementations as-is (documented bugs are intentional for regression testing in test_metrics.R). Production uses correct `compute_metrics()` in `R/harvest.R` and Python via C++ core.

**Manual Version Synchronization:**
- Issue: Package version lives in both `DESCRIPTION` (line 3) and `python/pyproject.toml` (line 10). No automation enforces sync.
- Files: `DESCRIPTION`, `python/pyproject.toml`
- Impact: Version mismatch creates confusion; Python `pip install` would install mismatched version if DESCRIPTION bumped but pyproject.toml not updated on release.
- Fix approach: CLAUDE.md documents this as manual requirement. Add pre-release checklist item or CI check to verify both files have identical `Version:`/`version =` values.

**Solver Formula Review Resistance:**
- Issue: Two verified-correct solver formulas have been questioned/attacked by reviewers 2+ times: (1) Chebyshev Mehrotra corrector's linear `y·Δs_aff` term (actually the `−Δs_aff·Δy_aff` cross-term, not a stray residual); (2) ORIS ALM Newton step `X̃(1−λ+μz)/(1+ρ)` (correct for un-normalized-KL generator, not missing `−ρ`).
- Files: `src/chebyshev.cpp`, `src/oris.cpp`
- Impact: Each reintroduction costs code review cycles; formulas are guarded by comments but reviewers search by formula shape, not derivation.
- Fix approach: Add per-formula ticket link and derivation reference to code comments. Keep CLAUDE.md section 4.7 as canonical "do NOT fix" list; link from solver files.

---

## Known Bugs

**Chebyshev Division by Zero (n_d):**
- Symptoms: If `st.n == 0` (edge case with zero observations), code at `src/chebyshev.cpp:77` (`T_flat[m] = Tgt[m] / n_d`) performs undefined division.
- Files: `src/chebyshev.cpp:56-77`
- Trigger: Call `harvest()` with zero-length input vectors
- Current mitigation: None (unchecked)
- Recommendations: Add guard at line 56: `if (st.n <= 0) { res.status = RK_ERR_BADARG; return res; }` before any division

**Chebyshev Zero-Size Matrix Allocation (nct_red):**
- Symptoms: If all margins have `cat_counts[k] < 2`, then `nct_red = 0`, allocating a 0×0 matrix at `src/chebyshev.cpp:175`. Code later assumes `nct_red > 0` for Hessian factorization.
- Files: `src/chebyshev.cpp:85-175`
- Trigger: All input margins have only one category
- Current mitigation: Overflow check at line 51 (`n_cats_total_with_na > 2048`) does not guard zero-cat-count case
- Recommendations: Add guard after line 85: validate `nct_red > 0` or handle single-category-only case explicitly

**Sinkhorn Redundant Empty Check:**
- Symptoms: `W_best` vector allocated unconditionally at `src/sinkhorn.cpp:84` with `ct.M_cell` elements. The `.empty()` check at line 204 is always false — dead code that will never trigger the fallback.
- Files: `src/sinkhorn.cpp:84, 204`
- Trigger: Best metric finalization path on any valid input
- Current mitigation: Fallback exists but is unreachable
- Recommendations: Remove `&& !W_best.empty()` from line 204 condition; rely only on `std::isfinite(best_metric_seen)`

**Greg Integer Overflow on Large K:**
- Symptoms: Accumulation `n_cats_total += st.cat_counts[k]` at `src/greg.cpp:40` can overflow if total categories exceed INT_MAX (unlikely but not validated).
- Files: `src/greg.cpp:40-57`
- Trigger: K > 64 with very large cat_counts (total > 2.1B)
- Current mitigation: Overflow check at line 42 (`n_cats_total > 2048`), but uses int not size_t; no underflow check
- Recommendations: Change to `size_t n_cats_total = 0;` or add explicit overflow check before accumulation loop

**SRAA Best-Iterate Selection (Re-introduced Bug):**
- Symptoms: At `kErrCheckInterval` iterations, code should use `select_metric(sraa_cfg.metric, cm)` on convergence metric, NOT the fast proxy `errRp`. Using `errRp` misses best iterate on problems with non-L2 convergence metrics.
- Files: `src/oris.cpp` (exact line not specified in CLAUDE.md but flagged as "bug has been re-introduced twice")
- Trigger: Any solve with SRAA acceleration and non-default convergence metric
- Current mitigation: CLAUDE.md documents as known re-introduction pattern
- Recommendations: Add assertion at best-iterate selection point; add test fixture for each convergence metric type

---

## Security Considerations

**C++ Stream I/O in Python Module:**
- Risk: Using `std::ofstream`/`iostream` in Python `_leafblower.so` (which static-links libstdc++) causes SIGSEGV in `std::codecvt do_unshift` when called from Python (leafblower-9nuo).
- Files: `src/oris_trajectory.cpp` (write_trajectory_csv function)
- Current mitigation: Code reviewed; now uses C stdio (`fopen`/`fprintf`/`fclose`); existing instances verified
- Recommendations: Add to CLAUDE.md linting checks; flag any new `#include <iostream>` or `#include <fstream>` in src/ files that are linked into the Python module

**Shadow .so Loading from ~/.local:**
- Risk: Running bare `python -m pytest` (without environment variable setup) can load stale `_leafblower.so` from `~/.local`, bypassing the newly-built module and causing parity test failures or memory corruption.
- Files: `CLAUDE.md` test instructions
- Current mitigation: Test invocation requires explicit environment prefix: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -m pytest`
- Recommendations: Add pre-test cleanup: `rm -rf ~/.local/lib/python*/site-packages/leafblower*` in CI; document in README

**Integer Overflow on Matrix Allocation:**
- Risk: Greg allocator at `src/greg.cpp:57` computes `n_cats_total² doubles` without runtime allocation success check. If `n_cats_total = 2048` (limit), allocation = 32MB; no validation that alloc succeeded.
- Files: `src/greg.cpp:57`
- Current mitigation: Overflow check at line 42 (`n_cats_total > 2048`)
- Recommendations: Add `if (N.capacity() < required_size) { res.status = RK_ERR_BADARG; return res; }` or rely on allocator exceptions with try/catch

---

## Performance Bottlenecks

**Newton-KL Converges to High-Error Fixed Point (K≥5 Severe-Skew):**
- Problem: Newton-KL solver converges to suboptimal fixed point at `gap=6.24e-2` on K≥20 severely-skewed problems (e.g., kk1204 fixture: K=20, target_skew=20, max_weight=3).
- Files: `src/newton_calib.cpp`, `python/leafblower/newton_*.py`
- Cause: Dual landscape basin floor appears intrinsic to Newton-KL formulation on ill-conditioned problems. Not a solver tuning issue.
- Improvement path: Epic-H routing: `harvest(method="auto", ...)` now selects ORIS (method="oris") instead of Newton-KL for K≥5 severe-skew. For K≥5 moderate-skew (max_T/min_T ≤ 5), Newton-KL still selected. Users can pin `method="newton_kl"` to force it.
- Metrics: Stepstone K=9 basin floor ~2.6e-4 (from 2.79e-4 master, 6.5% improvement); T2 failures documented as partial (Epic-Dβ). kk1204 K=20 gap=6.24e-2 vs master divergence.

**Large Source Files (Complexity Risk):**
- Problem: Four solver implementations exceed 900 lines each, making code review and maintenance harder.
- Files: `src/oris.cpp` (2246 lines), `src/r_bridge.cpp` (1274 lines), `src/newton_calib.cpp` (919 lines), `src/raking.cpp` (756 lines)
- Cause: Cell-table representation, multi-solver dispatch, diagnostics payload all in one file per solver family
- Improvement path: `oris.cpp` split into `oris_finalize.cpp` (339 lines) and `oris_trajectory.cpp` (hot loop remains in oris.cpp per CLAUDE.md TU split rules); r_bridge.cpp remains monolithic due to SEXP packing density. No split planned for newton_calib/raking unless new feature forces it.

**ORIS Cell-Compressed Speedup Not Exercised in All Tests:**
- Problem: ORIS achieves up to 1000× speedup on surveys with low tuple diversity by operating at cell-level instead of observation-level. Stepstone (K=9, n=200K) shows speedup but many test fixtures do not exercise this.
- Files: `src/cell_table.cpp`, `R/harvest.R` (cell compression integration)
- Impact: Performance regressions on real datasets (low-diversity K) could be missed if regression suite only uses high-diversity fixtures
- Recommendation: Add T1/T4 fixture with duplicate rows (high compress ratio) to stepstone suite

---

## Fragile Areas

**Cell Table Computation with Edge Cases:**
- Files: `src/cell_table.cpp:51-90`, `src/cell_table.hpp`
- Why fragile: Early return at line 51 (`if (K <= 0 || n <= 0) return -1;`) suppresses error message (returns `-1` but caller cannot diagnose K/n issue). Compare to line 90+ which sets `res.message`. Also: if all cells have `n_per_cell==0`, caller assumption of non-empty table breaks silently.
- Safe modification: Wrap `build_cell_table()` calls in input validation layer (`src/validation.cpp`); check K/n > 0 before calling. Add explicit zero-cell-count guards at caller sites (oris.cpp, chebyshev.cpp).
- Test coverage: `tests/test_cell_table_edge_cases.R` exists (leafblower-hkox.1 commit a4806c8) — validates n_per_cell==0 dispatch

**Weight Recovery Logic (Homotopy Levels):**
- Files: `src/oris.cpp`, `src/calib_dispatch.hpp` (finalize_weights, finalize_weights_buf functions)
- Why fragile: Order-sensitive: scale to n → THEN bounds_mode dispatch. Renormalizing AFTER water-fill silently breaks `bounds_mode="unit"` clamps. Also: struct field comment on `homotopy_levels_used` is wrong (returns 1 for n_levels=1, not 0).
- Safe modification: Add assertion: `if (bounds_mode == kBoundsUnit && renormalized_after_waterfill) { UNREACHABLE(); }`. Update struct comment.
- Test coverage: `test-bounds-order.R` validates order (leafblower-fmot.2 Epic-B)

**R-Bridge SEXP Packing (35-Argument .Call Signature):**
- Files: `src/r_bridge.cpp:1-50` (harvest C wrapper), `R/harvest.R` (roxygen2 entry)
- Why fragile: `.Call("harvest", w, n, K, grouping, targets, bounds, ...)` passes 35 positional args. Adding a single parameter requires:
  1. Update C signature + arg unpacking (r_bridge.cpp)
  2. Update roxygen2 docstring + body (harvest.R)
  3. Add R tests + Python parity test
  4. Update ABI tripwire (type comment on rk_params_t size in types.hpp per leafblower-8aex.3)
- Safe modification: Use feature flags (e.g., `LEAFBLOWER_NEWTON_TRACE` env var) for optional diagnostics; add to rk_result_t as optional fields, unpack only when flag set
- Test coverage: `tests/test-abi-tripwire.R` guards rk_params_t size (leafblower-8aex.3)

---

## Scaling Limits

**Category Count Hard Limit:**
- Current capacity: n_cats_total ≤ 2048 (enforced at `src/greg.cpp:42`, `src/chebyshev.cpp:51`)
- Limit: Beyond 2048 categories, matrix allocations (`N_red`, `N` Hessian) exceed reasonable heap budget (e.g., 2048² = 4M doubles = 32MB per solver invocation)
- Scaling path: Factored solvers (low-rank approximation) or sketching (random projection) for K > 32. Not planned for near term.

**Observation Count:**
- Current capacity: n ≤ 2B (fits in int, but cell_table uses int indexing)
- Limit: Cell table stores per-cell indices; with K=20 dense cells, index explosion on n > 100M
- Scaling path: Streaming or out-of-core cell table; not planned

**BLAS Thread Parallelization:**
- Current limit: Single-threaded BLAS enforced for deterministic tests (`OMP_NUM_THREADS=1`). Production code does not disable BLAS OpenMP by default.
- Limit: On shared-memory systems, uncontrolled BLAS parallelism + package's own OpenMP can oversubscribe
- Scaling path: Document recommendation for BLAS thread setting per system; consider thread-safe thread pool wrapper (not planned)

---

## Dependencies at Risk

**Stale `.venv` or Mismatched numpy/scikit-build:**
- Risk: Python build relies on scikit-build-core + pybind11 + CMake. If `.venv` has stale wheels or mismatched ABI (numpy < 1.21 on Python 3.12), build silently produces wrong `.so`
- Files: `python/pyproject.toml`, `python/CMakeLists.txt`, `python/.venv/`
- Impact: Import `leafblower` crashes or loads stale C++ symbols
- Migration plan: Add `uv pip install -e . --reinstall-package leafblower` (with `--reinstall-package` flag) to setup instructions; document `.venv/` is disposable
- Current mitigation: CLAUDE.md explicitly forbids bare `pip install`; requires `uv pip install -e .`

**OpenMP Optional (System Compiler Dependent):**
- Risk: `SystemRequirements` in DESCRIPTION lists "OpenMP (optional, for SIMD exp vectorisation)" but code compiles without it. If system compiler lacks OpenMP, exp() loops unvectorized, performance degrades by ~10-30% on large n.
- Files: `src/lbw_config.h`, configure script
- Impact: Silent performance regression on systems without OpenMP
- Migration plan: Add build summary to `R CMD INSTALL` output naming detected OpenMP flag; warn if absent

---

## Missing Critical Features

**Continuation Methods for Newton-KL Basin Escape:**
- Problem: Newton-KL converges to basin floor (gap=6.24e-2) on K≥5 severe-skew; no cold-start multi-start or continuation homotopy planned yet
- Blocks: Tier-1 quality improvement for ill-conditioned problems
- Status: Deferred to Epic-E (leafblower-zpkd closure)

**Chambolle-Pock and IPM Solvers Ruled Out:**
- Problem: Epic-J evaluated CP (Chambolle-Pock primal-dual) and IPM (Interior Point, Schur-complement). Both ruled out: CP kk1204=0.308, IPM=0.483 vs ieppa+sraa=0.0501 baseline. Research code deleted per policy.
- Blocks: No valid use case identified; Epic-K cancelled
- Status: Closed (leafblower-pcs9 Epic-K cancelled)

---

## Test Coverage Gaps

**Untested Optimization Paths (No -O Flag Set):**
- What's not tested: R package honors user's `-O` flag (e.g., from `~/.R/Makevars`); Python wheel build uses `cmake.build-type = "Release"` (implies `-O3`). If user has broken `-O` setting (e.g., `-O3 -ffast-math` on x86), behavior differs from test CI.
- Files: `configure`, `Makevars.in`, `python/pyproject.toml`
- Risk: Floating-point instability on user systems; rounding differences in tests not caught locally
- Priority: Medium — rare but high-impact if reported

**Newton-KL Suboptimal Fixed Points (Partial Acceptance):**
- What's not tested: How to detect/escape basin floor on new fixtures; current tests only validate convergence status, not quality ceiling
- Files: `tests/test-newton-kl-convergence.R`, `python/leafblower/test_newton_kl_parity.py`
- Risk: User introduces new problem type with same K/skew signature as kk1204; solver reports status=0 but quality is poor (gap > 1e-2)
- Priority: High — Epic-E will add Tier-1 metrics (margin_kl, DEFF) to results, enabling user to diagnose after-the-fact

**Cell Compression Not Exercised in Standard Stepstone:**
- What's not tested: ORIS cell-compression speedup on production-scale low-diversity data
- Files: `benchmarks/stepstone*.R`, `python/benchmarks/stepstone*.py`
- Risk: Performance regression on real datasets (e.g., survey with many duplicate responses) missed in CI
- Priority: Low — Epic-K (CP productionization) was cancelled, so no new solver to regress; stepstone coverage sufficient for current solvers

**Design Weights (d_i != 1) Edge Cases:**
- What's not tested: Extreme design weights (d_i ranging 0.001 to 1000) or all-zero design weights
- Files: `tests/test-design-weights.R`, `python/leafblower/test_design_weights_parity.py`
- Risk: Weight recovery with extreme d_i breaks silently; metric computation diverges from formula
- Priority: Medium — Epic-52jc added margin_kl and weight_kl to results (leafblower-52jc), but edge-case coverage not expanded

**Metric Divergence on Starved Categories (Test-Only):**
- What's not tested: User code calling legacy `stepstone_all_methods.R` metrics (for benchmark comparison) on data with starved categories; new `compute_metrics()` correctly returns Inf, but stale metrics return finite value
- Files: `benchmarks/study/common/test_metrics.R:130-152` (documents bug reproduction, doesn't fix legacy)
- Risk: Benchmark comparisons on user data with zero-weight cells silently underestimate error
- Priority: Low — only affects benchmark users; production harvest() uses correct C++ metrics via compute_metrics()

---

*Concerns audit: 2026-08-15*
