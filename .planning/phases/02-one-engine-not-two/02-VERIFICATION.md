---
phase: 02-one-engine-not-two
verified: 2026-08-15T10:15:59Z
status: passed
score: 5/5 must-haves verified
behavior_unverified: 0
overrides_applied: 0
---

# Phase 2: One Engine, Not Two — Verification Report

**Phase Goal:** An R user and a Python user get the same number because they ran the same
code, not because two independently maintained code paths happen to agree today.
**Verified:** 2026-08-15T10:15:59Z
**Status:** passed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths (ROADMAP Success Criteria, SC1-SC5)

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | Adding a solver/result field/routing rule requires editing ONE dispatch site; R bridge reaches solvers through the same path as the C ABI, not its own method-string chain | ✓ VERIFIED | `src/calib_dispatch.hpp::dispatch_solver()` has all 9 non-AUTO `case RK_ALG_*` arms (lines 753-1121); `src/r_bridge.cpp` has **zero** `strcmp(method_str, ...)` occurrences in live code (the sole hit at line 351 is `strcmp(...,"newton_kl")` for an unrelated `outer_max_iter` default, not solver dispatch) and calls `lbw::dispatch_solver()` at exactly 3 sites (AUTO primary L621, AUTO newton_kl fallback L678, unified explicit-method L739); `src/c_api.cpp` calls `lbw::dispatch_solver()` at all 10 of its algorithm branches. Regression-locked by `python/leafblower/test_single_dispatch_site.py::test_r_bridge_has_no_method_dispatch_chain` — **ran independently, PASSED**. |
| 2 | R and Python binaries built at the same optimization level, or the asymmetry is documented at the parity assertion as deliberate/bounded | ✓ VERIFIED | `CLAUDE.md:66` documents R sets no `-O` (CRAN `tools:::.check_make_vars` constraint) while `python/CMakeLists.txt:99` hard-sets `-O3`; both `python/CMakeLists.txt:105` and `python/leafblower/test_solver_parity.py` header cross-reference the decision. `leafblower-qzto` (P1) closed with this evidence. |
| 3 | Per-cell unit-mode water-fill exists once; a fix to `bounds_mode="unit"` cannot land in one solver family and miss the other | ✓ VERIFIED | `finalize_weights_buf` is defined exactly once, in `src/calib_dispatch.hpp:402`. `src/oris_finalize.cpp` explicitly documents the prior duplication was removed ("were duplicated here and in calib_dispatch.hpp::finalize_weights_buf. Delegate to the single source", L153) and calls the shared helper (L135, L161). All 7 bounded solvers (`oris_finalize.cpp`, `raking.cpp`, `chebyshev.cpp`, `greenkhorn.cpp`, `greg.cpp`, `logit_calib.cpp`, `sinkhorn.cpp`) delegate; `newton_kl`'s pre-existing, honestly-documented exclusion (no box-constrained inner step) is tracked on `leafblower-og7d.5`, not silently omitted. Regression-locked by `python/leafblower/test_finalize_weights_sync.py` (2 tests) — **ran independently, PASSED**. |
| 4 | A build-list divergence between `src/*.cpp` and `python/CMakeLists.txt:CORE_SOURCES` fails the test gate instead of surfacing as a link error | ✓ VERIFIED | `python/leafblower/test_core_sources_sync.py::test_core_sources_matches_src_glob` diffs the live `src/*.cpp` glob (minus R-only `r_bridge.cpp`) against `CORE_SOURCES` — **ran independently, PASSED** (17 == 17, zero drift today). |
| 5 | The full DoD gate is green after the change and the stepstone benchmark shows no regression | ✓ VERIFIED | Independently re-ran the full DoD gate from clean: `R CMD INSTALL --preclean .` succeeded; `Rscript devtools::test()` → **FAIL 0, WARN 141, SKIP 11, PASS 1839**; `uv pip install -e . --reinstall-package leafblower` succeeded; `pytest -q` (single-thread BLAS) → **228 passed, 0 failed**. Independently re-ran the stepstone gate (`LBW_BENCH_GATE=1 NOT_CRAN=true`, `filter="bench-gate"`): `kk1204 gate: status=0 iters=10 best_error=-7.376e-14 time=1.5s`, **FAIL 0, PASS 3, SKIP 2** — numbers match 02-08-SUMMARY.md's recorded baseline exactly (byte-identical `best_error`/`time`). |

**Score:** 5/5 truths verified (0 present, behavior-unverified)

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `src/calib_dispatch.hpp::lbw::dispatch_solver()` | Shared enum→solver dispatch table, all 9 non-AUTO algorithms | ✓ VERIFIED | 9 case arms confirmed (SINKHORN, GREG, GREENKHORN, LOGIT, CHEBYSHEV, RAKING, ORIS, ORIS_SOFT, NEWTON_KL) |
| `src/calib_dispatch.hpp::lbw::route_auto()` | Single AUTO routing decision, shared by both bridges | ✓ VERIFIED | Present (L1173); thresholds unchanged (0.9 compression, K≥5, skew>5.0 — verified by inspection against ROADMAP's stated formulas) |
| `src/calib_dispatch.hpp::lbw::kAlgNames` | Single enum-to-name table, holes at slots 2/7, static_assert guard | ✓ VERIFIED | 12-entry positional table (L59-73) with empty strings at slots 2/7; two `static_assert`s tie it to `RK_ALG_NEWTON_KL==11` and array length 12 |
| `python/leafblower/test_core_sources_sync.py` | SC4 regression guard | ✓ VERIFIED | Exists, PASSED under independent execution |
| `python/leafblower/test_finalize_weights_sync.py` | SC3 regression guard | ✓ VERIFIED | Exists, both tests PASSED under independent execution |
| `python/leafblower/test_single_dispatch_site.py` | SC1 regression guard | ✓ VERIFIED | Exists, PASSED under independent execution |
| `sizeof(rk_params_t)==264`, `sizeof(rk_result_t)==536` | ABI static_asserts untouched | ✓ VERIFIED | `src/leafblower.h:246,270` static_asserts present; R build (which compiles this header) succeeded, so both asserts held |

### Key Link Verification

| From | To | Via | Status | Details |
|------|-----|-----|--------|---------|
| `r_bridge.cpp` (AUTO branch) | `lbw::route_auto` → `lbw::dispatch_solver` | direct call | ✓ WIRED | L613 `route_auto(...)`, L621 `dispatch_solver(route.algorithm, ...)` |
| `r_bridge.cpp` (explicit-method branch) | `lbw::dispatch_solver` | direct call | ✓ WIRED | L739 single unified call, `dispatch_solver(alg, st, dres)` |
| `c_api.cpp` (all 9 explicit methods + AUTO) | `lbw::dispatch_solver` | direct call | ✓ WIRED | 10 call sites confirmed, each followed by `pack_dispatch_result_c`/`pack_dispatch_oris_extras_c` |
| `c_api.cpp` (AUTO resolution switch) | `lbw::route_auto` | direct call | ✓ WIRED | L268 |
| `test_single_dispatch_site.py` | `src/r_bridge.cpp` | text scan | ✓ WIRED | Reads file, asserts 0 strcmp(method_str) + ≤3 dispatch_solver call sites |
| `test_core_sources_sync.py` | `src/*.cpp` glob + `CORE_SOURCES` | text scan | ✓ WIRED | Set-equality assertion |
| `test_finalize_weights_sync.py` | `src/calib_dispatch.hpp::finalize_weights_buf` + 7 solver files | text scan | ✓ WIRED | Both tests exercised |
| `.coverage-thresholds.json::enforcement.command` | pytest | DoD gate wiring | ✓ WIRED | Command includes `pytest -q`, which collects the 3 new guard tests automatically (confirmed: no `-k`/path restriction in the enforcement command) |

### Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| SC1 guard passes | `pytest test_single_dispatch_site.py` | PASSED | ✓ PASS |
| SC3 guard passes | `pytest test_finalize_weights_sync.py` | 2 PASSED | ✓ PASS |
| SC4 guard passes | `pytest test_core_sources_sync.py` | PASSED | ✓ PASS |
| R build succeeds | `R CMD INSTALL --preclean .` | DONE (leafblower) | ✓ PASS |
| R full test suite green | `Rscript devtools::test()` | FAIL 0, PASS 1839 | ✓ PASS |
| Python full test suite green | `pytest -q` (single-thread BLAS) | 228 passed, 0 failed | ✓ PASS |
| Stepstone benchmark unregressed | `LBW_BENCH_GATE=1 filter="bench-gate"` | status=0 iters=10 best_error=-7.376e-14 time=1.5s, FAIL 0 | ✓ PASS |

### Probe Execution

No `scripts/*/tests/probe-*.sh` convention used by this project; the phase's own regression guards (SC1/SC3/SC4 pytest modules) and the `.coverage-thresholds.json` DoD command serve this role and were executed directly above.

### Requirements Coverage

| Requirement | Source Plans | Description | Status | Evidence |
|-------------|--------------|--------------|--------|----------|
| US-004 | 02-01 through 02-08 (all 8 plans) | Stable C API; both R `.Call()` bridge and Python pybind11 module route through the same shared dispatch table | ✓ SATISFIED | `REQUIREMENTS.md` already records this as Implemented citing Phase 2 plans 01-08 and `test_single_dispatch_site.py`; independently confirmed above (dispatch_solver/route_auto single-site, guard tests green, zero strcmp method chain) |

No orphaned requirements: `REQUIREMENTS.md`'s Phase 2 traceability table lists only US-004, `leafblower-rywn`, `leafblower-qzto`, and the water-fill/CORE_SOURCES items — all covered by the 8 plans and verified above.

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| — | — | No `TBD`/`FIXME`/`XXX`/`TODO`/`HACK`/`PLACEHOLDER` markers found in any phase-modified file (`calib_dispatch.hpp`, `c_api.cpp`, `r_bridge.cpp`, the 3 new test files) | — | None |
| `src/calib_dispatch.hpp:1120-1122` | `dispatch_solver`'s `default:` case | ⚠️ Warning (from 02-REVIEW.md, WR-02/WR-04) | Robustness gap, not a live bug | An unreachable-today `default:` leaves `out` at its default-constructed values (`status=RK_ERR_NOCONV`, `alg_used=RK_ALG_ORIS`), indistinguishable from "ORIS ran and failed" if `dispatch_solver` were ever called with `RK_ALG_AUTO` or a reserved slot (2/7). Both current callers pre-resolve AUTO before calling, so this is dead code today; flagged by the phase's own code review as a hardening opportunity, not a correctness defect present in the codebase. |
| `src/calib_dispatch.hpp:751-1121` | 9 near-identical ~25-line field-copy blocks in `dispatch_solver`'s case arms | ℹ️ Info (from 02-REVIEW.md, WR-01) | Maintainability | The migration moved solver-*selection* duplication (was: 2 bridges × 9 solvers) down to 1 site, but the field-*copy* boilerplate is still repeated 9× within that one site. Does not violate SC1's literal wording ("editing ONE dispatch site" — case arms are all in the one file/function) but is a legitimate future-refactor candidate the phase's own reviewer flagged. |
| `src/calib_dispatch.hpp` (7 of 9 case arms) | best_weights move without empty-guard | ℹ️ Info (from 02-REVIEW.md, WR-03) | Latent robustness gap | 7 of 9 case arms move `best_weights` unconditionally, relying on a comment-only claim ("never empty on any path") rather than a machine-checked invariant. Pre-existing pattern, not introduced by this phase's migration (the same asymmetry existed pre-migration); noted by the phase's own reviewer as a natural next hardening step. |

These three review findings (WR-01 through WR-04, one collapsed above as it's an "info" pairing) are all severity **Warning/Info** in the phase's own `02-REVIEW.md` (0 Critical), do not contradict any of the 5 success criteria, and do not block phase completion — recorded here for traceability, not as gaps.

### Human Verification Required

None. All 5 success criteria have machine-checked evidence (regression-guard tests independently re-run, full DoD gate independently re-run from a clean R build, stepstone benchmark independently re-run with numbers matching the recorded baseline). SC5's human-sign-off step (02-08-SUMMARY.md's Task 3, "D4" — user confirms unchanged user-visible surface) was already obtained by the user on 2026-08-15 per the SUMMARY's own recorded approval, and is corroborated here by the independently-passing full test suites (R: 0 FAIL/1839 PASS; Python: 228 passed/0 failed) which exercise exactly the R-only-field-presence and Python-result-dict-shape assertions that sign-off covered.

### Gaps Summary

None. All 5 ROADMAP success criteria for Phase 2 are independently verified against the live codebase (not SUMMARY.md narration): the shared dispatch table and AUTO-routing function exist and are the sole path both FFI bridges use; the R bridge's per-method `strcmp` chain is gone and enforced by a regression test; the `-O` asymmetry is documented; the unit-mode water-fill has one implementation with the `newton_kl` gap honestly tracked (not silently dropped); the `CORE_SOURCES`/`src/*.cpp` sync is a regression-gated test; and the full DoD gate plus stepstone benchmark are green with numbers matching the phase's own recorded baseline. Both beads tickets this phase targets (`leafblower-rywn` P0, `leafblower-qzto` P1) are closed with matching evidence. Three non-blocking maintainability findings from the phase's own code review are carried forward for future hardening, not as phase gaps.

---

_Verified: 2026-08-15T10:15:59Z_
_Verifier: Claude (gsd-verifier)_
