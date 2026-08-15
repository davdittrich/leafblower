---
phase: 01-verification-coverage-closed
verified: 2026-08-15T00:00:00Z
status: passed
score: 7/7 must-haves verified
behavior_unverified: 0
overrides_applied: 0
---

# Phase 1: Verification Coverage Closed — Verification Report

**Phase Goal:** The Definition-of-Done gate actually detects the two failure classes it
claims to cover — R↔Python numerical divergence, and weight-bound violations — before any
code that could cause them is touched.
**Verified:** 2026-08-15
**Status:** passed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | All nine non-AUTO `rk_algorithm_t` solvers are R-vs-Python weight-vector compared inside the blocking gate (SC1) | VERIFIED | `python/leafblower/test_parity_weights.py:82-85` parametrize list `["greenkhorn","logit","raking","oris","sinkhorn","newton_kl","chebyshev","greg","oris_soft"]` — exact match against `src/leafblower.h:42-52` enum (9 non-AUTO values, slots 2/7 correctly absent). Live run: 156 Python tests collected, 156 passed, 0 failed. |
| 2 | `raking`/`sinkhorn` covered by convergence-rule / `max_error` parity checks (SC2) | VERIFIED | `python/leafblower/test_solver_parity.py:269-330` — `test_raking_parity`, `test_sinkhorn_parity`, `test_raking_default_rule_parity`, `test_sinkhorn_default_rule_parity` all present and passing (part of the 156/156 green run). Each default-rule test asserts the resolved rule string on both bindings separately before weight parity. |
| 3 | No unexplained parity tolerance survives; `logit`'s special case is gone or documented (SC3) | VERIFIED | `python/leafblower/test_parity_weights.py:103-113` — single `tol = 1e-10` for all methods, `method == "logit"` conditional absent (`grep -c` = 0), adjacent comment records the measured 5.33e-15 value, date, cross-solver range, and the two mechanisms considered and ruled out. |
| 4 | Property-based test asserts `min(w)>=min_weight` and `max(w)<=max_weight` within 1e-10 across 50 datasets (KPI-02 / SC4) | VERIFIED | `tests/testthat/test-bound-property.R:106-207` — 50 `set.seed()` literals (17/16/17 stratified skewed-weights/sparse-cells/both), `bounds_mode="unit"`, unconditional assertion (no convergence precheck), non-vacuity witness `expect_gte(n_engaged, 40L)`. Live re-run: `[ FAIL 0 \| WARN 0 \| SKIP 0 \| PASS 808 ]`. **Live mutation check**: switching the tracer's `bounds_mode` from `"unit"` to `"cell"` reproducibly FAILS the test (`min(w)=0.0045 < 0.2`, `max(w)=128.3 > 5`, `n_bounds_clamped=0`) — proves the assertion is load-bearing, not tautological. File restored clean after the check. |
| 5 | `newton_kl`'s conflict with KPI-02's literal wording is adjudicated, pinned by an assertion, and ticketed, not silently dropped | VERIFIED | `tests/testthat/test-bound-property.R:231-255` pins the report-not-clamp contract (`status != RK_OK` AND `n_bounds_violated > 0` whenever bounds are violated). Decision (option-a) recorded in `01-03-SUMMARY.md` and as a `bd comment` on `leafblower-og7d.5` (verified live via `bd show`), which remains OPEN and correctly scheduled no earlier than Phase 2. |
| 6 | `testthat::test_dir("tests/testthat")` runs under the edition `CLAUDE.md` documents, with `DESCRIPTION` in agreement (SC5) | VERIFIED | `DESCRIPTION:22` `Config/testthat/edition: 3`. Live re-run: `EDITION= 3 FAILED= 0 ERROR= 0 PASSED= 1833 WARNING= 141 SKIPPED= 13` — matches `01-04-SUMMARY.md`'s recorded numbers exactly. `context()` deprecation cleanup (11 files) and the flip are separate commits (`da7bff9`, `d6c2c5c`, `dc57317`), each independently green. |
| 7 | The gate-collection gap is closed: the previously-orphaned parity file is now inside the blocking Python command (`leafblower-x7n8`) | VERIFIED | `python/leafblower/test_parity_weights.py` exists at the gate's collection root; `tests/test_parity_weights.py` (source path) confirmed gone. **Live mutation check**: temporarily changing `parents[2]` to `parents[1]` makes pytest collection ERROR (`AssertionError: R_HELPER not found ... REPO_ROOT arithmetic is wrong`), not a silent skip — proves the loud-failure property. File restored clean after the check. |

**Score:** 7/7 truths verified (0 present-but-behavior-unverified)

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `python/leafblower/test_parity_weights.py` | relocated, 9-solver parametrize, uniform tolerance | VERIFIED | Exists, wired into pytest collection root, content matches all claims |
| `tests/test_parity_weights.py` | no longer present | VERIFIED | Confirmed absent (`git mv` rename, `git diff -M` shows rename) |
| `python/leafblower/test_solver_parity.py` | raking/sinkhorn parity + default-rule lock | VERIFIED | 4 new test functions present, delegate to existing `_assert_parity` helper |
| `tests/testthat/test-bound-property.R` | 50-dataset property sweep + newton_kl exception | VERIFIED | 256 lines, all structural claims confirmed by grep + live run + live mutation |
| `DESCRIPTION` | `Config/testthat/edition: 3` | VERIFIED | Field present, live run confirms `EDITION=3` |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|----|--------|---------|
| `python/leafblower/test_parity_weights.py` `REPO_ROOT` | `<repo_root>/tests/parity/run_*_r.R` | `Path(__file__).resolve().parents[2]` | WIRED | Import-time assertion confirms resolution; live mutation proves loud failure on wrong arithmetic |
| pytest rootdir `python/` | `.coverage-thresholds.json` `enforcement.command` | Python collection | WIRED | 156 tests collected and run by the exact gate command |
| `harvest(bounds_mode="unit")` | `lbw::finalize_weights_buf` water-fill | returned weight vector | WIRED | Live mutation to `bounds_mode="cell"` demonstrates the enforced/unenforced code-path distinction the test depends on |
| `attr(res,"result")$n_bounds_clamped` | non-vacuity witness | test assertion | WIRED | `expect_gte(n_engaged, 40L)` present and reads the real clamp counter |
| `DESCRIPTION Config/testthat/edition` | `testthat::edition_get()` | live test run | WIRED | Confirmed `EDITION= 3` at runtime |

### Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| R build gate | `R CMD INSTALL --preclean .` | `DONE (leafblower)` | PASS |
| Full R suite, edition 3 | `Rscript -e test_dir(...)` | `EDITION=3 FAILED=0 ERROR=0 PASSED=1833 WARNING=141 SKIPPED=13` | PASS |
| Python collection count | `pytest --collect-only -q` | `156 tests collected` | PASS |
| Full Python suite | `pytest -q` | `156 passed, 0 failed` | PASS |
| REPO_ROOT wrong-arithmetic loud failure | `parents[2]`→`parents[1]`, `pytest --collect-only` | `AssertionError ... REPO_ROOT arithmetic is wrong`, collection ERROR (not skip) | PASS |
| Bound-enforcement mode-dependence | `bounds_mode="unit"`→`"cell"` in tracer, `test_file(...)` | 3 failures: min/max both violated, `n_bounds_clamped=0` | PASS |

Both mutation checks confirmed the working tree was restored byte-identical to `HEAD`
afterward (`git diff --stat` empty).

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|-------------|-------------|--------|----------|
| KPI-02 | 01-03 | Weight-bound enforcement — property test over 50 datasets | SATISFIED | `tests/testthat/test-bound-property.R`, live run green, non-vacuity witness present |
| leafblower-x7n8 | 01-01 | Gate-collection gap for the weight-parity file | SATISFIED (ticket still OPEN — see gaps below) | 156 tests collected inside the blocking gate |
| SC1–SC5 (ROADMAP) | 01-01..04 | See Observable Truths #1–4, #6 above | SATISFIED | — |

No orphaned requirements: `REQUIREMENTS.md`'s Traceability table maps only `KPI-02` to
Phase 1, and it is satisfied. `leafblower-og7d` (epic) and its `.1`–`.5` sub-tickets are the
plan-time-filed beads work items for SC1/SC2/SC3/KPI-02 and the newton_kl tension — all
implemented (see Anti-Patterns/Info note on ticket-closure hygiene below).

### Anti-Patterns Found

No `TBD`/`FIXME`/`XXX`/`TODO`/`HACK`/`PLACEHOLDER` markers found in any of the 17 files this
phase touched.

### Info / Non-Blocking Observations

1. **Ticket hygiene:** `leafblower-x7n8`, `leafblower-og7d`, and `leafblower-og7d.1`–`.4` are
   still `OPEN` in beads despite `01-01-SUMMARY.md`/`01-04-SUMMARY.md` describing them as
   "closable." The underlying work is verified complete and green; only the `bd close` step
   was not run. Does not block the phase's technical goal (the gate demonstrably detects both
   failure classes) but is a process-completeness gap against this project's session-close
   protocol ("Close finished work"). `leafblower-og7d.5` is correctly left OPEN (deliberately
   deferred to Phase 2).
2. **REQUIREMENTS.md staleness:** KPI-02's entry still reads "Not located" / `[ ]` in
   `.planning/REQUIREMENTS.md`, predating this phase's work. The document states it derives
   from planning-time extraction; refreshing it is normally a milestone/docs-update step, not
   a per-phase one. Flagged for the next docs-update pass, not a phase-goal blocker.
3. **Stale SUMMARY claim:** `01-01-SUMMARY.md` and `01-02-SUMMARY.md` both describe a
   "pre-existing, out-of-scope" failure in `python/leafblower/test_trajectory_csv_smoke.py`
   (CSV header mismatch). Live re-run shows this test passing (`1 passed`) and the actual
   `src/oris_trajectory.cpp` header (`"iter,errRp\n"`) matches the test's assertion exactly —
   the claimed divergence does not reproduce on the current tree. This is a narrative
   inaccuracy in two SUMMARYs (over-cautious, not over-claiming), with zero effect on the
   phase's actual gate health: the full Python suite is 156/156 green either way.
4. **Code review (`01-REVIEW.md`):** 0 critical, 3 warning, 2 info findings. WR-01 (blanket
   `suppressWarnings()` muting all warnings, not just the documented incidental one) is fixed
   and committed as `5a25671` — confirmed live (file diff inspected, full R suite re-run
   green with the fix in place). WR-02 (subprocess-doubling in two new tests, copies a
   pre-existing pattern) and WR-03 (no per-test timeout guard on the 400-call sweep) are
   correctly deferred as advisory per the task brief — neither affects correctness or gate
   detection capability.

### Human Verification Required

None. All truths for this phase are deterministic, test-layer assertions with no UI, visual,
or real-time-behavior component; all were verified by live command execution plus two
independent mutation-based load-bearing checks.

### Gaps Summary

No gaps. All 5 ROADMAP success criteria, KPI-02, and the plan-level must-haves (loud-failure
property, non-vacuity witness, newton_kl adjudication, fix-first-then-flip bisectability) are
verified against the live codebase, not just SUMMARY claims. Both failure classes the gate is
supposed to catch — R↔Python numerical divergence (via the 9-solver uniform-tolerance parity
matrix, now actually collected) and weight-bound violations (via the 50-dataset unconditional
property sweep) — were independently reproduced as catchable by deliberately breaking each
mechanism and observing a loud failure, then restoring the tree clean.

---

_Verified: 2026-08-15_
_Verifier: Claude (gsd-verifier)_
