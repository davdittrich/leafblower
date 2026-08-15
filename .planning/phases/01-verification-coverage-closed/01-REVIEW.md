---
phase: 01-verification-coverage-closed
reviewed: 2026-08-15T00:00:34Z
depth: standard
files_reviewed: 16
files_reviewed_list:
  - python/leafblower/test_parity_weights.py
  - python/leafblower/test_solver_parity.py
  - tests/testthat/test-bound-property.R
  - tests/testthat/test-alm-config-grouping.R
  - tests/testthat/test-cell-table.R
  - tests/testthat/test-compare.R
  - tests/testthat/test-oris-b12-fallback-best-reset.R
  - tests/testthat/test-oris-b13-best-error-honesty.R
  - tests/testthat/test-oris-dispatch.R
  - tests/testthat/test-oris-faithful.R
  - tests/testthat/test-oris-nonuniform-d.R
  - tests/testthat/test-oris.R
  - tests/testthat/test-oris-sraa-log-path.R
  - tests/testthat/test-oris-sraa.R
  - tests/testthat/test-autocollapse.R
  - tests/testthat/test-harvest-rval.R
findings:
  critical: 0
  warning: 3
  info: 2
  total: 5
status: issues_found
---

# Phase 01: Code Review Report

**Reviewed:** 2026-08-15T00:00:34Z
**Depth:** standard
**Files Reviewed:** 16 (+ `DESCRIPTION`, metadata-only)
**Status:** issues_found

## Summary

Phase 01 closes four verification gaps: weight-vector parity coverage for all
nine non-AUTO solvers (`python/leafblower/test_parity_weights.py`,
`test_solver_parity.py`), a KPI-02 bound-invariant property test across 50
stratified datasets and 8 solvers plus a documented `newton_kl` exception
(`test-bound-property.R`), and a testthat edition-2→3 flip (11 mechanical
`context()` removals + a `Config/testthat/edition: 3` DESCRIPTION field + a
genuine `expect_warning()` scoping fix in two files).

This is unusually well-evidenced test work: tolerances are justified with
measured numbers rather than guessed, the `newton_kl` "report not clamp"
exception is backed by real result-struct fields (`status`,
`n_bounds_violated`, `n_bounds_clamped` — verified against `r_bridge.cpp` and
sibling test files, not invented), the raking/sinkhorn default-rule lock was
verified with a load-bearing mutation test per 01-02-SUMMARY.md, and the
9-solver parametrize list in `test_parity_weights.py` matches
`rk_algorithm_t` in `leafblower.h` exactly (verified against source: AUTO=0,
ORIS=1, [2 removed], RAKING=3, SINKHORN=4, CHEBYSHEV=5, GREG=6, [7 reserved],
ORIS_SOFT=8, GREENKHORN=9, LOGIT=10, NEWTON_KL=11 → 9 non-AUTO entries).

No blocking correctness or security defects found. The `expect_warning()` /
`suppressWarnings()` fix (the specific pattern flagged for scrutiny by this
review's brief) is directionally correct and does not break the assertion it
wraps — verified by tracing testthat 3e's calling-handler order and cross-
checked against the plan's own documented probe methodology — but it is
broader than it needs to be, and that breadth has a real (if modest) cost:
see WR-01. Two smaller robustness/consistency gaps are also flagged.

## Warnings

### WR-01: `suppressWarnings()` wrapper mutes ALL non-target warnings, not just the documented incidental one

**File:** `tests/testthat/test-autocollapse.R:26-34, 55-63, 65-72`
**File:** `tests/testthat/test-harvest-rval.R:72-80`
**Issue:** The fix for testthat 3e's non-swallowing `expect_warning()` wraps
each call as `suppressWarnings(expect_warning(harvest(...), "<regexp>"))`.
This is correctly ordered (outer, not inner — verified this does not
consume the target warning before `expect_warning()`'s handler observes it,
matching the documented probe in 01-04-SUMMARY.md) and does not weaken the
positive assertion: `expect_warning()` still fails if `"became NA"` /
`"out.of.vocabulary|OOV|..."` never fires.

However, the outer `suppressWarnings()` silences *every* warning that
`harvest()` raises during the call, not just the one documented incidental
warning (sparse-category or fixed-point-convergence diagnostic) that
motivated the fix. Under testthat edition 3, an unexpected/unmatched warning
propagating out of a test is itself a signal — that is the entire behavior
change SC5/D-11 was investigated to explain. By blanket-suppressing, these
four call sites opt back out of that signal for any *future* new warning
`harvest()` might raise here (e.g. a genuinely new numerical-instability or
deprecation warning introduced by an unrelated later change) — it will now
pass silently instead of surfacing as edition 3's stricter warning-count
increase, which is exactly the mechanism that caught this same class of bug
(D-11) during this phase.
**Fix:** Narrow the suppression to the two specific known-incidental message
patterns instead of `suppressWarnings()`'s catch-all, e.g.:
```r
withCallingHandlers(
  expect_warning(harvest(...), "became NA"),
  warning = function(w) {
    if (grepl("sparse categor|fixed.point converg", conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  }
)
```
This preserves edition 3's protection against any *other* unexpected warning
at these four call sites while still absorbing the two documented incidental
ones. If the narrower fix is judged not worth the extra code for four call
sites, downgrade this to Info and record it as an accepted trade-off in
CLAUDE.md or the test file's own comment (the current comment explains why
suppression is safe for the *known* warning, but not that it also silences
*unknown* ones).

### WR-02: Newly-added `raking`/`sinkhorn` default-rule tests invoke `_run_r`/`_run_py` twice per test, doubling subprocess cost for no added coverage

**File:** `python/leafblower/test_solver_parity.py:284-330`
**Issue:** `test_raking_default_rule_parity` and `test_sinkhorn_default_rule_parity`
each call `_run_py(method, conv_py={})` once directly (to check the resolved
rule string) and then call `_assert_parity(method, conv_py={}, conv_r="list()")`,
which internally calls `_run_py` and `_run_r` again with the identical
arguments — i.e. Python `harvest()` runs twice and `Rscript` is subprocessed
twice per test, quadrupling the work for what is a deterministic,
side-effect-free computation (embedded fixture, fixed convergence spec).
This exactly mirrors the pre-existing `test_logit_default_rule_parity`
pattern (not introduced by this diff), so it is a faithful, intentional
copy — but copying a pre-existing inefficiency into two more tests compounds
it rather than fixing it. Out of `v1` performance scope per review rules, so
not escalated to Warning-blocking, but flagged as it is copy-paste growth of
an existing quality gap.
**Fix:** Cache `_run_py`/`_run_r` results and reuse them for both the
rule-string assertion and the weight-parity assertion, e.g. change
`_assert_parity` to optionally accept pre-computed `(w_py, ri_py)` /
`r_out` rather than re-invoking. Low priority; note only if touching this
file again.

### WR-03: `test-bound-property.R`'s 400-call sweep and `test_parity_weights.py`'s 3 subprocess-heavy tests have no explicit per-test timeout guard against R-side hangs

**File:** `tests/testthat/test-bound-property.R:174-207`
**Issue:** The 50-dataset × 8-solver sweep (`harvest(..., max_iterations = 1000L)`,
400 calls) has no iteration/time ceiling other than `max_iterations`; if a
future core regression causes one solver to spin near its iteration cap
across many of the 400 calls, this single `test_that()` becomes a slow,
hard-to-diagnose CI drag with no per-call bound to distinguish "one solver
is slow" from "the sweep is slow." This is a robustness/maintainability
observation, not a defect in the current code — `max_iterations` does bound
each individual call — and no evidence found that it currently misbehaves.
**Fix:** No action required now; if this test becomes a recurring source of
slow/flaky CI, consider splitting per-stratum or adding
`testthat::skip_on_ci()`-gated timing diagnostics. Downgraded to Info-
adjacent; kept as Warning only because a 400-call unconditional sweep in one
`test_that()` block is a maintainability smell worth a ticket, not because
it is presently broken.

## Info

### IN-01: `test_parity_weights.py` duplicates the single-thread-BLAS `os.environ.setdefault` guard already enforced by `python/conftest.py`

**File:** `python/leafblower/test_parity_weights.py:15-17`
**Issue:** The file sets `OMP_NUM_THREADS`/`OPENBLAS_NUM_THREADS`/`MKL_NUM_THREADS`
via `os.environ.setdefault(...)` before importing `numpy`, with a comment
explaining the ordering requirement. `python/conftest.py` (unchanged this
phase) already sets the identical three env vars at collection time, before
any test module (including this one) is imported by pytest, making this
file's copy redundant for any pytest-driven run. Not a bug — harmless
defense-in-depth, and correctly ordered even if it were the only guard — but
worth noting so a future edit doesn't assume conftest.py is the sole source
of truth or drops this file's copy assuming conftest.py alone is sufficient
in a non-pytest invocation context (e.g. running the file's functions
directly).
**Fix:** None required. Optionally add a one-line comment cross-referencing
`python/conftest.py` so the duplication reads as intentional rather than
copy-paste.

### IN-02: `test_solver_parity.py` module docstring's `raking`/`sinkhorn` max_error figures are unlabeled magic numbers with no re-verification hook

**File:** `python/leafblower/test_solver_parity.py:108-109`
**Issue:** The `_CONV_TOL` rationale comment records `raking max_error =
5.551e-17` and `sinkhorn max_error = 1.110e-16` as measured constants
justifying why the shared `0.01` tolerance is safe. These are comments, not
assertions — if a future core change silently shifts either solver's
converged residual upward (while staying under 0.01), nothing re-checks that
the comment's specific numbers still hold; the comment will silently go
stale. This matches the pre-existing pattern for `greg`/`newton_kl`/etc. in
the same block (not a new practice introduced by this diff), so it is
Info-level consistency, not a new defect.
**Fix:** None required now — matches established file convention. If
tightening later, consider a single parametrized "each converges under
`_CONV_TOL`" assertion is already present via `_assert_parity`'s precheck;
the comment's numbers are documentation only and don't need their own gate.

---

_Reviewed: 2026-08-15T00:00:34Z_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_
