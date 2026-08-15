---
gsd_state_version: 1.0
milestone: v1.0
milestone_name: milestone
current_phase: 05
current_phase_name: cran-pypi-release
status: executing
stopped_at: Completed 05-01-PLAN.md
last_updated: "2026-08-15T17:35:37.248Z"
last_activity: 2026-08-15
last_activity_desc: Plan 02-08 Task 3 resolved — user approved after independent re-verification; phase 2 closed
progress:
  total_phases: 5
  completed_phases: 4
  total_plans: 23
  completed_plans: 19
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-08-15)

**Core value:** Calibrated weights that are numerically correct, bound-respecting, and
identical from R and from Python.
**Current focus:** Phase 05 — cran-pypi-release

## Current Position

Phase: 05 (cran-pypi-release) — EXECUTING
Plan: 2 of 5
Status: Ready to execute
Last activity: 2026-08-15 — Phase 05 execution started
user-visible surface (harvest() field parity across all 9 solvers + AUTO in R, all 9
in Python with auto correctly R-only) after a first "approved" claim was rejected for
coming from the orchestrator's own run rather than the user. SC1-SC5 all have durable
evidence: SC1 via test_single_dispatch_site.py (RED/GREEN verified), SC5 via the full
DoD gate (R 0 FAIL/1833 PASS, Python 160 passed/0 failed) and stepstone benchmark
(byte-identical to the 02-07 baseline). leafblower-rywn (P0 dispatch-unification epic)
closed with DoD evidence.

Progress: [████████░░] 83% (phase 2)

**Brownfield.** The package is at v0.1.0 with eight shipped solvers and 1478 closed beads
tickets. 0% here measures the *remaining* work in this roadmap, not the product.

## Repository Facts

- Branch `master`, local-only, **NO git remote**. Complete = committed locally + gates green.
- HEAD `2b94e8e` "docs: map existing codebase".
- `bd` (beads) holds the live task queue: 1495 total, 1478 closed, **13 open** (1 P0, 5 P1,
  7 P2/P3). Roadmap phases sit ABOVE beads and reference existing ticket IDs.

- Working tree is dirty at bootstrap (`.wolf/` deletions, `.beads/issues.jsonl`,
  `.mcp.json`). Commit with an explicit pathspec, NEVER `git add -A`.

## Performance Metrics

**Velocity:**

- Total plans completed: 14
- Average duration: —
- Total execution time: —

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01 | 4 | - | - |
| 2 | 8 | - | - |
| 04 | 2 | - | - |

**Recent Trend:** No data yet.
**Per-Plan Metrics:**

| Plan | Duration | Tasks | Files |
|------|----------|-------|-------|
| Phase 02 P01 | 25min | 3 tasks | 3 files |
| Phase 02 P02 | ~35min | 3 tasks | 5 files |
| Phase 02 P03 | ~30min | 3 tasks | 3 files |
| Phase 2 P4 | ~35min | 2 tasks | 3 files |
| Phase 02 P05 | ~35min | 2 tasks | 3 files |
| Phase 02 P06 | ~20min | 2 tasks | 3 files |
| Phase 02 P07 | 50min | 3 tasks | 3 files |
| Phase 02 P08 | ~30min | 3 tasks | 1 file |
| Phase 03-honest-performance-gate P01 | ~55min | 3 tasks | 5 files |
| Phase 03-honest-performance-gate P02 | ~35min | 3 tasks | 3 files |
| Phase 03 P03 | 30min | 3 tasks | 2 files |
| Phase 03-honest-performance-gate P04 | ~30min | 3 tasks | 4 files |
| Phase 04 P01 | 25min | 2 tasks | 4 files |
| Phase 04 P02 | 15min | 2 tasks | 3 files |
| Phase 05 P01 | 19min | 2 tasks | 149 files |

## Accumulated Context

### Decisions

**There is no locked-decision layer** — the 44-doc corpus contains 0 ADRs. See PROJECT.md
§ Key Decisions for the two documents that act as decision records without being ADRs.

Carried into planning:

- Phase order is TDD-first: verification coverage (Phase 1) precedes the P0 dispatch
  unification (Phase 2) so the rewire lands against a net that can catch it.

- `leafblower-2ouc` (benchmark-study article) is tracked outside this roadmap; Phase 3
  reuses its `benchmarks/` infrastructure rather than duplicating it.

- [Phase ?]: Plan 02-01: DispatchResult dropped the plan's specified const CellTable& parameter — no solver needs it externally (each builds its own from CalibState); documented as a Rule-3 deviation.
- [Phase ?]: Plan 02-01: found a 6th superset-only field (aa_accepted_count, ORISResult) beyond the 5 named in leafblower-rywn/RESEARCH.md — filed on the ticket, not fixed; relevant to a future oris/oris_soft migration plan.
- [Phase ?]: Plan 02-02: SC3/SC4 closed as verification-only (D-03/D-05) — no solver code changed. newton_calib.cpp's exclusion from the SC3 shared-finalize sweep is deliberate and cross-referenced to the pre-existing leafblower-og7d.5 (no duplicate ticket filed — that ticket already fully records the same gap). SC2 closed as documentation-only per leafblower-qzto's own DoD wording (CRAN forecloses equalizing -O flags); leafblower-qzto closed.
- [Phase ?]: [Phase 02] Plan 02-03: greg, greenkhorn, logit migrated onto the shared dispatch table (RK_ALG_GREG/GREENKHORN/LOGIT case arms in lbw::dispatch_solver); ran full testthat suite instead of the plan's filter="greg" verify command, which matches zero files (testthat filter matches basenames only) — filed as a Rule 3 deviation, not a weaker check.
- [Phase ?]: [Phase 02] Plan 02-04: chebyshev and raking migrated onto the shared dispatch table (RK_ALG_CHEBYSHEV/RK_ALG_RAKING case arms). Chebyshev's oris warm-start now exists once inside dispatch_solver (no divergence found between the two prior per-bridge implementations). Raking is the first solver with a superset-only field (sraa_demoted); rk_result_t stays 536 bytes, Python's result dict unchanged.
- [Phase ?]: [Phase 02] Plan 02-05: oris, oris_soft migrated onto the shared dispatch table (RK_ALG_ORIS/RK_ALG_ORIS_SOFT case arms). Completed the plan's under-enumerated DispatchResult field list to the full oris-family diagnostic set actually exported to R (Rule 2 deviation); rk_result_t already carries all but aa_accepted_count, so a new c_api.cpp helper (pack_dispatch_oris_extras_c) packs them without touching the shared pack function used by the other 6 migrated solvers. No pack_oris_result vs pack_oris_result_c value divergence found. capacity_penalty/estimate_M_cell resolves identically on both bridges, at most once per explicit oris_soft solve.
- [Phase ?]: [Phase 02] Plan 02-06: newton_kl migrated onto the shared dispatch table (RK_ALG_NEWTON_KL case arm) -- the last named method on the string chain. c_api.cpp's explicit branch reuses plan 02-05's pack_dispatch_oris_extras_c (verified field-by-field byte-identical to the pre-migration reset) instead of a near-duplicate helper. All 9 non-AUTO dispatch_solver arms confirmed to assign stall_kind, verified to reach harvest.R's production accelerate heuristic unbroken. Only RK_ALG_AUTO routing remains unmigrated (plan 07).
- [Phase ?]: [Phase 02] Plan 02-07: lbw::route_auto() and lbw::kAlgNames added to calib_dispatch.hpp; R bridge's method-string chain fully collapsed to one lbw::dispatch_solver() call per explicit method (plus route_auto's up to two for AUTO). Two genuine, behaviorally-inert AUTO-implementation divergences found (st.oris_auto_selected set unconditionally vs. only-when-ORIS; accelerate not restored on c_api.cpp's fallback) -- fixed by adopting R's fixture-pinned behavior, tracked on leafblower-uqyf (closed same session). SC1 fully satisfied: one dispatch site for every solver and for AUTO routing.
- [Phase ?]: [Phase 02] Plan 02-08 (phase gate, complete): test_single_dispatch_site.py added -- SC1 now enforced by a pytest guard (RED/GREEN verified), not only satisfied by current source. Corrected the plan's stated "at most two dispatch_solver call sites" bound to the actual, architecturally-correct 3 (AUTO primary + AUTO fallback + the one unified explicit-method call), cross-checked against 02-07-SUMMARY.md's own D4 coverage entry. Full DoD gate proven green (R 0 FAIL/1833 PASS/141 WARN/13 SKIP; Python 160 passed/0 failed) and stepstone benchmark (kk1204) byte-identical to the 02-07 baseline (SC5 proven with numbers). leafblower-rywn closed with DoD evidence. Task 3's checkpoint:human-verify: a first "approved" message was rejected because it explicitly disclosed the orchestrator ran the verification itself, not the user -- no agent message substitutes for the user's own sign-off on a gate="blocking" checkpoint. The user then approved directly, after independent corroboration (test-newton-kl.R + test-cr-d5-auto-fallback-fields.R: 0 FAIL/29 PASS; test_solver_parity.py + test_parity_weights.py: 21 passed/0 failed). **Phase 2 (One Engine, Not Two) is complete: SC1-SC5 all have durable evidence.**
- [Phase ?]: Task 1: oris_soft convergence must use its own canonical marginal_kl/improvement rule, never a competitor's absolute/max_err/threshold shorthand — reintroduces the stopped-early confound.
- [Phase ?]: Task 2/3: testthat::test_dir(filter=...) requires NOT_CRAN=true (unlike devtools::test()) to actually run skip_on_cran()-gated assertions; new file paths inside such filtered runs must anchor on testthat::test_path(), not a bare relative string.
- [Phase ?]: 03-02: known_limit_k20_uniform uses the original 2026-04-23 kk1204 investigation's skewed target (0.3/0.175x4), not test-bench-gate.R's literal uniform 1/5 -- uniform targets on this exact data converge trivially (max_error ~4e-15), which cannot back the 'known limit is unachievable' claim. Fixture n/K/seed/max_weight/naming stay byte-identical to the in-repo test.
- [Phase ?]: 03-02: deff/n_eff measured for known_limit_k20_uniform (deff=2.44-2.74, n_eff~183k-205k) diverge ~3000-5700x from leafblower-ylsy's cited DEFF 8000-14000/n_eff 71-118 -- traced to that closure using a fifth, more-severe kk1204 skew variant (0.6/0.2/0.1/0.07/0.03 at n=1e6), confirming RESEARCH.md Pitfall 2's divergent-fixture-description finding.
- [Phase ?]: 03-02: full D-07 competitor set (survey, icarus, ReGenesees) measured on the medium class, all agreeing with oris_soft to within 5e-4 max_error; icarus requires method='logit' (not 'raking') to honour bounds; ReGenesees totals filled manually via pop.template() column-name parsing (no sampling frame available for fill.template()).
- [Phase ?]: D-06/D-10 paired framing: honest gate asserts wall_s<=0.5s AND n_eff>=60000 alongside existing bound/accuracy checks (03-03, developer-selected)
- [Phase ?]: README.md's headline claim excludes any autumn mention (D-05, user-emphatic); orienting sentence rephrased around DESCRIPTION's own Title/Description instead
- [Phase ?]: kk1204 block's fixture uses a uniform 1/5 target, not 03-02's skewed target -- does not exercise the documented degenerate case; this plan states the mismatch in-comment rather than changing the fixture (out of Task 3's scope)
- [Phase ?]: US-003/KPI-04 stay Traceability status Partial after Phase 3 (requirements.mark-complete correctly no-oped, not_found, 0 bytes written) -- verbose=1 clause unexercised and no single measurement clears the PRD's literal 1M+ rows AND 20+ margins target together
- [Phase ?]: D-04 resolved option-a: harvest(..., weights=w) now hard-stop()s naming design_weights=, not a warning
- [Phase ?]: [Phase 04] Plan 02: deleted docs/raking.md's §8.2/§12 ORIS/L-BFGS-B misattribution passages outright (no replacement text, per D-01: docs/methods/oris.md is the single authoritative description); annotated rk_algorithm_t slot 7 as removed GRAKE matching slot 2's convention; SC4 re-audit found zero new true-positive grake/lbfgsb/cp references. SC1/SC2/SC4 closed. R suite 0 FAIL/1837 PASS.
- [Phase ?]: [Phase 05] Plan 01: removed 11 tracked dev-artifact strays + 120-file leafblower.Rcheck/ tree from git; extended .Rbuildignore 23->58 lines (14 planned + 21 more found by the real R CMD build run, since R CMD build scans the working tree not the git index). Fixed a LaTeX-manual ERROR + non-ASCII WARNING at the root (Unicode math symbols in R/harvest.R roxygen comments + 4 sibling files), a NEWS.md unparseable NOTE (headers retitled to leafblower 0.1.0), and an unstated-test-deps WARNING (arrow/callr/jsonlite/withr added to Suggests). Final local check: 0 ERRORs, 1 WARNING, 3 NOTEs -- does not literally meet the plan's 0-warnings/<=1-NOTE bar; remaining findings traced to this machine's absent checkbashisms/tidy/V8 (Rule-3 excludes package-manager installs) plus a local R-Makeconf NOTE, both documented in cran-comments.md and the WINDOWS.md ledger, expected to close on the 05-03/05-04 CI matrix.

### Pending Todos

None captured yet.

### Blockers/Concerns

- **`leafblower-kk1.20.4` is a decision, not just work.** The composite gate
  "<30 s AND <1e-6" is structurally unachievable on K=20 uniform-random input. Phase 3
  cannot start planning until the REFRAME option is chosen (three options on the ticket).

- **No CI exists.** KPI-06 (Python 3.9–3.13 matrix) needs a pipeline or a documented manual
  matrix — decide during Phase 5 planning.

- `.wolf/buglog.json` is deleted in the working tree, so no recurring-bug history was
  available to the codebase audit.

## Deferred Items

| Category | Item | Status | Deferred At |
|----------|------|--------|-------------|
| *(none)* | | | |

## Session Continuity

Last session: 2026-08-15T17:35:37.241Z
Stopped at: Completed 05-01-PLAN.md
Next: Plan Phase 3 (Honest Performance Gate) — reuse leafblower-2ouc's benchmarks/ infrastructure per the carried-forward decision above; leafblower-kk1.20.4's REFRAME decision (30s/<1e-6 gate on kk1204) must be chosen before Phase 3 planning starts.
Resume file: None
