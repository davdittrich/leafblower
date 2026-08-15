---
gsd_state_version: 1.0
milestone: v1.0
milestone_name: milestone
current_phase: 2
current_phase_name: One Engine, Not Two
status: blocked
stopped_at: 02-08-PLAN.md Tasks 1-2 complete (SC1 guard test, SC5 DoD+benchmark proof); Task 3 checkpoint:human-verify (gate=blocking) awaiting human sign-off
last_updated: "2026-08-15T03:27:09Z"
last_activity: 2026-08-15
last_activity_desc: Plan 02-08 Tasks 1-2 executed — test_single_dispatch_site.py (SC1) added, full DoD gate + stepstone benchmark proven green (SC5); Task 3 human-verify checkpoint pending
progress:
  total_phases: 2
  completed_phases: 1
  total_plans: 12
  completed_plans: 11
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-08-15)

**Core value:** Calibrated weights that are numerically correct, bound-respecting, and
identical from R and from Python.
**Current focus:** Phase 2 — One Engine, Not Two (dispatch unification, solver-by-solver)

## Current Position

Phase: 2 of 5 (One Engine, Not Two)
Plan: 8 of 8 — Tasks 1-2 complete, Task 3 (human-verify checkpoint, gate=blocking) pending
Status: BLOCKED on checkpoint — phase gate awaiting human sign-off before phase 2 can close
Last activity: 2026-08-15 — Plan 02-08 Tasks 1-2 executed: test_single_dispatch_site.py
makes SC1 enforceable by test (RED/GREEN verified against a temporary probe); full DoD
gate proven green (R testthat 0 FAIL/1833 PASS/141 WARN/13 SKIP; Python pytest 160
passed/0 failed); stepstone benchmark (LBW_BENCH_GATE=1, kk1204) byte-identical to the
02-07 baseline (status=0 iters=10 best_error=-7.376e-14 time=1.5s) — SC5 proven with
numbers. leafblower-rywn (P0 dispatch-unification epic) closed with DoD evidence.
Task 3 requires a human to run harvest() for all 9 solvers + AUTO in R and Python and
confirm the user-visible surface (R-only fields present, Python fields unchanged) —
see 02-08-SUMMARY.md's "Next Phase Readiness" for the exact steps.

Progress: [█████████░] 96%

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

- Total plans completed: 4
- Average duration: —
- Total execution time: —

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01 | 4 | - | - |

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
- [Phase ?]: [Phase 02] Plan 02-08 (Tasks 1-2, phase gate): test_single_dispatch_site.py added -- SC1 now enforced by a pytest guard (RED/GREEN verified), not only satisfied by current source. Corrected the plan's stated "at most two dispatch_solver call sites" bound to the actual, architecturally-correct 3 (AUTO primary + AUTO fallback + the one unified explicit-method call), cross-checked against 02-07-SUMMARY.md's own D4 coverage entry. Full DoD gate proven green (R 0 FAIL/1833 PASS/141 WARN/13 SKIP; Python 160 passed/0 failed) and stepstone benchmark (kk1204) byte-identical to the 02-07 baseline (SC5 proven with numbers). leafblower-rywn closed with DoD evidence. **Task 3 (checkpoint:human-verify, gate=blocking) NOT resolved -- awaiting human sign-off before phase 2 closes.**

### Pending Todos

- Task 3 of 02-08-PLAN.md: human must run `harvest()` for all 9 solvers + AUTO in R and Python, confirm R-only fields (`n_projected_dims`, `lm_mu_final`, `sraa_demoted`, `convergence_stall_kind`) present in R and absent from Python's result dict, and judge whether the recorded stepstone delta is acceptable. See 02-08-SUMMARY.md.

### Blockers/Concerns

- **Phase 2 cannot close until 02-08 Task 3 resolves.** A human must respond "approved" (or describe what changed unexpectedly) to the checkpoint documented in 02-08-SUMMARY.md's "Next Phase Readiness" section.

- **`leafblower-kk1.20.4` is a decision, not just work.** The composite gate
  "<30 s AND <1e-6" is structurally unachievable on K=20 uniform-random input. Phase 3
  cannot start planning until the REFRAME option is chosen (three options on the ticket).

- **Phase 2 is the high-risk phase.** Unifying the two dispatch tables rewires the R path
  for all eight solvers under a no-LTO constraint and frozen ABI tripwires.

- **No CI exists.** KPI-06 (Python 3.9–3.13 matrix) needs a pipeline or a documented manual
  matrix — decide during Phase 5 planning.

- `.wolf/buglog.json` is deleted in the working tree, so no recurring-bug history was
  available to the codebase audit.

## Deferred Items

| Category | Item | Status | Deferred At |
|----------|------|--------|-------------|
| *(none)* | | | |

## Session Continuity

Last session: 2026-08-15T03:27:09Z
Stopped at: 02-08-PLAN.md Tasks 1-2 complete; Task 3 checkpoint:human-verify (gate=blocking) pending
Next: Resolve 02-08 Task 3 (human sign-off) — see 02-08-SUMMARY.md's "Next Phase Readiness" for exact verification steps. Phase 2 closes once approved.
Resume file: .planning/phases/02-one-engine-not-two/02-08-PLAN.md (Task 3)
