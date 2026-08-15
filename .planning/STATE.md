---
gsd_state_version: 1.0
milestone: v1.0
milestone_name: milestone
current_phase: 2
current_phase_name: One Engine, Not Two
status: executing
stopped_at: Completed 02-05-PLAN.md (oris/oris_soft dispatch migration)
last_updated: "2026-08-15T02:38:10.921Z"
last_activity: 2026-08-15
last_activity_desc: Plan 02-02 executed — build-list sync test, shared-finalize delegation test, -O asymmetry documented
progress:
  total_phases: 2
  completed_phases: 1
  total_plans: 12
  completed_plans: 9
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-08-15)

**Core value:** Calibrated weights that are numerically correct, bound-respecting, and
identical from R and from Python.
**Current focus:** Phase 2 — One Engine, Not Two (dispatch unification, solver-by-solver)

## Current Position

Phase: 2 of 5 (One Engine, Not Two)
Plan: 5 of 8 complete (SC4/SC3/SC2 guards)
Status: In progress — 7 solvers remain to migrate through the shared dispatch table (D-01)
Last activity: 2026-08-15 — Plan 02-02 executed: two new pytest regression guards
(python/leafblower/test_core_sources_sync.py for SC4, test_finalize_weights_sync.py for
SC3) plus SC2 (-O optimization-level asymmetry) documented in CLAUDE.md,
python/CMakeLists.txt, and test_solver_parity.py's header; leafblower-qzto closed.
Full DoD gate green (R 0 FAIL/1839 PASS, Python 159 passed/0 failed).

Progress: [████████░░] 75%

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

### Pending Todos

None captured yet.

### Blockers/Concerns

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

Last session: 2026-08-15T02:38:10.914Z
Stopped at: Completed 02-05-PLAN.md (oris/oris_soft dispatch migration)
Next: 02-03-PLAN.md — migrate greg, greenkhorn, logit onto the shared dispatch table (SC1).
Resume file: None
