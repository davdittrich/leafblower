---
gsd_state_version: 1.0
milestone: v1.0
milestone_name: milestone
current_phase: 2
current_phase_name: One Engine, Not Two
status: executing
stopped_at: Completed 02-01-PLAN.md (sinkhorn dispatch tracer)
last_updated: "2026-08-15T01:52:52.782Z"
last_activity: 2026-08-15
last_activity_desc: ROADMAP.md, REQUIREMENTS.md, PROJECT.md created from doc ingest
progress:
  total_phases: 2
  completed_phases: 1
  total_plans: 12
  completed_plans: 5
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-08-15)

**Core value:** Calibrated weights that are numerically correct, bound-respecting, and
identical from R and from Python.
**Current focus:** Phase 2 — One Engine, Not Two (dispatch unification, solver-by-solver)

## Current Position

Phase: 2 of 5 (One Engine, Not Two)
Plan: 01 of 8 complete (sinkhorn dispatch tracer)
Status: In progress — 7 solvers remain to migrate through the shared dispatch table (D-01)
Last activity: 2026-08-15 — Plan 02-01 executed: lbw::DispatchResult + lbw::dispatch_solver()
added to calib_dispatch.hpp; sinkhorn routed through it from both c_api.cpp and r_bridge.cpp.
Full DoD gate green (R 0 FAIL/1833 PASS, Python 156 passed/0 failed).

Progress: [████░░░░░░] 42%

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

Last session: 2026-08-15T01:52:52.776Z
Stopped at: Completed 02-01-PLAN.md (sinkhorn dispatch tracer)
Awaiting roadmap approval before `/gsd-plan-phase 1`.
Resume file: None
