# Phase 2: One Engine, Not Two - Context

**Gathered:** 2026-08-15
**Status:** Ready for planning

<domain>
## Phase Boundary

Unify R and Python's numerical dispatch paths so both languages call `rk_calibrate()` for
every solver. Today Python already does this (`python/leafblower/_bindings.cpp:194`); R's
`C_rk_calibrate` (`src/r_bridge.cpp:654-899`) runs its own `strcmp(method_str, ...)` chain
straight into the solver functions, bypassing `rk_calibrate()` entirely and hand-duplicating
AUTO routing, result packing, and (per SC3) potentially bound-enforcement logic. This phase
closes that gap — one dispatch path, not two hand-synced ones — plus the two smaller
consistency gates in ROADMAP SC2/SC4.

</domain>

<decisions>
## Implementation Decisions

### Dispatch migration strategy (leafblower-rywn, P0)
- **D-01:** Migrate solver-by-solver, incrementally — not one atomic cutover. Each task
  migrates one solver from `r_bridge.cpp`'s `strcmp` branch to routing through
  `rk_calibrate()`, verified (R↔Python parity green, full DoD gate green) before the next
  solver starts. The `strcmp` chain shrinks task by task until empty; AUTO routing
  (currently duplicated inline in `r_bridge.cpp:663-788` vs `c_api.cpp:292-397`) is
  consolidated as part of whichever task removes the last string-dispatched branch that
  depends on it. — **Reversibility:** costly — each migrated solver's commit changes a
  live R entry point (`C_rk_calibrate`); reverting a single solver's migration after later
  solvers have landed means re-adding one `strcmp` branch back into a chain that has
  otherwise moved on, not a clean `git revert`.
- **D-02:** No feature-flag or env-var fallback to the old dispatch path was requested — a
  bake-in safety net was explicitly not selected as a discussion area. Rollback is via git
  history (revert the offending solver's migration commit), consistent with the
  solver-by-solver granularity in D-01.

### SC3 — water-fill duplication scope
- **D-03:** Direct code read (this session) shows every solver (`oris`, `raking`,
  `chebyshev`, `greenkhorn`, `greg`, `logit_calib`, `sinkhorn`) already routes its
  per-cell `bounds_mode="unit"` finalization through the shared
  `finalize_weights`/`finalize_weights_buf` in `calib_dispatch.hpp` (per `CLAUDE.md`'s
  canonical-home convention) — `oris_finalize.cpp:150-159` even has an explicit
  "Delegate to the single source" comment documenting a prior consolidation. Raking's own
  `water_fill_cat` (`raking.cpp:168`) is different math — an inline per-margin IPF box
  projection integral to the algorithm, not the post-hoc per-cell finalize step SC3 is
  about. **Treat SC3 as already satisfied by the existing architecture.** Phase 2's SC3
  work is a verification task, not a code change: assert (via grep or a targeted test)
  that every solver's finalize path routes through the shared helper, so a future
  regression that reintroduces a local per-solver copy fails loudly instead of silently
  reintroducing the duplication SC3 exists to prevent. — **Reversibility:** reversible —
  a verification-only task with no behavioral change.
- **D-04:** If the researcher or planner finds a genuine remaining duplication during
  deeper investigation (this session's read was not exhaustive), it supersedes D-03 as
  the actual SC3 target — flag it explicitly rather than silently closing SC3.

### SC4 — build-list divergence gate
- **D-05:** Wire the `src/*.cpp` vs `python/CMakeLists.txt:CORE_SOURCES` sync check into
  the existing DoD gate (`.coverage-thresholds.json`'s `enforcement.command`) as a test
  assertion — R or Python, whichever the researcher finds cleaner to wire given the
  existing test infrastructure — so it runs automatically every time the DoD gate runs,
  the same way every other check in this repo does (no CI exists; the DoD gate command IS
  the CI). Not a separate manual script or Makefile target the user has to remember to
  run. — **Reversibility:** reversible — a new test assertion, no production code change.

### Claude's Discretion
- Whether the SC4 sync-check assertion lives in the R or Python test suite is left to the
  researcher/planner — pick whichever existing test infrastructure makes the check
  simplest to wire and least likely to rot.
- Exact task boundaries for the solver-by-solver migration order (D-01) — which solver
  goes first, how AUTO routing consolidation is sequenced relative to the 9 individual
  solver migrations — are a planning decision, not locked here.

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Dispatch unification (SC1)
- `src/r_bridge.cpp:654-899` — the `strcmp(method_str, ...)` dispatch chain to be
  replaced, including the duplicated inline AUTO routing at `:663-788`
- `src/c_api.cpp:414+` — `rk_calibrate()`'s enum-based dispatch, the target single path;
  AUTO routing at `:292-397`
- `src/c_api.cpp:458`, `:473`, `:505` — existing "mirrors r_bridge.cpp" comments pinned to
  specific line numbers in the other file; these go stale on migration and must be
  updated or removed as each solver moves
- `python/leafblower/_bindings.cpp:194` — Python's existing correct pattern (already
  calls `rk_calibrate()` directly) — the shape R migrates toward
- `.beads` ticket `leafblower-rywn` (P0) — full objective, reference data, and drift
  markers already captured

### Optimization-level asymmetry (SC2, separately ticketed)
- `.beads` ticket `leafblower-qzto` (P1) — already scoped as a documentation fix
  (CLAUDE.md's "no -O level" claim is R-only true; Python hard-sets `-O3`). Not
  discussed further here — the CRAN `tools:::.check_make_vars` constraint forecloses
  equalizing the flags, so "document as deliberate, bounded decision" is the only live
  option per SC2's own wording.
- `python/CMakeLists.txt:99`, `src/Makevars.in:12-15`, `configure:11-15` — the asymmetry's
  source locations

### Water-fill consolidation (SC3)
- `src/calib_dispatch.hpp` — canonical home for `finalize_weights`/`finalize_weights_buf`
  per `CLAUDE.md`
- `src/oris_finalize.cpp:150-159` — comment documenting a prior consolidation
  ("Delegate to the single source")
- `src/raking.cpp:164-168` — `water_fill_cat`, the algorithm-internal IPF box projection
  that is NOT the same duplication SC3 targets (see D-03)

### Build-list divergence gate (SC4)
- `python/CMakeLists.txt` (`CORE_SOURCES` list) vs `src/*.cpp` (R auto-globs; Python does
  not) — `CLAUDE.md`'s "Two build sites" footgun section documents this split
- `.coverage-thresholds.json` — `enforcement.command`, the DoD gate this check must join

### Project-level constraints (apply throughout)
- `CLAUDE.md` — no LTO (hot per-iteration code must stay with its caller across any TU
  boundary change); no `-O` in R's `PKG_CXXFLAGS`; ABI `static_assert` tripwires
  (`sizeof(rk_params_t) == 264`, `sizeof(rk_result_t) == 536`); enum slots 2 and 7
  permanently reserved; RAII-safe `Rf_error` unwinding convention
  (`src/r_bridge.cpp:270, 642-648`)
- `.planning/ROADMAP.md` Phase 2 section — full SC1-SC5, "Notes for planning"

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `src/calib_dispatch.hpp` — the shared dispatch-helper home; `finalize_weights`/
  `finalize_weights_buf` already prove the pattern works for bounds-mode finalization.
  The dispatch unification (SC1) should follow the same shape: pull the per-solver
  dispatch-and-pack logic that both `c_api.cpp` and `r_bridge.cpp` currently duplicate
  into one shared definition here, per `leafblower-rywn`'s stated objective.
- `python/leafblower/_bindings.cpp` — a working reference implementation of "call
  `rk_calibrate()`, don't reimplement dispatch" that R's migration should converge toward.

### Established Patterns
- Every solver already ends its finalize step by calling into `calib_dispatch.hpp`'s
  shared bounds-mode helper (see SC3 discussion, D-03) — this is the existing precedent
  for "shared helper, not per-solver copy" that SC1's dispatch unification should extend.
- RAII-safe error handling: `r_bridge.cpp` destroys all heap-backed locals before calling
  `Rf_error` (longjmp skips C++ dtors) — any migrated dispatch code touching this path
  must preserve that convention (`src/r_bridge.cpp:270, 642-648`).

### Integration Points
- `src/r_bridge.cpp::C_rk_calibrate()` is the sole R-side integration point being
  rewired — it currently extracts scalars, dispatches via `strcmp`, and packs results
  into SEXPs manually; post-migration it should extract scalars, call `rk_calibrate()`,
  and pack the returned `rk_result_t` into SEXPs (mirroring how `c_api.cpp` and
  `_bindings.cpp` both already consume `rk_calibrate()`'s output).

</code_context>

<specifics>
## Specific Ideas

No UI/visual specifics — this is a pure backend/numerical-core refactor. The concrete
detail captured is the migration granularity (D-01: solver-by-solver, not atomic) and the
SC4 enforcement mechanism (D-05: wired into the existing DoD gate, not a separate manual
script).

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope. The rollback/bake-in safety-net area (a
feature-flag fallback to the old dispatch during a bake-in period) was presented as a
discussion option but not selected; D-02 records that git-revert-based rollback at
solver-by-solver granularity is the accepted approach instead.

</deferred>

---

*Phase: 2-One Engine, Not Two*
*Context gathered: 2026-08-15*
