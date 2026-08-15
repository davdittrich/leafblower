# Phase 2: One Engine, Not Two - Research

**Researched:** 2026-08-15
**Domain:** C++17 calibration-weighting core, dual R (`.Call`) / Python (pybind11) FFI dispatch unification
**Confidence:** HIGH (all claims below are direct-read source verification against `master` on 2026-08-15; no web research was needed — this is a pure internal-refactor phase with zero new external dependencies)

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

**Dispatch migration strategy (leafblower-rywn, P0)**
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

**SC3 — water-fill duplication scope**
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

**SC4 — build-list divergence gate**
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

### Deferred Ideas (OUT OF SCOPE)
None — discussion stayed within phase scope. The rollback/bake-in safety-net area (a
feature-flag fallback to the old dispatch during a bake-in period) was presented as a
discussion option but not selected; D-02 records that git-revert-based rollback at
solver-by-solver granularity is the accepted approach instead.
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| US-004 (residual) | "A stable C API (`leafblower.h`) lets any language with C FFI call the engine without duplicating algorithm code — both the R `.Call()` bridge and the Python pybind11 module call the same `rk_calibrate()` symbol." Currently **Partial**: `src/r_bridge.cpp:654-899` dispatches on a method string straight into `lbw::<solver>_solve` and never calls `rk_calibrate`. | See **Architecture Patterns → Pattern 1** for the exact shape of the fix, and **Common Pitfalls → Pitfall 1** for the single biggest risk to closing this requirement without a silent regression (4 R-visible result fields have no home in `rk_result_t`, the literal C-ABI struct `rk_calibrate()` returns). |
</phase_requirements>

## Summary

This phase is a pure internal C++17 refactor with **zero new external dependencies** —
there is nothing to research in the "which library" sense. All research effort instead
went into direct source verification of the two dispatch paths named in CONTEXT.md,
because the single highest-risk fact for planning is not visible from the phase
description alone: **`rk_calibrate()`'s C-ABI result struct (`rk_result_t`) is missing
four fields that R's current `harvest()` result list exposes today and that live
`testthat` regression tests assert on** (`n_projected_dims`, `lm_mu_final`,
`sraa_demoted`, `convergence_stall_kind`). `stall_kind` is additionally consumed in
production by `R/harvest.R:702`, not just tests. This means SC1 ("the R bridge reaches
the solvers through the same path as the C ABI") **cannot** mean "R literally calls the
exported `rk_calibrate()` function and marshals `rk_result_t`" for the `newton_kl`,
`oris`, `oris_soft`, and `raking` migration tasks — that would silently drop
user-visible fields and break existing tests. The already-open ticket `leafblower-rywn`
independently arrived at the correct shape: a shared dispatch table **and a shared
result-extraction helper into a neutral internal struct**, from which `c_api.cpp`
marshals the ABI-frozen `rk_result_t` (narrower) and `r_bridge.cpp` marshals the SEXP
result list (fuller, unchanged field set). "Same path" = same dispatch table + same
neutral extraction helper in `calib_dispatch.hpp`, not literally the same exported
C symbol.

The second major finding is a full confirmation of D-03: exactly 7 of the 8 shipped
solvers (`oris`, `raking`, `chebyshev`, `greenkhorn`, `greg`, `logit_calib`, `sinkhorn`)
already route unit-mode bounds enforcement through the single shared
`finalize_weights`/`finalize_weights_buf` in `calib_dispatch.hpp`. The 8th,
**`newton_kl`, never calls it at all** — it has no per-obs water-fill and silently
ignores `bounds_mode="unit"` (only counts violations and falls back to `NOCONV` above a
5% violation threshold). This is a distinct, adjacent finding from SC3's stated
duplication-prevention scope (D-04's "genuine remaining duplication... supersedes D-03"
clause does not quite apply — there's no duplication here, there's an absence) and
should be flagged to the planner as an explicit scope decision, not silently folded into
or excluded from the SC3 verification task.

SC4 (`src/*.cpp` vs `CORE_SOURCES` sync) currently has **zero drift** — all 17
non-`r_bridge.cpp` files in `src/` are already listed in `python/CMakeLists.txt`'s
`CORE_SOURCES`. This is a pure regression-prevention gate: a new pytest assertion
comparing the two file lists, added once, never touched again.

**Primary recommendation:** Sequence D-01's solver-by-solver migration so that the
low-risk solvers with no extra R-only fields (`sinkhorn`, `greg`, `greenkhorn`, `logit`,
`chebyshev`) go first to prove the shared-table pattern, and treat `oris`/`oris_soft`
(needs `sraa_demoted`) and `newton_kl` (needs `n_projected_dims`, `lm_mu_final`,
`stall_kind`) and `raking` (needs `sraa_demoted`, `stall_kind`) as the tasks that must
also extend the neutral internal struct — file that extension as `rywn`'s own DoD item
already requires ("Any pre-existing field drift found in step 1 is filed as its OWN
ticket, not fixed here"). Do **not** touch `rk_result_t` or the frozen 264/536-byte ABI
sizes for this — the neutral struct is a new, unconstrained internal type; `rk_result_t`
stays exactly as-is and simply doesn't receive the 4 extra fields (matching Python's
current, already-accepted behavior).

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Solver dispatch (string/enum → solver call) | C++ core (`calib_dispatch.hpp`) | — | Both R and C-ABI callers must resolve to the same table; today it's duplicated in `r_bridge.cpp` and `c_api.cpp` — this phase consolidates it into the shared header, per `CLAUDE.md`'s canonical-home rule. |
| Result-field extraction (`res.base.*` → output) | C++ core (`calib_dispatch.hpp`, new neutral struct) | R bridge (SEXP marshal) / C API (`rk_result_t` marshal) | The extraction logic (which `res.base.*`/solver-specific fields exist per algorithm) is shared; the two *marshaling* targets differ in width (SEXP list has 4 fields `rk_result_t` lacks) and are each owned by their respective FFI boundary. |
| Per-obs bounds enforcement (`bounds_mode="unit"` water-fill) | C++ core (`calib_dispatch.hpp::finalize_weights_buf`) | — | Already the single shared implementation for 7/8 solvers (SC3, verified). `newton_kl` is the one solver that bypasses this tier entirely (see Pitfall 2) — a capability gap, not a tier misassignment. |
| Build-manifest sync (`src/*.cpp` vs `CORE_SOURCES`) | Python test suite (`python/leafblower/test_*.py`) | — | The divergence risk is Python-side only (R auto-globs; Python's CMake does not), so the check should live next to the artifact it protects — see D-05 discretion. |
| Optimization-level documentation (SC2) | Docs (`CLAUDE.md`) | — | Not a code-tier concern at all — CRAN's `tools:::.check_make_vars` forecloses equalizing flags; the only live option is documenting the asymmetry as deliberate (already scoped in `leafblower-qzto`, separately ticketed, not part of this research's code-level findings). |

## Standard Stack

Not applicable — this phase adds no new libraries, packages, or frameworks. It is a
pure internal refactor of existing C++17 code within `src/`, consumed by the existing R
and Python FFI layers. No `npm install` / `pip install` / `cargo add` of any kind is
implicated.

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Shared dispatch table + neutral result struct in `calib_dispatch.hpp` (header, `inline` functions) | A new `.cpp` translation unit (`calib_dispatch.cpp`) | Rejected per `leafblower-rywn`'s own guard ("If a new `src/*.cpp` is created despite the guard above, it [needs justification]") and `CLAUDE.md`'s no-LTO constraint: cross-TU calls don't inline, and while dispatch itself is cold (once-per-solve, so TU-boundary cost is irrelevant to the no-LTO rule), a new `.cpp` file adds a second `CORE_SOURCES` entry to keep in sync — the header-only approach `finalize_weights`/`finalize_weights_buf` already use is simpler and has zero new build-list surface. |
| Literal `rk_calibrate()` C-ABI call from R bridge | Shared dispatch table + wider neutral struct, narrowed differently per caller | The literal call would silently drop `n_projected_dims`/`lm_mu_final`/`sraa_demoted`/`stall_kind` from R's user-visible result (verified: `tests/testthat/test-newton-kl.R:25-29`, `test-newton-tsvd-projection.R:17-45`, `test-cr-d5-auto-fallback-fields.R:35`, and `R/harvest.R:702` production consumer of `stall_kind`) — a regression, not a refactor. The neutral-struct approach satisfies SC1's "same path" language (shared dispatch + shared extraction) without narrowing R's existing output. |

## Package Legitimacy Audit

Not applicable — no external packages are installed, upgraded, or introduced by this
phase. `Package Legitimacy Gate` protocol skipped; nothing to check against a registry.

## Architecture Patterns

### System Architecture Diagram (current state → target state)

```
CURRENT (two independent dispatch paths):

  R harvest()                          Python harvest()
       |                                     |
       v                                     v
  C_rk_calibrate()                    _bindings.cpp pybind11 wrapper
  (r_bridge.cpp:234-1136)                    |
       |                                     v
       | strcmp(method_str, ...)      rk_calibrate() [c_api.cpp:256]
       | chain (:654-899), own AUTO         |
       | routing copy (:663-788)      alg = enum switch [c_api.cpp:282-334]
       |                                     |
       v                                     v
  lbw::<solver>_solve(st)  <----+---->  lbw::<solver>_solve(st)
  (oris/raking/sinkhorn/          |     (same solver functions,
   greg/greenkhorn/logit/         |      reached via the enum path)
   newton_kl/chebyshev)           |
       |                          |
       v                          v
  res.base.* (CalibResult,   res.base.* -> pack_solver_result()/
  types.hpp) manually copied  pack_oris_result_c()/pack_newton_result_c()
  into ~49 res_* locals            (c_api.cpp:57-173)
       |                          |
       v                          v
  SEXP result list (49 named   rk_result_t (45 fields, leafblower.h;
  elements, r_bridge.cpp:      536-byte ABI-frozen struct)
  1017-1127) -- INCLUDES            |
  4 fields rk_result_t lacks        v
       |                       Python result_dict (mirrors rk_result_t
       v                       fields exactly, _bindings.cpp:203-247;
  R user's harvest() output    does NOT have the 4 R-only fields either)


TARGET (one shared dispatch + extraction site):

  R harvest()                          Python harvest()
       |                                     |
       v                                     v
  C_rk_calibrate()                    _bindings.cpp pybind11 wrapper
  (r_bridge.cpp, slimmed)                    |
       |                                     v
       +-----------+          +------------> rk_calibrate() [c_api.cpp]
                   v          v
        calib_dispatch.hpp: shared {enum<->string<->solver} table
        + shared result-extraction helper -> NEUTRAL internal struct
        (superset: includes n_projected_dims/lm_mu_final/sraa_demoted/
         stall_kind, since these already exist on res.base.*/solver-
         specific result types -- no new solver-side computation needed)
                   |          |
                   v          v
        R marshals NEUTRAL   c_api.cpp marshals NEUTRAL struct's
        struct -> SEXP        SUBSET -> rk_result_t (ABI-frozen,
        (49 elements,          536 bytes, unchanged) -> Python
        unchanged R surface)   result_dict (unchanged)
```

**Why the "target" R path does not shrink to zero extra logic:** `rk_calibrate()`'s
signature is `int rk_calibrate(int n, int K, double* weights, ..., const rk_params_t*
params, rk_result_t* result)` — fixed by `leafblower.h`, ABI-frozen, ships in the public
header. R cannot get `n_projected_dims` etc. out of a literal call to this function
because `rk_result_t` doesn't carry them ([VERIFIED: src/leafblower.h:112-161] — the
full `rk_result_t` field list has no `n_projected_dims`/`lm_mu_final`/`sraa_demoted`/
`stall_kind` member; `[VERIFIED: src/c_api.cpp:150-151]` — explicit comment: `"(n_projected_dims
/ lm_mu_final are NewtonCalibResult-only; rk_result_t does not carry them, so nothing
newton-specific to copy into the C result here.)"`). SC1's "same path" is satisfied by
routing through the **same dispatch table and extraction logic** in `calib_dispatch.hpp`
— not by literally calling the narrower exported C symbol from R.

### Recommended Project Structure

No new files/folders. The dispatch table and neutral extraction helper belong in the
existing canonical-home header:

```
src/
├── calib_dispatch.hpp   # ADD: {enum, string, solver-fn} table + neutral result
│                        #      extraction helper, alongside the existing
│                        #      finalize_weights/finalize_weights_buf (already
│                        #      the proven pattern for "shared helper, not
│                        #      per-solver copy" in this codebase)
├── c_api.cpp            # SLIM: rk_calibrate() calls the shared table, marshals
│                        #       neutral struct -> rk_result_t (narrower)
├── r_bridge.cpp          # SLIM: C_rk_calibrate() calls the shared table, marshals
│                        #       neutral struct -> SEXP list (unchanged 49 fields)
└── leafblower.h          # UNCHANGED: rk_params_t (264B), rk_result_t (536B),
                           #            rk_algorithm_t enum — no ABI edits in this phase
```

### Pattern 1: Shared dispatch table + neutral result struct (the SC1 fix)

**What:** One `{rk_algorithm_t, method_string, solver-invocation}` table in
`calib_dispatch.hpp`, consumed by both `c_api.cpp`'s `rk_calibrate()` and
`r_bridge.cpp`'s `C_rk_calibrate()`. Paired with one result-extraction function that
reads every field off `res.base.*` (and the solver-specific result type, e.g.
`ORISResult`/`NewtonCalibResult`) into a neutral struct — a strict superset of
`rk_result_t`'s field set, since it must also carry `n_projected_dims`, `lm_mu_final`,
`sraa_demoted`, `stall_kind`.

**When to use:** This is the target architecture for `leafblower-rywn`'s whole DoD, not
a per-task pattern — every solver migration task (D-01) converges toward it.

**Example (existing precedent this migration should mirror — verified working code):**
```cpp
// Source: src/calib_dispatch.hpp:359-433 (existing, unmodified) — this is the
// ALREADY-PROVEN shape for "shared helper both R and C-API callers invoke",
// applied today to bounds finalization. The dispatch-table work in this phase
// is architecturally the same move, applied to solver selection + result packing.
inline void finalize_weights_buf(double* w, int n, const CalibState& st,
                                 const CellTable& ct,
                                 int& n_bounds_violated, int& n_bounds_clamped) {
    // ... single implementation, called by oris_finalize.cpp, raking.cpp,
    // chebyshev.cpp, greenkhorn.cpp, greg.cpp, logit_calib.cpp, sinkhorn.cpp
}
inline void finalize_weights(CalibState& st, const CellTable& ct,
                             int& n_bounds_violated, int& n_bounds_clamped) {
    finalize_weights_buf(st.weights, st.n, st, ct, n_bounds_violated, n_bounds_clamped);
}
```

### Pattern 2: kAlgNames positional table with reserved-slot holes

**What:** Both `r_bridge.cpp` and `c_api.cpp` maintain a 12-element `const char*[]`
indexed by `rk_algorithm_t`, with empty/placeholder strings at indices 2 and 7 (the
permanently-reserved removed-`lbfgsb` and removed-`grake` slots).

**When to use:** Any migration that touches these tables (or creates a third,
consolidated one) must keep exactly 12 positional entries and the two documented holes.

**Example — the two EXISTING tables that must either be unified or kept in exact sync:**
```cpp
// Source: src/r_bridge.cpp:943-957 (verified verbatim)
static const char* kAlgNames[] = {
    "",           // 0 = RK_ALG_AUTO
    "oris",       // 1 = RK_ALG_ORIS
    "",           // 2 = (removed lbfgsb slot)
    "raking",     // 3 = RK_ALG_RAKING
    "sinkhorn",   // 4 = RK_ALG_SINKHORN
    "chebyshev",  // 5 = RK_ALG_CHEBYSHEV
    "greg",       // 6 = RK_ALG_GREG
    "",           // 7 = deprecated GRAKE
    "oris_soft",  // 8 = RK_ALG_ORIS_SOFT
    "greenkhorn", // 9 = RK_ALG_GREENKHORN
    "logit",      // 10 = RK_ALG_LOGIT
    "newton_kl",  // 11 = RK_ALG_NEWTON_KL
};
static const int kAlgNamesLen = 12;
static_assert(RK_ALG_NEWTON_KL == 11, "kAlgNames table needs update on enum change");
```
```cpp
// Source: src/c_api.cpp:30-47 (verified verbatim) — DIFFERENT display strings
// for the SAME enum ("ORIS" vs "oris", "(reserved)"/"(gap)" vs "") — exactly
// the kind of drift SC1's consolidation should eliminate.
static constexpr const char* kAlgNames[] = {
    "auto", "ORIS", "(reserved)", "raking", "sinkhorn", "chebyshev", "greg",
    "(gap)", "ORIS-soft", "greenkhorn", "logit", "newton_kl"
};
static_assert(sizeof(kAlgNames)/sizeof(kAlgNames[0]) == 12,
    "kAlgNames must cover all 12 enum slots 0..11");
```

### Pattern 3: `kAlgMap` string→enum lookup already deduplication-ready

**What:** `r_bridge.cpp` already has exactly ONE string→`rk_algorithm_t` map
(`kAlgMap`, an `unordered_map<string_view, rk_algorithm_t>`), used once to resolve
`method_sexp` into `p.algorithm` before the `strcmp` dispatch chain runs. This part of
SC1's "the string↔enum mapping must exist exactly once" (per `leafblower-rywn`'s Format
guard) is **already satisfied** — the remaining work is the dispatch/execution chain
below it, not the string resolution itself.

**Example (verified verbatim, already correct — do not duplicate this into the shared table; reuse it):**
```cpp
// Source: src/r_bridge.cpp:94-105
const std::unordered_map<std::string_view, rk_algorithm_t> kAlgMap = {
    {"oris",       RK_ALG_ORIS},
    {"oris_soft",  RK_ALG_ORIS_SOFT},
    {"raking",     RK_ALG_RAKING},
    {"greg",       RK_ALG_GREG},
    {"chebyshev",  RK_ALG_CHEBYSHEV},
    {"sinkhorn",   RK_ALG_SINKHORN},
    {"auto",       RK_ALG_AUTO},
    {"greenkhorn", RK_ALG_GREENKHORN},
    {"logit",      RK_ALG_LOGIT},
    {"newton_kl",  RK_ALG_NEWTON_KL},
};
```

### Anti-Patterns to Avoid
- **Literal `rk_calibrate()` call from R for the field-bearing solvers:** As shown in
  Pattern 1, this silently drops 4 R-visible, test-guarded fields for `newton_kl`
  (`n_projected_dims`, `lm_mu_final`), `oris`/`oris_soft`/`raking` (`sraa_demoted`), and
  every solver that can emit `RK_ERR_STALL` (`stall_kind`, consumed in production by
  `R/harvest.R:702`). Route through the shared dispatch + neutral-extraction layer
  instead.
- **Renumbering or resizing `rk_algorithm_t`, `rk_params_t`, or `rk_result_t` to "fix"
  the field gap:** Forbidden by `CLAUDE.md` ("enum slots 2 and 7 stay holes") and the
  phase's own hard constraint ("the ABI `static_assert` tripwires (264 / 536 bytes) must
  be re-measured if any struct changes"). The neutral struct is a NEW, unconstrained
  internal type — it does not touch `rk_result_t`'s frozen layout.
- **Moving hot per-iteration solver code into `calib_dispatch.hpp` or a new TU as part
  of "consolidation":** The dispatch/extraction logic runs once per `harvest()` call
  (cold); nothing here should touch the per-iteration inner loops inside
  `oris.cpp`/`raking.cpp`/etc. `CLAUDE.md`'s no-LTO rule means a hot loop moved across a
  TU boundary loses inlining — this phase's scope is entirely the once-per-solve
  dispatch-and-pack wrapper, never the solvers' internals.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Per-obs bounds water-fill for `bounds_mode="unit"` | A new per-solver water-fill copy | `calib_dispatch.hpp::finalize_weights`/`finalize_weights_buf` (already exists, already used by 7/8 solvers) | This is SC3 itself — the shared helper already exists; the task is verifying/enforcing its use, not building a new one. |
| Result-field extraction per FFI boundary | Two independent field-copy blocks (current state) | One shared extraction into a neutral struct in `calib_dispatch.hpp`, narrowed per-caller at the marshal step | Prevents the exact drift SC1 exists to close — a field added to one path and forgotten on the other. |

**Key insight:** Every "don't hand-roll" item in this phase already has its shared
implementation living in `calib_dispatch.hpp` (or, for the extraction helper, should be
added there) — the pattern is proven by `finalize_weights`, not novel.

## Common Pitfalls

### Pitfall 1: `rk_result_t` cannot carry 4 fields R currently exposes to users
**What goes wrong:** A naive migration of `newton_kl`, `oris`, `oris_soft`, or `raking`
that routes R through the literal `rk_calibrate()` C-ABI function (rather than a shared
dispatch table with a wider internal struct) silently drops
`n_projected_dims`/`lm_mu_final` (newton_kl-only), `sraa_demoted` (oris/oris_soft/raking
only), and `stall_kind` (any solver reaching `RK_ERR_STALL`) from R's `harvest()` result.
**Why it happens:** `rk_result_t` (`leafblower.h:112-161`) is ABI-frozen at 536 bytes
(`static_assert(sizeof(rk_result_t) == EXPECTED_RK_RESULT_BYTES)`,
`leafblower.h:245-247`) and was never extended to carry these 4 fields — Python's
binding has the identical gap today (verified: `_bindings.cpp:194-247` builds its
`result_dict` from exactly `rk_result_t`'s fields, nothing more) and nobody has noticed
because Python users never had these fields to begin with. R's richer SEXP list
(`r_bridge.cpp:1017-1127`, 49 named elements) was built by hand-copying off
`res.base.*`/solver-specific result structs directly — a path this migration is
removing.
**How to avoid:** Route the field-bearing solvers through a **new, unconstrained
internal struct** in `calib_dispatch.hpp` that both callers extract from; `c_api.cpp`
marshals a *subset* of it into `rk_result_t` (unchanged ABI); `r_bridge.cpp` marshals
the *full* struct into the SEXP list (unchanged R-visible surface). Per
`leafblower-rywn`'s own DoD ("Any pre-existing field drift found in step 1 is filed as
its OWN ticket, not fixed here"), this gap should be named explicitly in the plan as
existing-drift-to-preserve, not silently regressed.
**Warning signs:** `tests/testthat/test-newton-kl.R:25-29`,
`test-newton-tsvd-projection.R:17-45`, `test-newton-kl-tsvd-ratio.R:40-41`, and
`test-cr-d5-auto-fallback-fields.R:35` all assert directly on these fields — any of
these going red (or silently returning `NULL`/`NA`) during the `newton_kl`/`oris`/
`raking` migration tasks is this pitfall firing. `R/harvest.R:694-702` is the
**production** consumer of `stall_kind` (`calib_result$convergence_stall_kind`,
labeled "leafblower-8eod: use solver-emitted convergence_stall_kind... replaces
accelerate_bool heuristic in harvest.R") — a regression here degrades a real user-facing
diagnostic silently (falls back to `NA`/wrong heuristic), not just a test failure.

### Pitfall 2: `newton_kl` never enforces `bounds_mode="unit"` — an SC3-adjacent gap, not SC3 itself
**What goes wrong:** A user calling `harvest(..., method="newton_kl", bounds_mode="unit")`
gets weights that may violate `[min_weight, max_weight]` per-observation — `newton_kl`
only *counts* violations (`newton_calib.cpp:869`, `res.n_bounds_violated = n_violated`)
and falls back to `RK_ERR_NOCONV` above a 5% violation fraction; it never calls
`finalize_weights`/`finalize_weights_buf` at all (confirmed: zero matches for
`finalize_weights` in `newton_calib.cpp`, versus one call each in all 7 other solver
`.cpp` files), so there is no per-cell water-fill redistribution for this solver.
**Why it happens:** `newton_kl` is a smooth-dual Newton method with no box-constrained
inner step (unlike ORIS/raking's cell-table water-fill machinery) — this may be a
deliberate algorithmic choice (bounds enforcement isn't structurally available inside a
KKT/dual-Newton step the way it is inside an IPF/BCD sweep) rather than an oversight,
but nothing in `validation.cpp` rejects or warns on the `newton_kl` + `bounds_mode="unit"`
combination — it is accepted silently and simply does less than the other 7 solvers.
**How to avoid:** This is **not** what D-03 verified SC3 to be about (D-03 named exactly
the 7 solvers that share the finalize helper) and is **not** duplication (D-04's trigger
condition), so it does not automatically supersede D-03's "SC3 = verification only"
disposition. Flag this explicitly to the user/planner as a scope decision: either (a)
leave it out of Phase 2 entirely (file a standalone ticket, since it's an `newton_kl`
capability gap unrelated to the *duplication* SC3 targets), or (b) fold a
documentation-only note into the SC3 verification task ("`newton_kl` is excluded from
this invariant by design; `bounds_mode="unit"` has no effect for this method"). Silently
ignoring it (treating SC3's verification-passing as proof `bounds_mode="unit"` works
uniformly across all 8 solvers) would be wrong — say so explicitly per CONTEXT.md's own
instruction not to silently close gaps.
**Warning signs:** Any SC3 verification test written as "grep every solver `.cpp` for a
`finalize_weights` call" will correctly find 7/8 and must decide what to do about the
8th — don't let it silently pass by only checking the 7 solvers D-03 named.

### Pitfall 3: The `estimate_M_cell` AUTO-routing cache is currently duplicated per dispatch site, not just per solver
**What goes wrong:** AUTO routing computes `M_cell_est` (compression-ratio estimate) at
up to two call sites within the *same* `C_rk_calibrate()` invocation today — once for
`oris_soft`'s `capacity_penalty` auto-resolution (`r_bridge.cpp:479-499`, cached in
`m_cell_est_cache`) and again inside the `"auto"` branch (`r_bridge.cpp:671-675`, reusing
the cache if already computed). `c_api.cpp`'s AUTO path (`c_api.cpp:299`) computes it
independently a third time with no cache-sharing across the R/C-API boundary (each
process has its own call, this is expected) but WITHIN `c_api.cpp` there is only one
site (`oris_soft`'s capacity resolution at `c_api.cpp:530` and AUTO's routing at
`c_api.cpp:299` are two more potentially-redundant computations on the same inputs).
**Why it happens:** Historical incremental feature addition (`oris_soft`'s
`capacity_penalty` auto-resolution and AUTO's `target_skew`/`M_cell` routing were added
at different times) rather than a single up-front computation.
**How to avoid:** When consolidating into the shared dispatch table, decide once whether
`M_cell_est` is computed eagerly (before the table dispatches) or lazily-and-cached (the
current `r_bridge.cpp` pattern) — and apply the SAME policy on both the R and C-API
paths, rather than porting R's cache-with-two-guarded-sites pattern into `c_api.cpp`
verbatim (which has its own, textually different, redundant computation today). This is
a candidate for `leafblower-rywn`'s step-1 "field present on one path, absent on the
other" inventory, though it's a computation-redundancy finding rather than a field gap.
**Warning signs:** A profiling regression on `AUTO` + `oris_soft`-adjacent paths after
consolidation, or (more likely) simply an inconsistency in whether the migrated shared
table computes `M_cell_est` once or twice per call.

### Pitfall 4: RAII-unwind discipline in `r_bridge.cpp` must survive the refactor
**What goes wrong:** `Rf_error()` performs a `longjmp` that skips C++ destructors
(R-exts §5.5). `C_rk_calibrate` currently guarantees every heap-backed local
(`std::vector`s, `std::string`s) is explicitly destroyed via `std::vector<T>().swap(...)`
/ `std::string().swap(...)` *before* any `Rf_error()` call, on both the error path
(`r_bridge.cpp:914-938`) and the success path just before SEXP allocation
(`r_bridge.cpp:983-994`, since `Rf_allocVector`/`Rf_mkChar` can themselves OOM-trigger
`Rf_error`). A refactor that introduces new heap-backed locals into the shared dispatch
call (e.g. a neutral result struct with `std::vector<double> best_weights` embedded)
must extend this same swap-before-`Rf_error` discipline to the new locals, or a new
leak-on-error path is introduced.
**Why it happens:** This is exactly the kind of easy-to-miss discipline that a
mechanical "extract shared code" refactor can silently drop if the neutral struct isn't
audited against every `Rf_error()` call site.
**How to avoid:** After defining the neutral struct, re-walk `r_bridge.cpp:270,
642-648` and both swap-release blocks (`:914-938` error path, `:983-994` success path)
and confirm every new heap-backed member of the neutral struct is included in both.
**Warning signs:** `CLAUDE.md` calls this out explicitly as a hard constraint ("the
RAII-safe `Rf_error` unwinding convention... must be preserved") — treat any migration
task touching `r_bridge.cpp` as required to re-verify this, not just re-run tests (a
leak-on-error path won't fail a test unless the test specifically forces an error and
checks for a leak, which none currently do per a grep of `tests/testthat` for
`Rf_error`-adjacent leak assertions — this is a manual-review item, not a
test-gate-catchable one).

## Code Examples

### Before: the two divergent dispatch chains (what SC1 replaces)

```cpp
// Source: src/r_bridge.cpp:654-663 (verified verbatim) — string dispatch,
// bypasses rk_calibrate() entirely
if (strcmp(method_str, "raking") == 0) {
    auto res = lbw::raking_solve(st);
    res_status     = res.base.status;
    res_iterations = res.base.iterations;
    res_max_error  = res.base.max_error;
    res_alg_used   = (int)RK_ALG_RAKING;
    pack_solver_result(res);
    res_sraa_demoted = res.sraa_demoted ? 1 : 0;
    res_best_weights = std::move(res.base.best_weights);
} else if (strcmp(method_str, "auto") == 0) {
    // ... 124 more lines of AUTO routing duplicated from c_api.cpp:292-397
```

```cpp
// Source: src/c_api.cpp:414-420 (verified verbatim) — enum dispatch, the
// target shape SC1's shared table generalizes
if (alg == RK_ALG_RAKING) {
    // Classical raking: IPF + Dykstra box + Dykstra hyperplane (renamed from ORIS)
    auto res = lbw::raking_solve(st);
    status = res.base.status;
    iterations = res.base.iterations;
    max_error = res.base.max_error;
    used = RK_ALG_RAKING;
    // ... 20 more lines of manual field copy, NO res.sraa_demoted anywhere
```

### Existing regression tests that pin the 4 at-risk fields (must stay green)

```r
# Source: tests/testthat/test-newton-kl.R:25-29 (verified verbatim)
# WH-d: lm_mu_final must surface in R result via r_bridge SEXP-pack
expect_true("lm_mu_final" %in% names(res),
  label="WH-d: lm_mu_final must surface in R result")
expect_true(is.finite(res$lm_mu_final),
  label="WH-d: lm_mu_final must be finite on converging fixture")
```

```r
# Source: R/harvest.R:694-702 (verified verbatim) — PRODUCTION consumer, not a test
# leafblower-8eod: use solver-emitted convergence_stall_kind (set at RK_ERR_STALL
# ... oris.cpp fires for both SRAA (accelerate=TRUE, stall_kind=1) and plain-BCD
# (accelerate=FALSE, stall_kind=2) — NOT bijective with the user flag; required
# route (a). stall_kind=0 -> NA (no stall).
sk <- calib_result$convergence_stall_kind
```

### SC4 sync-check target shape (new pytest assertion, D-05)

```python
# NEW FILE, e.g. python/leafblower/test_core_sources_sync.py
# Compares src/*.cpp (minus r_bridge.cpp, which is R-only) against
# python/CMakeLists.txt's CORE_SOURCES list. Verified baseline (2026-08-15):
# both sides currently list exactly the same 17 files -- see the grep below,
# reproduced from this research session:
#   src/*.cpp (18 files) minus r_bridge.cpp = 17 files:
#     calib_linalg.cpp calib_validate.cpp c_api.cpp cell_table.cpp
#     chebyshev.cpp design_effect.cpp greenkhorn.cpp greg.cpp logit.cpp
#     logit_calib.cpp newton_calib.cpp oris.cpp oris_finalize.cpp
#     oris_trajectory.cpp raking.cpp sinkhorn.cpp validation.cpp
#   python/CMakeLists.txt CORE_SOURCES (python/CMakeLists.txt:60-77):
#     ../src/c_api.cpp ../src/calib_linalg.cpp ../src/calib_validate.cpp
#     ../src/cell_table.cpp ../src/chebyshev.cpp ../src/greenkhorn.cpp
#     ../src/greg.cpp ../src/oris.cpp ../src/oris_trajectory.cpp
#     ../src/oris_finalize.cpp ../src/logit.cpp ../src/logit_calib.cpp
#     ../src/newton_calib.cpp ../src/raking.cpp ../src/sinkhorn.cpp
#     ../src/design_effect.cpp ../src/validation.cpp
# -> zero drift today; this test is pure regression prevention.
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

def test_core_sources_matches_src_glob():
    src_files = {p.name for p in (REPO_ROOT / "src").glob("*.cpp")} - {"r_bridge.cpp"}
    cmake_text = (REPO_ROOT / "python" / "CMakeLists.txt").read_text()
    listed = set(re.findall(r"\.\./src/(\w+\.cpp)", cmake_text))
    assert src_files == listed, (
        f"src/*.cpp vs CORE_SOURCES drift: "
        f"missing from CORE_SOURCES={src_files - listed}, "
        f"extra in CORE_SOURCES={listed - src_files}"
    )
```

## State of the Art

Not applicable in the "library/framework evolution" sense — this phase touches no
external tooling. The one relevant "state of the art" fact is internal: `leafblower-rywn`
(filed 2026-08-14) already independently designed the shared-table-plus-neutral-struct
approach this research arrived at (see its Step-by-Step Logic item 4: "Define one
shared result-extraction helper that reads `res.base.*` into a neutral struct. Both
callers then marshal that struct into their own output format (SEXP vs `rk_result_t`)."),
confirming this is not a novel research recommendation but convergent with existing
project planning already captured in beads.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | `newton_kl`'s lack of `bounds_mode="unit"` enforcement is an existing, accepted design gap rather than an unfiled bug — no beads ticket or doc explicitly says so, this research only confirmed the absence of enforcement code and the absence of a validation-layer rejection. | Common Pitfalls → Pitfall 2 | If wrong (i.e., it actually is an unnoticed regression), the planner should file it as its own defect ticket rather than treat it as pre-existing accepted scope; either way Phase 2 should not silently paper over it in the SC3 verification task. |
| A2 | The SC4 sync-check is best placed in the Python test suite (`python/leafblower/test_*.py`) rather than R's `testthat`. | Code Examples → SC4 sync-check | Low risk — D-05 explicitly leaves this to researcher/planner discretion; if the planner prefers R-side (e.g. because R's `R CMD INSTALL` step runs first in the DoD gate command and would fail faster), that's an equally valid choice with no functional difference, since the check only reads text files and needs no compiled artifact from either side. |

## Open Questions (RESOLVED)

1. **RESOLVED by plan 02-01** — Exact shape of the neutral internal struct (field types, ownership of `best_weights`)
   - What we know: It must be a strict superset of `rk_result_t`'s fields (verified list:
     `leafblower.h:112-161`) plus the 4 R-only fields (verified:
     `r_bridge.cpp:558-568,589,637,661,714-715,730,757-758,782,842-843`) plus the
     obs-level `best_weights` vector (currently `res_best_weights`, a
     `std::vector<double>` — note `c_api.cpp` already normalizes this into the flat
     `weights` output array for `greenkhorn`/`logit` at `c_api.cpp:459-462,474-477`, so
     R's separate tracking of `res_best_weights` may become redundant once R also copies
     from the output `weights` buffer the same way — worth checking whether R's
     `res_best_weights`-into-`weights` copy at `r_bridge.cpp:974-981` becomes
     unnecessary once dispatch is shared, since the shared dispatch could do this copy
     once for both callers).
   - What's unclear: Whether to define this as a new named struct in `calib_dispatch.hpp`
     or as an out-parameter bundle threaded through the existing per-solver `*Result`
     types (avoiding a new struct definition entirely, at the cost of the extraction
     helper needing per-solver-type overloads/templates, similar to the existing
     `pack_solver_result`/`pack_oris_result`/`pack_newton_result_c` template pattern
     already in both files).
   - Recommendation: This is a planning/implementation-detail decision, not a research
     gap — either shape satisfies SC1 as long as it's ONE definition consumed by both
     callers. Leave to the planner/implementer; flag in the plan that
     `leafblower-rywn`'s own Step 1 ("Build the inventory FIRST... tabulate... every
     res.base.* field copied) on BOTH paths") is the mechanism for finalizing this
     shape, and should run before the struct is defined.

2. **RESOLVED by Phase 1** — Whether `oris_soft`'s missing entry in the R↔Python weight-parity matrix
   (`tests/test_parity_weights.py:73`, per REQUIREMENTS.md's un-ticketed concerns list)
   should block or merely accompany the SC1 migration.**
   - What we know: REQUIREMENTS.md flags this as "the one shipped solver absent from the
     R↔Python weight-parity matrix... High priority per the concerns audit" but it is not
     named in CONTEXT.md's decisions, ROADMAP's SC1-SC5, or either beads ticket
     (`rywn`/`qzto`) for Phase 2.
   - What's unclear: Whether Phase 2's SC1 migration of `oris_soft` should add the
     missing parity test as part of the same task (since the migration touches
     `oris_soft`'s dispatch anyway) or stay strictly scoped to dispatch unification and
     leave the parity-matrix gap to whichever phase REQUIREMENTS.md's traceability table
     assigns it (not shown in the excerpt read this session as explicitly phase-mapped).
   - Recommendation: Out of this research's scope to decide — surfacing it here so the
     planner can make an explicit in/out-of-scope call rather than the gap going
     unmentioned.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| R | R build/test gate (`R CMD INSTALL --preclean .`, testthat) | Yes | 4.6.1 | — |
| cmake | Python build (`scikit-build-core` invokes it) | Yes | 4.4.2 | — |
| Python (project `.venv`, uv-managed) | Python build/test gate | Yes | 3.14 (`.venv/bin/python3.14`) | — |
| LAPACK | Linked by `python/CMakeLists.txt:93-94` (`find_package(LAPACK REQUIRED)`); R side links via R's own BLAS/LAPACK | Not independently re-verified this session (existing DoD gate already exercises this successfully per STATE.md's Phase 1 completion) | — | — |

No missing dependencies identified — this phase requires no new tools beyond what the
existing DoD gate (`R CMD INSTALL` + testthat + `uv pip install -e .` + pytest) already
exercises successfully per Phase 1's completion.

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | R: testthat edition 3 (`DESCRIPTION:22`, verified: `Config/testthat/edition: 3`); Python: pytest |
| Config file | `DESCRIPTION` (R); `python/pyproject.toml` + `python/conftest.py` (Python) |
| Quick run command | R: `Rscript -e 'library(testthat); library(leafblower); test_dir("tests/testthat", filter="<name>")'`; Python: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -m pytest -k <name> -q` |
| Full suite command | `.coverage-thresholds.json`'s `enforcement.command` (verified verbatim): `R CMD INSTALL --preclean . && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript -e 'library(testthat); library(leafblower); out <- test_dir("tests/testthat", stop_on_failure=TRUE)' && cd python && uv pip install -e . --reinstall-package leafblower && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -m pytest -q` |

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| US-004 (SC1: one dispatch site) | R and Python route through the same shared dispatch table; result fields identical to pre-migration baseline for every solver | regression (existing suites) | Full DoD gate (`.coverage-thresholds.json`) — both `tests/testthat/*.R` (100 files) and `python/leafblower/test_*.py` (20+ files, incl. `test_parity_weights.py`, `test_solver_parity.py`) | ✅ (all exist; per D-01, run per-solver-migration, not just at phase end) |
| US-004 (SC1, field-gap risk) | The 4 R-only fields (`n_projected_dims`, `lm_mu_final`, `sraa_demoted`, `convergence_stall_kind`) survive migration unchanged | regression (existing suites) | `tests/testthat/test-newton-kl.R`, `test-newton-tsvd-projection.R`, `test-newton-kl-tsvd-ratio.R`, `test-cr-d5-auto-fallback-fields.R` | ✅ — these are the exact tripwires; no new test needed, just don't break them |
| SC3 (water-fill single source) | Every solver's unit-mode finalize routes through `calib_dispatch.hpp::finalize_weights[_buf]` | new verification test | grep-based or reflection-based assertion (e.g. a script/test asserting each solver `.cpp` calls the shared helper) | ❌ Wave 0 — new test, per D-03's "verification task, not a code change" |
| SC4 (build-list sync) | `src/*.cpp` (minus `r_bridge.cpp`) == `python/CMakeLists.txt` `CORE_SOURCES` | new regression-prevention test | `python/leafblower/test_core_sources_sync.py` (see Code Examples) or R-side equivalent | ❌ Wave 0 — new test, per D-05 |

### Sampling Rate
- **Per task commit (per D-01, per solver migration):** the targeted solver's own R
  testthat file(s) + `python/leafblower/test_solver_parity.py` /
  `test_parity_weights.py` for that method.
- **Per wave merge:** Full DoD gate (`.coverage-thresholds.json`'s `enforcement.command`).
- **Phase gate:** Full DoD gate green, plus `LBW_BENCH_GATE=1` stepstone benchmark run
  on every commit touching a TU boundary (per the phase's own "Notes for planning" —
  though this migration is dispatch-only and should not create new TU boundaries per
  Pattern 1/Alternatives Considered, so this is a should-not-trigger check, not an
  expected-to-trigger one).

### Wave 0 Gaps
- [ ] `python/leafblower/test_core_sources_sync.py` (or R-side equivalent per D-05) —
  covers SC4.
- [ ] A new SC3 verification test/script asserting all 7 (or 8, pending Pitfall 2's
  scope decision) solvers route through `finalize_weights`/`finalize_weights_buf` —
  covers SC3.
- [ ] `leafblower-rywn`'s own Step 1 field inventory (a `bd comment`, not a test file)
  should run before any solver migration task starts, per its own Step-by-Step Logic.

## Security Domain

Not applicable in the network/auth/injection sense — this package has no network
surface, no authentication, no user-facing string interpolation into a query/shell/HTML
context. The one relevant ASVS-adjacent control already exists and is unaffected by this
phase's scope:

| ASVS Category | Applies | Standard Control |
|---------------|---------|-----------------|
| V5 Input Validation | yes (pre-existing, unchanged) | `lbw::validate_calibrate_inputs` (`validation.cpp`), called identically from both `r_bridge.cpp:444-456` and `c_api.cpp:336` before any solver runs — this phase's dispatch consolidation does not touch validation, only what happens after it passes. |

No STRIDE-relevant threat patterns apply to a single-process, single-user numerical
library with no I/O beyond R/Python-supplied in-memory arrays.

## Sources

### Primary (HIGH confidence — direct source read, this session, 2026-08-15)
- `src/r_bridge.cpp` (full dispatch chain, SEXP construction, RAII-unwind blocks)
- `src/c_api.cpp` (full `rk_calibrate()` body, `pack_solver_result`/`pack_oris_result_c`/
  `pack_newton_result_c` templates)
- `src/leafblower.h` (`rk_algorithm_t`, `rk_params_t`, `rk_result_t`, ABI `static_assert`s)
- `src/calib_dispatch.hpp` (`finalize_weights`/`finalize_weights_buf`/`regate_unit_status`)
- `src/oris_finalize.cpp`, `src/raking.cpp`, `src/newton_calib.cpp` (per-solver
  finalize/bounds behavior)
- `python/leafblower/_bindings.cpp` (Python's `rk_calibrate()`-based reference pattern)
- `python/CMakeLists.txt` (`CORE_SOURCES` list, `-O3` flag)
- `R/harvest.R` (`stall_kind` production consumer)
- `tests/testthat/test-newton-kl.R`, `test-newton-tsvd-projection.R`,
  `test-newton-kl-tsvd-ratio.R`, `test-cr-d5-auto-fallback-fields.R`
- `.beads` tickets `leafblower-rywn`, `leafblower-qzto` (`bd show`)
- `.coverage-thresholds.json`, `DESCRIPTION`, `.planning/REQUIREMENTS.md`,
  `.planning/STATE.md`, `.planning/phases/02-one-engine-not-two/02-CONTEXT.md`

### Secondary / Tertiary
None used — this phase required no web research; all findings are direct codebase
verification.

## Metadata

**Confidence breakdown:**
- Standard stack: N/A — no new dependencies
- Architecture: HIGH — every claim traces to a quoted, line-numbered source read this
  session
- Pitfalls: HIGH — the field-gap finding (Pitfall 1) is cross-verified via source
  (`leafblower.h`, `c_api.cpp`'s own comment), a live test suite, and a production R
  consumer; the `newton_kl` bounds gap (Pitfall 2) is verified via a negative grep
  (absence of `finalize_weights` calls) confirmed against all 8 solver files

**Research date:** 2026-08-15
**Valid until:** No expiry driver — this is a point-in-time codebase state snapshot, not
a library/API research finding. Re-verify line numbers if `master` moves significantly
before planning executes (git HEAD was `2b94e8e` at CONTEXT.md gathering; confirm
unchanged before relying on exact line citations).
