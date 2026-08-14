# Removed solver slots — supersession record

**Date:** 2026-08-14
**Status:** Accepted
**Type:** SPEC (supersession record)

## Purpose

Three solvers were specified in this corpus and are absent from the live
`rk_algorithm_t` enum in `src/leafblower.h`. No earlier spec records their
removal, so a reader resolving the corpus by date concludes each is still live.
This document is the missing supersession record. It supersedes every earlier
spec's treatment of `grake`, `lbfgsb`, and `cp` as live methods.

Each disposition below is traced to the commit that made it.

## `grake` — slot 7, removed

**Disposition:** Specified, partially built, removed before any release. Never
shipped. Slot 7 is permanently reserved and MUST NOT be reused.

**Evidence:** commit `9a67891` (2026-04-30), *"remove(grake): drop
RK_ALG_GRAKE=7 enum value — no release, no ABI constraint"*. Commit body:
"Pre-release codebase; no ABI compatibility obligation. Remove the deprecated
sentinel entirely rather than leaving dead enum values." The commit touches
`src/leafblower.h` (one deletion) and `.beads/issues.jsonl` only.

**Superseded claims.** These specs treat `grake` as live and are overridden by
this record:

- `2026-04-25-calibration-solvers-design.md` — §3 `#define RK_ALG_GRAKE 7`;
  Goal 4 `method='grake'`; §7 algorithm; §11 acceptance criterion A4
  (match `survey::calibrate(epsilon=1e-10)` within 1%); files `src/grake.cpp`,
  `src/grake.hpp`. **All withdrawn.** No acceptance criterion for `grake`
  carries forward; `src/grake.cpp` and `src/grake.hpp` do not exist.
- `2026-04-27-sinkhorn-a1-fix-ieppa-admm-method.md` — `case RK_ALG_GRAKE:` in
  `select_solver_objective`, `src/grake.hpp` in files-changed. **Withdrawn.**
- `2026-04-29-greenkhorn-solver.md` — the R lookup
  `alg_names <- c("", "ieppa", "lbfgsb", "raking", "sinkhorn", "chebyshev",
  "greg", "grake", "ieppa_soft", "greenkhorn", "logit")` pins `"grake"` at
  index 7 and `"lbfgsb"` at index 2. **Both entries are stale.** The live
  lookup has holes at 2 and 7.
- `2026-04-29-chebyshev-greg-fix.md` — Out of Scope: "grake separate fix
  (inherits chebyshev improvements automatically)". **Moot.**

**Do not confuse the solver with the metric.** `grake_norm` is a live
convergence metric and is NOT affected by this removal. It is computed and
returned today at `src/greg.cpp:153`, `src/logit_calib.cpp:558`, and packed
into results at `src/c_api.cpp:74`, `:110`, `:426`. Grepping `grake` in the
tree returns these live hits; they are the metric, not the removed solver.

## `lbfgsb` — slot 2, removed

**Disposition:** Removed. Slot 2 is permanently reserved and MUST NOT be
reused.

**Evidence:** annotated in the enum itself at `src/leafblower.h:44` —
`/* 2 = removed (was RK_ALG_LBFGSB) */`. Documentation drift from this removal
was purged in commit `7fa211c`, *"docs(cr-g): purge L-BFGS-B fiction; fix
convergence/param/header doc drift"*. A related header comment at
`src/leafblower.h:65` records "(L-BFGS-B reference removed)".

**Superseded claims.** `2026-04-29-greenkhorn-solver.md`'s `alg_names` lookup
lists `"lbfgsb"` at index 2 (see above). **Withdrawn.**

## `cp` — slot 12, never landed

**Disposition:** Ported, then reverted the same day. `RK_ALG_CP` is not in the
enum and slot 12 was never occupied.

**Evidence:** commit `00a3f10` (2026-05-03) *"feat(cp): K-1 port
research/cp_calib to src/ + wire harvest(method='cp') (Epic-K)"*, reverted by
`3fac1d6` (2026-05-03). `grep -c 'RK_ALG_CP' src/leafblower.h` returns 0.

**Superseded claims.** `2026-05-02-epic-k-cp-productionization-design.md`
specifies `method="cp"` / `RK_ALG_CP = 12` in full. **Not live.** Treat the
spec as a withdrawn proposal, not as pending work, unless Epic-K is
deliberately revived — in which case slot 12 becomes available and this record
should be amended.

## Live enum (authoritative)

As of 2026-08-14, `src/leafblower.h:40-53`:

| Value | Name | Status |
|---|---|---|
| 0 | `RK_ALG_AUTO` | live |
| 1 | `RK_ALG_ORIS` | live (renamed from `ieppa`) |
| 2 | — | **reserved** (was `RK_ALG_LBFGSB`) |
| 3 | `RK_ALG_RAKING` | live |
| 4 | `RK_ALG_SINKHORN` | live |
| 5 | `RK_ALG_CHEBYSHEV` | live |
| 6 | `RK_ALG_GREG` | live |
| 7 | — | **reserved** (was `RK_ALG_GRAKE`) |
| 8 | `RK_ALG_ORIS_SOFT` | live |
| 9 | `RK_ALG_GREENKHORN` | live |
| 10 | `RK_ALG_LOGIT` | live |
| 11 | `RK_ALG_NEWTON_KL` | live |

Eight solvers, matching `docs/methods/00-overview.md`. Enum values are frozen
(commit `77d0614`, *"enum values frozen"*); slots 2 and 7 stay holes.

Note: one spec in the corpus gives `RK_ALG_RAKING` as 2. That is wrong —
raking is 3, and 2 is the reserved `lbfgsb` hole.
