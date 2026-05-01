# Newton-KL Target Homotopy — Result

**Date:** 2026-05-01
**Epic:** `leafblower-91u7` (Epic-B; NEEDS_HOMOTOPY follow-up to epic 5k08)
**Plan:** `docs/superpowers/plans/2026-05-01-newton-kl-homotopy-plan.md` (rev 2)
**Verdict:** **BLOCKED** — target homotopy regressed T2 from 2.8e-4 (LM baseline) to 2.29e-3 (8× worse). Implementation NOT merged to master. Follow-up: `Epic-C` (IEPPA warm-start).

## What landed on master

| Commit | Change | Status |
|---|---|---|
| `0dde818` | TH-0: spec amendment — "Target homotopy" subsection | Shipped |
| `5355f69` | TH-1a: refactor iter loop into `run_newton_inner` lambda + `n_homotopy_levels_used` field | Shipped (no behavior change; T2 still 2.8e-4) |

## What did NOT land

TH-1b-impl (target homotopy outer ε loop, NEWS.md update) was implemented in a
worktree (commit `7b8310e` on `worktree-agent-ace1c893bdacf6b3b`, plus an
uncommitted fix to relax break-on-NOCONV). Tests confirmed an 8× regression on
T2 stepstone K=9 (max_err: 2.8e-4 baseline → 2.29e-3 with homotopy). The
worktree was discarded; master remains at TH-1a state.

## Why homotopy regressed

The plan rev 2 §Context (corrected from rev 1) framed homotopy as a
**warm-start basin trick** rather than a rank fix:

> "rank deficiency in `H` persists at every ε, but Newton starting near a warm
> start with smaller `||λ_init − λ*||` stays inside the basin where its
> quadratic model is locally accurate."

This argument **assumes** that `λ*(T_eps)` for small ε is close to `λ*(T_0)`.
For the stepstone fixture this assumption fails: per-margin uniform `T_eps`
on overlapping margins (e.g., `rk_age10_gender`, `rk_gender_time`,
`rk_i_loc_time_gender`) produces intermediate problems whose dual optimum is
**structurally** offset from the original problem's optimum — not just
quantitatively close. Warm-starting Newton from `λ*(T_eps→0+)` lands it in a
basin distinct from the basin around `λ*(T_0)`. The Newton step toward
`λ*(T_0)` from the wrong basin then over-shoots, producing the 2.29e-3 result
(8× worse than the 2.8e-4 we got starting from `λ=0`).

Diagnostic: with the relaxed break-on-NOCONV (continue past intermediate ε
failures), total iterations rose to 75 (out of 6 × 20 = 120 budget) and
`status=1` (NOCONV) at final ε=0.0. The chain ran, but Newton at final ε=0
hits the same `consecutive_failed=3` line-search-failure path that the LM-only
solver hits. The warm start did not help. It also did not hurt convergence
quality at ε=0 specifically — the 2.29e-3 result is what Newton reaches from
the warm-started λ before hitting the same drift wall.

## Key empirical numbers

| Configuration | T2 (stepstone K=9) max_err | T2 status |
|---|---|---|
| Master pre-homotopy (LM only, TH-1a refactor) | **2.8e-4** | 1 (NOCONV) |
| Worktree homotopy + strict break-on-failure | 2.29e-3 | 1 |
| Worktree homotopy + relaxed break-on-NOCONV | 2.29e-3 | 1 |

Identical result for strict vs relaxed break confirms: the regression is
*structural* (warm-start basin mismatch), not driven by intermediate-ε
status handling.

## Why Epic-C is the right follow-up

IEPPA (iterative proportional fitting / coordinate descent in log space) is
known to converge robustly on overlapping-margin fixtures (stepstone
included; that's the original ieppa+sraa baseline that hits 1.13e-4 on the
same problem). IEPPA's coordinate-by-coordinate updates cannot overshoot
across margin boundaries — the algorithm always lands `λ` inside the
correct basin.

Epic-C runs IEPPA for 5-10 iterations FIRST, hands off `λ` to the LM-Newton
inner. Newton then does the fast quadratic polishing from inside the basin.
This is the standard pattern for combining first-order coordinate methods
(robust, slow tail) with second-order methods (fast, basin-sensitive).

## Plan rev 2 reviewer concerns that proved relevant

The Feasibility reviewer (F1, prior plan-review-gate output) explicitly
flagged:

> "[T_eps blends targets but] does not change which obs share `(j_k, j_{k'})`
> patterns, so `S_{k1j1,k2j2}/Z` retains the same support holes."

Plan rev 2 acknowledged this and reframed the mechanism as a basin warm-start
trick — but the rev 2 mechanism still assumed the warm start would land in
the *same basin* as `λ*(T_0)`. Empirically, on stepstone, it does not.

Reviewer-Feasibility also flagged (F2):
> "Final 0.01 → 0 jump is the riskiest step and the plan acknowledges it
> least."

Rev 2 tightened to 6-level `{0.5, 0.1, 0.03, 0.01, 0.003, 0}`. Even with
`0.003 → 0` (333× ratio), the basin mismatch persists.

In hindsight, the rev 2 mechanism rewrite was correctly diagnosing the
*absence* of a rank fix, but it overestimated how reliably warm starts on
this problem class would land in the right basin. The reviewer's residual
"prove the mechanism, not just assert it" comment turned out to identify
the actual failure mode.

## Files of record

- `src/newton_calib.cpp` (master HEAD `5355f69` — TH-1a refactor; no homotopy)
- `src/newton_calib.hpp` (master HEAD has `lm_mu_final`, `n_homotopy_levels_used` fields)
- `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md` (TH-0 spec amendment landed)
- `docs/superpowers/plans/2026-05-01-newton-kl-homotopy-plan.md` (plan rev 2 — supersedes rev 1; documents target-homotopy mechanism)
- Tickets: `leafblower-91u7` Epic-B (closing BLOCKED); `leafblower-<TBD>` Epic-C (filed as follow-up)

## Outstanding scope (carried into Epic-C)

- Run IEPPA inner for 5-10 iterations at the original `T`; hand `λ` to LM-Newton.
- Test T2 (stepstone), T8 (K=20 severe-skew unit) PASSES at <1e-4 / <1e-3.
- Benchmark kk1204 severe-skew under IEPPA-warm-start + LM-Newton.
- AUTO routing safety guard only on Epic-C's ESCAPE verdict (not before).
- `lm_mu_final` and `n_homotopy_levels_used` R-side surfacing deferred until a method actually needs it (LM-1c can re-open if Epic-C reuses the homotopy field; otherwise the field stays unused).
