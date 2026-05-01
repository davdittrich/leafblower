# Newton-KL IEPPA Warm-Start — Result

**Date:** 2026-05-02
**Epic:** `leafblower-usg8` (Epic-C, BLOCKED follow-up to Epic-B target homotopy BLOCKED)
**Plan:** `docs/superpowers/plans/2026-05-01-newton-kl-ieppa-warmstart-plan.md` (rev 3, gate-approved iter 2)
**Verdict:** **BLOCKED** — warm-start regresses T2 stepstone K=9 from 2.79e-4 (cold) to 4.39e-4 (any K_warm ≥ 1). Hypothesis falsified empirically.

## What landed on master

| Commit | Change | Status |
|---|---|---|
| `32fcee6` | WI-0: spec amendment — "IEPPA warm-start" subsection | Shipped (will need erratum note) |
| `ecd3dec` | WI-0b: basin-overlap kill-switch script + CSV | Shipped |
| `af22726` | WI-1: `lf_capture` parameter on `ieppa_solve` + RAII guard + K_warm sweep | Shipped (harmless additive API; default-nullptr keeps existing callers bit-identical) |

## What did NOT land

WI-2 (warm-start wiring), WI-1c (r_bridge surfacing), WI-3 (tests), WI-4 (verify), WI-5a (bench), WI-5b (verdict). The WI-2 worktree was discarded after empirical regression measured.

## Empirical finding

K_warm sweep on stepstone K=9 fixture (n=200000), worktree commit `7b8310e`-equivalent:

| K_warm | max_err   | iters | Verdict             |
|--------|-----------|-------|---------------------|
| 0 (cold)  | 2.79e-4 | 7  | master baseline preserved |
| 1      | 4.39e-4   | 1     | regression          |
| 2      | 4.39e-4   | 1     | regression          |
| 4      | 4.39e-4   | 1     | regression          |
| 8      | 4.39e-4   | 1     | regression          |

K_warm ≥ 1: warm-started Newton runs exactly 1 iteration before best-iterate-restore freezes λ at the warm-start value. Newton's first step would push λ further (rank-deficient drift), so best-iterate restores to the warm-start. Result: warm-start max_err exactly matches IEPPA's K_warm-iter plateau (4.39e-4 — see WI-1 K_warm sweep).

## Why warm-start fails

The kill-switch (WI-0b) verdict was PROCEED with `max_log_ratio = 1.50` — IEPPA and cold-start Newton land on demonstrably distinct basins. Plan §Mechanism interpreted this as evidence that IEPPA-warm-started Newton would converge to a *better* basin than cold start.

Empirically, the opposite happens. The two basins are:

1. **Cold-start Newton basin** (master, K_warm=0): λ trajectory `0 → λ*_cold` with quadratic descent through gap=0.03 → 0.001 → 0.00014 (basin floor) → drift up; best-iterate restores to gap=0.00014 → max_err=2.79e-4.

2. **IEPPA basin** (K_warm=8): λ trajectory `0 →_IEPPA λ*_ieppa` where λ*_ieppa achieves max_err=4.39e-4. Newton from λ*_ieppa: gradient near zero in IEPPA-coordinate-descent direction but NOT zero in Newton's full-Newton direction. First Newton step amplifies the orthogonal component, λ drifts, best-iterate restores back to λ*_ieppa. max_err=4.39e-4.

The plan's mechanism rewrite assumed `λ*_ieppa` was *inside* Newton's quadratic basin. The empirical finding: `λ*_ieppa` is NOT inside Newton's quadratic basin — it's at a different fixed point (IEPPA fixed point) which is *outside* Newton's quadratic-convergence radius for this problem. Newton's first step from λ*_ieppa goes to a worse point, immediately rejecting the warm-start.

## Why the plan-review-gate did not catch this

Plan rev 3 passed gate-iter 2 (Feasibility/Completeness/Scope all PASS). The Feasibility reviewer verified:
- File paths and line refs.
- Conversion formula `λ_{k,j} = lf_{k,j} - lf_{k,0}` and the NA-bucket distinction.
- All 9 tickets exist with correct dep chain.

What the reviewer COULD NOT verify (and that no static analysis can verify) was whether the IEPPA fixed point sits inside Newton's quadratic basin. That is an empirical property of the dual landscape on a specific fixture. The kill-switch (WI-0b) was designed as the empirical falsification gate; its `max_log_ratio=1.5` outcome was misinterpreted as evidence FOR the hypothesis when it is actually a NECESSARY-but-not-sufficient condition (basins differ → warm-start *might* land Newton differently; says nothing about *whether the new landing zone is better*).

A more discriminating kill-switch would have run a small-K_warm warm-started Newton and checked whether the result improves over cold-start. That experiment requires the WI-1 + WI-2 wiring, which is exactly the work that BLOCKED the epic.

## What this means for the gate

Stepstone K=9 max_err remains at 2.79e-4 (master, T2 fails by 13%). On kk1204 K=20 severe-skew, master Newton-KL still drifts (status=1, gap saturated) and IEPPA-alone is the only working method (per the spec table 3.7s/2.4e-14 — though the WI-0b numbers suggest the spec table may have used different fixture parameters; this discrepancy was raised as F-R4 in plan rev 2 review).

The 1e-4 gate on stepstone is unreachable from any Newton-from-warm-start configuration tested. Reaching <1e-4 requires either:
1. **Different starting strategy** — e.g., random restart with lowest-gap-pick across multiple seeds (avoids basin-floor lock-in).
2. **Different solver** — e.g., trust-region Newton with explicit rank-deficient-direction projection.
3. **Tightening the gate** — accept 2.8e-4 as the achievable floor and amend T2.
4. **Routing K=9 overlapping-margin problems to a different method** — IEPPA-alone via AUTO routing change (per plan §F-R4: IEPPA misses gate by 13% at 1.13e-4 ON THE OTHER FIXTURE; on the WI-0b fixture, IEPPA reaches 4.39e-4 which is the same as warm-start, so this option doesn't help here either).

None of (1)–(4) is in scope for Epic-C. Filing follow-ups is conditional on user direction.

## Carryover scope

WI-1 (`lf_capture` parameter on `ieppa_solve`) is on master and is harmless. The `n_warmstart_iters_used` field on `NewtonCalibResult` (added in earlier work for Epic-B) remains on master — never set by master code (warm-start wiring not landed), always reads as 0. The `n_homotopy_levels_used` field is the still-stale name; rename was scheduled for WI-2 but did not land.

Suggest a small cleanup ticket if user wants:
- Either delete the unused `n_warmstart_iters_used` / `n_homotopy_levels_used` field (no consumers).
- Or note in NEWS.md that the field is deprecated / experimental.

## Files of record

- `src/newton_calib.cpp` (master HEAD `af22726` — TH-1a refactor + LM/LSE/best-iterate; no warm-start)
- `src/newton_calib.hpp` (master HEAD has `lm_mu_final`, `n_homotopy_levels_used` — stale field)
- `src/ieppa.cpp` / `src/ieppa.hpp` (master HEAD `af22726` — `ieppa_solve(state, lf_capture=nullptr)` shim)
- `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md` (TH-0 + WI-0 spec amendments landed; WI-0 amendment overstates warm-start utility — needs erratum)
- `docs/superpowers/plans/2026-05-01-newton-kl-ieppa-warmstart-plan.md` (Epic-C plan rev 3)
- `benchmarks/basin_overlap_kill_switch.R` + `benchmarks/results/basin_overlap_killswitch.csv` (WI-0b)
- `benchmarks/k_warm_sweep.R` + `benchmarks/results/k_warm_sweep.csv` (WI-1)
- Beads: `leafblower-usg8` (epic, closing BLOCKED); `usg8.1`/`.2`/`.8` (WI-0/WI-1/WI-0b shipped); `usg8.3`-`.7`/`.9` (closed BLOCKED).

## Lessons

1. The WI-0b kill-switch's `max_log_ratio` test is necessary but not sufficient. A real warm-start utility test requires running both pipelines and comparing FINAL Newton output, not just IEPPA output vs cold-Newton output.
2. The Epic-B + Epic-C sequence both attacked the basin-of-attraction problem from a different direction (target homotopy / IEPPA warm-start). Both failed for the same underlying reason: stepstone's overlapping-margin Hessian creates a basin-floor structure where the LM-Newton can converge to gap=1.4e-4 from cold start, and ANY non-cold start lands λ at a fixed point that is NOT inside the LM-Newton quadratic basin.
3. The kk1204 severe-skew K=20 case is a different pathology (Newton-KL drifts unboundedly from cold; IEPPA-alone converges). On stepstone, Newton-KL has bounded behavior but a high basin floor; on kk1204, Newton-KL has unbounded behavior. These are TWO different problems and warm-start does not solve EITHER.
