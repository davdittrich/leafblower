# Newton-KL LM-Damping Convergence Fix — Result

**Date:** 2026-05-01
**Epic:** `leafblower-5k08`
**Plan:** `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md` (LM section, amended via LM-0) + `docs/superpowers/plans/2026-05-01-newton-kl-lm-plan.md` (rev 2).
**Verdict:** **NEEDS_HOMOTOPY** — Newton-KL converges, but max_err on stepstone K=9 (T2) plateaus at 2.8e-4, missing the <1e-4 gate. Drafted follow-up: `leafblower-91u7` (Epic-B: target homotopy).

## What shipped

| Commit | Change |
|---|---|
| `8d0d68b` | LM-0: spec amendment — "Levenberg-Marquardt damping" subsection. |
| `88033d1` | LM-1: scale-invariant LM `(H + μ·diag(H))·δ = G` with Marquardt gain-ratio adaptive μ; three-way Armijo (FULL_ACCEPT / BACKTRACKED / FAILED) with 3-retry cap; `lm_mu_final` field on `NewtonCalibResult`. |
| `e4276cc` | LSE stabilization (the LM-1 ticket assumed it shipped, but it was prior-session uncommitted) + best-iterate fallback + `δᵀHδ` non-negative clamp. |

## Why NEEDS_HOMOTOPY (not GATE_MET)

Stepstone K=9 (T2) trace under final implementation:

| iter | gap | g | accepted | α | μ |
|---|---|---|---|---|---|
| 0 | 3.1e-2 | 1.22e1 | yes | 0.5 | 1e-3 |
| 1 | 1.5e-2 | 1.21e1 | yes | 1.0 | 1e-3 |
| 2 | 7.7e-3 | 1.21e1 | yes | 1.0 | 3.3e-4 |
| 3 | 2.8e-3 | 1.18e1 | yes | 1.0 | 1.1e-4 |
| 4 | 9.6e-4 | 1.11e1 | yes | 1.0 | 3.7e-5 |
| 5 | 3.0e-4 | 8.88 | yes | 1.0 | 1.2e-5 |
| **6** | **1.4e-4** | **1.85** | yes | 1.0 | 4.1e-6 |
| 7 | 1.6e-4 | -1.99e1 | yes | 1.0 | 1.4e-6 |
| 8+ | ~1.6e-4 | → -∞ | yes (drift) | 1.0 | shrinks |

**Iter 6 is the natural plateau.** From iter 7 onward, `g` keeps decreasing (Armijo accepts because the step *is* in a descent direction) but the dual gradient stops shrinking. The mechanism: rank-deficient `H_pre` (Schur-complemented covariance with many empty cells in stepstone's overlapping margins) gives Newton a direction along which `g(λ) → −∞`. The dual is unbounded along this ray *numerically* even though the primal calibration optimum is finite.

Best-iterate fallback rescues this: λ is restored to iter 6's value before recovery. Final `max_err = 2.8e-4` corresponds to gap ≈ 1.4e-4 with reference-category amplification (ref-cat error = sum of non-ref errors → factor of nj−1 ≈ 2 for stepstone's 12-cat margins).

## Why LM alone cannot reach gap < 1e-4

The plateau at gap ≈ 1.4e-4 is *intrinsic* to the dual landscape with the current sample / target combination. LM damping, three-way Armijo, gain-ratio scheduling, and best-iterate fallback together stabilize convergence and prevent the catastrophic max_err = 0.988 we saw under the bandaid baseline — but they do not change the convergence basin.

To break the plateau, the **landscape itself** must be deformed. Options:

1. **Target homotopy** (Epic-B): smooth targets toward uniform, anneal to original. Each intermediate problem has a strictly-feasible interior optimum where Newton converges fully. Final ε=0 step starts from a near-optimal warm start, avoiding the rank-deficient region entirely.
2. **IEPPA warm-start** (Epic-C, alternative): run 5–10 IEPPA iters first (coordinate descent, no overshoot risk), hand λ to Newton for polishing.

The plan's verdict tree picks (1) for NEEDS_HOMOTOPY since it's a self-contained Newton-KL extension with no cross-method coupling. (2) is a fallback if homotopy doesn't ship.

## Comparison vs prior baselines

| Implementation | T2 (stepstone K=9) max_err | T2 status |
|---|---|---|
| Plain Newton (no LSE, no LM) | NaN on K≥9 | crashes / status=0 |
| LSE only (in-session, never committed) | < 1e-4 (passed at one point) | 0 |
| LSE + bandaid (trust=1, pivot=1e-4) | 3.85e-4 | 1 |
| LSE + bandaid (trust=1, pivot=1e-12) | 3.11e-2 | 1 |
| LM-1 (this work, no LSE, no best-iter) | 9.88e-1 | 1 |
| **LM-1 + LSE + best-iter (final)** | **2.79e-4** | **1** |

The final implementation is the strongest yet on the convergence-quality side, but still 2.8× the gate. Direction of progress is correct; further reduction requires homotopy.

## Files of record

- `src/newton_calib.cpp` (final implementation)
- `src/newton_calib.hpp` (`lm_mu_final` field)
- `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md` (spec, amended)
- `docs/superpowers/plans/2026-05-01-newton-kl-lm-plan.md` (plan rev 2)
- Beads: `leafblower-5k08` (epic, closing as NEEDS_HOMOTOPY); `leafblower-91u7` (Epic-B follow-up).

## Outstanding scope (carried into Epic-B)

- LM-1c: surface `lm_mu_final` in `src/r_bridge.cpp` so R callers see it (small, can land before or after Epic-B).
- LM-1b: T6 / T6b μ-schedule unit tests — deferred since they need `lm_mu_final` exposed in R.
- AUTO routing safety guard: defer to Epic-B verdict (if Epic-B succeeds, no guard needed; if Epic-B fails, route K=20 severe to ieppa+sraa via target-skew detection).
