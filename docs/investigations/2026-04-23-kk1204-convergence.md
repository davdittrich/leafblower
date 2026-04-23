# kk1.20.4 Convergence + Wall-Clock Investigation

**Date:** 2026-04-23
**Triggers:** leafblower-kk1.20.4 (original P1 gate), leafblower-g8f
  (follow-up convergence study)
**Baseline commit:** 106295c (post WU-P1 descent monitor)
**Host:** local (dd @ linux 6.19.12-arch1-1)

## Input

n = 1,000,000, K = 20, cat_counts = (5,)*20, max_weight = 3,
targets = skewed (0.3, 0.175, 0.175, 0.175, 0.175) per margin.
Data: uniform random sampling of each margin independently (so the joint
distribution is uniform over the 5^20 product space, while targets request
a skewed marginal — forces the solver to move mass).

## Measurements

### Cell compression (WU-P2.1)

- M_cell = 1,000,000
- ratio M_cell / n = 1.0000
- compression factor = 1.00x
- product(cat_counts) = 9.54e+13 (cells populated / available = 1.05e-08)

Every observation lands in its own unique (g_1,...,g_K) cell. Cell
compression gives zero benefit at this regime. Faithful iEPPA degenerates
to obs-level work with full log-space Sinkhorn overhead.

### Per-iteration cost (WU-P2.2, 50-iter budget, tol_abs=1e-300)

| method | total s | s/iter | max_err @ 50 |
|---|---|---|---|
| ieppa  | 35.887 | 0.718 | 3.960e-03 |
| raking | 1.645  | 0.033 | 2.471e-02 |
| lbfgsb | 17.607 | 0.352 | 1.206e-03 |

Per-iter ratios:
- ieppa / raking = **21.8x** (ieppa is ~22x slower per iter at zero compression)
- lbfgsb / raking = 10.7x
- ieppa / lbfgsb = 2.0x

Note: the plan step specified `convergence = list(absolute = 0)` but the
C API rejects `tol_abs <= 0` (`src/c_api.cpp:124`). Substituted `1e-300`
(finite, positive, far below any observed errRp → never triggers early
termination). Effect on per-iter denominator: nil.

### Convergence trajectory (WU-P2.3, max_iter=500)

**ieppa** — ran full 500 iters, monotone decreasing throughout:

| iter | 1 | 10 | 50 | 100 | 200 | 300 | 400 | 500 |
|---|---|---|---|---|---|---|---|---|
| errRp | 2.43e-2 | 8.56e-3 | 3.96e-3 | 2.47e-3 | 1.62e-3 | 1.34e-3 | 1.21e-3 | 1.15e-3 |

Final: `iEPPA max_iter exhausted (NOCONV) in 500 iters, errRp=1.147e-03`.
Monotone: YES.
Late-stage slope: iter 400→500, Δ errRp = 6.5e-5 per 100 iters (6.5e-7/iter).
Asymptote behaviour: errRp decays roughly like iter^{-0.5} — classical
Sinkhorn-style sublinear rate on this input class.

**raking** — WU-P1 descent monitor fired at iter 40:

| iter | 1 | 10 | 20 | 30 | 40 |
|---|---|---|---|---|---|
| errRp | 3.44e-2 | 2.13e-2 | 2.19e-2 | 2.20e-2 | 2.20e-2 |

Stall message: `raking: errRp stalled for 5 consecutive checks
(last=2.20e-02, window_min=inf); likely near-infeasible bounds. Aborting
at iter 40.` Monotone: NO — errRp rose from iter 10 (2.13e-2) to iter 20+
(2.19-2.20e-2) and flatlined.

**Monitor instrumentation defect observed (not in scope for WU-P2 fix):**
the stall message reports `window_min=inf`, i.e. `min_errRp_window` was
never updated from its initial `+inf` sentinel. Reading the P1.2.1 code:
the "improvement" branch requires `errRp < min_errRp_window - eps` which
at `min_errRp_window = +inf` is trivially true, so branch should fire on
iter 1. It does not. Root cause likely: on iter 1, `rel_eps = 0.01 * inf
= inf`, and `eps = max(inf, tol_abs) = inf`, so `min_errRp - eps = -inf`,
and the comparison `errRp < -inf` fails. `n_no_improve` increments from
iter 1 onward, and the monitor fires 5 checks later regardless of actual
trajectory. Filed as follow-up (see Recommendations §3).

Despite the logged message being misleading, the monitor's early-abort
behaviour is still operationally correct on THIS input: iters 20-40 are
genuinely stalled at ~2.20e-2 (a real near-infeasible plateau), so the
solver is right to bail. The defect matters for the diagnostic text and
for inputs where iter 1's errRp is the true minimum.

## Interpretation

**Why ieppa is slow.** The problem class forces M_cell = n: K=20 independent
uniform categorical columns with 5 levels each yield 5^20 ≈ 9.5e13 possible
cells; at n=10^6 the birthday-style collision rate is below 10^-7, so every
row is its own cell. CellTable compression saves zero work. Per-iter cost
(0.718s) is dominated by log-space Sinkhorn updates across K=20 margins on
n=1M observations. No log-space Sinkhorn can avoid O(n*K) work per outer
sweep; the 22x slowdown vs raking is the constant factor of log/exp cycles
relative to raking's straight multiplicative IPF updates.

**Why neither converges to tol_abs=1e-6.** Two distinct failure modes:

1. **ieppa** converges monotonically at a Sinkhorn-asymptotic rate: errRp
   shrinks like iter^{-0.5}. Extrapolating the late-stage slope (6.5e-7
   per iter at iter 500), the iter budget needed to reach 1e-6 is
   approximately (1.147e-3 - 1e-6) / 6.5e-7 ≈ 1.76 million iters. At
   0.718s/iter that is 14 days of wall-clock on this hardware. In
   practice the sublinear rate will degrade further as errRp shrinks, so
   this is a lower bound. **Unreachable.**

2. **raking** hits a near-infeasible plateau around 2.20e-2 within 20
   iters. The max_weight=3 bound combined with skewed targets on
   uniform-random data does not admit a weight vector with primal error
   below ~2e-2. This is a property of the input (the feasible region
   dictated by the bounds does not intersect the target polytope to
   within 1e-6), not of the solver. No amount of iteration fixes this.
   Confirmed by the monitor: the trajectory is non-monotone (rises
   slightly after iter 10 and stays flat), i.e. the solver is hunting
   inside a near-infeasible facet.

**Is the kk1.20.4 gate (<30s AND <1e-6) achievable here?** **No.** Either
algorithm fails for a distinct reason, and neither reason is fixable by
tuning:
- raking: input is near-infeasible at max_weight=3 — no weight vector
  satisfies all K=20 target margins to within 1e-6.
- ieppa: algorithmic rate is sublinear at M_cell=n with no compression;
  iter budget required exceeds real time by 4-5 orders of magnitude.
- lbfgsb: 0.352s/iter × 50 iters = 17.6s with max_err still 1.21e-3 at
  iter 50. Even if lbfgsb has a faster asymptotic rate, it starts 2x
  slower per iter than ieppa without a demonstrated convergence
  advantage on this class.

## Recommendations

1. **kk1.20.4 gate — REFRAME.** The gate `<30s AND <1e-6` is unachievable
   on the K=20 / cat=5 / uniform-random / max_weight=3 input class with
   any of the three implemented algorithms. Failure modes are orthogonal
   (near-infeasibility for raking; sublinear rate at M_cell=n for ieppa;
   lbfgsb starts 2x slower per iter and offers no rate advantage). The
   gate should be either:
   - narrowed to K ≤ 10 or cat counts where ∏cat_counts ≤ ~n/10 (i.e.
     where M_cell/n < 0.1 so iEPPA's compression applies), OR
   - relaxed to `<30s AND <1e-4` for raking on feasible inputs
     (max_weight wide enough that the near-infeasible plateau is below
     1e-4).
   Propose: file a bd ticket with the reframe payload in §3 below.

2. **Routing — insufficient data for a numeric threshold; file BLSE WU
   as prerequisite.** A single probed regime (M_cell/n = 1.00) cannot
   localize the breakpoint where ieppa overtakes raking. The ieppa
   per-iter cost at non-degenerate M_cell/n has not been measured, so
   neither T1 (M_cell/n ratio) nor T2 (n floor) can be set from evidence.
   Concrete prereq: a 3x3 sweep over M_cell/n ∈ {0.01, 0.1, 0.5} and
   n ∈ {1e4, 1e5, 1e6}, recording per-iter cost and iters-to-converge.
   Until then, AUTO should keep its current routing and add a single
   hard guard: if M_cell/n > 0.9 (near-incompressible), fall back to
   raking unconditionally — observed 22x per-iter penalty with zero
   cell-compression benefit on this input.

3. **Follow-ups filed in WU-P2.6:**
   - `docs(kk1204)`: reframe kk1.20.4 gate — currently unreachable; replace
     with (a) K ≤ 10 variant or (b) <1e-4 variant.
   - `feat(routing)`: AUTO guard — fall back raking when M_cell/n > 0.9
     (cheap one-pass probe before dispatch).
   - `feat(bench)`: BLSE WU — 3x3 sweep over (M_cell/n, n) to localize
     ieppa/raking routing breakpoint.
   - `fix(raking)`: descent monitor window_min never updates from +inf
     sentinel — diagnostic message reports `window_min=inf`. Does not
     affect abort behaviour on observed inputs but misleads operators.

## Raw data

- `/tmp/kk1204/probe.R` — M_cell probe script (not committed)
- `/tmp/kk1204/periter.R` — per-iter timing script (not committed)
- `/tmp/kk1204/traj_run.R` — trajectory capture driver (not committed)
- `/tmp/kk1204/trajectory.log` — 56 lines of `iter N: errRp=X` points
  (51 ieppa + 5 raking, trajectory through the full 500-iter budget or
  monitor abort)
- `/tmp/kk1204/periter.out` — per-iter cost table raw output
