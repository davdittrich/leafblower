# ORIS Iterate-Change θ₂ for Box-Constrained Over-Relaxation — Design (v2)

**Date:** 2026-06-02 (rev 2 — plan-review-gate iter 1: cross-language default/Python wiring, §3.7)
**Epic:** leafblower-e18t (reopened) — Phase-2 v2 remediation
**Predecessor:** e18t Phase-2 closed **NO-SHIP** — the free-subspace **marginal-residual** θ₂
estimator (e18t.3) stalled on bounded stepstone mw=5 (500-iter budget vs 140 fixed).
This is mj1p.2's failure mode reborn.

## Root cause (3-way verified)

The spectral theory `ω_opt(I) = 2/(1+√(1−ρ(M_II)²))`, `M_II` = free-block principal
submatrix, is **correct** (Lehmann 2022). The runtime **estimator** chose the wrong
observable.

1. **Algebra.** e18t.6's `R2_free = Σ_j (freeColSum_j − qf_j)²` with `qf_j = q_j −
   clamped_mass_j` is **identically** the global marginal residual `Σ_j (S_lin_j − q_j)²`
   — the clamped mass cancels into the free target. So `oris.cpp`'s `errF` *is* the global
   residual. (This is the "S_lin_free ≡ S_lin" identity the e18t.3 code-quality reviewer
   flagged BLOCKING; it was waved through — it was the bug.)
2. **SageMath** (probes in `/tmp/e18t_rootcause/`, to be formalized + committed by e18t.7).
   On an infeasible-after-clamp margin (column target unreachable), the marginal residual
   **plateaus at a nonzero floor from iteration 0** (lag-1 ratio → 1), while the free
   iterate-change `‖ΔX_free‖²` **decays to 0** at the true free-block rate. theta2→1 drives
   `ω = 2/(1+√0) = 2.0` (capped 1.8) = maximum over-relaxation on an already-converged
   block → stall.
3. **NotebookLM** (Lehmann SOR notebook `1e3036a1`, 191 sources). Confirmed: the iterate
   increment `Δx_k = x_k − x_{k−1}` forms a Krylov sequence under the iteration matrix, so
   its norm ratio is a **power iteration** converging to ρ(M_II) **regardless of global
   feasibility**. Lehmann/Soma-Uschmajew use residual ratios *only because their
   unconstrained problems are strictly feasible*. Classical SOR (Hageman & Young 1981)
   uses increment-based estimation precisely to avoid dependence on the residual offset.

**Why e18t.6's GO was false:** its synthetic verification used a single-clamp,
**feasible-after-clamp** toy where the marginal residual → 0 at the true rate. It never
exercised the bounds-active-infeasible regime (stepstone mw=5) — the exact regime the epic
exists to fix.

## Mechanism

> **Mechanism:** Free-coordinate iterate-change power-iteration estimator.
> Replace the mode-2 marginal-residual θ₂ with the lag-recovered ratio of
> `‖ΔX_free‖²` (free cells only). **Single global ω**, applied to all margins.

> **Forbidden:** marginal/constraint-residual ratios as the θ₂ observable; per-margin ω
> (use one global dominant-eigenvalue estimate); ω ceiling > 1.8; ω ≥ 2; touching the SRAA
> path, ALM/`oris_soft`, or the solver fixed point. No new INPUT ABI field; reuse
> `omega_mode_id=2`.

> **Audit:** SageMath artifact derives the observable and the cadence recovery and gates the
> implementation (GO). Orchestrator independently re-runs every ship-gate benchmark number.

## 1. The estimator

Asymptotically `e_free(k) = X_free(k) − X*_free ≈ ρ(M_II)^k v`. The increment
`ΔX_free(k) = X_free(k) − X_free(k−1)` obeys the **same** linear recurrence, so
`‖ΔX_free(k)‖² / ‖ΔX_free(k−1)‖² → ρ(M_II)²` — a power iteration that needs **no fixed
point** and is **feasibility-agnostic**.

### 1.1 Cadence recovery (block-root, primary)

`oris.cpp` adapts ω inside the `iter % kErrCheckInterval == 0` gate (`I = 10`). Measuring
`ΔX_free` over an interval of `I` sweeps gives a block ratio `→ ρ^(2I)`; recover the
per-sweep rate by the `I`-th root:

```
every check m (every I sweeps), in the flat path, mode 2, post-burnin:
  S_dX(m)  = Σ_{c : !is_pinned[c]} (X[c] − X_snapshot[c])²       # free-cell change since last check
  ratio    = S_dX(m) / S_dX(m−1)                                  # → ρ^(2I)
  theta2   = ratio ^ (1.0 / I)                                    # → ρ(M_II)²
  ω_global = 2 / (1 + √(1 − clamp(theta2, 0, 1−1e-9)))            # ceiling kSorProdCeiling = 1.8
  for k in 0..K-1: sor_omega[k] = ω_global
  X_snapshot ← X
```

This corrects a **latent e18t.3 cadence bug**: e18t.3 treated `errF_k/errF_prev` (a lag-`I`
ratio = ρ^(2I)) as ρ² directly, systematically under-reading ρ. The `^(1/I)` root fixes it.

SageMath (`verify_cadence`): block-root recovers ρ²=0.0236 as 0.0199 (ω 1.0050 vs true
1.0060) — slightly conservative (errs toward ω=1, never toward oscillation = **safe**).

### 1.2 Lag-1 per-sweep (fallback, pre-registered)

If the §4 ship gate shows block-root is too conservative and loses the T2 win, switch to
single-sweep lag-1 measurement (adapt ω every iteration, one extra O(M_cell) pass):
`theta2 = ‖ΔX_free(k)‖² / ‖ΔX_free(k−1)‖²` directly → ρ². SageMath: 0.0252 vs 0.0236 target
(ω 1.0064) — more accurate. This is a **bounded** fallback, decided by data at e18t.9, not a
free hand.

## 2. State (solve-local, scalar)

Collapses e18t.3's per-`k` arrays to scalars (the estimate is one global ρ):

| name | type | role |
| :--- | :--- | :--- |
| `X_snapshot` | `vector<double>(M_cell)` | free-cell state at last check (init = X at first post-burnin check) |
| `S_dX_prev` | `double` (init `+inf`) | lag-1 block iterate-change |
| `sor_theta2_ema` | `double` (init `-1`) | EMA of theta2 (reuse `kSorEmaAlpha=0.2`) |
| `sor_consec_up` | `int` | consecutive `S_dX` increases (gate 10) |
| `sor_cooldown_left` | `int` | remaining cooldown sweeps (gate 10) |
| `sor_cooldown_trips` | `int` | total trips this solve (gate 10 latch) |
| `sor_latched` | `bool` | permanent ω=1 latch |
| `is_pinned` | reuse e18t.3 `vector<char>(M_cell)` | active set |

No `CalibSorCfg` field. No per-`k` vectors for the v2 path. `is_pinned` clamp-loop logic
(linear `oris.cpp:1393`, log path) unchanged from e18t.3.

## 3. Guard envelope (10 gates, global)

Applied at each post-burnin check, in order, before writing `ω_global`:

1. **Scope:** `alm_active` → block not reached (mode-2 oris flat only).
2. **Free-set:** `n_free == 0` → `ω=1`, reset EMA. (`n_free` = Σ `!is_pinned[c]`.)
3. **Warm-up:** `iter < sor_burnin` OR `S_dX_prev` not yet set → `ω=1`.
4. **Finiteness/mass:** `!isfinite(S_dX)` OR `S_dX_prev < kResidFloor (1e-12)` OR free mass
   `< kMinSafeTotalWeight (1e-100)` → `ω=1`, reset EMA.
5. **(folded into 4)** — single denominator is `S_dX_prev`.
6. **Increment-grew:** `ratio ≥ 1` → hard-reset `sor_theta2_ema=0`, `ω=1`, `sor_consec_up++`.
7. **EMA:** `theta2 = clamp(ratio^(1/I), 0, 1−1e-9)`; `ema = α·theta2 + (1−α)·ema`.
8. **Formula + ceiling:** `ω_global = omega_from_theta2(ema, kSorProdCeiling=1.8)`.
9. **Oscillation damp:** on `S_dX` trend sign-flip (was-decreasing → now-increasing),
   `ω_global = max(omega_min, ω_global · kSorOscillationDamp(0.7))`; advance lag + reset EMA.
10. **Monotone latch:** `sor_consec_up ≥ 3` → `ω=1`, `cooldown_left = kSorCooldown(5)`,
    `cooldown_trips++`, reset EMA; `cooldown_trips ≥ kSorLatchTrips(3)` → `sor_latched=true`
    (permanent ω=1). During cooldown: `ω=1`, decrement, reset EMA.

Lag update at end of each check: `S_dX_prev ← S_dX`.

## 3.7 Cross-language default + Python wiring (rev 2 — plan-review-gate iter 1)

**Pre-existing divergence (surfaced by the gate).** `omega_mode_id`'s C default is **2**
(`src/types.hpp:74`, `src/c_api.cpp:114`). R `parse_sor` (`R/harvest.R`) defaults to **1L**
since e18t.5's NO-SHIP. The Python binding is asymmetric: `_parse_sor` returns a 6-tuple
(`enabled,auto,omega_init,omega_min,omega_fixed,burnin`) — **`omega_mode_id` is dropped** —
and `_bindings.cpp` never forwards it. So a Python `harvest(..., sor=list(auto=TRUE))`
silently runs the C default (mode-2) while the same R call runs mode-1. e18t.4's Python
mode-2 test passes only because the C default coincides with the dict value it (inertly)
passes. This is a real bug e18t.5 left; v2 must not inherit it.

**Two corrections, both in scope (no new ABI — `sor_omega_mode_id` already exists in
`rk_params_t`; this wires an existing C field that R already forwards):**

1. **Wire `omega_mode_id` through Python (e18t.10).** `_parse_sor` 6-tuple → 7-tuple
   including `omega_mode_id` (passthrough the dict key, default to the C default);
   `_bindings.cpp` forwards it to `p.sor_omega_mode_id`. After this, Python can select any
   mode explicitly, exactly like R — required for deterministic cross-language parity.
2. **Unify the default decision (e18t.9).** The §4 ship decision flips `omega_mode_id`'s
   default in **all three sites together**: `R/harvest.R` `parse_sor`, `src/types.hpp:74`,
   `src/c_api.cpp:114`. SHIP → 2 everywhere; NO-SHIP → 1 everywhere. No site is left behind.

## 4. Ship gate (re-run; reuse e18t fixtures)

SHIP only if **ALL** pass (TIE = NO-SHIP), per-fixture numbers re-verified by orchestrator:

1. **Unconstrained win:** iterate-change (mode 2) iters `<` fixed (mode 1) on the T2
   slow-unconstrained fixture (`oris_shipgate_fixture.R`, seed 20260531). Pre-registered
   fixed baseline = 350. (e18t.5 marginal mode-2 = 280; v2 must stay `< 350`.)
2. **No bounded regression:** mode-2 iters on stepstone mw=5 `≤ 300` (pre-registered fixed
   baseline). This is the condition e18t.3 **failed** (500). The fix must pass it.
3. **No NOCONV flip:** no fixture in `oris_fixed_omega_baseline.json` (e18t.2) that converged
   (status 0/5) under fixed-ω flips to NOCONV (status 1) under v2.
4. **Wall-clock not regressed** (iters proxy acceptable).

On SHIP → `omega_mode_id` default = 2 in **all three sites** (§3.7). On NO-SHIP → default = 1
everywhere, v2 opt-in (selectable in R *and* Python after e18t.10), commit the documented
negative. Either way the corrected derivation stands.

**Python parity (resolves gate iter-1 Completeness blocks).** e18t.8's "Python parity" check
= R(`omega_mode_id=2L`) vs Python weights at rtol=1e-6 (during e18t.8 the C default is still
2, so Python's SOR-enabled path hits mode-2; the e18t.4 test path is reused). After e18t.10,
the parity test forces mode-2 **explicitly on both sides**, making it deterministic regardless
of the eventual default — required because e18t.9 may flip the C default to 1 (NO-SHIP), after
which only explicit selection can exercise mode-2 in Python.

## 5. Observability (reuse e18t.4)

`sor_omega_mean`, `sor_any_latched`, `sor_n_pinned_fb`, `sor_n_warmup_fb`, `sor_n_conv_fb`,
`sor_n_resid_grew` (rename intent → increment-grew), `sor_n_monotone_cd` — already in the
result struct, wired R+Python. v2 increments them at the global gate sites. No new ABI.

## 6. Tickets (reopen e18t)

- **e18t.7** — Corrected SageMath derivation: iterate-change observable + cadence recovery +
  feasible/infeasible verification + NotebookLM corroboration. GO/NO-GO gate. (Probes exist
  in `/tmp/e18t_rootcause/`; formalize + commit under `docs/superpowers/derivations/`.)
- **e18t.8** — Implement: replace mode-2 marginal estimator with global iterate-change in
  `oris.cpp` flat path. `R CMD INSTALL` clean; `devtools::test()` 0 fail; Python parity
  rtol=1e-6 on the T2 mode-2 fixture (C default still 2 at this point — Python SOR-enabled
  hits mode-2; reuse the e18t.4 test path).
- **e18t.10** — Wire `omega_mode_id` through the Python binding (§3.7): `_parse_sor`
  6-tuple→7-tuple, `_bindings.cpp` forwards `p.sor_omega_mode_id`. No new ABI (existing C
  field). Adds a Python explicit-mode parity test. Depends on e18t.8.
- **e18t.9** — Re-run ship gate (§4), decide `omega_mode_id` default, **flip it in all three
  sites** (`R/harvest.R`, `src/types.hpp`, `src/c_api.cpp`), finalize docs/tests/NEWS in one
  commit cycle. Depends on e18t.10. Orchestrator re-verifies all numbers.

## 7. Alternatives considered

- **Renormalize the marginal residual by free mass** — rejected: still a residual, still
  plateaus on infeasible margins (the floor is in the constraint, not the normalization).
- **Per-margin ω from per-margin iterate-change** — rejected: SageMath shows all per-margin
  ratios converge to the *same* dominant ρ (the free block has one dominant eigenvalue); a
  single global estimate is the correct power-iteration object and simpler.
- **Generic LCP solver (Siconos/PATH)** — rejected: wrong layer; changes the fixed point.
- **Keep mode-2 marginal, document NO-SHIP permanently** — rejected by user; the root cause
  is understood and the fix is a single well-motivated observable swap.

## 8. Out of scope

SRAA adapt path; `oris_soft`/ALM (`alm_active`); any change to the solver fixed point, the
clamp/active-set rule, or the convergence criterion; new INPUT ABI.
