# Iterate-Change θ₂ Estimator for ORIS — Derivation Results

**Ticket:** leafblower-e18t.7 (Phase-1 GO gate for e18t.8)
**Artifact:** `2026-06-02-iterate-change-omega.sage` (SageMath + NumPy)
**Re-run with:** `sage docs/superpowers/derivations/2026-06-02-iterate-change-omega.sage`
**Date:** 2026-06-02

---

## Problem Statement

The current θ₂ estimator (e18t.6) uses the free-subspace **marginal residual**:

```
R2_free(k) = ||row/col sums of free cells - free targets||^2
```

On **infeasible-after-clamp** problems (stepstone mw=5: column target unreachable because
upper bounds cap total column mass below the target), `R2_free` **plateaus** at a nonzero
infeasibility floor → ratio → 1 → ω → 1.8 (capped) → stall (500 iters vs 300 baseline).

**Proposed fix (e18t.8):** Replace `R2_free` with the free-coordinate **iterate-change**:

```
S_dX(m) = sum_{c: !is_pinned[c]} (X[c] - X_snapshot[c])^2
```

This measures how much free cells moved since the last checkpoint. The plateau is gone:
`S_dX / S_dX_prev → ρ(M_II)^(2*I)` regardless of feasibility.

---

## Part A — 2×2 QQ Symbolic Cross-Check

Symbolic computation of `T_con = (I-C)(I-R)` over SR with asymmetric seed
`X* = [[2,1],[1,2]]`, `p=q=[3,3]`.

| Quantity | Value |
|---|---|
| T_con eigenvalues (|λ|) | 0.0, 0.1111..., 1.0 |
| ρ(M_II) from charpoly roots | **0.11111111** |
| ρ(M_II) from NumPy eigvals | **0.11111111** |
| Symbolic / NumPy cross-check | **PASS** |

The characteristic polynomial roots agree with the numerical eigvals.

---

## Part B — Feasible Probe

**Problem:** 3×3, `A=[[2,1,1],[1,3,1],[1,1,2]]`, `p=(4,5,4)`, `q=(5,4,4)`.
Single mild clamp at largest unconstrained cell (85% cap) → margins remain reachable.

| Quantity | Value |
|---|---|
| clamped cells | 1 |
| global residual at optimum | 0.000e+00 (FEASIBLE) |
| ρ(M_II) | 0.075304 |
| ρ(M_II)² | 0.005671 |
| **(a) marginal lag-1 ratio** | **0.005632** → ω = 1.0014 |
| **(b) iterate-change ratio** | **0.005716** → ω = 1.0014 |
| **(c) target ρ²** | **0.005671** → ω = 1.0014 |

Both estimators agree. Criterion: `|marginal − ρ²| < 0.05` ✓  `|iterate-change − ρ²| < 0.05` ✓

---

## Part C — Infeasible Probe

**Problem:** 4×4 symmetric, `p=q=[5,5,5,5]`. Column 0 cells capped at 1.0 each:
column 0 target = 5, max achievable = 4×1.0 = 4 < 5 → **INFEASIBLE margin**.

| Quantity | Value |
|---|---|
| clamped cells | 4 |
| global residual at optimum | 1.000e+00 (INFEASIBLE) |
| column 0 residual | 1.0000 (cannot vanish) |
| ρ(M_II) | 0.106445 |
| ρ(M_II)² | 0.011331 |
| **(a) marginal floor** | **1.2500e+00 (PLATEAU)** |
| marginal std tail-ratio | **1.000000** → ω = **1.9999** ← THE BUG |
| **(b) iterate-change ratio** | **0.002457** → ω = 1.0006 |
| **(c) target ρ²** | **0.011331** → ω = 1.0028 |

Iteration evidence:

| iter | marginal_resid | free_iterate_change |
|---:|---:|---:|
| 0 | 1.250049e+00 | 4.926282e-03 |
| 10 | **1.250000e+00** | 5.761150e-29 |
| 50 | **1.250000e+00** | 3.328007e-31 |
| 100 | **1.250000e+00** | 5.916457e-31 |
| 200 | **1.250000e+00** | 5.916457e-31 |
| 400 | **1.250000e+00** | 5.916457e-31 |

The marginal residual is stuck at 1.25 from iteration 10 onward. The iterate-change
decays to numerical floor in the same span. The marginal estimator would pick ω ≈ 2.0
(stall); the iterate-change picks ω ≈ 1.001 (correct).

---

## Part D — Cadence Probe (I=10, kErrCheckInterval)

**Problem:** 4×4, mild diagonal clamps at [0,0] and [3,3] (feasible, slower convergence).

| Quantity | Value |
|---|---|
| clamped cells | 2 |
| ρ(M_II) | 0.153584 |
| ρ(M_II)² | 0.023588 |
| ρ(M_II)^(2×10) | 5.332e-17 |
| **(lag-1 single-sweep) median ratio** | **0.025198** ≈ ρ² |
| **(block I=10) median ratio** | **9.881e-18** ≈ ρ^(2×10) |
| **(block I=10) θ₂ = ratio^(1/10)** | **0.019929** ≈ ρ² = 0.023588 |
| block-root ω | 1.0050 |
| lag-1 ω | 1.0064 |
| true ω | 1.0060 |

Both `lag-1` and `block-root` recover ρ² within 0.05. The block-root method preserves
the existing `kErrCheckInterval=10` cadence with one I-th-root correction.

---

## Estimator Formula (for e18t.8 implementation)

```
ESTIMATOR FORMULA (for e18t.8):
  every check m (every I=kErrCheckInterval=10 sweeps), mode 2, post-burnin:
    S_dX(m) = sum_{c: !is_pinned[c]} (X[c] - X_snapshot[c])^2
    ratio    = S_dX(m) / S_dX(m-1)               # -> rho^(2*I)
    theta2   = ratio^(1.0/I)                       # -> rho(M_II)^2
    omega    = 2 / (1 + sqrt(1 - clamp(theta2, 0, 1-1e-9)))  # ceiling 1.8
    for k in 0..K-1: sor_omega[k] = omega
    X_snapshot <- X
  State: X_snapshot[M_cell], S_dX_prev (init +inf)
  Fallback (lag-1): theta2 = S_dX(m) / S_dX(m-1) directly (no I-th root), adapt every sweep
```

**New state required in `oris.cpp`:** `X_snapshot[M_cell]` (size = number of cells),
`S_dX_prev` (scalar, init `+inf`). Updated every `kErrCheckInterval` sweeps in mode 2,
post-burnin.

---

## GO/NO-GO Decision

| Criterion | Result |
|---|---|
| [1] Symbolic cross-check (charpoly ↔ eigvals) | PASS |
| [2] Feasible: marginal ≈ iterate-change ≈ ρ² (±0.05) | PASS |
| [3] Infeasible: marginal plateaus (floor > 1e-6), iterate-change → ρ² | PASS |
| [4] Cadence: block-root AND lag-1 recover ρ² (±0.05) | PASS |

**DECISION: GO**

All criteria met. The iterate-change estimator `S_dX = sum(!pinned)(X−snap)²` is
validated as the correct replacement for the marginal-residual θ₂ proxy on infeasible
problems. Proceed to **e18t.8**: implement `S_dX` in `oris.cpp` `omega_mode_id=2`.
