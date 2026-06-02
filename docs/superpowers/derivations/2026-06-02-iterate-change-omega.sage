#!/usr/bin/env sage
# -*- coding: utf-8 -*-
r"""
Iterate-change theta2 estimator for ORIS omega_mode_id=2 (spectral SOR).

Ticket: leafblower-e18t.7 — GO/NO-GO gate for the e18t.8 iterate-change fix.

PROBLEM STATEMENT
=================
The current theta2 estimator (e18t.6) uses the free-subspace MARGINAL RESIDUAL:
  R2_free(k) = ||row/col sums of free cells - free targets||^2

On INFEASIBLE-after-clamp problems (e.g. stepstone mw=5: column target unreachable
because upper bounds cap total column mass below the target), R2_free PLATEAUS at a
nonzero infeasibility floor → ratio → 1 → omega → 1.8 (capped) → stall.

PROPOSED FIX (e18t.8)
======================
Replace R2_free with the FREE-COORDINATE ITERATE-CHANGE:
  S_dX(m) = sum_{c: !is_pinned[c]} (X[c] - X_snapshot[c])^2

This quantity measures HOW MUCH the free cells moved since the last checkpoint.
At the constrained fixed point, free cells stop moving → S_dX → 0 geometrically.
The plateau is gone: S_dX / S_dX_prev → rho(M_II)^(2*I) regardless of feasibility.

WHAT THIS SCRIPT PROVES
========================
1. 2x2 QQ symbolic cross-check: extract rho(M_II) analytically from T_con = (I-C)(I-R).
2. Feasible probe: marginal ≈ iterate-change ≈ rho^2 (both work; baseline confirmation).
3. Infeasible probe: marginal PLATEAUS (ratio→1), iterate-change → rho^2 (fix works).
4. Cadence probe: block-I=10 and lag-1 both recover rho^2 via I-th root / direct.
5. GO/NO-GO decision.

Run with:
    sage docs/superpowers/derivations/2026-06-02-iterate-change-omega.sage
"""

import numpy as np
np.set_printoptions(precision=6, suppress=True, linewidth=140)

# =============================================================================
# PART A — 2x2 QQ SYMBOLIC CROSS-CHECK
# =============================================================================
# Build T_con = (I-C)(I-R) symbolically on SR over QQ for a 2x2 free block.
# Verify that the smaller root of the characteristic polynomial = rho(M_II).

print("=" * 78)
print("PART A — 2x2 QQ SYMBOLIC CROSS-CHECK")
print("=" * 78)

x11, x12, x21, x22, p1s, p2s, q1s, q2s = var(
    'x11 x12 x21 x22 p1 p2 q1 q2', domain='real')

# 2x2 free block (no clamped cells), free targets = full p, q.
# Free cells indexed (0,0)=x11, (0,1)=x12, (1,0)=x21, (1,1)=x22
# R[a,b] = X*[row_a, col_b] / p_row_a  if row_a == row_b
# C[a,b] = X*[row_a, col_b] / q_col_b  if col_b == col_a
R_sym = matrix(SR, [
    [x11/p1s, x12/p1s, 0,       0      ],
    [x11/p1s, x12/p1s, 0,       0      ],
    [0,       0,       x21/p2s, x22/p2s],
    [0,       0,       x21/p2s, x22/p2s],
])
C_sym = matrix(SR, [
    [x11/q1s, 0,       x21/q1s, 0      ],
    [0,       x12/q2s, 0,       x22/q2s],
    [x11/q1s, 0,       x21/q1s, 0      ],
    [0,       x12/q2s, 0,       x22/q2s],
])
I4 = identity_matrix(SR, 4)
T_sym = (I4 - C_sym) * (I4 - R_sym)

# Substitute concrete values: balanced 2x2, p=q=[3,3], X*=[[1.5,1.5],[1.5,1.5]]
# -> R and C are both rank-1 projection matrices; rho(T_con) = 0 (trivial case).
# Use an asymmetric case: X* = [[2,1],[1,2]], p=[3,3], q=[3,3]
subs = {x11:2, x12:1, x21:1, x22:2, p1s:3, p2s:3, q1s:3, q2s:3}
T_num_sym = T_sym.subs(subs)
cp = T_num_sym.charpoly()
rts = sorted([abs(complex(r).real) for r in cp.roots(ring=CC, multiplicities=False)])
print(f"T_con eigenvalues (|λ|, sorted): {[float(r) for r in rts]}")
# The characteristic polynomial has roots; pick the subdominant non-unit root.
non_unit = [r for r in rts if abs(r - 1.0) > 1e-6 and r > 1e-12]
rho_sym = max(non_unit) if non_unit else 0.0
# Cross-check via numpy eigvals
T_np = np.array([[float(T_num_sym[i,j]) for j in range(4)] for i in range(4)], dtype=float)
evs_np = sorted([abs(v) for v in np.linalg.eigvals(T_np)], reverse=True)
rho_np = next((v for v in evs_np if v < 1 - 1e-6), evs_np[-1])
print(f"Symbolic rho(M_II) from charpoly roots = {float(rho_sym):.8f}")
print(f"NumPy   rho(M_II) from eigvals         = {rho_np:.8f}")
match_sym = abs(float(rho_sym) - rho_np) < 1e-6
print(f"Symbolic / numpy cross-check: {'PASS' if match_sym else 'FAIL'}")
print()

# =============================================================================
# SHARED UTILITIES
# =============================================================================

def settle(A, p, q, L, U, iters=8000):
    """IPF to constrained optimum. Clamps held at BOUND U_ij (not precomputed)."""
    m, n = A.shape
    X = np.clip(A.copy(), L, U)
    cl = (X >= U - 1e-12) | (X <= L + 1e-12)
    for _ in range(iters):
        for i in range(m):
            f = ~cl[i]; frz = X[i, ~f].sum(); t = p[i] - frz; s = X[i, f].sum()
            if s > 0 and t > 0: X[i, f] *= t / s
        X = np.clip(X, L, U); cl = (X >= U - 1e-12) | (X <= L + 1e-12)
        for j in range(n):
            f = ~cl[:, j]; frz = X[~f, j].sum(); t = q[j] - frz; s = X[f, j].sum()
            if s > 0 and t > 0: X[f, j] *= t / s
        X = np.clip(X, L, U); cl = (X >= U - 1e-12) | (X <= L + 1e-12)
    return X, cl


def build_T_con(Xstar, pf, qf, free_idx):
    """Linearized free-block operator T_con = (I-C)(I-R) with free denominators."""
    F = len(free_idx)
    R = np.zeros((F, F)); C = np.zeros((F, F))
    for a, (i, j) in enumerate(free_idx):
        for b, (i2, j2) in enumerate(free_idx):
            if i == i2: R[a, b] = Xstar[i2, j2] / pf[i]
            if j == j2: C[a, b] = Xstar[i2, j2] / qf[j]
    I = np.eye(F)
    return (I - C) @ (I - R)


def rho_M_II(A, p, q, L, U):
    """Spectral radius of the constrained free-block operator."""
    Xs, cl = settle(A, p, q, L, U)
    fr = ~cl
    pf = p - np.where(cl, Xs, 0.).sum(1)
    qf = q - np.where(cl, Xs, 0.).sum(0)
    fi = [(i, j) for i in range(A.shape[0]) for j in range(A.shape[1]) if fr[i, j]]
    T = build_T_con(Xs, pf, qf, fi)
    evs = sorted([abs(v) for v in np.linalg.eigvals(T)], reverse=True)
    return next((v for v in evs if v < 1 - 1e-9), evs[-1])


def sweep1(X, cl, p, q, L, U):
    """One IPF row+col sweep in-place. Clamps at BOUND."""
    m, n = X.shape
    for i in range(m):
        f = ~cl[i]; frz = X[i, ~f].sum(); t = p[i] - frz; s = X[i, f].sum()
        if s > 0 and t > 0: X[i, f] *= t / s
    X[:] = np.clip(X, L, U); cl[:] = (X >= U - 1e-12) | (X <= L + 1e-12)
    for j in range(n):
        f = ~cl[:, j]; frz = X[~f, j].sum(); t = q[j] - frz; s = X[f, j].sum()
        if s > 0 and t > 0: X[f, j] *= t / s
    X[:] = np.clip(X, L, U); cl[:] = (X >= U - 1e-12) | (X <= L + 1e-12)


def tail_ratio(seq, lo=1e-22, hi=1e-2):
    """Median lag-1 ratio, excluding floor and large-amplitude entries."""
    rr = [seq[k+1] / seq[k]
          for k in range(len(seq) - 1)
          if lo < seq[k] < hi and seq[k+1] > 0]
    return float(np.median(rr)) if rr else float('nan')


def omega_of(t):
    t = min(max(t, 0.0), 1.0 - 1e-9)
    return 2.0 / (1.0 + (1.0 - t) ** 0.5)


# =============================================================================
# PART B — FEASIBLE PROBE
# =============================================================================
# Single mild clamp: free cells CAN reach the targets. Both marginal and
# iterate-change should converge to rho^2.

print("=" * 78)
print("PART B — FEASIBLE PROBE (single mild clamp, margins reachable)")
print("=" * 78)

A_F = np.array([[2., 1., 1.], [1., 3., 1.], [1., 1., 2.]])
p_F = np.array([4., 5., 4.]); q_F = np.array([5., 4., 4.])
# Unconstrained optimum; clamp largest cell to 85% -> still feasible
Xu = A_F.copy()
for _ in range(5000):
    Xu *= (p_F / Xu.sum(1))[:, None]; Xu *= (q_F / Xu.sum(0))[None, :]
ij_F = np.unravel_index(np.argmax(Xu), Xu.shape)
L_F = np.zeros_like(A_F); U_F = np.full_like(A_F, np.inf); U_F[ij_F] = 0.85 * Xu[ij_F]

rho_F = rho_M_II(A_F, p_F, q_F, L_F, U_F)
Xs_F, cl_F = settle(A_F, p_F, q_F, L_F, U_F)
fr_F = ~cl_F
global_resid_F = max(abs(Xs_F.sum(1) - p_F).max(), abs(Xs_F.sum(0) - q_F).max())
print(f"clamped cells: {int(cl_F.sum())} | global residual at optimum: {global_resid_F:.3e}")
print(f"FEASIBLE: {global_resid_F < 1e-6}")
print(f"rho(M_II) = {rho_F:.6f}   rho^2 = {rho_F**2:.6f}")

rng = np.random.default_rng(int(20260602))
X = Xs_F.copy(); X[fr_F] *= np.exp(0.02 * rng.standard_normal(int(fr_F.sum())))
X = np.clip(X, L_F, U_F); cl = cl_F.copy(); prev = X.copy()
marg_F, dX_F = [], []
for k in range(600):
    sweep1(X, cl, p_F, q_F, L_F, U_F)
    # marginal residual (free targets)
    pf = p_F - np.where(cl_F, X, 0.).sum(1); qf = q_F - np.where(cl_F, X, 0.).sum(0)
    frs = np.where(fr_F, X, 0.).sum(1); fcs = np.where(fr_F, X, 0.).sum(0)
    marg_F.append(float(((frs - pf)**2).sum() + ((fcs - qf)**2).sum()))
    # iterate-change on free cells
    d = X - prev; dX_F.append(float((d[fr_F]**2).sum())); prev = X.copy()

marg_F = np.array(marg_F); dX_F = np.array(dX_F)
mr_F = tail_ratio(marg_F); ic_F = tail_ratio(dX_F)
print(f"(a) marginal lag-1 ratio    = {mr_F:.6f}  -> omega = {omega_of(mr_F):.4f}")
print(f"(b) iterate-change ratio    = {ic_F:.6f}  -> omega = {omega_of(ic_F):.4f}")
print(f"(c) rho(M_II)^2 (target)    = {rho_F**2:.6f}  -> omega = {omega_of(rho_F**2):.4f}")
feasible_marginal_ratio = mr_F
feasible_iterate_ratio  = ic_F
rho_M_II_val = rho_F
print()

# =============================================================================
# PART C — INFEASIBLE PROBE
# =============================================================================
# Column 0 target = 5, upper bound 1.0/cell, max achievable = 4*1.0 = 4 < 5.
# Marginal residual PLATEAUS (ratio -> 1, the bug).
# Iterate-change decays geometrically at rho^2.

print("=" * 78)
print("PART C — INFEASIBLE PROBE (column 0 target unreachable; stepstone regime)")
print("=" * 78)

A_I = np.array([[2., 1., 1., 1.],
                [1., 2., 1., 1.],
                [1., 1., 2., 1.],
                [1., 1., 1., 2.]])
p_I = np.array([5., 5., 5., 5.]); q_I = np.array([5., 5., 5., 5.])
L_I = np.zeros_like(A_I); U_I = np.full_like(A_I, np.inf); U_I[:, 0] = 1.0

rho_I = rho_M_II(A_I, p_I, q_I, L_I, U_I)
Xs_I, cl_I = settle(A_I, p_I, q_I, L_I, U_I)
fr_I = ~cl_I
global_resid_I = max(abs(Xs_I.sum(1) - p_I).max(), abs(Xs_I.sum(0) - q_I).max())
print(f"clamped cells: {int(cl_I.sum())} | global residual at optimum: {global_resid_I:.3e}")
print(f"INFEASIBLE: {global_resid_I > 1e-6}")
print(f"rho(M_II) = {rho_I:.6f}   rho^2 = {rho_I**2:.6f}")
print(f"column 0 residual at optimum: {abs(Xs_I.sum(0) - q_I)[0]:.4f} (nonzero = infeasible margin)")

rng2 = np.random.default_rng(int(20260602))
X = Xs_I.copy(); X[fr_I] *= np.exp(0.02 * rng2.standard_normal(int(fr_I.sum())))
X = np.clip(X, L_I, U_I); cl = cl_I.copy(); prev = X.copy()
marg_I, dX_I = [], []
for k in range(600):
    sweep1(X, cl, p_I, q_I, L_I, U_I)
    pf = p_I - np.where(cl_I, X, 0.).sum(1); qf = q_I - np.where(cl_I, X, 0.).sum(0)
    frs = np.where(fr_I, X, 0.).sum(1); fcs = np.where(fr_I, X, 0.).sum(0)
    marg_I.append(float(((frs - pf)**2).sum() + ((fcs - qf)**2).sum()))
    d = X - prev; dX_I.append(float((d[fr_I]**2).sum())); prev = X.copy()

marg_I = np.array(marg_I); dX_I = np.array(dX_I)
# For marginal: filter above the plateau floor to show the plateau ratio
marg_floor = np.median(marg_I[-50:])
mr_I_raw = tail_ratio(marg_I, lo=marg_floor * 1.001, hi=marg_I[0])
# Also compute with standard filter to show it's near 1
mr_I_std = tail_ratio(marg_I, lo=1e-20, hi=marg_I[0])
ic_I = tail_ratio(dX_I)

mr_I_raw_omega = f"{omega_of(mr_I_raw):.4f}" if not (mr_I_raw != mr_I_raw) else "N/A (no valid pairs)"
print(f"(a) marginal floor:            {marg_floor:.4e}   (PLATEAU — nonzero infeasibility)")
print(f"    marginal tail-ratio (raw)  = {mr_I_raw}  -> omega = {mr_I_raw_omega}  [ratio at floor=NaN → no valid pairs]")
print(f"    marginal tail-ratio (std)  = {mr_I_std:.6f}  -> omega = {omega_of(mr_I_std):.4f}  [the bug: near 2.0]")
print(f"(b) iterate-change ratio       = {ic_I:.6f}  -> omega = {omega_of(ic_I):.4f}")
print(f"(c) rho(M_II)^2 (target)       = {rho_I**2:.6f}  -> omega = {omega_of(rho_I**2):.4f}")
infeasible_marginal_plateaus = (marg_floor > 1e-6)
infeasible_iterate_recovers  = (abs(ic_I - rho_I**2) < 0.05)
print()
print("Iteration table (marginal vs iterate-change):")
print("  iter   marginal_resid   free_iterate_change")
for i in [0, 10, 50, 100, 200, 400, 599]:
    print(f"  {i:4d}   {marg_I[i]:.6e}   {dX_I[i]:.6e}")
print()

# =============================================================================
# PART D — CADENCE PROBE (I=10, matching kErrCheckInterval)
# =============================================================================
# Show that measuring over I=10 sweeps gives block ratio -> rho^(2I).
# I-th root recovers rho^2. Lag-1 single-sweep also works.
#
# Uses a 4x4 problem with MILD diagonal clamps (feasible, slow convergence)
# so the iterate-change series has enough range to measure before hitting floor.
# The infeasible problem converges too fast (~10 sweeps) for block cadence.

print("=" * 78)
print("PART D — CADENCE PROBE (I=10 block vs lag-1 single-sweep; mild-clamp 4x4)")
print("=" * 78)

A_C = np.array([[2., 1., 1., 1.], [1., 2., 1., 1.],
                [1., 1., 2., 1.], [1., 1., 1., 2.]])
p_C = np.array([5., 5., 5., 5.]); q_C = np.array([5., 5., 5., 5.])
# Mild upper bounds on diagonal corners only -> feasible, slow convergence
L_C = np.zeros_like(A_C); U_C = np.full_like(A_C, np.inf)
U_C[0, 0] = 1.2; U_C[3, 3] = 1.2

rho_C = rho_M_II(A_C, p_C, q_C, L_C, U_C)
Xs_C, cl_C = settle(A_C, p_C, q_C, L_C, U_C)
fr_C = ~cl_C
print(f"clamped cells: {int(cl_C.sum())} | rho(M_II)={rho_C:.6f}  rho^2={rho_C**2:.6f}")

rng3 = np.random.default_rng(int(20260531))  # same seed as original cadence probe
X = Xs_C.copy(); X[fr_C] *= np.exp(0.02 * rng3.standard_normal(int(fr_C.sum())))
X = np.clip(X, L_C, U_C); cl = cl_C.copy()
I_int = 10; prev = X.copy(); snap = X.copy()
lag1_seq = []; block_seq = []
for k in range(1, 401):
    sweep1(X, cl, p_C, q_C, L_C, U_C)
    d = X - prev; lag1_seq.append(float((d[fr_C]**2).sum())); prev = X.copy()
    if k % I_int == 0:
        db = X - snap; block_seq.append(float((db[fr_C]**2).sum())); snap = X.copy()

lag1_seq = np.array(lag1_seq); block_seq = np.array(block_seq)
r_lag1  = tail_ratio(lag1_seq)
r_block = tail_ratio(block_seq, lo=1e-24, hi=block_seq[0] * 2.0)
theta2_block = r_block ** (1.0 / I_int)
rho2_target_C = rho_C ** 2

print(f"(lag-1 single-sweep) median ratio    = {r_lag1:.6f}  ≈ rho^2  (direct)")
print(f"(block I=10)         median ratio    = {r_block:.4e}  ≈ rho^(2*10)={rho_C**(2*I_int):.3e}")
print(f"(block I=10) theta2 = ratio^(1/10)  = {theta2_block:.6f}  ≈ rho^2 = {rho2_target_C:.6f}")
print(f"omega: lag1->{omega_of(r_lag1):.4f}  block-root->{omega_of(theta2_block):.4f}  "
      f"true->{omega_of(rho2_target_C):.4f}")
block_root_ok = abs(theta2_block - rho2_target_C) < 0.05
lag1_ok       = abs(r_lag1 - rho2_target_C) < 0.05
print(f"block-root recovers rho^2: {'YES' if block_root_ok else 'NO'}")
print(f"lag-1 recovers rho^2:      {'YES' if lag1_ok else 'NO'}")
print()

# =============================================================================
# PART E — ESTIMATOR FORMULA (for e18t.8)
# =============================================================================

print("=" * 78)
print("ESTIMATOR FORMULA (for e18t.8):")
print("  every check m (every I=kErrCheckInterval=10 sweeps), mode 2, post-burnin:")
print("    S_dX(m) = sum_{c: !is_pinned[c]} (X[c] - X_snapshot[c])^2")
print("    ratio    = S_dX(m) / S_dX(m-1)               # -> rho^(2*I)")
print("    theta2   = ratio^(1.0/I)                       # -> rho(M_II)^2")
print("    omega    = 2 / (1 + sqrt(1 - clamp(theta2, 0, 1-1e-9)))  # ceiling 1.8")
print("    for k in 0..K-1: sor_omega[k] = omega")
print("    X_snapshot <- X")
print("  State: X_snapshot[M_cell], S_dX_prev (init +inf)")
print("  Fallback (lag-1): theta2 = S_dX(m) / S_dX(m-1) directly (no I-th root), adapt every sweep")
print("=" * 78)
print()

# =============================================================================
# PART F — GO/NO-GO DECISION
# =============================================================================

print("=" * 78)
print("GO/NO-GO DECISION")
print("=" * 78)

crit1 = match_sym
crit2 = abs(mr_F - rho_F**2) < 0.05 and abs(ic_F - rho_F**2) < 0.05
crit3 = infeasible_marginal_plateaus and infeasible_iterate_recovers
crit4 = block_root_ok and lag1_ok

print(f"[1] Symbolic cross-check passes:              {'PASS' if crit1 else 'FAIL'}")
print(f"[2] Feasible: marginal ≈ iterate-change ≈ ρ²: {'PASS' if crit2 else 'FAIL'}")
print(f"    marginal={mr_F:.4f}  iterate={ic_F:.4f}  rho^2={rho_F**2:.4f}")
print(f"[3] Infeasible: marginal plateaus, iterate→ρ²:{'PASS' if crit3 else 'FAIL'}")
print(f"    marginal_floor={marg_floor:.2e}  iterate_ratio={ic_I:.4f}  rho^2={rho_I**2:.4f}")
print(f"[4] Cadence (block+lag-1) both recover ρ²:   {'PASS' if crit4 else 'FAIL'}")

all_pass = crit1 and crit2 and crit3 and crit4
decision = "GO" if all_pass else "NO-GO"
print()
print(f"DECISION: {decision}")
if all_pass:
    print("  All criteria met. The iterate-change estimator is validated.")
    print("  Proceed to e18t.8: implement S_dX in oris.cpp omega_mode_id=2.")
else:
    print("  One or more criteria failed. Do NOT proceed to e18t.8 until resolved.")
print("=" * 78)
