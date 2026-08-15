#!/usr/bin/env python3
"""Honest Performance Gate -- leafblower_greenkhorn/leafblower_sinkhorn vs. POT

Phase: 03-honest-performance-gate, Plan 07 (G-03-1 gap closure, continued)

`docs/methods/greenkhorn.md` and `docs/methods/sinkhorn.md` both name Python
Optimal Transport (POT, [flamary2021pot]) as "the primary open-source
vehicle" for these two methods -- available as `ot.bregman.greenkhorn` /
`ot.sinkhorn`. POT solves the classical TWO-marginal transport problem
(a transport matrix P, cost matrix M, scalar regularization `reg`, NO bounds
argument); leafblower's greenkhorn/sinkhorn solve the K-margin (K >= 1)
box-bounded survey-calibration problem. Per `docs/methods/oris.md`'s
"Relation to Sinkhorn and Greenkhorn" section (Axis 1): "In the pure
2-marginal case Sinkhorn and cyclic IPF coincide." This script therefore
builds the ONE mathematically faithful comparison available: a dedicated
K=2-margin fixture, run through both engines.

WHY the M=-log(prior), reg=1 construction is faithful, not invented:
POT's entropic-OT objective uses kernel `exp(-M/reg)`. Setting
`M = -log(K_prior)` and `reg = 1` recovers `exp(-M/reg) = K_prior` exactly
-- i.e. POT's kernel becomes leafblower's own `X_init[c]` prior (the
observed cross-tabulation count matrix, leafblower's implicit prior when
`start_weights` is left at its Python default of `None`). POT's transport
plan then minimizes the SAME `min_X KL(X || X_init)` s.t. margins objective
leafblower's greenkhorn/sinkhorn solve for K=2 -- the entropic-OT/
KL-projection equivalence documented in Peyre & Cuturi (2019), already
cited as [peyrecuturi2019computational] in docs/methods/sinkhorn.md, and
the "Sinkhorn and cyclic IPF coincide" statement in docs/methods/oris.md
cited above.

Both arms therefore run effectively UNBOUNDED (POT cannot enforce a box;
leafblower's max_weight is set beyond what the fixture can reach) -- a
caveat stated explicitly in every CSV row's `note`, never silently
generalized to leafblower's normal bounded, K>2 workload.
"""
import os

# --- Determinism guard (CLAUDE.md single-thread BLAS protocol, ENFORCED) ---
# Must run BEFORE any numpy/scipy-backed import: numpy/scipy read BLAS-thread
# env vars at first use, and `import ot` transitively imports numpy. Mirrors
# benchmarks/oris_soft_vs_competitors.R's require_single_thread_blas() guard:
# refuse to run rather than silently measure under multi-thread BLAS. Does
# NOT call os.environ[...] = "1" itself -- the guard only has teeth if the
# caller's environment is what is checked.
def require_single_thread_blas() -> None:
    bad = [v for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
           if os.environ.get(v) != "1"]
    if bad:
        raise RuntimeError(
            "refusing to measure: {} must be set to \"1\" for a reproducible "
            "single-thread BLAS run (CLAUDE.md protocol); got {{{}}}".format(
                ", ".join(bad),
                ", ".join(f"{v}={os.environ.get(v) or '<unset>'}" for v in bad),
            )
        )


require_single_thread_blas()

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import ot  # noqa: E402
from leafblower import harvest  # noqa: E402

# --- Fixture: K=2-margin, n=10000, 4 levels per margin ---
# Skewed, non-uniform per-margin target (both columns) -- an exactly-uniform
# target would make both sides converge trivially and prove nothing,
# mirroring oris_soft_vs_competitors.R's own documented reasoning for why
# its medium fixture uses a skewed target rather than uniform.
SEED = 307
N = 10_000
LEVELS = ["a", "b", "c", "d"]

rng = np.random.default_rng(SEED)
df = pd.DataFrame({
    "m1": pd.Categorical(rng.choice(LEVELS, size=N), categories=LEVELS),
    "m2": pd.Categorical(rng.choice(LEVELS, size=N), categories=LEVELS),
})

p_skew = np.array([0.40, 0.30, 0.20, 0.10])
p_skew = p_skew / p_skew.sum()
tgt_m1 = pd.Series(p_skew, index=LEVELS)
tgt_m2 = pd.Series(p_skew, index=LEVELS)
targets = {
    "m1": dict(zip(LEVELS, tgt_m1)),
    "m2": dict(zip(LEVELS, tgt_m2)),
}

# Observed cross-tabulation COUNT matrix -- leafblower's implicit X_init[c]
# prior when start_weights is left at its Python default (None).
K_prior = pd.crosstab(df["m1"], df["m2"]).reindex(
    index=LEVELS, columns=LEVELS, fill_value=0
)

# Structurally-empty cells legitimately produce +inf here, matching
# leafblower's own documented "X_init[c] = 0 cells are skipped" semantics
# (docs/methods/oris.md "How leafblower deviates" table) -- not a numerical
# bug to suppress.
with np.errstate(divide="ignore"):
    M = -np.log(K_prior.to_numpy(dtype=np.float64))

# Degeneracy guard: every row/column of K_prior must have at least one
# nonzero entry, else the fixture (and its M=-log(prior) construction) is
# degenerate. n=10000 uniform draws over 4x4 makes this astronomically
# unlikely to fire; it is a defensive check, not a workaround for observed
# failures.
row_sums = K_prior.sum(axis=1)
col_sums = K_prior.sum(axis=0)
empty_rows = row_sums[row_sums == 0].index.tolist()
empty_cols = col_sums[col_sums == 0].index.tolist()
if empty_rows or empty_cols:
    raise ValueError(
        f"degenerate fixture: empty K_prior row(s) {empty_rows}, "
        f"empty col(s) {empty_cols} -- regenerate with a different seed"
    )

# a, b: both scaled to the SAME total mass N as the observed data, matching
# leafblower's own T_kj = target_proportion * n convention (read from
# benchmarks/oris_soft_vs_competitors.R).
a = tgt_m1.to_numpy(dtype=np.float64) * N
b = tgt_m2.to_numpy(dtype=np.float64) * N
