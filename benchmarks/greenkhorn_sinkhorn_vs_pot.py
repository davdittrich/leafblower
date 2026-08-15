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

# ============================================================================
# Task 2: run leafblower's greenkhorn/sinkhorn vs POT's greenkhorn/sinkhorn
# ============================================================================
import importlib.metadata  # noqa: E402
import platform  # noqa: E402
import time  # noqa: E402
from leafblower import design_effect, effective_sample_size  # noqa: E402

CAVEAT = (
    "POT has no bounds mechanism; both sides run effectively unbounded "
    "(max_weight set beyond what the fixture can reach) to isolate the pure "
    "marginal-fit computation -- NOT a test of leafblower's normal bounded, "
    "K>2 workload."
)

SCHEMA_COLUMNS = [
    "input_class", "n", "n_margins", "n_categories", "m_cell", "m_cell_over_n",
    "max_weight", "arm", "wall_s", "max_error", "max_w", "min_w", "deff", "n_eff",
    "iterations", "ok", "note",
]

n_categories = len(LEVELS)
K_prior_arr = K_prior.to_numpy(dtype=np.float64)
m_cell = int((K_prior_arr > 0).sum())
tgt_m1_arr = tgt_m1.to_numpy(dtype=np.float64)
tgt_m2_arr = tgt_m2.to_numpy(dtype=np.float64)
m1_codes = df["m1"].cat.codes.to_numpy()
m2_codes = df["m2"].cat.codes.to_numpy()


def timed_median(fn, iterations=2):
    """Python analogue of bench::mark(iterations=2)'s median."""
    times = []
    result = None
    for _ in range(iterations):
        t0 = time.perf_counter()
        result = fn()
        times.append(time.perf_counter() - t0)
    return result, float(np.median(times))


def margin_max_error_2d(X_cell, tgt_m1_v, tgt_m2_v):
    """Shared, arm-independent accuracy metric applied identically to both
    engines' recovered 4x4 cell-mass table. Never reads a package's
    self-reported convergence status (mirrors oris_soft_vs_competitors.R's
    margin_max_error())."""
    total = X_cell.sum()
    row_props = X_cell.sum(axis=1) / total
    col_props = X_cell.sum(axis=0) / total
    return float(max(np.max(np.abs(row_props - tgt_m1_v)),
                      np.max(np.abs(col_props - tgt_m2_v))))


def cell_table_from_weights(w):
    """leafblower side: aggregate the returned per-observation weight vector
    into the same 4x4 cell shape via a weighted crosstab."""
    ct = pd.crosstab(df["m1"], df["m2"], values=pd.Series(w, index=df.index),
                      aggfunc="sum")
    return ct.reindex(index=LEVELS, columns=LEVELS, fill_value=0.0).to_numpy(dtype=np.float64)


def implied_weights_from_cell_table(T_cell):
    """POT side has no per-observation identity, only a solved 4x4 cell-mass
    table. Under the K=2 IPF equivalence this script exists to demonstrate,
    observations sharing a cell are exchangeable for this margin-only
    objective, so each observation in cell (i,j) is assigned the implied
    weight T[i,j] / K_prior[i,j] -- the same cell-uniform-weight structure
    leafblower's own (unbounded) unit-mode water-fill produces. K_prior[i,j]
    is guaranteed > 0 for every (i,j) actually occupied by an observation
    (the fixture's degeneracy guard above); POT's kernel is exactly 0 on any
    K_prior[i,j] == 0 cell (exp(-M/reg) = 0 there), so T is 0 there too and
    no such observation exists to index."""
    ratio = np.where(K_prior_arr > 0, T_cell / np.where(K_prior_arr > 0, K_prior_arr, 1.0), 0.0)
    return ratio[m1_codes, m2_codes]


def lb_call(method):
    return harvest(df, targets, method=method, max_weight=float(N),
                   bounds_mode="unit", convergence=None, attach_weights=False)


def pot_call(method):
    if method == "greenkhorn":
        return ot.bregman.greenkhorn(a, b, M, reg=1.0)
    return ot.sinkhorn(a, b, M, reg=1.0, method="sinkhorn")


rows = []

print(f"\n=== k2_margin_pot_equiv n={N} K=2 nj={n_categories} m_cell={m_cell} "
      f"(m_cell/n={m_cell / N:.4f}) ===")

for method in ("greenkhorn", "sinkhorn"):
    res, wall_s = timed_median(lambda m=method: lb_call(m))
    w = res["weights"]
    result_dict = res["result"]
    status = result_dict["status"]
    iterations = result_dict["iterations"]
    X_cell = cell_table_from_weights(w)
    max_error = margin_max_error_2d(X_cell, tgt_m1_arr, tgt_m2_arr)
    ok = bool(status in (0, 5))
    note = (
        f"{CAVEAT} leafblower method='{method}', max_weight={N} (effectively "
        f"unbounded); status={status}, iterations={iterations}"
    )
    rows.append({
        "input_class": "k2_margin_pot_equiv", "n": N, "n_margins": 2,
        "n_categories": n_categories, "m_cell": m_cell, "m_cell_over_n": m_cell / N,
        "max_weight": float(N), "arm": f"leafblower_{method}", "wall_s": wall_s,
        "max_error": max_error, "max_w": float(np.max(w)), "min_w": float(np.min(w)),
        "deff": design_effect(w), "n_eff": effective_sample_size(w),
        "iterations": iterations, "ok": ok, "note": note,
    })
    print(f"  {'leafblower_' + method:<22} wall={wall_s:7.4f}s status={status} "
          f"max_err={max_error:.3e} max_w={np.max(w):.3f} n_eff={effective_sample_size(w):.1f}")

for method in ("greenkhorn", "sinkhorn"):
    T, wall_s = timed_median(lambda m=method: pot_call(m))
    max_error = margin_max_error_2d(T, tgt_m1_arr, tgt_m2_arr)
    w_implied = implied_weights_from_cell_table(T)
    ok = bool(np.all(np.isfinite(T)))
    pot_fn = "ot.bregman.greenkhorn(a,b,M,reg=1.0)" if method == "greenkhorn" \
        else "ot.sinkhorn(a,b,M,reg=1.0,method='sinkhorn')"
    note = (
        f"{CAVEAT} {pot_fn}; max_w/min_w/deff/n_eff derived from the implied "
        f"per-observation weight T[i,j]/K_prior[i,j] within each cell (K=2 IPF "
        f"equivalence: observations sharing a cell are exchangeable under this "
        f"margin-only objective)"
    )
    rows.append({
        "input_class": "k2_margin_pot_equiv", "n": N, "n_margins": 2,
        "n_categories": n_categories, "m_cell": m_cell, "m_cell_over_n": m_cell / N,
        "max_weight": np.nan, "arm": f"pot_{method}", "wall_s": wall_s,
        "max_error": max_error, "max_w": float(np.max(w_implied)),
        "min_w": float(np.min(w_implied)), "deff": design_effect(w_implied),
        "n_eff": effective_sample_size(w_implied), "iterations": np.nan, "ok": ok,
        "note": note,
    })
    print(f"  {'pot_' + method:<22} wall={wall_s:7.4f}s max_err={max_error:.3e} "
          f"max_w={np.max(w_implied):.3f} n_eff={effective_sample_size(w_implied):.1f}")

results = pd.DataFrame(rows, columns=SCHEMA_COLUMNS)
os.makedirs("benchmarks/results", exist_ok=True)
results.to_csv("benchmarks/results/greenkhorn_sinkhorn_vs_pot.csv", index=False)
print("\nWrote benchmarks/results/greenkhorn_sinkhorn_vs_pot.csv")

# --- Machine and provenance capture (mirrors oris_soft_vs_competitors.R) ---
def _cpu_model():
    try:
        if os.path.exists("/proc/cpuinfo"):
            with open("/proc/cpuinfo") as f:
                for line in f:
                    if line.startswith("model name"):
                        return line.split(":", 1)[1].strip()
        return None
    except OSError:
        return None


cpu_model = _cpu_model()
env_lines = [
    f"Python version: {platform.python_version()}",
    f"Platform: {platform.platform()}",
    f"POT: {ot.__version__}",
    f"leafblower: {importlib.metadata.version('leafblower')}",
    f"OMP_NUM_THREADS: {os.environ.get('OMP_NUM_THREADS', '')}",
    f"OPENBLAS_NUM_THREADS: {os.environ.get('OPENBLAS_NUM_THREADS', '')}",
    f"MKL_NUM_THREADS: {os.environ.get('MKL_NUM_THREADS', '')}",
    f"CPU model: {cpu_model if cpu_model else 'NA (non-Linux host)'}",
]
with open("benchmarks/results/greenkhorn_sinkhorn_vs_pot_env.txt", "w") as f:
    f.write("\n".join(env_lines) + "\n")
print("Wrote benchmarks/results/greenkhorn_sinkhorn_vs_pot_env.txt")
