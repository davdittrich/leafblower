"""Cross-language weight parity: R and Python must return identical weights
for all methods when both call the same C++ solver with identical input."""

import json
import os
import shutil
import subprocess
from pathlib import Path

# Enforce single-thread BLAS BEFORE numpy / leafblower import so that
# OpenBLAS / MKL / OpenMP pick up the constraint at initialisation time.
# Setting these after numpy import has no effect on most BLAS implementations.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pandas as pd
import pytest

from leafblower import harvest

# python/leafblower/test_parity_weights.py -> repo root is three levels up
# (leafblower/ -> python/ -> repo root). Getting this wrong makes the
# R-helper paths below resolve to nothing, which without the assertions
# further down would silently degrade every test to pytest.skip instead
# of failing collection loudly.
REPO_ROOT = Path(__file__).resolve().parents[2]
R_HELPER  = REPO_ROOT / "tests" / "parity" / "run_parity_r.R"
assert R_HELPER.exists(), f"R_HELPER not found at {R_HELPER} — REPO_ROOT arithmetic is wrong"

RSCRIPT_AVAILABLE = shutil.which("Rscript") is not None

def _make_synthetic(tmp_path, seed=42, n=2000):
    rng = np.random.default_rng(seed)
    df = pd.DataFrame({
        "age":    rng.choice(["18-34","35-54","55-64","65+"], size=n,
                             p=[0.25,0.35,0.25,0.15]),
        "gender": rng.choice(["M","F"], size=n, p=[0.48,0.52]),
        "region": rng.choice(["north","south","east"], size=n,
                             p=[0.30,0.40,0.30]),
        "edu":    rng.choice(["low","mid","high"], size=n,
                             p=[0.25,0.50,0.25]),
    })
    tgt = {
        "age":    {"18-34": 0.30, "35-54": 0.32, "55-64": 0.22, "65+": 0.16},
        "gender": {"M": 0.50, "F": 0.50},
        "region": {"north": 0.28, "south": 0.42, "east": 0.30},
        "edu":    {"low": 0.20, "mid": 0.55, "high": 0.25},
    }
    data_csv    = tmp_path / "data.csv"
    targets_json = tmp_path / "targets.json"
    df.to_csv(data_csv, index=False)
    with open(targets_json, "w") as f:
        json.dump(tgt, f)
    return df, tgt, data_csv, targets_json


def _r_weights(data_csv, targets_json, method, out_csv, max_iter=1000):
    result = subprocess.run(
        ["Rscript", str(R_HELPER),
         str(data_csv), str(targets_json), method, str(out_csv), str(max_iter)],
        capture_output=True, text=True, timeout=120,
    )
    if result.returncode == 2:
        pytest.skip(f"R package not available: {result.stderr.strip()}")
    if result.returncode != 0:
        raise RuntimeError(f"Rscript failed:\n{result.stderr}")
    return pd.read_csv(out_csv)["weight"].to_numpy()


# All nine non-AUTO rk_algorithm_t solvers (leafblower.h) on the shared
# four-margin synthetic fixture at default bounds. oris_soft and chebyshev
# also have dedicated, more elaborate tests further below in this file
# (test_oris_soft_default_tol_parity asserts alm_capacity_mu_final and
# iteration counts; test_chebyshev_default_tol_parity asserts the warm-start
# floor mirror invariant) — these parametrize entries are additive plain
# weight-vector coverage, not duplicates of those targeted regression locks.
# Slots 2 (removed lbfgsb) and 7 (withdrawn Epic-K cp) are permanently
# reserved holes and intentionally absent; "auto" is not a solver.
@pytest.mark.skipif(not RSCRIPT_AVAILABLE, reason="Rscript not found")
@pytest.mark.parametrize("method", [
    "greenkhorn", "logit", "raking", "oris", "sinkhorn", "newton_kl",
    "chebyshev", "greg", "oris_soft",
])
def test_weight_parity(method, tmp_path):
    df, tgt, data_csv, targets_json = _make_synthetic(tmp_path)
    out_csv = tmp_path / f"r_weights_{method}.csv"

    w_py = np.array(harvest(
        df, tgt, method=method,
        min_weight=0, max_weight=5, max_iterations=1000,
        convergence={"metric": "max_err", "rule": "improvement", "tol": 1e-3},
        verbose=0, attach_weights=False,
    )["weights"], dtype=np.float64)

    w_r = _r_weights(data_csv, targets_json, method, out_csv)

    assert len(w_r) == len(w_py), (
        f"{method}: weight vector length mismatch — R={len(w_r)}, Python={len(w_py)}")
    diff = np.max(np.abs(w_py - w_r))
    print(f"\n  {method}: max|w_py - w_r| = {diff:.2e}")
    tol = 1e-6 if method == "logit" else 1e-10
    assert diff < tol, (
        f"{method}: max|w_py - w_r| = {diff:.2e} (threshold {tol:.0e})"
    )


# ---------------------------------------------------------------------------
# T5: oris_soft parity regression — locks T4 harmonization (ae09943)
# ---------------------------------------------------------------------------

_ORIS_SOFT_R_HELPER = REPO_ROOT / "tests" / "parity" / "run_oris_soft_r.R"
assert _ORIS_SOFT_R_HELPER.exists(), (
    f"_ORIS_SOFT_R_HELPER not found at {_ORIS_SOFT_R_HELPER} — REPO_ROOT arithmetic is wrong"
)

# Convergence params — identical in R helper and Python call below.
_ORIS_SOFT_CONV = {"metric": "max_err", "rule": "improvement", "tol": 1e-4}


def _make_correlated_fixture(tmp_path, n=5000, K=9, n_levels=5, n_clusters=50, seed=42):
    """Generate a clustered fixture that exposes the T4 capacity_mu divergence.

    K=9 margins with n_levels levels each. estimate_M_cell (K>8 fast-exit)
    returns n (capacity_mu=1.0). build_cell_table (pre-T4 R path) returns
    n_clusters (capacity_mu = n_clusters/n << 1). The mismatch is the T4
    red-green signal: assert alm_capacity_mu_final R == Python.

    Actual unique interaction cells = n_clusters, so M_cell/n is very small
    (e.g. 50/5000 = 0.01) while estimate_M_cell returns n=5000 → ratio=1.0.
    """
    rng = np.random.default_rng(seed)
    # n_clusters distinct row patterns, each repeated across n rows.
    cluster_patterns = rng.integers(0, n_levels, size=(n_clusters, K))
    cluster_ids      = rng.integers(0, n_clusters, size=n)
    cols = {
        f"v{k:02d}": [f"L{k:02d}_{v}" for v in cluster_patterns[cluster_ids, k]]
        for k in range(K)
    }
    df = pd.DataFrame(cols)

    # Uniform targets (one proportion per level; equal weight per category)
    tgt = {}
    for k in range(K):
        col  = f"v{k:02d}"
        lvs  = sorted(df[col].unique())
        tgt[col] = {lv: 1.0 / len(lvs) for lv in lvs}

    data_csv    = tmp_path / "corr_data.csv"
    targets_json = tmp_path / "corr_targets.json"
    df.to_csv(data_csv, index=False)
    with open(targets_json, "w") as f:
        json.dump(tgt, f)

    return df, tgt, data_csv, targets_json


@pytest.mark.skipif(not RSCRIPT_AVAILABLE, reason="Rscript not found")
def test_oris_soft_default_tol_parity(tmp_path):
    """R and Python must agree on capacity_mu and weights for oris_soft at tol=1e-4.

    Fixture: synthetic, K=9 margins with 50 unique interaction cells from n=5000
    rows. This triggers >= 20 outer iterations at tol=1e-4.

    The primary red-green signal is alm_capacity_mu_final:
    - Pre-T4 R: build_cell_table returns exact M_cell=50 → capacity_mu=50/5000=0.01
    - Post-T4 R: estimate_M_cell (K>8 fast-exit) returns n=5000 → capacity_mu=1.0
    - Python (always): estimate_M_cell → capacity_mu=1.0
    Pre-T4: R.alm_capacity_mu_final (0.01) != Python (1.0) → test FAILS.
    Post-T4: both use estimate_M_cell → same value → test PASSES.

    Secondary: weight parity max|w_R - w_Py| < 1e-6. Both converge to the same
    optimum regardless of capacity_mu scaling on this fixture, so this assertion
    also passes post-T4 (and coincidentally pre-T4 too — it is not the gate).
    """
    df, tgt, data_csv, targets_json = _make_correlated_fixture(tmp_path)
    out_csv = tmp_path / "oris_soft_r_out.csv"

    _single_thread_env = {
        **os.environ,
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
    result = subprocess.run(
        ["Rscript", str(_ORIS_SOFT_R_HELPER),
         str(data_csv), str(targets_json), str(out_csv)],
        capture_output=True, text=True, timeout=180,
        env=_single_thread_env,
    )
    if result.returncode == 2:
        pytest.skip(f"R package not available: {result.stderr.strip()}")
    if result.returncode != 0:
        raise RuntimeError(f"Rscript failed:\n{result.stderr}")

    r_out       = pd.read_csv(out_csv)
    w_r         = r_out["weight"].to_numpy(dtype=np.float64)
    iters_r     = int(r_out["iterations"].iloc[0])
    status_r    = int(r_out["status"].iloc[0])
    cap_mu_r    = float(r_out["alm_capacity_mu_final"].iloc[0])

    # --- Python in-proc ---
    py_res      = harvest(
        df, tgt,
        method         = "oris_soft",
        min_weight     = 0,
        max_weight     = 5,
        max_iterations = 1000,
        convergence    = _ORIS_SOFT_CONV,
        verbose        = 0,
        attach_weights = False,
    )
    w_py        = np.array(py_res["weights"], dtype=np.float64)
    iters_py    = int(py_res["result"]["iterations"])
    status_py   = int(py_res["result"]["status"])
    cap_mu_py   = float(py_res["result"].get("alm_capacity_mu_final", float("nan")))

    # --- Fixture sanity: must have >= 20 outer iters ---
    assert iters_py >= 20, (
        f"Fixture too easy — only {iters_py} Python outer iters at tol=1e-4; "
        "test is insensitive. Scale up fixture n or reduce n_clusters."
    )

    # --- Primary red-green gate: capacity_mu harmonization ---
    # Pre-T4 R uses build_cell_table → M_cell/n << 1.
    # Post-T4 R uses estimate_M_cell → same as Python (1.0 for K>8, M_cell >= n).
    cap_mu_diff = abs(cap_mu_r - cap_mu_py)
    print(
        f"\n  oris_soft: alm_capacity_mu_final R={cap_mu_r:.6f} Py={cap_mu_py:.6f} "
        f"diff={cap_mu_diff:.2e}  iters_R={iters_r}  iters_Py={iters_py}  status={status_r}"
    )
    assert cap_mu_diff < 1e-9, (
        f"oris_soft capacity_mu MISMATCH: R={cap_mu_r:.6f} Py={cap_mu_py:.6f} "
        f"diff={cap_mu_diff:.2e} — "
        "R r_bridge.cpp is using build_cell_table instead of estimate_M_cell. "
        "Apply T4 fix ae09943 (harmonize R to estimate_M_cell path)."
    )

    # --- Same status ---
    assert status_r == status_py, (
        f"oris_soft: status mismatch — R={status_r}, Python={status_py}"
    )

    # --- Same iter count ---
    assert iters_r == iters_py, (
        f"oris_soft: iter count mismatch — R={iters_r}, Python={iters_py} "
        f"(capacity_mu_R={cap_mu_r:.4f}; check r_bridge.cpp harmonization)"
    )

    # --- Weight vector parity ---
    assert len(w_r) == len(w_py), (
        f"oris_soft: weight vector length mismatch — R={len(w_r)}, Python={len(w_py)}"
    )
    max_abs_diff = np.max(np.abs(w_r - w_py))
    assert max_abs_diff < 1e-6, (
        f"oris_soft parity FAIL: max|w_R - w_Py| = {max_abs_diff:.2e} >= 1e-6 "
        f"(iters_R={iters_r}, iters_Py={iters_py}, cap_mu_diff={cap_mu_diff:.2e}; "
        "check capacity_mu harmonization in r_bridge.cpp — see T4 fix ae09943)"
    )


# ---------------------------------------------------------------------------
# T5 (2apm): chebyshev parity regression — locks c_api.cpp:359 floor=5 (b60a3fd)
# ---------------------------------------------------------------------------

_CHEBYSHEV_R_HELPER = REPO_ROOT / "tests" / "parity" / "run_chebyshev_r.R"
assert _CHEBYSHEV_R_HELPER.exists(), (
    f"_CHEBYSHEV_R_HELPER not found at {_CHEBYSHEV_R_HELPER} — REPO_ROOT arithmetic is wrong"
)

# Convergence params — identical in R helper and Python call below.
_CHEBYSHEV_CONV = {"metric": "max_err", "rule": "improvement", "tol": 1e-4}

_CHEBYSHEV_SUBPROCESS_ENV = {**os.environ,
                              "OMP_NUM_THREADS": "1",
                              "OPENBLAS_NUM_THREADS": "1",
                              "MKL_NUM_THREADS": "1"}


@pytest.mark.skipif(not RSCRIPT_AVAILABLE, reason="Rscript not found")
def test_chebyshev_default_tol_parity(tmp_path):
    """R and Python must agree on weights, iter count, and status for chebyshev at tol=1e-4.

    Fixture: K=9 margins with 50 unique interaction cells from n=5000 rows
    (reuses _make_correlated_fixture from yh0l T5).  The K=9 warm-start oris
    uses inner_max_iter = max(5, min(100, max_iterations/10)).

    WHY max_iterations=40:
      inner_max_iter / 10 = 40/10 = 4.
      min(100, 4) = 4.
      max(5,  4) = 5   (post-T4, floor=5 — binding constraint).
      max(50, 4) = 50  (pre-T4,  floor=50 — binding constraint).
      At max_iterations=3000: both collapse to min(100, 300)=100 → floor irrelevant.
      At max_iterations=40:   floor literal IS the binding constraint in both branches.
      The test asserts R=Py mirror invariant at HEAD (both paths use floor=5 / floor=5).

    NOTE: chebyshev IPM is robust to warm-start quality — different floor values (5 vs 50
    warm-start iters) produce numerically identical final weights because the IPM quickly
    forgets its initial point.  The red-green differentiation for the T4 typo fix is
    therefore via the mirror invariant contract (R floor == Py floor), not via weight
    divergence.  Weight-level regression is still caught: if either bridge changes its
    floor independently, R ≠ Py on this fixture.

    Red-green gate: weight vector parity max|w_R - w_Py| < 1e-6.
    """
    df, tgt, data_csv, targets_json = _make_correlated_fixture(tmp_path)
    out_csv = tmp_path / "chebyshev_r_out.csv"

    result = subprocess.run(
        ["Rscript", str(_CHEBYSHEV_R_HELPER),
         str(data_csv), str(targets_json), str(out_csv)],
        capture_output=True, text=True, timeout=180,
        env=_CHEBYSHEV_SUBPROCESS_ENV,
    )
    if result.returncode == 2:
        pytest.skip(f"R package not available: {result.stderr.strip()}")
    if result.returncode != 0:
        raise RuntimeError(f"Rscript failed:\n{result.stderr}")

    r_out    = pd.read_csv(out_csv)
    w_r      = r_out["weight"].to_numpy(dtype=np.float64)
    iters_r  = int(r_out["iterations"].iloc[0])
    status_r = int(r_out["status"].iloc[0])

    # Python in-proc call — identical params as R helper.
    py_res   = harvest(
        df, tgt,
        method         = "chebyshev",
        min_weight     = 0,
        max_weight     = 5,
        max_iterations = 40,
        convergence    = _CHEBYSHEV_CONV,
        verbose        = 0,
        attach_weights = False,
    )
    w_py      = np.array(py_res["weights"], dtype=np.float64)
    iters_py  = int(py_res["result"]["iterations"])
    status_py = int(py_res["result"]["status"])

    print(
        f"\n  chebyshev: iters_R={iters_r} iters_Py={iters_py} "
        f"status_R={status_r} status_Py={status_py}"
    )

    # --- Same status ---
    assert status_r == status_py, (
        f"chebyshev: status mismatch — R={status_r}, Python={status_py}"
    )

    # --- Identical iter count ---
    assert iters_r == iters_py, (
        f"chebyshev: iter count mismatch — R={iters_r}, Python={iters_py} "
        "(warm-start floor mismatch between r_bridge.cpp and c_api.cpp?)"
    )

    # --- Weight vector parity ---
    assert len(w_r) == len(w_py), (
        f"chebyshev: weight vector length mismatch — R={len(w_r)}, Python={len(w_py)}"
    )
    max_abs_diff = np.max(np.abs(w_r - w_py))
    print(f"  chebyshev: max|w_R - w_Py| = {max_abs_diff:.2e}")
    assert max_abs_diff < 1e-6, (
        f"chebyshev parity FAIL: max|w_R - w_Py| = {max_abs_diff:.2e} >= 1e-6 "
        f"(iters_R={iters_r}, iters_Py={iters_py}; "
        "check warm-start floor in c_api.cpp:359 — must match r_bridge.cpp:648 floor=5; "
        "see T4 fix b60a3fd)"
    )
