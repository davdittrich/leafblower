"""CR-A5 (mxcl.5): cross-solver returned-weights ≡ reported-metrics gate (Python mirror).

Mirror of tests/testthat/test-returned-weights-invariant.R. For every calibration
method, independently recompute max_error / kl / chi2 from the RETURNED weights
using a fresh NumPy implementation (never the package's own metric code) and assert
they match the reported result-dict fields.

The three Blocking defects in epic CR-A (Cholesky triangle mismatch, chebyshev
warm-start pairing, ORIS SRAA stale-iterate) all returned weights that did not
realize the reported metrics. This gate closes that gap on the Python binding too.

This is a SELF-CONSISTENCY invariant (returned weights vs the solver's own report),
not an R↔Python parity test, so the fixture is generated independently here.

Formulas mirror compute_cell_metrics (src/calib_dispatch.hpp:239), kMetricEps=1e-10:
  S_p    = Σ_{i∈cat(k,j)} w_i / W        (achieved proportion, W = Σw)
  maxerr = max_{k,j} |S_p − T_kj|                              (errRp, MAX)
  kl     = max_k Σ_j [T>0] T·log((T+ε)/(S_p+ε))               (reverse-KL, MAX)
  chi2   = Σ_{k,j: pop>ε} (obs − pop)² / pop,  obs=Σw, pop=T·W (Pearson)

Per-method tolerances match the R gate (see its header for the empirical sweep).
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import warnings

import numpy as np
import pandas as pd
import pytest

from leafblower import harvest

_EPS = 1e-10

# atol_maxerr, atol_kl, atol_chi2 per method (empirical sweep × safety headroom).
_TOL = {
    "raking":     (1e-8, 1e-8, 1e-5),
    "sinkhorn":   (1e-8, 1e-8, 1e-5),
    "greenkhorn": (1e-5, 1e-8, 1e-3),   # entropic: report-vs-final drift, Σw≠n
    "logit":      (1e-5, 1e-8, 1e-5),   # best-iterate/final-Newton drift
    "greg":       (1e-8, 1e-8, 1e-5),
    "chebyshev":  (1e-8, 1e-8, 1e-5),   # CR-A2 landed → GREEN
    "newton_kl":  (1e-8, 1e-8, 1e-5),   # BLAS reassociation
    "oris":       (1e-8, 1e-8, 1e-5),
    "oris_soft":  (1e-8, 1e-8, 1e-5),
    # "auto" is R-only (R harvest dispatches; Python API rejects it) — omitted.
}


def _fixture():
    rng = np.random.default_rng(99)
    n = 600
    df = pd.DataFrame({
        "x": pd.Categorical(rng.choice(["a", "b", "c"], size=n)),
        "y": pd.Categorical(rng.choice(["p", "q"],      size=n)),
    })
    target = {"x": {"a": 1/3, "b": 1/3, "c": 1/3}, "y": {"p": 0.5, "q": 0.5}}
    return df, target


def _rederive(w, df, target):
    """Fresh, package-independent re-derivation of the three cell metrics."""
    W = float(np.sum(w))
    maxerr = 0.0
    kl = 0.0
    chi2 = 0.0
    for v, tv in target.items():
        cats = list(tv.keys())
        col = df[v].astype(str).to_numpy()
        obs = np.array([w[col == c].sum() for c in cats], dtype=np.float64)
        T = np.array([tv[c] for c in cats], dtype=np.float64)
        Sp = obs / W
        maxerr = max(maxerr, float(np.max(np.abs(Sp - T))))
        kl_k = float(np.sum(np.where(T > 0, T * np.log((T + _EPS) / (Sp + _EPS)), 0.0)))
        kl = max(kl, kl_k)
        pop = T * W
        chi2 += float(np.sum(np.where(pop > _EPS, (obs - pop) ** 2 / pop, 0.0)))
    return maxerr, kl, chi2


@pytest.mark.parametrize("method", list(_TOL.keys()))
def test_returned_weights_realize_reported_metrics(method):
    df, target = _fixture()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = harvest(df, target, method=method, max_iterations=1000,
                      attach_weights=False)
    w = np.asarray(out["weights"], dtype=np.float64)
    r = out["result"]
    assert r["status"] == 0, f"status for {method}"

    maxerr, kl, chi2 = _rederive(w, df, target)
    t = _TOL[method]
    assert abs(maxerr - r["max_error"]) < t[0], f"max_error mismatch for {method}"
    assert abs(kl     - r["kl"])        < t[1], f"kl mismatch for {method}"
    assert abs(chi2   - r["chi2"])      < t[2], f"chi2 mismatch for {method}"
