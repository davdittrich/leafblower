"""R↔Python parity tests for greg, newton_kl, logit, chebyshev, greenkhorn.

Fixture generated from R: set.seed(42); sample(c("x","y","z"), 300, replace=TRUE)
and sample(c("p","q"), 300, replace=TRUE).  Embedded as constants so both sides
run on identical data without filesystem side-effects.

Pattern mirrors test_harvest_na_parity.py:
  1. Build Python result with explicit convergence spec.
  2. Build R result via subprocess Rscript with the SAME explicit convergence spec.
  3. Precheck both converged (max_error < CONV_TOL) before asserting np.allclose.
  4. np.allclose(w_r, w_py, rtol=1e-6, atol=0.0).

Explicit convergence spec (rule="improvement", tol=0.001) is passed on both
sides to avoid the R/Python default-spec mismatch: R logit defaults to
rule="threshold" (12 iters) while Python defaults to rule="improvement" (22
iters).  Same spec → same stopping point → byte-identical weights from the
shared C++ core.

greenkhorn/sinkhorn are entropic solvers — sum(w) != n by construction.
Parity assertion is still valid: both sides call the same C++ core.

greg is a one-shot linear solver (1 iter regardless); max_error ~0.005 is the
chi2-scaled residual, not a non-convergence; _CONV_TOL=0.01 covers it.
"""

import json
import subprocess
import warnings

import numpy as np
import pandas as pd
import pytest

from leafblower import harvest

# ---------------------------------------------------------------------------
# Shared fixture — generated from R set.seed(42) to guarantee identical data
# ---------------------------------------------------------------------------
_A = [
    "x","x","x","x","y","y","y","x","z","z","x","x","y","y","y","z","z","x","x","z",
    "x","z","x","x","y","z","y","x","y","y","z","z","y","y","y","y","x","x","y","y",
    "z","z","x","x","y","y","y","y","y","z","y","x","y","z","y","y","y","x","x","x",
    "z","x","y","x","x","x","z","z","x","x","y","x","y","y","y","x","y","x","y","z",
    "x","x","z","y","y","z","x","x","y","z","y","x","x","x","x","x","y","z","x","y",
    "x","z","y","x","x","y","x","x","z","y","z","y","y","x","z","x","x","x","z","z",
    "x","y","y","z","y","z","z","y","x","x","y","x","y","z","y","y","y","y","y","z",
    "y","z","y","y","z","y","y","x","x","y","y","x","z","x","y","y","y","y","z","z",
    "y","x","y","x","y","x","x","x","z","y","y","y","x","x","z","y","y","z","z","x",
    "z","y","x","z","z","y","z","x","x","x","y","x","z","y","z","x","x","x","x","x",
    "y","z","z","z","x","x","x","y","z","z","x","z","y","z","y","z","z","x","x","x",
    "z","y","x","z","z","x","x","x","y","x","x","y","z","x","y","y","y","y","z","y",
    "z","y","x","x","z","z","x","y","x","y","x","y","x","z","y","z","y","y","y","z",
    "x","z","x","z","y","x","y","y","y","x","z","y","y","x","z","z","y","z","x","z",
    "z","x","y","y","x","y","z","x","y","z","x","x","y","z","x","y","y","z","z","y",
]
_B = [
    "q","p","p","q","p","p","p","q","p","q","p","q","p","p","q","q","q","p","q","q",
    "q","p","q","p","q","p","p","q","q","p","p","q","p","q","q","q","q","q","p","p",
    "q","p","p","q","q","p","p","q","q","q","q","q","p","p","p","p","q","p","p","q",
    "q","p","p","q","q","p","q","q","p","p","p","q","p","p","q","p","q","q","q","p",
    "p","q","p","q","p","q","q","p","p","q","p","p","q","q","q","q","q","p","q","p",
    "p","p","q","q","q","p","q","p","q","q","p","q","q","p","p","p","q","q","p","q",
    "p","q","p","p","q","q","q","p","p","p","q","p","q","q","p","p","p","q","p","p",
    "p","q","p","p","p","q","q","q","q","p","q","p","p","q","p","p","q","q","p","p",
    "p","q","q","p","q","p","q","p","p","q","p","q","q","p","p","q","p","q","q","q",
    "q","q","p","p","p","q","p","p","q","q","q","p","q","p","p","q","p","p","p","q",
    "p","q","q","p","p","p","q","q","q","p","p","p","p","p","q","q","q","q","q","p",
    "q","q","p","p","q","q","q","q","p","p","q","q","p","p","p","q","p","q","p","p",
    "p","p","q","p","q","q","p","q","q","p","p","p","q","p","q","p","q","p","p","p",
    "q","p","q","q","p","p","q","q","q","p","q","p","q","q","q","p","q","p","p","p",
    "p","q","p","q","q","q","p","q","q","p","q","p","q","p","q","p","p","q","p","p",
]

_TARGETS = {
    "a": {"x": 1/3, "y": 1/3, "z": 1/3},
    "b": {"p": 0.5, "q": 0.5},
}

# Explicit convergence spec applied identically on both sides to ensure
# the same stopping point from the shared C++ core.
_CONV_PY = {"rule": "improvement", "tol": 0.001}
_CONV_R  = 'list(rule="improvement", tol=0.001)'

# Precheck tolerance: max_error must be below this before asserting parity.
# greg is a 1-iter linear solver with max_error ~0.005 (chi2 residual, not
# a convergence failure) — 0.01 is the appropriate upper bound.
_CONV_TOL = 0.01


def _make_py_df():
    return pd.DataFrame({
        "a": pd.Categorical(_A, categories=["x", "y", "z"]),
        "b": pd.Categorical(_B, categories=["p", "q"]),
    })


_FIXTURE_JSON = json.dumps({"a": _A, "b": _B})

_R_TEMPLATE = (
    'fixture <- jsonlite::fromJSON(\'{fixture_json}\'); '
    'df <- data.frame('
    '  a = factor(fixture$a, levels=c("x","y","z")),'
    '  b = factor(fixture$b, levels=c("p","q")),'
    '  stringsAsFactors=FALSE'
    '); '
    'tgt <- list(a=c(x=1/3,y=1/3,z=1/3), b=c(p=0.5,q=0.5)); '
    'res <- suppressWarnings(leafblower::harvest(df, tgt, method="{method}",'
    '  max_iterations=500L, convergence={conv_r})); '
    'ri <- attr(res,"result"); '
    'cat(jsonlite::toJSON(list('
    '  max_error=ri$max_error,'
    '  weights=as.numeric(res$weights)'
    '), auto_unbox=TRUE, digits=15))'
)


def _run_r(method: str) -> dict:
    """Run harvest in R via subprocess; return {max_error, weights}."""
    r_script = _R_TEMPLATE.format(
        fixture_json=_FIXTURE_JSON,
        method=method,
        conv_r=_CONV_R,
    )
    proc = subprocess.run(
        ["Rscript", "-e", r_script],
        capture_output=True, text=True, timeout=90,
    )
    assert proc.returncode == 0, f"Rscript failed for method={method}:\n{proc.stderr}"
    return json.loads(proc.stdout.strip())


def _run_py(method: str):
    """Run harvest in Python; return (weights_array, result_dict)."""
    df = _make_py_df()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = harvest(df, _TARGETS, method=method, max_iterations=500,
                      convergence=_CONV_PY, attach_weights=True)
    return res["weights"].to_numpy(), res.attrs.get("result", {})


def _assert_parity(method: str):
    w_py, ri_py = _run_py(method)
    assert ri_py.get("max_error", 1.0) < _CONV_TOL, (
        f"Python {method} did not converge: max_error={ri_py.get('max_error')}"
    )
    r_out = _run_r(method)
    assert r_out["max_error"] < _CONV_TOL, (
        f"R {method} did not converge: max_error={r_out['max_error']}"
    )
    w_r = np.array(r_out["weights"])
    assert len(w_py) == len(w_r), (
        f"{method}: length mismatch Python={len(w_py)} R={len(w_r)}"
    )
    assert np.allclose(w_r, w_py, rtol=1e-6, atol=0.0), (
        f"{method} R↔Python mismatch: max|Δw|={np.max(np.abs(w_r - w_py)):.3e}"
    )


# ---------------------------------------------------------------------------
# Per-solver tests
# ---------------------------------------------------------------------------

def test_greg_parity():
    """greg: R↔Python weights match to rtol=1e-6.

    greg is a one-shot linear solver (1 iter); max_error ~0.005 is the
    chi2-scaled residual. _CONV_TOL=0.01 covers it.
    """
    _assert_parity("greg")


def test_newton_kl_parity():
    """newton_kl: R↔Python weights match to rtol=1e-6."""
    _assert_parity("newton_kl")


def test_logit_parity():
    """logit: R↔Python weights match to rtol=1e-6.

    R logit defaults to rule="threshold"; Python defaults to rule="improvement".
    Explicit rule="improvement" on both sides forces the same stopping point.
    """
    _assert_parity("logit")


def test_chebyshev_parity():
    """chebyshev: R↔Python weights match to rtol=1e-6."""
    _assert_parity("chebyshev")


def test_greenkhorn_parity():
    """greenkhorn: R↔Python weights match to rtol=1e-6.

    greenkhorn is an entropic solver — sum(w) != n by construction.
    Parity assertion is still valid: both sides call the same C++ core.
    max_error precheck uses the marginal-error metric, not sum-of-weights.
    """
    _assert_parity("greenkhorn")
