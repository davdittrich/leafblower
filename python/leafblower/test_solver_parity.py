"""R↔Python parity tests for greg, newton_kl, logit, chebyshev, greenkhorn.

Fixture generated from R: set.seed(42); sample(c("x","y","z"), 300, replace=TRUE)
and sample(c("p","q"), 300, replace=TRUE).  Embedded as constants so both sides
run on identical data without filesystem side-effects.

Pattern mirrors test_harvest_na_parity.py:
  1. Build Python result with explicit convergence spec.
  2. Build R result via subprocess Rscript with the SAME explicit convergence spec.
  3. Precheck both converged (max_error < CONV_TOL) before asserting np.allclose.
  4. np.allclose(w_r, w_py, rtol=1e-6, atol=0.0).

Most solver cases pass an explicit convergence spec (rule="improvement",
tol=0.001) on both sides so the stopping point is fixed by the spec, not the
per-method default.  Same spec → same stopping point → byte-identical weights
from the shared C++ core.

The logit case (test_logit_default_rule_parity) deliberately passes NO explicit
rule (convergence=list()/{}) to LOCK the per-method DEFAULT-rule resolution
across bindings.  Both R parse_convergence() and Python _parse_convergence()
resolve method="logit" + empty convergence to (metric=max_err, rule=improvement,
pct_tol=0.001) — integer codes metric=0, rule=1.  See leafblower-6uhm.2.

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
    '  rule=ri$convergence_used$rule,'
    '  weights=as.numeric(res$weights)'
    '), auto_unbox=TRUE, digits=15))'
)


def _run_r(method: str, conv_r: str = _CONV_R) -> dict:
    """Run harvest in R via subprocess; return {max_error, weights}.

    conv_r is the literal R `convergence=` expression (e.g. 'list()' to
    exercise the per-method default rule).
    """
    r_script = _R_TEMPLATE.format(
        fixture_json=_FIXTURE_JSON,
        method=method,
        conv_r=conv_r,
    )
    proc = subprocess.run(
        ["Rscript", "-e", r_script],
        capture_output=True, text=True, timeout=90,
    )
    assert proc.returncode == 0, f"Rscript failed for method={method}:\n{proc.stderr}"
    return json.loads(proc.stdout.strip())


def _run_py(method: str, conv_py: dict = _CONV_PY):
    """Run harvest in Python; return (weights_array, result_dict).

    conv_py is the Python `convergence=` dict (e.g. {} to exercise the
    per-method default rule).
    """
    df = _make_py_df()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = harvest(df, _TARGETS, method=method, max_iterations=500,
                      convergence=conv_py, attach_weights=True)
    return res["weights"].to_numpy(), res.attrs.get("result", {})


def _assert_parity(method: str, conv_py: dict = _CONV_PY, conv_r: str = _CONV_R):
    w_py, ri_py = _run_py(method, conv_py)
    assert ri_py.get("max_error", 1.0) < _CONV_TOL, (
        f"Python {method} did not converge: max_error={ri_py.get('max_error')}"
    )
    r_out = _run_r(method, conv_r)
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
    """logit: R↔Python weights match to rtol=1e-6 (explicit improvement spec)."""
    _assert_parity("logit")


def test_logit_default_rule_parity():
    """logit with NO explicit rule: R↔Python weights match to rtol=1e-6.

    Locks the per-method DEFAULT-rule resolution: harvest(..., method="logit",
    convergence=list()/{}) must resolve to the SAME (metric, rule, tol) on both
    bindings. Both resolve metric=max_err (int 0), rule=improvement (int 1),
    pct_tol=0.001 — so the shared C++ core stops at the identical iterate and
    returns byte-identical weights over the shared fixture.

    This is the regression guard for leafblower-6uhm.2: an earlier comment in
    this file claimed R logit defaulted to rule="threshold" while Python
    defaulted to rule="improvement". That mismatch never existed on a shared
    fixture (the apparent divergence was RNG in non-shared test data); this test
    pins the default-rule alignment so any future drift fails loudly.

    Weight parity alone cannot catch a rule divergence on a fixture that
    converges to machine precision (threshold and improvement stop at the same
    iterate), so we ALSO assert the resolved rule string on BOTH bindings.
    """
    # Resolved-rule lock: both bindings must default logit to "improvement".
    _, ri_py = _run_py("logit", conv_py={})
    assert ri_py.get("convergence_used", {}).get("rule") == "improvement", (
        f"Python logit default rule = "
        f"{ri_py.get('convergence_used', {}).get('rule')!r}, expected 'improvement'"
    )
    r_out = _run_r("logit", conv_r="list()")
    assert r_out.get("rule") == "improvement", (
        f"R logit default rule = {r_out.get('rule')!r}, expected 'improvement'"
    )
    # Weight parity over the shared fixture.
    _assert_parity("logit", conv_py={}, conv_r="list()")


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
