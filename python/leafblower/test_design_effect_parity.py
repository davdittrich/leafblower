"""R-parity test for design_effect Python wrapper."""

import json
import subprocess

import numpy as np
import pandas as pd
import pytest

from leafblower import design_effect


def test_r_python_parity_kish_1arg():
    """1-arg Kish: Python and R formula produce identical float."""
    w = np.array([1.0, 2.0, 3.0, 4.0])
    d_py = design_effect(w)
    expected = len(w) * np.sum(w**2) / np.sum(w)**2
    assert d_py == pytest.approx(expected, rel=1e-12)


def test_r_python_parity_4arg_same_c_entry():
    """4-arg: Python and R both call rk_design_effect → byte-identical deff_H."""
    np.random.seed(2024)
    n = 50
    regions = np.random.choice(["N", "S", "E", "W"], n)
    data = pd.DataFrame({"region": regions})
    target = {"region": {"N": 0.25, "S": 0.25, "E": 0.25, "W": 0.25}}
    w = np.ones(n)
    y = 10.0 + 3.0 * (regions == "N") - 2.0 * (regions == "S") + np.random.randn(n)

    d_py = design_effect(w, outcome=y, data=data, target=target)

    # Call same fixture from R
    fixture = {"weights": w.tolist(), "outcome": y.tolist(), "region": regions.tolist()}
    with open("/tmp/de_fixture.json", "w") as f:
        json.dump(fixture, f)

    r_script = (
        'f <- jsonlite::fromJSON("/tmp/de_fixture.json"); '
        'data <- data.frame(region = f$region, stringsAsFactors = FALSE); '
        'target <- list(region = c(N = 0.25, S = 0.25, E = 0.25, W = 0.25)); '
        'd <- leafblower::design_effect(f$weights, outcome = f$outcome, data = data, target = target); '
        'cat(sprintf("%.20f", d))'
    )
    result = subprocess.run(["Rscript", "-e", r_script], capture_output=True, text=True, timeout=30)
    d_r = float(result.stdout.strip())
    assert d_py == pytest.approx(d_r, rel=1e-10), (
        f"Python deff_H={d_py:.15f} R deff_H={d_r:.15f} diff={abs(d_py-d_r):.2e}"
    )
