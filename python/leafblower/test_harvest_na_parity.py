"""R-parity test for harvest(add_na_proportion=True) — tu15."""

import json
import subprocess

import numpy as np
import pandas as pd

from leafblower import harvest


def _make_fixture(seed: int = 1, n: int = 500, na_frac: float = 0.20):
    rng = np.random.default_rng(seed)
    levels = ["a", "b", "c"]
    g_raw = rng.choice(levels, size=n).tolist()
    na_idx = rng.choice(n, size=round(n * na_frac), replace=False).tolist()
    for i in na_idx:
        g_raw[i] = None  # NA
    return g_raw


def test_harvest_na_parity_r_python():
    """Python and R harvest(add_na_proportion=True) produce rtol=1e-6 matching weights."""
    g_raw = _make_fixture(seed=1, n=500, na_frac=0.20)
    target = {"g": {"a": 0.4, "b": 0.35, "c": 0.25}}

    # --- Python side ---
    g_str = [None if v is None else str(v) for v in g_raw]
    df_py = pd.DataFrame({"g": pd.Categorical(
        g_str, categories=["a", "b", "c"]
    )})
    res_py = harvest(
        df_py,
        targets=target,
        method="raking",
        add_na_proportion=True,
        max_iterations=500,
        attach_weights=True,
        weight_column="weights",
    )
    r_result = res_py.attrs.get("result", {})
    assert r_result.get("status", -1) == 0, (
        f"Python harvest did not converge: result={r_result}"
    )
    w_py = res_py["weights"].to_numpy()

    # --- R side via subprocess ---
    fixture = {"g": [("NA" if v is None else v) for v in g_raw]}
    with open("/tmp/harvest_na_fixture.json", "w") as fh:
        json.dump(fixture, fh)

    r_script = (
        'f <- jsonlite::fromJSON("/tmp/harvest_na_fixture.json"); '
        'g_raw <- f$g; '
        'g_raw[g_raw == "NA"] <- NA; '
        'g <- factor(g_raw, levels = c("a", "b", "c")); '
        'df <- data.frame(g = g, stringsAsFactors = FALSE); '
        'target <- list(g = c(a = 0.4, b = 0.35, c = 0.25)); '
        'res <- leafblower::harvest(df, target, method = "raking", '
        '  add_na_proportion = TRUE, max_iterations = 500L); '
        'conv <- (attr(res, "result")$status == 0L); '
        'w <- res$weights; '
        'out <- list(converged = conv, weights = as.numeric(w)); '
        'cat(jsonlite::toJSON(out, auto_unbox = TRUE, digits = 15))'
    )
    proc = subprocess.run(
        ["Rscript", "-e", r_script],
        capture_output=True, text=True, timeout=60
    )
    assert proc.returncode == 0, f"Rscript failed:\n{proc.stderr}"
    r_out = json.loads(proc.stdout.strip())

    assert r_out["converged"], (
        f"R harvest did not converge. stderr:\n{proc.stderr}"
    )
    w_r = np.array(r_out["weights"])

    assert len(w_py) == len(w_r), (
        f"Weight vector length mismatch: Python={len(w_py)}, R={len(w_r)}"
    )
    assert np.allclose(w_r, w_py, rtol=1e-6), (
        f"R↔Python weight mismatch: max|Δw|={np.max(np.abs(w_r - w_py)):.3e}, "
        f"rtol=1e-6"
    )
