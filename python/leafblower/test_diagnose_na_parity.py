"""R-parity test for diagnose_weights NA-bin handling (NABIN.3 / leafblower-4ihf.3).

Verifies the injected add_na_proportion "NA" bin is counted identically in R
(R/diagnose_weights.R) and Python (_harvest.py diagnose_weights), with an
all-observations denominator so shares sum to 1.
"""

import json
import os
import subprocess

import numpy as np
import pandas as pd

from leafblower import harvest, diagnose_weights

_REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), os.pardir, os.pardir)
)


def _make_fixture(seed: int = 1, n: int = 500, na_frac: float = 0.20):
    rng = np.random.default_rng(seed)
    levels = ["a", "b", "c"]
    g_raw = rng.choice(levels, size=n).tolist()
    na_idx = rng.choice(n, size=round(n * na_frac), replace=False).tolist()
    for i in na_idx:
        g_raw[i] = None
    return g_raw


def test_diagnose_weights_na_bin_r_python_parity():
    """R and Python diagnose_weights agree on the injected 'NA' bin (rtol=1e-6)."""
    na_frac = 0.20
    g_raw = _make_fixture(seed=1, n=500, na_frac=na_frac)
    target = {"g": {"a": 0.4, "b": 0.35, "c": 0.25}}

    # --- Python side ---
    df_py = pd.DataFrame({"g": pd.Categorical(g_raw, categories=["a", "b", "c"])})
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

    # diagnose_weights uses the RENORMALISED target that harvest built (existing
    # levels scaled by (1 - na_frac), plus an explicit "NA" bin). Rebuild it
    # here identically (harvest does not expose it on the result).
    base = {"a": 0.4, "b": 0.35, "c": 0.25}
    tgt_used = {k: v * (1.0 - na_frac) for k, v in base.items()}
    tgt_used["NA"] = na_frac
    diag_py = diagnose_weights(df_py, {"g": tgt_used}, w_py)
    # Index by level for deterministic comparison.
    diag_py = diag_py.set_index("level").sort_index()

    # --- R side via subprocess (load_all picks up this worktree's source) ---
    fixture = {
        "g": [("NA" if v is None else v) for v in g_raw],
        "weights": w_py.tolist(),
        "target_used": {k: float(v) for k, v in tgt_used.items()},
    }
    fx_path = "/tmp/diagnose_na_fixture.json"
    with open(fx_path, "w") as fh:
        json.dump(fixture, fh)

    r_script = (
        f'devtools::load_all("{_REPO_ROOT}", quiet = TRUE); '
        f'f <- jsonlite::fromJSON("{fx_path}"); '
        'g_raw <- f$g; g_raw[g_raw == "NA"] <- NA; '
        'g <- factor(g_raw, levels = c("a", "b", "c")); '
        'df <- data.frame(g = g, stringsAsFactors = FALSE); '
        # Use the exact renormalised target Python used (named numeric vector).
        'tu <- unlist(f$target_used); '
        'target <- list(g = tu); '
        'w <- as.numeric(f$weights); '
        'd <- diagnose_weights(df, target, w); '
        'out <- list('
        '  level = as.character(d$level), '
        '  prop_original = as.numeric(d$prop_original), '
        '  prop_weighted = as.numeric(d$prop_weighted)); '
        'cat(jsonlite::toJSON(out, auto_unbox = FALSE, digits = 15))'
    )
    proc = subprocess.run(
        ["Rscript", "-e", r_script],
        capture_output=True, text=True, timeout=120,
    )
    assert proc.returncode == 0, f"Rscript failed:\n{proc.stderr}"
    r_out = json.loads(proc.stdout.strip())
    diag_r = pd.DataFrame(
        {
            "level": [str(x) for x in r_out["level"]],
            "prop_original": r_out["prop_original"],
            "prop_weighted": r_out["prop_weighted"],
        }
    ).set_index("level").sort_index()

    # Align levels.
    assert set(diag_r.index) == set(diag_py.index), (
        f"Level mismatch: R={sorted(diag_r.index)}, Py={sorted(diag_py.index)}"
    )

    # --- The NA-bin must be nonzero on BOTH sides (regression guard). ---
    assert diag_r.loc["NA", "prop_original"] > 0.0
    assert diag_py.loc["NA", "prop_original"] > 0.0
    # NA-bin observed share == na_frac (uniform observed count / n).
    assert abs(diag_py.loc["NA", "prop_original"] - na_frac) < 1e-9

    # --- Full R<->Python parity (rtol=1e-6, atol=0). ---
    for col in ("prop_original", "prop_weighted"):
        assert np.allclose(
            diag_r[col].to_numpy(), diag_py[col].to_numpy(),
            rtol=1e-6, atol=0.0,
        ), (
            f"R<->Python {col} mismatch:\n"
            f"R={diag_r[col].to_dict()}\nPy={diag_py[col].to_dict()}"
        )
