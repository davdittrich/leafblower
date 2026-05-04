"""Cross-language weight parity: R and Python must return identical weights
for all methods when both call the same C++ solver with identical input."""

import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# Python leafblower is installed or in python/ subdir
try:
    from leafblower import harvest
except ImportError:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "python"))
    from leafblower import harvest

REPO_ROOT = Path(__file__).resolve().parent.parent
R_HELPER  = REPO_ROOT / "tests" / "parity" / "run_parity_r.R"

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
    if result.returncode != 0:
        raise RuntimeError(f"Rscript failed:\n{result.stderr}")
    return pd.read_csv(out_csv)["weight"].to_numpy()


@pytest.mark.skipif(not RSCRIPT_AVAILABLE, reason="Rscript not found")
@pytest.mark.parametrize("method", [
    "greenkhorn", "logit", "raking", "ieppa", "sinkhorn", "newton_kl",
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

    diff = np.max(np.abs(w_py - w_r))
    print(f"\n  {method}: max|w_py - w_r| = {diff:.2e}")
    assert diff < 1e-10, (
        f"{method}: weight vectors differ by {diff:.2e} "
        f"(expected < 1e-10 for shared C++ solver)"
    )
