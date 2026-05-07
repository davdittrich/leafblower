#!/usr/bin/env python3
"""
benchmarks/yh0l/run_traj_py.py
T2: run ieppa_soft with verbose=2 on stepstone fulldata; capture log to traj_Py.log
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 must be set by caller (env vars).
"""

import json
import os
import sys
from pathlib import Path

# Force single-thread BLAS for deterministic divergence trace
# Must be set BEFORE numpy/leafblower import — BLAS thread count locks at import time.
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

ROOT = Path(__file__).resolve().parent.parent.parent
DATA_PATH    = ROOT / "benchmarks" / "stepstone_fulldata_bench_data.parquet"
TARGETS_PATH = ROOT / "benchmarks" / "stepstone_fulldata_bench_targets.json"
LOG_PATH     = ROOT / "benchmarks" / "yh0l" / "traj_Py.log"

sys.path.insert(0, str(ROOT / "python"))

import pandas as pd
from leafblower import harvest


def main():
    print(f"[run_traj_py] ROOT={ROOT}")
    print(f"[run_traj_py] DATA={DATA_PATH}")

    # ── Load data ────────────────────────────────────────────────────────────
    df = pd.read_parquet(DATA_PATH)
    with open(TARGETS_PATH) as f:
        targets_raw = json.load(f)
    targets = {}
    for k, v in targets_raw.items():
        total = sum(v.values())
        targets[k] = {lv: prop / total for lv, prop in v.items()}
    print(f"[run_traj_py] n={len(df):,} margins={len(targets)}")

    # ── Capture verbose=2 stdout to log ─────────────────────────────────────
    print(f"[run_traj_py] writing verbose=2 log to {LOG_PATH}")
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)

    import datetime
    with open(LOG_PATH, "w") as log:
        log.write(f"## run_traj_py.py start: {datetime.datetime.now()}\n")
        log.write(f"## OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS', '(unset)')} "
                  f"OPENBLAS_NUM_THREADS={os.environ.get('OPENBLAS_NUM_THREADS', '(unset)')}\n")

        # Redirect stdout to log while harvest runs
        old_stdout = sys.stdout
        sys.stdout = log
        try:
            res = harvest(
                df,
                targets,
                method         = "ieppa_soft",
                convergence    = {"tol": 1e-4},
                max_weight     = 5,
                max_iterations = 3000,
                attach_weights = False,
                verbose        = 2,
            )
            sys.stdout = old_stdout
        except Exception as exc:
            sys.stdout = old_stdout
            log.write(f"## ERROR: {exc}\n")
            log.write(f"## TERMINATE status=ERROR iterations=NA max_error=NA\n")
            log.write(f"## run_traj_py.py end: {datetime.datetime.now()}\n")
            print(f"[run_traj_py] ERROR: {exc}")
            return

        status_map = {0: "converged", 1: "no_conv", 2: "infeasible",
                      3: "bad_arg",   4: "budget",   5: "stall"}
        # attach_weights=False → {"weights": arr, "result": result_dict}
        rd          = res.get("result", res) if isinstance(res, dict) else {}
        st          = rd.get("status", -1)
        iterations  = rd.get("iterations", -1)
        max_error   = rd.get("max_error", float("nan"))
        log.write(f"## TERMINATE status={st} ({status_map.get(st, 'unknown')}) "
                  f"iterations={iterations} max_error={max_error:.6e}\n")
        log.write(f"## run_traj_py.py end: {datetime.datetime.now()}\n")

    print(f"[run_traj_py] done. Log: {LOG_PATH}")


if __name__ == "__main__":
    main()
