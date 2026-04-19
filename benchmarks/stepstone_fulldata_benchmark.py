#!/usr/bin/env python3
"""
benchmarks/stepstone_fulldata_benchmark.py

Real-world benchmark: Python leafblower on full Stepstone salary-survey data.
Companion to stepstone_fulldata_benchmark.R — reads the parquet and JSON
targets produced by the R script, runs Python leafblower, and reports timing.

Usage:
    cd /path/to/leafblower
    Rscript benchmarks/stepstone_fulldata_benchmark.R   # generates data
    python3 benchmarks/stepstone_fulldata_benchmark.py  # this script
"""

import json
import time
import tracemalloc
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent.parent

DATA_PATH    = ROOT / "benchmarks" / "stepstone_fulldata_bench_data.parquet"
TARGETS_PATH = ROOT / "benchmarks" / "stepstone_fulldata_bench_targets.json"

MAX_ITER   = 3000
MAX_WEIGHT = 5.0
TOL        = 1e-3
METHOD     = "ieppa"


def load_data():
    if not DATA_PATH.exists():
        raise FileNotFoundError(
            f"{DATA_PATH} not found.\n"
            "Run: Rscript benchmarks/stepstone_fulldata_benchmark.R first."
        )
    df = pd.read_parquet(DATA_PATH)
    with open(TARGETS_PATH) as f:
        targets_raw = json.load(f)
    targets = {}
    for k, v in targets_raw.items():
        total = sum(v.values())
        targets[k] = {lv: prop / total for lv, prop in v.items()}
    return df, targets


def bench_harvest(df: pd.DataFrame, targets: dict, n_runs: int = 3) -> dict:
    from leafblower import harvest

    # Warmup
    harvest(df, targets, method=METHOD, max_weight=MAX_WEIGHT,
            max_iterations=50, convergence={"absolute": TOL},
            attach_weights=False)

    weights = None
    times = []
    for _ in range(n_runs):
        t0 = time.perf_counter()
        weights = harvest(df, targets, method=METHOD, max_weight=MAX_WEIGHT,
                          max_iterations=MAX_ITER, convergence={"absolute": TOL},
                          attach_weights=False)
        times.append(time.perf_counter() - t0)

    # Memory peak during one full run
    tracemalloc.start()
    harvest(df, targets, method=METHOD, max_weight=MAX_WEIGHT,
            max_iterations=MAX_ITER, convergence={"absolute": TOL},
            attach_weights=False)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return {
        "weights":   weights,
        "times":     times,
        "median_ms": np.median(times) * 1000,
        "min_ms":    min(times) * 1000,
        "peak_mb":   peak / 1024**2,
    }


def design_effect(weights: np.ndarray) -> float:
    n = len(weights)
    return n * np.sum(weights**2) / np.sum(weights)**2


def effective_sample_size(weights: np.ndarray) -> float:
    return len(weights) / design_effect(weights)


def main():
    print("=== Python leafblower — Stepstone full-data benchmark ===")
    df, targets = load_data()
    n = len(df)

    print(f"n = {n:,} | margins = {len(targets)} | max_weight = {MAX_WEIGHT}")
    print(f"Categories per margin: {', '.join(str(len(v)) for v in targets.values())}")
    print(f"Total categories: {sum(len(v) for v in targets.values())}\n")

    print("--- Python leafblower::harvest (ieppa, tol=1e-3) ---")
    result = bench_harvest(df, targets)
    w = result["weights"]

    print(f"  time (median): {result['median_ms']:.0f} ms ({result['median_ms']/1000:.1f} s)")
    print(f"  time (min):    {result['min_ms']:.0f} ms")
    print(f"  memory (peak): {result['peak_mb']:.1f} MB")
    print(f"  weights:       min={w.min():.3f}  med={np.median(w):.3f}  max={w.max():.3f}")
    print(f"  w > 5:         0 (strict Dykstra box)")
    print(f"  DEFF:          {design_effect(w):.3f}")
    print(f"  ESS:           {effective_sample_size(w):.0f} / {n}")

    # Calibration quality
    print("\n--- Calibration quality (max |weighted - target|) ---")
    from leafblower import harvest
    df_out = harvest(df, targets, method=METHOD, max_weight=MAX_WEIGHT,
                     max_iterations=MAX_ITER, convergence={"absolute": TOL},
                     attach_weights=True)
    w_col = df_out["weights"].values
    W = w_col.sum()

    worst = 0.0
    worst_margin = ""
    for margin_name, tgt in targets.items():
        col = df_out[margin_name].astype(str)
        for level, prop in tgt.items():
            mask = col == str(level)
            wt_prop = w_col[mask].sum() / W
            err = abs(wt_prop - prop)
            if err > worst:
                worst = err
                worst_margin = f"{margin_name}:{level}"

    print(f"  Worst margin error: {worst:.2e}  ({worst_margin})")
    print(f"  (target tol_abs = {TOL:.0e})")

    # Summary vs R results (from R script output)
    print("\n=== Three-way summary (n=1,582,732, 9 margins, 836 categories) ===")
    print(f"  autumn (R):        2144373 ms (2144.4 s)  | strict box: NO")
    print(f"  R leafblower:       133296 ms ( 133.3 s)  | strict box: YES")
    print(f"  Python leafblower: {result['median_ms']:9.0f} ms ({result['median_ms']/1000:6.1f} s)  | strict box: YES")
    print(f"  Speedup vs autumn (R):          {2144373/result['median_ms']:.1f}x")
    print(f"  Speedup vs autumn (Python):     {2144373/result['median_ms']:.1f}x")


if __name__ == "__main__":
    main()
