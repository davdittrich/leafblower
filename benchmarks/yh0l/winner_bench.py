#!/usr/bin/env python3
"""
benchmarks/yh0l/winner_bench.py
T3: multi-DGP winner determination — Python side
Run BEFORE winner_bench.R (memory: gotcha-openblas-thread-oversubscription).
Outputs: benchmarks/yh0l/results_Py_<dgp>.csv + weights feather per DGP/rep.

BLAS threads must be set BEFORE numpy import — set here unconditionally.
"""

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import json, math, sys, time, csv, statistics
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.feather as feather

ROOT    = Path(__file__).resolve().parent.parent.parent
YDIR    = ROOT / "benchmarks" / "yh0l"
WDIR    = YDIR / "weights_Py"
WDIR.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(ROOT / "python"))
from leafblower import harvest

NREPS  = 3
MAX_IT = 3000
MAX_WT = 5.0
TOL    = 1e-4

# ── DGP definitions ───────────────────────────────────────────────────────────
DGPS = {
    "fulldata": {
        "data":    ROOT / "benchmarks" / "stepstone_fulldata_bench_data.parquet",
        "targets": ROOT / "benchmarks" / "stepstone_fulldata_bench_targets.json",
    },
    "medium": {
        "data":    ROOT / "benchmarks" / "stepstone_bench_data.parquet",
        "targets": ROOT / "benchmarks" / "stepstone_bench_targets.json",
    },
    "small_hicard": {
        "data":    YDIR / "small_hicard_data.parquet",
        "targets": YDIR / "small_hicard_targets.json",
    },
}

_STATUS = {0: "converged", 1: "no_conv", 2: "infeasible",
           3: "bad_arg", 4: "budget", 5: "stall"}


def compute_metrics(w_arr, df, tgt):
    """Tier-1 metrics: margin_kl, weight_kl, max_err, chi2."""
    w = pd.Series(w_arr, index=df.index)
    W = float(w.sum())
    max_err, mkl, chi2 = 0.0, 0.0, 0.0
    for v, targets in tgt.items():
        if v not in df.columns:
            continue
        s_map = w.groupby(df[v]).sum() / W
        for lv, T in targets.items():
            S = float(s_map.get(str(lv), 0.0))
            max_err = max(max_err, abs(S - T))
            if T > 0 and S > 0:
                mkl += T * math.log(T / S)
            chi2 += (S - T) ** 2 / max(T, 1e-12)
    wm  = W / len(w_arr)
    pos = w_arr > 0
    wkl = float(np.sum(w_arr[pos] * np.log(w_arr[pos] / wm)) / W)
    return {"margin_kl": mkl, "max_err": max_err, "chi2": chi2, "weight_kl": wkl}


def load_dgp(dgp_cfg):
    df = pd.read_parquet(dgp_cfg["data"])
    with open(dgp_cfg["targets"]) as f:
        tgt_raw = json.load(f)
    tgt = {k: {lv: v / sum(d.values()) for lv, v in d.items()}
           for k, d in tgt_raw.items()}
    # Keep only target columns (same filter as R)
    df = df[[c for c in tgt if c in df.columns]]
    return df, tgt


def run_dgp(dgp_name, dgp_cfg):
    print(f"\n=== Py DGP: {dgp_name} ===")
    df, tgt = load_dgp(dgp_cfg)
    n = len(df)
    print(f"  n={n:,} margins={len(tgt)}")

    rows = []
    for rep in range(1, NREPS + 1):
        print(f"  rep {rep}/{NREPS} ... ", end="", flush=True)
        t0 = time.perf_counter()
        res = harvest(
            df,
            tgt,
            method         = "ieppa_soft",
            convergence    = {"tol": TOL},
            max_weight     = MAX_WT,
            max_iterations = MAX_IT,
            attach_weights = False,
            verbose        = 0,
        )
        t1      = time.perf_counter()
        wall_ms = (t1 - t0) * 1000.0

        rd      = res.get("result", res)
        st      = rd.get("status", -1)
        iters   = rd.get("iterations", -1)
        me      = rd.get("max_error", float("nan"))
        print(f"status={st}({_STATUS.get(st,'?')}) iter={iters} max_err={me:.3e} wall={wall_ms:.0f}ms")

        w_arr = np.asarray(res["weights"], dtype=np.float64)

        # Save weights
        w_path = WDIR / f"{dgp_name}_rep{rep}.feather"
        feather.write_feather(pd.DataFrame({"w": w_arr}), str(w_path))

        m = compute_metrics(w_arr, df, tgt)
        rows.append({
            "dgp":        dgp_name,
            "rep":        rep,
            "n":          n,
            "margin_kl":  m["margin_kl"],
            "weight_kl":  m["weight_kl"],
            "max_err":    m["max_err"],
            "chi2":       m["chi2"],
            "wall_ms":    wall_ms,
            "status":     st,
            "iterations": iters,
            "lang":       "Py",
        })
    return rows


def main():
    all_rows = []
    for dgp_name, dgp_cfg in DGPS.items():
        if not dgp_cfg["data"].exists():
            print(f"[SKIP] {dgp_name}: data not found at {dgp_cfg['data']}")
            continue
        try:
            rows = run_dgp(dgp_name, dgp_cfg)
            all_rows.extend(rows)
            out_path = YDIR / f"results_Py_{dgp_name}.csv"
            with open(out_path, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=rows[0].keys())
                w.writeheader()
                w.writerows(rows)
            print(f"  -> saved {out_path}")
        except Exception as e:
            print(f"[ERROR] DGP {dgp_name}: {e}")

    if all_rows:
        out_path = YDIR / "results_Py_all.csv"
        with open(out_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=all_rows[0].keys())
            w.writeheader()
            w.writerows(all_rows)
        print(f"\nAll Py results saved to {out_path}")
        print("\nSummary:")
        for r in all_rows:
            print(f"  {r['dgp']} rep={r['rep']} margin_kl={r['margin_kl']:.6f} "
                  f"wall={r['wall_ms']:.0f}ms iter={r['iterations']}")


if __name__ == "__main__":
    main()
