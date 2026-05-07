#!/usr/bin/env python3
# benchmarks/parity_bench.py
# Python side of R-Python parity benchmark on stepstone fulldata.
# All methods, defaults (per-method natural metric, mirrored from R harvest),
# tol=1e-4. Saves per-method weights to feather for cross-language comparison.

import sys, time, json, math, csv, statistics
from pathlib import Path
import re

import pandas as pd
import numpy as np
import pyarrow.feather as feather

repo = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(repo / "python"))
from leafblower import harvest

DATA    = repo / "benchmarks/stepstone_fulldata_bench_data.parquet"
TARGETS = repo / "benchmarks/stepstone_fulldata_bench_targets.json"
WDIR    = repo / "benchmarks/results/parity/weights_Py"
OUT_CSV = repo / "benchmarks/results/parity/parity_bench_Py.csv"
MAX_IT  = 3000
MAX_WT  = 5.0
TOL     = 1e-4
NREPS   = 3
DEF     = {"tol": TOL}   # defaults: per-method natural metric (harvest.py switch)

WDIR.mkdir(parents=True, exist_ok=True)

_STATUS = {0: "converged", 1: "no_conv", 2: "infeasible", 3: "bad_arg",
           4: "budget",    5: "stall"}

print("Loading data...")
df  = pd.read_parquet(DATA)
with open(TARGETS) as f:
    tgt_raw = json.load(f)
tgt = {k: {lv: v / sum(d.values()) for lv, v in d.items()} for k, d in tgt_raw.items()}
df  = df[[c for c in tgt if c in df.columns]]
n = len(df)
print(f"n = {n:,} | margins = {len(tgt)} | max_iter = {MAX_IT} | tol = {TOL:.0e}\n")


def compute_metrics(w_arr, df, tgt):
    w = pd.Series(w_arr, index=df.index)
    W = float(w.sum())
    max_err, mkl, chi2 = 0.0, 0.0, 0.0
    for v, targets in tgt.items():
        if v not in df.columns:
            continue
        s_map = w.groupby(df[v]).sum() / W
        for lv, T in targets.items():
            S = float(s_map.get(lv, 0.0))
            max_err = max(max_err, abs(S - T))
            if T > 0 and S > 0:
                mkl += T * math.log(T / S)
            chi2 += (S - T) ** 2 / max(T, 1e-12)
    wm = W / len(w_arr)
    pos = w_arr > 0
    wkl = float(np.sum(w_arr[pos] * np.log(w_arr[pos] / wm)) / W)
    return {"max_err": max_err, "marginal_kl": mkl, "kl": wkl, "chi2": chi2}


def slug(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", s)


def run(label, conv, **kwargs):
    times = []
    res_last = None
    for _ in range(NREPS):
        t0 = time.perf_counter()
        try:
            res = harvest(df, tgt, min_weight=0.0, max_weight=MAX_WT,
                          max_iterations=MAX_IT, convergence=conv,
                          verbose=0, attach_weights=False, **kwargs)
        except Exception as e:
            print(f"  ERROR: {e}")
            return None
        times.append((time.perf_counter() - t0) * 1000)
        res_last = res
    assert res_last is not None
    ri   = res_last["result"]
    w_np = np.asarray(res_last["weights"], dtype=np.float64)
    m    = compute_metrics(w_np, df, tgt)

    fpath = WDIR / f"{slug(label)}.feather"
    feather.write_feather(pd.DataFrame({"w": w_np}), fpath)

    return {
        "method":      label,
        "median_ms":   round(statistics.median(times)),
        "min_ms":      round(min(times)),
        "max_ms":      round(max(times)),
        "status":      _STATUS.get(ri["status"], str(ri["status"])),
        "iterations":  ri["iterations"],
        "algorithm":   ri.get("algorithm_name", kwargs.get("method", "?")),
        "max_err":     m["max_err"],
        "marginal_kl": m["marginal_kl"],
        "kl":          m["kl"],
        "chi2":        m["chi2"],
    }


cfg = [
    ("ieppa",              {"method": "ieppa"}),
    ("ieppa_accel",        {"method": "ieppa",      "accelerate": True}),
    ("ieppa_greedy",       {"method": "ieppa",      "scheduler": "greedy"}),
    ("ieppa_greedy_accel", {"method": "ieppa",      "scheduler": "greedy", "accelerate": True}),
    ("ieppa_soft",         {"method": "ieppa_soft"}),
    ("ieppa_soft_accel",   {"method": "ieppa_soft", "accelerate": True}),
    ("raking",             {"method": "raking"}),
    ("greenkhorn",         {"method": "greenkhorn"}),
    ("greenkhorn_greedy",  {"method": "greenkhorn", "scheduler": "greedy"}),
    ("greenkhorn_accel",   {"method": "greenkhorn", "accelerate": True}),
    ("sinkhorn",           {"method": "sinkhorn"}),
    ("chebyshev",          {"method": "chebyshev"}),
    ("greg",               {"method": "greg"}),
    ("logit",              {"method": "logit"}),
    ("newton_kl",          {"method": "newton_kl"}),
]

results = []
for i, (label, kwargs) in enumerate(cfg, 1):
    print(f"[{i:2d}/{len(cfg):2d}] {label:<22}", end="", flush=True)
    r = run(label, DEF, **kwargs)
    if r is not None:
        results.append(r)
        print(f"  {r['median_ms']:7d}ms  {r['status']:<10}  iters={r['iterations']}  "
              f"mkl={r['marginal_kl']:.4g}  max_err={r['max_err']:.4g}")

hdr = ["method","median_ms","min_ms","max_ms","status","iterations","algorithm",
       "max_err","marginal_kl","kl","chi2"]
OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
with open(OUT_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=hdr)
    w.writeheader()
    w.writerows(results)
print(f"\nSaved metrics: {OUT_CSV}\nSaved weights: {WDIR}/*.feather")
