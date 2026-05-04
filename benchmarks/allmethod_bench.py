#!/usr/bin/env python3
# benchmarks/allmethod_bench.py
# Python/R parity benchmark: all methods + accel/greedy on stepstone fulldata.
# Uses improvement rule (stall detection) matching allmethod_bench.R.
# Tier-1 metrics: marginal_kl + wall_time.

import sys, time, json, math, csv
from pathlib import Path

import pandas as pd
import numpy as np

repo = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(repo / "python"))
from leafblower import harvest

DATA    = repo / "benchmarks/stepstone_fulldata_bench_data.parquet"
TARGETS = repo / "benchmarks/stepstone_fulldata_bench_targets.json"
MAX_IT  = 3000
MAX_WT  = 5.0
TOL     = 1e-3
OUTFILE = repo / "benchmarks/results/allmethod_bench_py.csv"

# Natural convergence per method — improvement rule (stall detection)
KL   = {"metric": "marginal_kl", "rule": "improvement", "tol": TOL}
MERR = {"metric": "max_err",     "rule": "improvement", "tol": TOL}
WKL  = {"metric": "kl",         "rule": "improvement", "tol": TOL}
CHI  = {"metric": "chi2",       "rule": "improvement", "tol": TOL}

_STATUS = {0: "converged", 1: "no_conv", 2: "infeasible", 3: "bad_arg",
           4: "budget",    5: "stall"}

print("Loading data...")
df  = pd.read_parquet(DATA)
with open(TARGETS) as f:
    tgt_raw = json.load(f)
tgt = {k: {lv: v / sum(d.values()) for lv, v in d.items()} for k, d in tgt_raw.items()}
df  = df[[c for c in tgt if c in df.columns]]
print(f"n = {len(df):,} | margins = {len(tgt)} | max_iter = {MAX_IT} | tol = {TOL:.0e}\n")


def post_hoc_mkl(w_arr, df, tgt):
    """Compute obs-level marginal_kl from returned weights via groupby."""
    w = pd.Series(w_arr, index=df.index)
    W = w.sum()
    mkl, max_err = 0.0, 0.0
    for v, targets in tgt.items():
        if v not in df.columns: continue
        s_map = w.groupby(df[v]).sum() / W
        for lv, T in targets.items():
            S = float(s_map.get(lv, 0.0))
            max_err = max(max_err, abs(S - T))
            if T > 0 and S > 0:
                mkl += T * math.log(T / S)
    return mkl, max_err


def run(label, conv, **kwargs):
    t0 = time.perf_counter()
    try:
        res = harvest(df, tgt, min_weight=0.0, max_weight=MAX_WT,
                      max_iterations=MAX_IT, convergence=conv,
                      verbose=0, attach_weights=False, **kwargs)
    except Exception as e:
        elapsed = (time.perf_counter() - t0) * 1000
        print(f"  ERROR: {e}")
        return {"method": label, "loss_fn": conv.get("metric","?"),
                "time_ms": round(elapsed), "status": "error",
                "iterations": -1, "algorithm": "error",
                "marginal_kl": float("nan"), "max_err": float("nan")}
    elapsed = (time.perf_counter() - t0) * 1000
    ri   = res["result"]
    w_np = np.array(res["weights"])
    mkl, merr = post_hoc_mkl(w_np, df, tgt)
    return {
        "method":      label,
        "loss_fn":     conv.get("metric", "?"),
        "time_ms":     round(elapsed),
        "status":      _STATUS.get(ri["status"], str(ri["status"])),
        "iterations":  ri["iterations"],
        "algorithm":   ri.get("algorithm_name", kwargs.get("method", "?")),
        "marginal_kl": round(mkl, 6),
        "max_err":     round(merr, 6),
    }


cfg = [
    # ieppa family — natural: marginal_kl
    ("ieppa",                  KL,   {"method": "ieppa"}),
    ("ieppa + accel",          KL,   {"method": "ieppa",      "accelerate": True}),
    ("ieppa + greedy",         KL,   {"method": "ieppa",      "scheduler": "greedy"}),
    ("ieppa + greedy+accel",   KL,   {"method": "ieppa",      "scheduler": "greedy", "accelerate": True}),
    ("ieppa_soft (auto cp)",   KL,   {"method": "ieppa_soft"}),
    ("ieppa_soft + accel",     KL,   {"method": "ieppa_soft", "accelerate": True}),
    # raking family
    ("raking",                 MERR, {"method": "raking"}),
    ("greenkhorn",             MERR, {"method": "greenkhorn"}),
    ("greenkhorn + greedy",    MERR, {"method": "greenkhorn", "scheduler": "greedy"}),
    ("greenkhorn + accel",     MERR, {"method": "greenkhorn", "accelerate": True}),
    # KL/dual family
    ("sinkhorn",               WKL,  {"method": "sinkhorn"}),
    ("newton_kl",              WKL,  {"method": "newton_kl"}),
    # logit
    ("logit",                  MERR, {"method": "logit"}),
    # chi2 family (R benchmark: chebyshev ~17s, greg ~3s)
    ("chebyshev",              CHI,  {"method": "chebyshev"}),
    ("greg",                   CHI,  {"method": "greg"}),
    ## lbfgsb removed — strictly dominated by newton_kl (O(n) vs O(M_cell) per iter)
]

results = []
for i, (label, conv, kwargs) in enumerate(cfg, 1):
    print(f"[{i:2d}/{len(cfg):2d}] {label:<36}", end="", flush=True)
    r = run(label, conv, **kwargs)
    results.append(r)
    mkl = f"{r['marginal_kl']:.5f}" if math.isfinite(r['marginal_kl']) else "  NaN"
    print(f"  {r['time_ms']:7d}ms  {r['status']:<10}  mkl={mkl}  iters={r['iterations']}")

hdr = ["method","loss_fn","time_ms","status","iterations","algorithm","marginal_kl","max_err"]
print("\n=== Python Results (sorted by time_ms) ===")
col_w = [36, 12, 8, 11, 10, 14, 13, 10]
fmt = "".join(f"{{:<{w}}}" for w in col_w)
print(fmt.format(*hdr))
print("-" * sum(col_w))
for r in sorted(results, key=lambda x: x["time_ms"]):
    mkl = f"{r['marginal_kl']:.5f}" if math.isfinite(r.get('marginal_kl', float('nan'))) else "NaN"
    merr = f"{r['max_err']:.5f}" if math.isfinite(r.get('max_err', float('nan'))) else "NaN"
    print(fmt.format(r["method"], r["loss_fn"], r["time_ms"], r["status"],
                     r["iterations"], r["algorithm"], mkl, merr))

OUTFILE.parent.mkdir(exist_ok=True)
with open(OUTFILE, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=hdr)
    w.writeheader()
    w.writerows(results)
print(f"\nSaved: {OUTFILE}")
