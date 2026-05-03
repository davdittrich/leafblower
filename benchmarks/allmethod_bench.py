#!/usr/bin/env python3
# benchmarks/allmethod_bench.py
# Full-method comparison on stepstone fulldata (1.58M rows, 9 margins).
# MAX_ITER=3000, absolute tol=1e-3 — mirrors allmethod_bench.R settings.

import sys, time, json, csv, math
from pathlib import Path

import pandas as pd

repo = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(repo / "python"))
from leafblower import harvest

DATA    = repo / "benchmarks/stepstone_fulldata_bench_data.parquet"
TARGETS = repo / "benchmarks/stepstone_fulldata_bench_targets.json"
MAX_IT  = 3000
MAX_WT  = 5.0
TOL     = 1e-3
OUTFILE = repo / "benchmarks/results/allmethod_bench_py.csv"

_STATUS = {0: "converged", 1: "no_conv", 2: "infeasible", 3: "bad_arg",
           4: "budget",    5: "stall"}
_ALG    = {0: "auto", 1: "ieppa", 2: "lbfgsb", 3: "raking", 4: "sinkhorn",
           5: "chebyshev", 6: "greg", 8: "ieppa_soft", 9: "greenkhorn",
           10: "logit", 11: "newton_kl"}

print("Loading data...")
df  = pd.read_parquet(DATA)
with open(TARGETS) as f:
    tgt_raw = json.load(f)
# normalize targets (match R behavior)
tgt = {k: {lv: v / sum(d.values()) for lv, v in d.items()} for k, d in tgt_raw.items()}
# drop non-margin columns
df  = df[[c for c in tgt]]
print(f"n = {len(df):,} | margins = {len(tgt)} | max_iter = {MAX_IT} | tol = {TOL:.0e}\n")


def fmt(v, sig=4):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "NaN"
    return f"{v:.{sig}g}"


def run(label, **kwargs):
    t0 = time.perf_counter()
    try:
        res = harvest(
            df, tgt,
            min_weight    = 0.0,
            max_weight    = MAX_WT,
            max_iterations = MAX_IT,
            convergence   = {"absolute": TOL},
            verbose       = 0,
            attach_weights = False,
            **kwargs,
        )
    except Exception as e:
        elapsed = (time.perf_counter() - t0) * 1000
        print(f"  ERROR: {e}")
        return {"method": label, "time_ms": round(elapsed), "status": "error",
                "iterations": -1, "algorithm": "error",
                "max_error": float("nan"), "mean_error": float("nan"),
                "kl": float("nan"), "marginal_kl": float("nan"),
                "chi2": float("nan"), "grake_norm": float("nan"),
                "l1_weight": float("nan"), "solver_obj": float("nan")}
    elapsed = (time.perf_counter() - t0) * 1000
    ri = res["result"]
    cu = ri.get("convergence_used", {})
    return {
        "method":     label,
        "time_ms":    round(elapsed),
        "status":     _STATUS.get(ri["status"], str(ri["status"])),
        "iterations": ri["iterations"],
        "algorithm":  _ALG.get(ri["algorithm_used"], str(ri["algorithm_used"])),
        "max_error":  ri.get("max_error",  float("nan")),
        "mean_error": ri.get("mean_error", float("nan")),
        "kl":         ri.get("kl",         float("nan")),
        "marginal_kl": ri.get("marginal_kl", float("nan")),
        "chi2":       ri.get("chi2",        float("nan")),
        "grake_norm": ri.get("grake_norm",  float("nan")),
        "l1_weight":  ri.get("l1_weight",   float("nan")),
        "solver_obj": cu.get("objective",   float("nan")),
    }


cfg = [
    ("ieppa",                  {"method": "ieppa"}),
    ("ieppa + accel",          {"method": "ieppa",      "accelerate": True}),
    ("ieppa + greedy",         {"method": "ieppa",      "scheduler": "greedy"}),
    ("ieppa + greedy + accel", {"method": "ieppa",      "scheduler": "greedy", "accelerate": True}),
    ("raking",                 {"method": "raking"}),
    ("greenkhorn",             {"method": "greenkhorn"}),
    ("greenkhorn + greedy",    {"method": "greenkhorn", "scheduler": "greedy"}),
    ("greenkhorn + accel",     {"method": "greenkhorn", "accelerate": True}),
    ("sinkhorn",               {"method": "sinkhorn"}),
    ("chebyshev",              {"method": "chebyshev"}),
    ("greg",                   {"method": "greg"}),
    ("logit",                  {"method": "logit"}),
    ("newton_kl",              {"method": "newton_kl"}),
    ("lbfgsb",                 {"method": "lbfgsb"}),
    ("ieppa_soft (cp=1)",      {"method": "ieppa_soft", "capacity_penalty": 1.0}),
]

results = []
for i, (label, kwargs) in enumerate(cfg, 1):
    print(f"[{i:2d}/{len(cfg):2d}] {label:<38}", end="", flush=True)
    r = run(label, **kwargs)
    results.append(r)
    print(f"  {r['time_ms']:6d}ms  {r['status']:<10}  max_err={fmt(r['max_error'])}  iters={r['iterations']}")

hdr = ["method","time_ms","status","iterations","algorithm",
       "max_error","mean_error","kl","marginal_kl","chi2","grake_norm","l1_weight","solver_obj"]

print("\n=== Results ===")
col_w = [38, 8, 11, 10, 13] + [12]*8
fmt_row = "".join(f"{{:<{w}}}" for w in col_w)
print(fmt_row.format(*hdr))
print("-" * sum(col_w))
for r in results:
    print(fmt_row.format(
        r["method"], r["time_ms"], r["status"], r["iterations"], r["algorithm"],
        fmt(r["max_error"]), fmt(r["mean_error"]), fmt(r["kl"]),
        fmt(r["marginal_kl"]), fmt(r["chi2"]), fmt(r["grake_norm"]),
        fmt(r["l1_weight"]), fmt(r["solver_obj"])))

OUTFILE.parent.mkdir(exist_ok=True)
with open(OUTFILE, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=hdr)
    w.writeheader()
    w.writerows(results)
print(f"\nSaved: {OUTFILE}")
