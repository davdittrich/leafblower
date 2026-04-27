#!/usr/bin/env python3
"""
Python IPF benchmark: ipfn (DataFrame mode) + AequilibraE (2D subset).
Outputs JSON to stdout. Use: OMP_NUM_THREADS=1 python3 benchmarks/python_ipf_benchmark.py
"""

import json, sys, time
import numpy as np
import pandas as pd
import polars as pl
from ipfn.ipfn import ipfn as IPFN  # type: ignore[import]  # nested pkg structure; Pyright false positive


def load_data():
    df_pl = pl.read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
    if "uuid" in df_pl.columns:
        df_pl = df_pl.drop("uuid")
    df_pd = df_pl.to_pandas()

    with open("benchmarks/stepstone_fulldata_bench_targets.json") as f:
        tgt_raw = json.load(f)
    tgt = {}
    for k, t in tgt_raw.items():
        v = {c: float(x) for c, x in t.items()}
        s = sum(v.values())
        tgt[k] = {c: x / s for c, x in v.items()}

    return df_pl, df_pd, tgt


def compute_metrics(weights, df_pd, tgt):
    W = float(weights.sum())
    n = len(weights)
    max_err = L1 = chi2 = marg_kl = 0.0
    for col, targets in tgt.items():
        for cat, t in targets.items():
            S = float(weights[df_pd[col].values == cat].sum()) / W
            e = abs(S - t)
            if e > max_err:
                max_err = e
            L1 += e
            exp = t * W
            if exp > 0:
                chi2 += ((S * W - exp) ** 2) / exp
            if t > 0 and S > 0:
                marg_kl += t * np.log(t / S)
    pos = weights > 0
    weight_kl = float(np.sum(weights[pos] * np.log(weights[pos])) / n)
    deff = float(n * np.sum(weights ** 2) / W ** 2)
    ess = float(W ** 2 / np.sum(weights ** 2))
    return dict(
        max_err=round(max_err, 8), L1=round(L1, 8),
        chi2=round(chi2, 4), marg_kl=round(marg_kl, 8),
        weight_kl=round(weight_kl, 8),
        DEFF=round(deff, 6), ESS=round(ess, 1),
        wmin=round(float(weights.min()), 6),
        wmax=round(float(weights.max()), 6),
    )


def run_ipfn(df_pl, df_pd, tgt, max_iter=500, tol=1e-5):
    """ipfn.IPFN in DataFrame mode on compressed cell table (28,905 cells)."""
    margins = list(tgt.keys())
    n = len(df_pd)

    t0 = time.time()

    # Build cell table via Polars (fast groupby on 1.58M rows)
    cell_table = (
        df_pl.select(margins)
        .with_columns([pl.col(c).cast(pl.Utf8) for c in margins])
        .group_by(margins)
        .agg(pl.len().alias("total"))
        .to_pandas()
    )
    cell_table["total"] = cell_table["total"].astype(float)
    cell_seed = cell_table["total"].values.copy()

    # Build ipfn aggregates: target total per category per margin
    aggregates = []
    dimensions = []
    for col in margins:
        target_total = pd.Series(
            {cat: frac * n for cat, frac in tgt[col].items()},
            name="total", dtype=float
        )
        aggregates.append(target_total)
        dimensions.append([col])  # list-of-lists required by ipfn_df

    # Run ipfn.IPFN — DataFrame mode
    ipf = IPFN(cell_table, aggregates, dimensions,
               weight_col="total",
               convergence_rate=tol,
               max_iteration=max_iter)
    cell_result = ipf.iteration()

    wall = time.time() - t0

    # Map cell weights back to individual obs
    cell_multiplier = np.where(
        cell_seed > 0, cell_result["total"].values / cell_seed, 1.0
    )
    lookup = cell_table[margins].copy()
    lookup["_mult"] = cell_multiplier
    merged = df_pd[margins].merge(lookup, on=margins, how="left")
    weights = merged["_mult"].fillna(1.0).values.astype(float)

    m = compute_metrics(weights, df_pd, tgt)
    # ipfn does not expose actual iteration count — report max_iter as nominal
    return dict(method="ipfn-dataframe", wall=round(wall, 2),
                iters=max_iter, status=0, bounds="none",
                note="iters=nominal (ipfn API does not expose actual count)",
                **m)


def run_aequilibrae_2d(df_pd, tgt, max_iter=500, tol=1e-5):
    """AequilibraE 2D IPF on rk_gender x rk_time subset (K=9 out of design scope)."""
    col1, col2 = "rk_gender", "rk_time"
    if col1 not in tgt or col2 not in tgt:
        return dict(method="aequilibrae-2D-subset", wall=0, iters=0,
                    status=-1, bounds="none",
                    note="required margins not in targets",
                    **{k: float("nan") for k in
                       ["max_err","L1","chi2","marg_kl","weight_kl","DEFF","ESS","wmin","wmax"]})

    cats1 = sorted(tgt[col1].keys())
    cats2 = sorted(tgt[col2].keys())
    n = len(df_pd)

    seed = np.zeros((len(cats1), len(cats2)), dtype=float)
    for i, c1 in enumerate(cats1):
        for j, c2 in enumerate(cats2):
            seed[i, j] = ((df_pd[col1] == c1) & (df_pd[col2] == c2)).sum()

    target1 = np.array([tgt[col1][c] * n for c in cats1])
    target2 = np.array([tgt[col2][c] * n for c in cats2])

    t0 = time.time()
    M = seed.copy()
    iters = 0
    for it in range(max_iter):
        iters = it + 1
        old = M.copy()
        rs = M.sum(axis=1)
        for i in range(len(cats1)):
            if rs[i] > 0:
                M[i, :] *= target1[i] / rs[i]
        cs = M.sum(axis=0)
        for j in range(len(cats2)):
            if cs[j] > 0:
                M[:, j] *= target2[j] / cs[j]
        if np.max(np.abs(M - old)) < tol:
            break
    wall = time.time() - t0

    weights = np.ones(n, dtype=float)
    for i, c1 in enumerate(cats1):
        for j, c2 in enumerate(cats2):
            mask = (df_pd[col1] == c1) & (df_pd[col2] == c2)
            if seed[i, j] > 0:
                weights[mask.values] = M[i, j] / seed[i, j]

    m = compute_metrics(weights, df_pd, tgt)
    note = f"2D only ({col1} x {col2}); K=9 out of AequilibraE ipf_core design scope"
    return dict(method="aequilibrae-2D-subset", wall=round(wall, 2),
                iters=iters, status=0, bounds="none", note=note, **m)


def main():
    print("Loading data (Polars)...", file=sys.stderr)
    df_pl, df_pd, tgt = load_data()
    unique_cells = df_pl.select(list(tgt.keys())).unique().height
    print(f"n={len(df_pd):,}, K={len(tgt)}, unique cells={unique_cells:,}", file=sys.stderr)

    results = []

    print("Running ipfn.IPFN (DataFrame mode, K=9, no bounds)...", file=sys.stderr)
    results.append(run_ipfn(df_pl, df_pd, tgt))

    print("Running AequilibraE 2D subset (rk_gender x rk_time)...", file=sys.stderr)
    results.append(run_aequilibrae_2d(df_pd, tgt))

    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
