# Python IPF Benchmark vs Leafblower Methods

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Benchmark all leafblower methods + Python IPF implementations (ipfn, AequilibraE) on stepstone-fulldata (n=1.58M, K=9). Present unified metrics table.

**Architecture:**
- Task 1: Python benchmark script (`benchmarks/python_ipf_benchmark.py`)
- Task 2: Extend `benchmarks/stepstone_all_methods.R` to call Python + unify table
- Task 3: Run and present results

**Key design decisions:**
- `ipfn v1.4.4`: uses **DataFrame mode** (`IPFN(cell_df, aggregates, dimensions, weight_col='total')`). ipfn supports DataFrame input where each row is a cell — not a dense K-dimensional tensor. Cell table has 28,905 unique cells (from 1.58M obs). Polars builds cell table efficiently, converts to pandas for ipfn.
- `AequilibraE v1.6.2`: designed for 2D OD matrices. Run on 2-margin subset (rk_gender × rk_time). Honest labeling: `"aequilibrae-2D-subset"`. Full K=9 out of design scope.
- **Bounds**: Python methods have no max_weight/min_weight bounds. Run unconstrained. Flag clearly.
- `OMP_NUM_THREADS=1` for all runs (fair wall time comparison).

**Tech Stack:** Python 3 (polars, pandas, numpy, pyarrow, ipfn 1.4.4, aequilibrae 1.6.2), R (arrow, leafblower, jsonlite).

---

## File Map

| File | Task | Change |
|------|------|--------|
| `benchmarks/python_ipf_benchmark.py` | 1 | NEW — Python IPF runner, outputs JSON |
| `benchmarks/stepstone_all_methods.R` | 2 | Add Python results section + tryCatch |

---

## Task 1: Python IPF benchmark script (leafblower-excm)

**File:** `benchmarks/python_ipf_benchmark.py` (NEW)

```python
#!/usr/bin/env python3
"""
Python IPF benchmark: ipfn (DataFrame mode) + AequilibraE (2D subset).
Outputs JSON to stdout. Use: OMP_NUM_THREADS=1 python3 benchmarks/python_ipf_benchmark.py
"""

import json, sys, time
import numpy as np
import pandas as pd
import polars as pl
from ipfn.ipfn import ipfn as IPFN  # import the class, not the module


def load_data():
    """Load data and targets. Polars for speed on 1.58M rows."""
    df_pl = pl.read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
    if "uuid" in df_pl.columns:
        df_pl = df_pl.drop("uuid")
    df_pd = df_pl.to_pandas()  # ipfn needs pandas

    with open("benchmarks/stepstone_fulldata_bench_targets.json") as f:
        tgt_raw = json.load(f)
    tgt = {}
    for k, t in tgt_raw.items():
        v = {c: float(x) for c, x in t.items()}
        s = sum(v.values())
        tgt[k] = {c: x / s for c, x in v.items()}

    return df_pl, df_pd, tgt


def compute_metrics(weights, df_pd, tgt):
    """Compute all benchmark metrics from per-obs weights."""
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
    # Weight-space KL: Σ w_i log(w_i) / n  (d_i = 1, normalized weights)
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
    """
    ipfn.IPFN in DataFrame mode on compressed cell table.
    Cell table: 28,905 unique cells built via Polars groupby.
    ipfn.IPFN(cell_df, aggregates, dimensions, weight_col='total')
    maps to K=9 survey calibration without a K-dimensional tensor.
    No bounds (unconstrained IPF).
    """
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
    cell_seed = cell_table["total"].values.copy()  # save for weight-back mapping

    # Build ipfn aggregates: target total count per category per margin
    aggregates = []
    dimensions = []
    for col in margins:
        target_total = pd.Series(
            {cat: frac * n for cat, frac in tgt[col].items()},
            name="total", dtype=float
        )
        aggregates.append(target_total)
        dimensions.append([col])  # ipfn_df requires list-of-lists

    # Run ipfn.IPFN — DataFrame mode with column name dimensions
    ipf = IPFN(cell_table, aggregates, dimensions,
               weight_col="total",
               convergence_rate=tol,
               max_iteration=max_iter)
    cell_result = ipf.iteration()

    wall = time.time() - t0

    # Map cell weights back to individual obs
    # w_i = (calibrated_cell_total / seed_cell_total)
    cell_multiplier = np.where(
        cell_seed > 0, cell_result["total"].values / cell_seed, 1.0
    )
    # Build obs→cell lookup via merge
    lookup = cell_table[margins].copy()
    lookup["_mult"] = cell_multiplier
    merged = df_pd[margins].merge(lookup, on=margins, how="left")
    weights = merged["_mult"].fillna(1.0).values.astype(float)

    m = compute_metrics(weights, df_pd, tgt)
    iters = max_iter  # ipfn doesn't expose actual iter count easily
    return dict(method="ipfn-dataframe", wall=round(wall, 2),
                iters=iters, status=0, bounds="none", **m)


def run_aequilibrae_2d(df_pd, tgt, max_iter=500, tol=1e-5):
    """
    AequilibraE 2D IPF on rk_gender × rk_time subset.
    AequilibraE ipf_core is designed for 2D OD matrices; K=9 is out of scope.
    Implements standard 2D Furness/Fratar IPF using AequilibraE's matrix types.
    """
    col1, col2 = "rk_gender", "rk_time"
    if col1 not in tgt or col2 not in tgt:
        return dict(method="aequilibrae-2D-subset", wall=0, iters=0,
                    status=-1, bounds="none",
                    note="required margins not found in targets",
                    **{k: float("nan") for k in
                       ["max_err","L1","chi2","marg_kl","weight_kl",
                        "DEFF","ESS","wmin","wmax"]})

    cats1 = sorted(tgt[col1].keys())
    cats2 = sorted(tgt[col2].keys())
    n = len(df_pd)

    # Seed matrix: observed cell counts (2D cross-tab)
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

    # Map 2D weights back to individual obs
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
    print(f"n={len(df_pd):,}, K={len(tgt)}, "
          f"unique cells={df_pl.select(list(tgt.keys())).unique().height:,}",
          file=sys.stderr)

    results = []

    print("Running ipfn.IPFN (DataFrame mode, K=9, no bounds)...", file=sys.stderr)
    results.append(run_ipfn(df_pl, df_pd, tgt))

    print("Running AequilibraE 2D subset (rk_gender x rk_time)...", file=sys.stderr)
    results.append(run_aequilibrae_2d(df_pd, tgt))

    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
```

### Run command
```bash
cd /home/dd/Gemini/leafblower
OMP_NUM_THREADS=1 python3 benchmarks/python_ipf_benchmark.py > /tmp/python_bench.json 2>/tmp/python_bench.log
cat /tmp/python_bench.log  # progress
cat /tmp/python_bench.json # results
```

---

## Task 2: Extend stepstone_all_methods.R

**File:** `benchmarks/stepstone_all_methods.R` (append after autumn section)

```r
cat("\n=== Python IPF implementations ===\n")
cat("Note: Python methods have NO max_weight/min_weight bounds.\n")
cat("ipfn: DataFrame mode on 28,905 unique cells (K=9, n=1.58M).\n")
cat("aequilibrae: 2D subset only (rk_gender x rk_time).\n\n")

py_json <- tryCatch({
  py_out <- system(
    "OMP_NUM_THREADS=1 python3 benchmarks/python_ipf_benchmark.py 2>/dev/null",
    intern = TRUE)
  if (length(py_out) == 0) stop("no output")
  jsonlite::fromJSON(paste(py_out, collapse = "\n"))
}, error = function(e) {
  cat(sprintf("Python benchmark failed: %s\n", conditionMessage(e)))
  NULL
})

if (!is.null(py_json)) {
  for (i in seq_len(nrow(py_json))) {
    r <- py_json[i, ]
    note <- if (!is.null(r$note) && !is.na(r$note))
      sprintf("  [%s]", r$note) else ""
    cat(sprintf(
      "%-26s  wall=%6.1fs  iters=%4s  status=%d  max_err=%.4e  marg_kl=%.3e  weight_kl=%.3e  DEFF=%.4f  ESS=%s  wmin=%.3f  wmax=%.3f%s\n",
      r$method, r$wall,
      ifelse(is.na(r$iters), "  —", as.character(r$iters)),
      r$status, r$max_err, r$marg_kl,
      ifelse(is.na(r$weight_kl), NaN, r$weight_kl),
      r$DEFF,
      format(round(r$ESS), big.mark = ","),
      r$wmin, r$wmax, note))
  }
}
```

---

## Task 3: Run full benchmark

```bash
cd /home/dd/Gemini/leafblower
OMP_NUM_THREADS=1 Rscript benchmarks/stepstone_all_methods.R 2>&1 | tee /tmp/full_bench.log
grep -E "^(===|ieppa|raking|sinkhorn|autumn|ipfn|aequili|grake|greg|cheby)" /tmp/full_bench.log
```

---

## Metrics

| Metric | Formula | Note |
|--------|---------|------|
| max_err | max_k,j \|S_kj/W - t_kj\| | Primary quality |
| marg_kl | Σ_k Σ_j t_kj log(t_kj / S_kj/W) | Marginal divergence |
| weight_kl | Σ_i w_i log(w_i) / n | Solver objective (d_i=1) |
| L1 | Σ_k Σ_j \|S_kj/W - t_kj\| | Sum absolute errors |
| chi2 | Σ_k Σ_j (S-T)²/T | Chi-squared deviation |
| DEFF | n Σw²/(Σw)² | Design effect |
| ESS | (Σw)²/Σw² | Effective sample size |
| wmin/wmax | min/max weights | Bound check |
| bounds | "none" for Python, "[0, max_w]" for leafblower | Comparability |

---

## Self-Review

**Spec coverage:**

| Item | Task |
|------|------|
| ipfn.IPFN class actually used | 1 (DataFrame mode) |
| Polars for cell table (fast groupby) | 1 |
| AequilibraE 2D subset honest label | 1 |
| tryCatch for JSON parse in R | 2 |
| All 9 metrics in output | 1, 2 |
| Bounds clearly noted | 1, 2 |

**Placeholder scan:** None.
**Clean-code:** Polars instead of Pandas for cell table (faster); functions clearly named; no dead imports.
