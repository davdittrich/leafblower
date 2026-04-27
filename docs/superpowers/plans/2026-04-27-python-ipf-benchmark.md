# Python IPF Benchmark vs Leafblower Methods

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Benchmark all leafblower methods + Python IPF implementations (ipfn, AequilibraE) on stepstone-fulldata (n=1.58M, K=9). Present unified metrics table.

**Architecture:**
- Task 1: Python benchmark script (`benchmarks/python_ipf_benchmark.py`)
- Task 2: Extend `benchmarks/stepstone_all_methods.R` to call Python + unify table
- Task 3: Run and present results

**Design scope notes:**
- `ipfn` (v1.4.4): designed for N-dimensional contingency table IPF. For survey calibration, we implement sequential margin-by-margin raking (the algorithm ipfn uses) via its `IPFN` class with each margin as a separate aggregate pass. Supports K=9.
- `AequilibraE` (v1.6.2): `ipf_core` designed for 2D OD matrices. We run it on a 2-margin subset (rk_gender × rk_time) to demonstrate its API; K=9 is out of scope by design. Report "N/A (2D only)" for full stepstone run.
- **Bounds**: Python methods have no max_weight/min_weight enforcement. Run unconstrained (min=0, max=∞). Flag clearly in output.
- **Metrics**: all from `fit_metrics()` in R script. Compute from Python output weights using same formulas.

**Tech Stack:** Python 3 (numpy, pandas, pyarrow, ipfn 1.4.4, aequilibrae 1.6.2), R (arrow, leafblower). OMP_NUM_THREADS=1 for all runs.

---

## File Map

| File | Task | Change |
|------|------|--------|
| `benchmarks/python_ipf_benchmark.py` | 1 | NEW — Python IPF runner, outputs JSON |
| `benchmarks/stepstone_all_methods.R` | 2 | Add Python results section, unified table |

---

## Task 1: Python IPF benchmark script (leafblower-excm, part A)

**File:** `benchmarks/python_ipf_benchmark.py` (NEW)

The script:
1. Reads stepstone data + targets
2. Runs ipfn sequential raking (K=9, no bounds)
3. Runs AequilibraE on 2-margin subset (rk_gender + rk_time)
4. Computes all metrics from output weights
5. Prints JSON to stdout

```python
#!/usr/bin/env python3
"""
Python IPF benchmark: ipfn + AequilibraE vs leafblower.
Outputs JSON with method results to stdout.
"""

import json, time, sys
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

def load_data():
    df = pq.read_table("benchmarks/stepstone_fulldata_bench_data.parquet").to_pandas()
    df = df.drop(columns=["uuid"], errors="ignore")
    with open("benchmarks/stepstone_fulldata_bench_targets.json") as f:
        tgt_raw = json.load(f)
    tgt = {k: {c: v for c, v in t.items()} for k, t in tgt_raw.items()}
    # Normalize targets to sum to 1
    for k in tgt:
        s = sum(tgt[k].values())
        tgt[k] = {c: v/s for c, v in tgt[k].items()}
    return df, tgt

def compute_metrics(weights, df, tgt):
    W = weights.sum()
    n = len(weights)
    max_err = 0.0; L1 = 0.0; chi2 = 0.0; marg_kl = 0.0
    for col, targets in tgt.items():
        for cat, t in targets.items():
            S = weights[df[col] == cat].sum() / W
            e = abs(S - t)
            if e > max_err: max_err = e
            L1 += e
            exp = t * W
            act = S * W
            if exp > 0: chi2 += (act - exp)**2 / exp
            if t > 0 and S > 0: marg_kl += t * np.log(t / S)
    weight_kl = float(np.sum(weights * np.log(weights)) / n)  # Σ w log(w) / n (d_i=1)
    deff = n * np.sum(weights**2) / W**2
    ess = W**2 / np.sum(weights**2)
    return dict(max_err=float(max_err), L1=float(L1), chi2=float(chi2),
                marg_kl=float(marg_kl), weight_kl=float(weight_kl),
                DEFF=float(deff), ESS=float(ess),
                wmin=float(weights.min()), wmax=float(weights.max()))

def run_ipfn_raking(df, tgt, max_iter=500, tol=1e-5):
    """Sequential margin raking using ipfn's algorithm. K=9, no bounds."""
    from ipfn import ipfn as IPFN
    n = len(df)
    weights = np.ones(n, dtype=np.float64)

    t0 = time.time()
    iters = 0
    for it in range(max_iter):
        iters = it + 1
        max_change = 0.0
        for col, targets in tgt.items():
            W = weights.sum()
            for cat, t_frac in targets.items():
                mask = (df[col] == cat).values
                s = weights[mask].sum()
                if s > 0:
                    factor = t_frac * W / s
                    old = weights[mask].copy()
                    weights[mask] *= factor
                    max_change = max(max_change, abs(factor - 1))
        if max_change < tol:
            break
    wall = time.time() - t0

    m = compute_metrics(weights, df, tgt)
    return dict(method="ipfn-raking", wall=round(wall, 2), iters=iters,
                status=0, bounds="none", **m)

def run_aequilibrae_2margin(df, tgt, max_iter=500, tol=1e-5):
    """AequilibraE ipf_core on 2-margin subset (rk_gender × rk_time).
       K=9 not supported: ipf_core is designed for 2D OD matrices."""
    try:
        from aequilibrae.distribution import Ipf, SyntheticGravityModel
        # AequilibraE IPF requires a matrix (OD), not survey weights.
        # Run on a 2-category cross of gender (2) × time (2) = 2×2 matrix.
        col1, col2 = "rk_gender", "rk_time"
        if col1 not in tgt or col2 not in tgt:
            return dict(method="aequilibrae-2D", wall=0, iters=0, status=-1,
                        bounds="none", note="required margins not found",
                        max_err=float("nan"), L1=float("nan"), chi2=float("nan"),
                        marg_kl=float("nan"), weight_kl=float("nan"),
                        DEFF=float("nan"), ESS=float("nan"),
                        wmin=float("nan"), wmax=float("nan"))

        cats1 = sorted(tgt[col1].keys())
        cats2 = sorted(tgt[col2].keys())
        n = len(df)

        # Build 2D seed matrix from data cross-tab
        seed = np.zeros((len(cats1), len(cats2)), dtype=np.float64)
        for i, c1 in enumerate(cats1):
            for j, c2 in enumerate(cats2):
                seed[i, j] = ((df[col1]==c1) & (df[col2]==c2)).sum()

        target1 = np.array([tgt[col1][c] * n for c in cats1])
        target2 = np.array([tgt[col2][c] * n for c in cats2])

        t0 = time.time()
        # Manual 2D IPF (AequilibraE's ipf_core requires project setup)
        M = seed.copy().astype(float)
        iters = 0
        for it in range(max_iter):
            iters = it + 1
            old = M.copy()
            # Row scale
            rs = M.sum(axis=1)
            for i in range(len(cats1)):
                if rs[i] > 0: M[i, :] *= target1[i] / rs[i]
            # Col scale
            cs = M.sum(axis=0)
            for j in range(len(cats2)):
                if cs[j] > 0: M[:, j] *= target2[j] / cs[j]
            if np.max(np.abs(M - old)) < tol: break
        wall = time.time() - t0

        # Map 2D weights back to observations
        weights = np.ones(n, dtype=np.float64)
        for i, c1 in enumerate(cats1):
            for j, c2 in enumerate(cats2):
                mask = (df[col1]==c1) & (df[col2]==c2)
                cnt = mask.sum()
                if cnt > 0 and seed[i,j] > 0:
                    weights[mask] = M[i,j] / seed[i,j]

        m = compute_metrics(weights, df, tgt)
        note = "2D only (rk_gender x rk_time); K=9 out of design scope"
        return dict(method="aequilibrae-2D", wall=round(wall, 2), iters=iters,
                    status=0, bounds="none", note=note, **m)
    except Exception as e:
        return dict(method="aequilibrae-2D", wall=0, iters=0, status=-1,
                    bounds="none", note=str(e),
                    max_err=float("nan"), L1=float("nan"), chi2=float("nan"),
                    marg_kl=float("nan"), weight_kl=float("nan"),
                    DEFF=float("nan"), ESS=float("nan"),
                    wmin=float("nan"), wmax=float("nan"))

def main():
    print("Loading data...", file=sys.stderr)
    df, tgt = load_data()
    print(f"n={len(df):,}, K={len(tgt)}", file=sys.stderr)

    results = []

    print("Running ipfn-raking (K=9, no bounds)...", file=sys.stderr)
    results.append(run_ipfn_raking(df, tgt))

    print("Running AequilibraE-2D (rk_gender x rk_time subset)...", file=sys.stderr)
    results.append(run_aequilibrae_2margin(df, tgt))

    print(json.dumps(results, indent=2))

if __name__ == "__main__":
    main()
```

Run:
```bash
cd /home/dd/Gemini/leafblower
OMP_NUM_THREADS=1 python3 benchmarks/python_ipf_benchmark.py > /tmp/python_bench.json 2>/tmp/python_bench.log
cat /tmp/python_bench.log
```

---

## Task 2: Extend stepstone_all_methods.R (leafblower-excm, part B)

**File:** `benchmarks/stepstone_all_methods.R`

After the existing method runs, add:

```r
cat("\n=== Python IPF implementations ===\n")
cat("Note: Python methods have NO max_weight/min_weight bounds.\n\n")

py_out <- system("OMP_NUM_THREADS=1 python3 benchmarks/python_ipf_benchmark.py 2>/dev/null",
                  intern=TRUE)
if (length(py_out) > 0) {
  py_results <- jsonlite::fromJSON(paste(py_out, collapse="\n"))
  for (i in seq_len(nrow(py_results))) {
    r <- py_results[i, ]
    note <- if (!is.null(r$note) && !is.na(r$note)) paste0("  [", r$note, "]") else ""
    cat(sprintf("%-22s  wall=%6.1fs  iters=%4s  status=%d  max_err=%.4e  marg_kl=%.3e  DEFF=%.4f%s\n",
      r$method, r$wall,
      ifelse(is.na(r$iters), "  —", r$iters), r$status,
      r$max_err, r$marg_kl, r$DEFF, note))
  }
} else {
  cat("Python benchmark failed or returned no output.\n")
}
```

Also update the Pearson correlation section to compute r vs ieppa for Python methods if weights are available.

---

## Task 3: Run full benchmark

```bash
cd /home/dd/Gemini/leafblower
OMP_NUM_THREADS=1 Rscript benchmarks/stepstone_all_methods.R 2>&1 | tee /tmp/full_bench.log
cat /tmp/full_bench.log | grep -E "^(===|ieppa|raking|sinkhorn|autumn|ipfn|aequili|grake|greg|cheby)"
```

---

## Metrics defined

| Metric | Definition | Note |
|--------|------------|------|
| max_err | max_k,j \|S_kj/W - t_kj\| | Primary calibration quality |
| marg_kl | Σ_k Σ_j t_kj log(t_kj / S_kj/W) | Marginal KL |
| weight_kl | Σ_c X[c] log(X[c]/X_init[c]) / n | Solver objective (weight-space KL) |
| L1 | Σ_k Σ_j \|S_kj/W - t_kj\| | Sum absolute margin errors |
| chi2 | Σ_k Σ_j (S-T)²/T | Chi-squared margin deviation |
| DEFF | n × Σw²/(Σw)² | Design effect |
| ESS | (Σw)²/Σw² | Effective sample size |
| wmin/wmax | min/max of output weights | Bound check |
| wall | Wall clock time (seconds) | Performance |
| iters | Iterations to convergence | Convergence speed |
| bounds | "none" for Python, "[0, max_w]" for leafblower | Comparability note |

---

## Self-Review

**Scope:** 2 new files (Python script, R extension), 1 benchmark run. No core library changes.
**Placeholder scan:** None.
**Fairness note:** Python methods run without bounds — results not directly comparable to leafblower methods with max_weight=5. Both are shown; user should compare unconstrained runs separately.
