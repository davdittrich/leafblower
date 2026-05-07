#!/usr/bin/env python3
"""
benchmarks/yh0l/make_small_hicard.py
T3: synthesize small high-cardinality DGP (n≈5000, K=10 margins, ≥20 outer iters at tol=1e-4).

Design:
  - n=5000 observations
  - 10 margins, each with 10-15 categories = 100-150 total cells
  - Biased sample (5-15% relative bias) to force ≥20 iters
  - No NA bins (simplest path — same as stepstone DGPs)
  - All categories present in sample

Outputs:
  benchmarks/yh0l/small_hicard_data.parquet
  benchmarks/yh0l/small_hicard_targets.json
"""

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import json, sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
YDIR = ROOT / "benchmarks" / "yh0l"
YDIR.mkdir(parents=True, exist_ok=True)

DATA_PATH    = YDIR / "small_hicard_data.parquet"
TARGETS_PATH = YDIR / "small_hicard_targets.json"

RNG = np.random.default_rng(42)

N   = 5000
K   = 10   # number of margins
# Categories per margin: 10, 12, 11, 13, 10, 12, 11, 14, 10, 11 = 114 total
CAT_SIZES = [10, 12, 11, 13, 10, 12, 11, 14, 10, 11]
assert len(CAT_SIZES) == K

def make_margin(k, n_cats):
    """True population proportions (roughly uniform with small variation)."""
    alpha = np.ones(n_cats) * 2.0  # Dirichlet concentration > 1 → not too spiky
    props = RNG.dirichlet(alpha)
    return props


def bias_sample_probs(true_props, bias_strength=0.12):
    """
    Biased sample distribution: multiply each cell prob by exp(b) where
    b ~ Uniform(-bias_strength, bias_strength). Renormalize.
    Ensures ≥5% bias in some margins → forces solver to work.
    """
    bias = np.exp(RNG.uniform(-bias_strength, bias_strength, size=len(true_props)))
    biased = true_props * bias
    return biased / biased.sum()


def main():
    margins = {}
    targets_dict = {}

    for k in range(K):
        vname   = f"v{k:02d}"
        n_cats  = CAT_SIZES[k]
        levels  = [f"L{k:02d}_{j:02d}" for j in range(n_cats)]

        # True population proportions = target
        true_props  = make_margin(k, n_cats)
        # Sample distribution = biased
        samp_props  = bias_sample_probs(true_props, bias_strength=0.15)

        margins[vname] = {
            "levels":      levels,
            "true_props":  true_props,
            "samp_props":  samp_props,
        }
        targets_dict[vname] = {lv: float(p) for lv, p in zip(levels, true_props)}

    # Generate independent draws per margin (cross-tab not required for ieppa_soft)
    data = {}
    for vname, m in margins.items():
        codes = RNG.choice(len(m["levels"]), size=N, p=m["samp_props"])
        data[vname] = [m["levels"][c] for c in codes]

    df = pd.DataFrame(data)
    print(f"Small hicard DGP: n={N} K={K} total_cats={sum(CAT_SIZES)}")

    # Verify all levels present in sample
    for vname, m in margins.items():
        present = set(df[vname].unique())
        missing = set(m["levels"]) - present
        if missing:
            print(f"  WARNING: {vname} missing levels: {missing} — filling with rare obs")
            for i, lv in enumerate(missing):
                df.loc[i, vname] = lv  # force at least 1 obs per level

    # Check bias magnitude
    print("\nBias check (max |sample_prop - target_prop| per margin):")
    for vname, m in margins.items():
        max_bias = float(np.max(np.abs(m["samp_props"] - m["true_props"])))
        print(f"  {vname}: max_bias={max_bias:.4f}")

    # Save
    df.to_parquet(DATA_PATH)
    with open(TARGETS_PATH, "w") as f:
        json.dump(targets_dict, f, indent=2)

    print(f"\nSaved:")
    print(f"  {DATA_PATH}")
    print(f"  {TARGETS_PATH}")

    # Quick sanity: run ieppa_soft and check iter count
    print("\nSanity check: running ieppa_soft to verify ≥20 outer iters...")
    sys.path.insert(0, str(ROOT / "python"))
    from leafblower import harvest
    tgt = {k: {lv: p for lv, p in d.items()} for k, d in targets_dict.items()}
    res = harvest(df, tgt, method="ieppa_soft",
                  convergence={"tol": 1e-4}, max_weight=5.0, max_iterations=3000,
                  verbose=0, attach_weights=False)
    rd    = res.get("result", res)
    iters = rd.get("iterations", -1)
    st    = rd.get("status", -1)
    me    = rd.get("max_error", float("nan"))
    print(f"  status={st} iterations={iters} max_error={me:.3e}")
    if iters < 20:
        print(f"  WARNING: only {iters} iters — DGP may not be challenging enough")
    else:
        print(f"  OK: {iters} >= 20 outer iters confirmed")


if __name__ == "__main__":
    main()
