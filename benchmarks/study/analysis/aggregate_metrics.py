# STUDY-BRANCH-ONLY-DO-NOT-MERGE
"""benchmarks/study/analysis/aggregate_metrics.py -- WU-12a (leafblower-2ouc.40)

Walks results/runs.parquet + results/weights/*.parquet and computes one
quality-metrics row per (solver, problem, thread, build) cell that has REAL
weights, plus pairwise RQ5 agreement rows -- REUSING common/metrics.py's
compute_metrics()/native_divergence()/agreement() verbatim (no new metric
formulas; see common/metrics.py for the formula definitions and rationale).

Writes:
    benchmarks/study/results/metrics.parquet
    benchmarks/study/results/agreement.parquet

Analysis/read-only: does NOT edit any leafblower-package file, or any
benchmarks/study/common/ or benchmarks/study/spec/ file.

Run (repo root, single-thread BLAS -- see /home/dd/Gemini/leafblower/CLAUDE.md):
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        python/.venv/bin/python benchmarks/study/analysis/aggregate_metrics.py
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import json  # noqa: E402
import math  # noqa: E402
import sys  # noqa: E402
from itertools import combinations  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any  # noqa: E402

import numpy as np  # noqa: E402
import polars as pl  # noqa: E402

_THIS_DIR = Path(__file__).resolve().parent                # benchmarks/study/analysis
STUDY_DIR = _THIS_DIR.parent                                # benchmarks/study
_COMMON = str(STUDY_DIR / "common")
if _COMMON not in sys.path:
    sys.path.insert(0, _COMMON)

import instance_family  # noqa: E402
import metrics  # noqa: E402
import problem_io  # noqa: E402

instance_family.install_gen_resolver()  # gen:instance_family data_refs, BEFORE any load_problem_spec

SPEC_DIR = STUDY_DIR / "spec"
RESULTS_DIR = STUDY_DIR / "results"
WEIGHTS_DIR = RESULTS_DIR / "weights"
REGISTRY_PATH = STUDY_DIR / "registry.json"
RUNS_PATH = RESULTS_DIR / "runs.parquet"
APISTRAT_PARQUET = RESULTS_DIR / "_problem_data" / "canonical_survey_apistrat.parquet"
APISTRAT_SPEC = SPEC_DIR / "canonical_survey_apistrat.json"
METRICS_OUT = RESULTS_DIR / "metrics.parquet"
AGREEMENT_OUT = RESULTS_DIR / "agreement.parquet"

STRICTLY_CONVEX = metrics.KL_STRICTLY_CONVEX_FAMILIES  # {"kl", "chi2", "logit"} (frozen in metrics.py)


# ---- registry -----------------------------------------------------------

def load_registry() -> dict[str, dict[str, Any]]:
    with open(REGISTRY_PATH) as f:
        return json.load(f)["solvers"]


def family_primary_for(entry: dict[str, Any]) -> str:
    """Singleton families -> that family. Multi-family dispatch solvers
    (gecal/icarus/weightit/regenesees/laeken) -> 'dispatch' (plan §12a step 3)."""
    families = entry["families"]
    return families[0] if len(families) == 1 else "dispatch"


def base_solver_name(solver: str) -> str:
    """Strip the leafblower lang suffix (_r/_py); competitor names pass through."""
    if solver.endswith("_r"):
        return solver[:-2]
    if solver.endswith("_py"):
        return solver[:-3]
    return solver


# ---- problem spec resolution ---------------------------------------------

def spec_path_for(problem_id: str) -> Path | None:
    root = SPEC_DIR / f"{problem_id}.json"
    if root.exists():
        return root
    inst = SPEC_DIR / "instance_family" / f"{problem_id}.json"
    if inst.exists():
        return inst
    return None


def _load_apistrat_problem() -> dict[str, Any]:
    """pkg:survey::apistrat is unresolvable in Python (R-only canonical
    dataset, problem_io.py docstring) -- ticket 12d materialized an R sidecar
    parquet (stype + pw); synthesize the standardized problem object from it
    + the spec JSON's targets/bounds/tol (same normalization problem_io.py
    itself applies, mirrored here since we bypass load_problem_spec)."""
    with open(APISTRAT_SPEC) as f:
        spec = json.load(f)
    df = pl.read_parquet(APISTRAT_PARQUET)
    groups = {"stype": df["stype"].to_numpy()}
    design_weights = df["pw"].to_numpy().astype(np.float64)
    targets_raw = spec["targets"]["stype"]
    total = float(sum(targets_raw.values()))
    targets = {"stype": {k: v / total for k, v in targets_raw.items()}}
    bounds_spec = spec["bounds"]
    bounds = dict(
        min=0.0 if bounds_spec.get("min") is None else float(bounds_spec["min"]),
        max=math.inf if bounds_spec.get("max") is None else float(bounds_spec["max"]),
    )
    return dict(groups=groups, targets=targets, design_weights=design_weights, bounds=bounds)


def load_problem(problem_id: str) -> dict[str, Any] | None:
    if problem_id == "canonical_survey_apistrat":
        return _load_apistrat_problem()
    spec_path = spec_path_for(problem_id)
    if spec_path is None:
        return None
    raw = problem_io.load_problem_spec(spec_path)
    groups = {m: raw["data"][m].to_numpy() for m in raw["margins"]}
    return dict(groups=groups, targets=raw["targets"], design_weights=raw["design_weights"], bounds=raw["bounds"])


def bounds_for_metrics(bounds: dict[str, float]) -> dict[str, float]:
    """compute_metrics' bound_violation reads bounds['L']/['U'] (default
    0.0/inf); native_divergence's logit branch reads bounds['min']/['max']
    (KeyError if absent). One dict satisfying both call sites."""
    lo, hi = bounds["min"], bounds["max"]
    return dict(L=lo, U=hi, min=lo, max=hi)


# ---- weight loading (real-weight gate) ------------------------------------

def load_real_weight(weights_ref: str | None) -> np.ndarray | None:
    """REAL iff the file exists, height>1, and ALL entries finite (no NaN/±inf)."""
    if not weights_ref:
        return None
    path = WEIGHTS_DIR / os.path.basename(weights_ref)
    if not path.exists():
        return None
    df = pl.read_parquet(path)
    if df.height <= 1:
        return None
    w = df["weight"].to_numpy().astype(np.float64)
    if not np.all(np.isfinite(w)):
        return None
    return w


# ---- per-problem cell table -------------------------------------------

def build_cells(sub: pl.DataFrame) -> pl.DataFrame:
    return (
        sub.sort(["solver", "thread", "build", "rep"])
        .group_by(["solver", "thread", "build"], maintain_order=True)
        .agg(
            pl.col("weights_ref").first().alias("weights_ref"),
            pl.col("status").first().alias("status"),
            pl.len().alias("rep_count"),
            pl.col("wall_time_s").median().alias("wall_time_s_median"),
            pl.col("wall_time_s").min().alias("wall_time_s_min"),
            pl.col("wall_time_s").max().alias("wall_time_s_max"),
            pl.col("peak_rss_bytes").max().alias("peak_rss_bytes"),
            pl.col("iterations").first().alias("iterations"),
        )
    )


# ---- main aggregation ------------------------------------------------

def main() -> None:
    registry = load_registry()
    df = pl.read_parquet(RUNS_PATH)
    problem_ids = sorted(df.select("problem").unique().to_series().to_list())

    metrics_rows: list[dict[str, Any]] = []
    agreement_rows: list[dict[str, Any]] = []
    problem_failures: list[tuple[str, str]] = []
    cell_warnings: list[tuple[str, str, str]] = []

    for pid in problem_ids:
        try:
            problem = load_problem(pid)
        except Exception as e:  # gen:/pkg: resolution errors, malformed spec, etc.
            problem_failures.append((pid, f"{type(e).__name__}: {e}"))
            continue
        if problem is None:
            problem_failures.append((pid, "no spec file found under spec/ or spec/instance_family/"))
            continue

        bnds = bounds_for_metrics(problem["bounds"])
        n_expected = len(next(iter(problem["groups"].values())))

        sub = df.filter(pl.col("problem") == pid)
        cells = build_cells(sub)

        vectors_by_family: dict[str, list[tuple[str, str, int, np.ndarray]]] = {}
        minimax_records: list[tuple[str, str, int, float]] = []

        for cellrow in cells.iter_rows(named=True):
            solver = cellrow["solver"]
            entry = registry.get(solver)
            if entry is None:
                cell_warnings.append((pid, solver, "solver not found in registry.json"))
                continue
            family_primary = family_primary_for(entry)
            method_class = entry["method_class"]

            w = load_real_weight(cellrow["weights_ref"])
            if w is None:
                continue  # no real weights for this cell -- not a failure, just excluded
            if w.shape[0] != n_expected:
                cell_warnings.append(
                    (pid, solver, f"weight length {w.shape[0]} != problem n {n_expected}")
                )
                continue

            row = metrics.compute_metrics(
                w, problem["groups"], problem["targets"], d=problem["design_weights"],
                bounds=bnds, family=family_primary,
            )
            nd = metrics.native_divergence(
                w, problem["design_weights"], family_primary, bnds,
                groups=problem["groups"] if family_primary == "minimax" else None,
                targets=problem["targets"] if family_primary == "minimax" else None,
            )
            row.update(nd)
            row.update(
                solver=solver, problem=pid, thread=cellrow["thread"], build=cellrow["build"],
                status=cellrow["status"], rep_count=cellrow["rep_count"],
                wall_time_s_median=cellrow["wall_time_s_median"],
                wall_time_s_min=cellrow["wall_time_s_min"], wall_time_s_max=cellrow["wall_time_s_max"],
                peak_rss_bytes=cellrow["peak_rss_bytes"], iterations=cellrow["iterations"],
                method_class=method_class, family=family_primary,
            )
            metrics_rows.append(row)

            if family_primary in STRICTLY_CONVEX:
                vectors_by_family.setdefault(family_primary, []).append((solver, cellrow["build"], cellrow["thread"], w))
            elif family_primary == "minimax":
                minimax_records.append((solver, cellrow["build"], cellrow["thread"], row["margin_linf"]))

        # RQ5 agreement, bounded per-problem: pairwise over this problem's
        # already-loaded real vectors, then dropped before the next problem.
        for family, lst in vectors_by_family.items():
            for (sa, ba, ta, wa), (sb, bb, tb, wb) in combinations(lst, 2):
                agr = metrics.agreement(wa, wb, family)
                agreement_rows.append(
                    dict(
                        solver_a=sa, build_a=ba, thread_a=ta, solver_b=sb, build_b=bb, thread_b=tb, problem=pid, family=family,
                        mode=agr["mode"], pearson=agr["pearson"], spearman=agr["spearman"],
                        max_abs_diff=agr["max_abs_diff"], cosine=agr["cosine"],
                        degenerate=bool(not (np.isfinite(agr["pearson"]) and np.isfinite(agr["spearman"]))),
                        obj_val_a=None, obj_val_b=None, rel_gap=None,
                    )
                )
        for (sa, ba, ta, oa), (sb, bb, tb, ob) in combinations(minimax_records, 2):
            agr = metrics.agreement(None, None, "minimax", obj_val=oa, obj_val_ref=ob)
            rel_gap = agr["abs_diff"] / max(abs(oa), abs(ob), 1e-12)
            agreement_rows.append(
                dict(
                    solver_a=sa, build_a=ba, thread_a=ta, solver_b=sb, build_b=bb, thread_b=tb, problem=pid, family="minimax",
                    mode=agr["mode"], pearson=None, spearman=None, max_abs_diff=None, cosine=None,
                    degenerate=False,
                    obj_val_a=oa, obj_val_b=ob, rel_gap=rel_gap,
                )
            )

        del vectors_by_family, minimax_records  # bounded per-problem streaming: free before next problem

    # infer_schema_length=None: full-scan dtype inference. Both row lists are
    # family-heterogeneous -- native_divergence() emits family-dependent metric
    # keys, and agreement rows null the non-applicable columns per family
    # (convex: obj_val_*/rel_gap None; minimax: pearson/spearman/max_abs_diff/
    # cosine None). Default 100-row inference locks a null dtype from all-None
    # leading rows, then fails to append a later f64 (ComputeError).
    metrics_df = pl.DataFrame(metrics_rows, infer_schema_length=None)
    agreement_df = pl.DataFrame(agreement_rows, infer_schema_length=None)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    metrics_df.write_parquet(METRICS_OUT)
    agreement_df.write_parquet(AGREEMENT_OUT)

    print(f"n_metric_rows={metrics_df.height}")
    print(f"n_agreement_pairs={agreement_df.height}")
    if metrics_df.height:
        print("rows per family:")
        print(metrics_df.group_by("family").len().sort("family"))
    if problem_failures:
        print(f"problem_failures ({len(problem_failures)}):")
        for pid, reason in problem_failures:
            print(f"  {pid}: {reason}")
    if cell_warnings:
        print(f"cell_warnings ({len(cell_warnings)}):")
        for pid, solver, reason in cell_warnings[:50]:
            print(f"  {pid}/{solver}: {reason}")


if __name__ == "__main__":
    main()
