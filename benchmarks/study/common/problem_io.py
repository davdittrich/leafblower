"""problem_io.py — leafblower benchmark study problem-spec loader.

Resolves spec/*.json problem specs (data_ref origins file:/pkg:/gen:/inline)
into a standardized problem object. Design of record: docs/benchmark/DESIGN.md
§3; schema: benchmarks/study/spec/schema.json; ticket leafblower-2ouc.3 (WU-2).

All file paths in specs (data_ref 'file:', targets_ref) are relative to the
repo root, matching the existing benchmarks/*.py convention
(pl.read_parquet("benchmarks/...")) and CLAUDE.md's Build & Test invocation
(run with cwd = repo root).

Standardized problem object (dict), returned by load_problem_spec():
    id                  str
    data                pandas.DataFrame, len n; margin columns cast to str
    design_weights      np.ndarray[float64], shape (n,) -- d_i (defaults to
                        all-ones when omitted)
    margins             list[str], columns in `data`
    targets             dict[str, dict[str, float]] -- T_kj, each margin's
                        category values normalized to sum to 1
    bounds              dict(min=float, max=float); null sides resolved to
                        min=0.0 / max=float('inf')
    tol                 float
    objective_families  list[str]
    K                   int, == len(margins)
"""

from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd

_VALID_ORIGIN_RE = re.compile(r"^(file:.+|pkg:[^:]+::[^:]+|gen:.+|inline)$")

# Python-side registry for pkg: dataset resolvers, keyed by "<package>::<dataset>".
# Empty by default: the canonical R packages (survey, sampling, anesrake,
# ebal, optweight, balance) are R-only and their pkg: specs are single-arm
# home-turf problems by design (DESIGN.md §3) -- NOT resolving in Python is
# expected, not a bug. Python-native pkg: sources register here (WU-6 scope).
_PKG_LOADERS: dict[str, Callable[[], pd.DataFrame]] = {}


def load_problem_spec(spec_path: str | Path) -> dict[str, Any]:
    spec_path = Path(spec_path)
    if not spec_path.exists():
        raise FileNotFoundError(f"problem_io: spec file not found: {spec_path}")
    with open(spec_path) as f:
        spec = json.load(f)
    _validate_spec(spec)

    data = _resolve_data_ref(spec["data_ref"], spec.get("data"))
    n = len(data)

    design_weights = _resolve_design_weights(spec.get("design_weights"), data, n)

    targets_raw = _resolve_targets(spec.get("targets_ref"), spec.get("targets"))
    targets = {k: _normalize_target(v) for k, v in targets_raw.items()}

    margins = list(spec["margins"])
    missing = [m for m in margins if m not in data.columns]
    if missing:
        raise ValueError(
            f"problem_io: margin column(s) not found in resolved data: {missing}"
        )
    if set(targets.keys()) != set(margins):
        raise ValueError(
            "problem_io: targets_ref/targets keys must match `margins` exactly. "
            f"margins={margins} targets={list(targets.keys())}"
        )
    # String-categorical invariant (leafblower-2ouc.37): margins MUST be string-
    # coded categories. A numeric-coded margin coerces INCONSISTENTLY across arms --
    # R factor() levels render whole numbers as "1", Python astype(str) renders
    # float64 as "1.0" -- silently diverging from the JSON target keys and breaking
    # every ==-keyed site (structural_infeasible_cats, margin_stats, ref_convex).
    # load_problem_spec is the single chokepoint all consumers share, so reject here
    # (before the astype(str) coercion, which would mask the raw dtype). CONTENT
    # check (not dtype alone): an object column holding Python ints passes str() but
    # breaks ==-keyed sites. No-op on the current string-categorical study set.
    for m in margins:
        vals = data[m].dropna().unique()
        if not all(isinstance(v, str) for v in vals):
            raise ValueError(
                f"problem_io: margin column {m!r} must be string-categorical, got dtype "
                f"{data[m].dtype} with non-string values -- benchmark margins MUST be "
                f"string-coded categories (numeric margins diverge R vs Python; "
                f"leafblower-2ouc.37). Add explicit canonical string coercion before "
                f"introducing a numeric margin."
            )
    for m in margins:
        data[m] = data[m].astype(str)

    K = int(spec["K"])
    if K != len(margins):
        raise ValueError(f"problem_io: spec['K'] ({K}) != len(margins) ({len(margins)})")

    bounds_spec = spec["bounds"]
    bounds = dict(
        min=0.0 if bounds_spec.get("min") is None else float(bounds_spec["min"]),
        max=float("inf") if bounds_spec.get("max") is None else float(bounds_spec["max"]),
    )

    return dict(
        id=spec["id"],
        data=data,
        design_weights=design_weights,
        margins=margins,
        targets=targets,
        bounds=bounds,
        tol=float(spec["tol"]),
        objective_families=list(spec["objective_families"]),
        K=K,
    )


def _normalize_target(t: dict[str, float]) -> dict[str, float]:
    s = float(sum(t.values()))
    return {k: v / s for k, v in t.items()}


def _validate_spec(spec: dict[str, Any]) -> None:
    required = ["id", "data_ref", "margins", "bounds", "tol", "objective_families", "K"]
    missing = [k for k in required if k not in spec]
    if missing:
        raise ValueError(f"problem_io: spec missing required field(s): {missing}")
    data_ref = spec["data_ref"]
    if not _VALID_ORIGIN_RE.match(data_ref):
        raise ValueError(
            f"problem_io: invalid data_ref origin: '{data_ref}' "
            "(expected file:/pkg:/gen:/inline)"
        )
    if data_ref == "inline" and spec.get("data") is None:
        raise ValueError("problem_io: data_ref='inline' requires a `data` field in the spec")
    has_ref = spec.get("targets_ref") is not None
    has_inline = spec.get("targets") is not None
    if has_ref == has_inline:
        raise ValueError("problem_io: spec must set exactly one of targets_ref / targets")


def _resolve_data_ref(data_ref: str, inline_data: list[dict[str, Any]] | None) -> pd.DataFrame:
    if data_ref.startswith("file:"):
        path = data_ref[len("file:") :]
        if not os.path.exists(path):
            raise FileNotFoundError(f"problem_io: file not found: {path}")
        ext = Path(path).suffix.lower()
        if ext == ".parquet":
            df = pd.read_parquet(path)
        elif ext == ".csv":
            df = pd.read_csv(path)
        else:
            raise ValueError(f"problem_io: unsupported file: extension '{ext}' (path={path})")
        if "uuid" in df.columns:
            df = df.drop(columns=["uuid"])
        return df.reset_index(drop=True)

    if data_ref.startswith("pkg:"):
        pkg_spec = data_ref[len("pkg:") :]
        if pkg_spec not in _PKG_LOADERS:
            raise NotImplementedError(
                f"problem_io: no Python resolver registered for data_ref='{data_ref}'. "
                "R-only canonical datasets (survey/sampling/anesrake/ebal/optweight/"
                "balance) are single-arm home-turf specs by design (DESIGN.md §3); "
                "Python-native pkg: sources register in problem_io._PKG_LOADERS (WU-6 scope)."
            )
        return _PKG_LOADERS[pkg_spec]()

    if data_ref.startswith("gen:"):
        recipe_id = data_ref[len("gen:") :]
        raise NotImplementedError(
            "problem_io: gen: data_ref origin (parametric synthetic instance "
            f"generator) is WU-3 scope, not yet implemented (recipe_id='{recipe_id}')"
        )

    if data_ref == "inline":
        if inline_data is None:
            raise ValueError("problem_io: data_ref='inline' requires a `data` field in the spec")
        return pd.DataFrame(inline_data)

    raise ValueError(
        f"problem_io: unrecognised data_ref origin: '{data_ref}' "
        "(expected file:/pkg:/gen:/inline)"
    )


def _resolve_design_weights(design_weights: Any, data: pd.DataFrame, n: int) -> np.ndarray:
    if design_weights is None:
        return np.ones(n, dtype=np.float64)
    if isinstance(design_weights, str):
        if not design_weights.startswith("column:"):
            raise ValueError(
                f"problem_io: unrecognised design_weights string: '{design_weights}' "
                "(expected 'column:<name>')"
            )
        col = design_weights[len("column:") :]
        if col not in data.columns:
            raise ValueError(f"problem_io: design_weights column '{col}' not found in data")
        return data[col].to_numpy(dtype=np.float64)
    if isinstance(design_weights, list):
        arr = np.asarray(design_weights, dtype=np.float64)
        if arr.shape[0] != n:
            raise ValueError(
                f"problem_io: inline design_weights length {arr.shape[0]} != n {n}"
            )
        return arr
    raise ValueError(f"problem_io: unrecognised design_weights spec: {design_weights!r}")


def _resolve_targets(
    targets_ref: str | None, inline_targets: dict[str, dict[str, float]] | None
) -> dict[str, dict[str, float]]:
    has_ref = targets_ref is not None
    has_inline = inline_targets is not None
    if has_ref == has_inline:
        raise ValueError("problem_io: spec must set exactly one of targets_ref / targets")
    if has_ref:
        if not os.path.exists(targets_ref):
            raise FileNotFoundError(f"problem_io: targets_ref file not found: {targets_ref}")
        with open(targets_ref) as f:
            return json.load(f)
    return inline_targets
