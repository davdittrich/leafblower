#!/usr/bin/env python3
"""test_instance_family.py -- tests for benchmarks/study/common/instance_family.py
(WU-3, leafblower-2ouc.4).

Usage:
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \\
    python/.venv/bin/python benchmarks/study/common/test_instance_family.py [roundtrip_r_dump_path]

When roundtrip_r_dump_path is given, independently regenerates the same
small instances (by params, not by re-reading R's matrices) and diffs
against R's dump for exact R<->Python parity (DoD item 4).
"""

from __future__ import annotations

import json
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import math
import re
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import instance_family as ifam  # noqa: E402
import problem_io  # noqa: E402

failures = 0


def check(desc: str, cond: bool) -> None:
    global failures
    if cond:
        print(f"  PASS: {desc}")
    else:
        print(f"  FAIL: {desc}")
        failures += 1


def near(a: float, b: float, tol: float = 1e-9) -> bool:
    return math.isclose(a, b, abs_tol=tol, rel_tol=0.0)


print("== build_instance_grid() ==")
grid = ifam.build_instance_grid()
check("30 <= n_instances <= 100", 30 <= len(grid) <= 100)
print(f"  (n_instances = {len(grid)})")

n_vals = [g["n"] for g in grid]
K_vals = [g["K"] for g in grid]
card_vals = [g["card"] for g in grid]
skew_vals = [g["skew"] for g in grid]
cond_vals = [g["cond"] for g in grid]
infeas_vals = [g["infeas"] for g in grid]

check("n axis covers all 4 frozen levels", set(n_vals) == set(ifam.IF_N_LEVELS))
check("K axis covers all 3 frozen levels", set(K_vals) == set(ifam.IF_K_LEVELS))
check("card axis covers low/medium/high", set(card_vals) == set(ifam.IF_CARD_LEVELS))
check("skew axis covers none/moderate/extreme", set(skew_vals) == set(ifam.IF_SKEW_LEVELS))
check("cond axis covers well/moderate/ill", set(cond_vals) == set(ifam.IF_COND_LEVELS))
check("infeas axis covers loose/moderate/tight", set(infeas_vals) == set(ifam.IF_INFEAS_MAX))
check("only a handful (<=4) of instances touch n=1.58M", sum(1 for n in n_vals if n == 1580000) <= 4)
check("no duplicate instance ids", len({ifam._instance_id(g) for g in grid}) == len(grid))

print("== determinism: _generate_categories is repeatable ==")
mat1, _ = ifam._generate_categories(200, 4, 8, 1, 0.5, 0)
mat2, _ = ifam._generate_categories(200, 4, 8, 1, 0.5, 0)
check("identical matrices across repeated calls", mat1 == mat2)
mat3, _ = ifam._generate_categories(200, 4, 8, 1, 0.5, 1)
check("different seed => different matrix", mat1 != mat3)
check("category values within [1,C]", all(1 <= v <= 8 for row in mat1 for v in row))

print("== generate_instance_family_specs() writes schema-conformant specs ==")
with tempfile.TemporaryDirectory() as tmp_dir:
    res = ifam.generate_instance_family_specs(Path(tmp_dir) / "instance_family_spec_test")
    check("n written ids == grid length", len(res["ids"]) == len(grid))
    sample_ids = [res["ids"][0], res["ids"][len(res["ids"]) // 2], res["ids"][-1]]
    for id_ in sample_ids:
        with open(Path(res["out_dir"]) / f"{id_}.json") as f:
            spec = json.load(f)
        required = ["id", "data_ref", "margins", "bounds", "tol", "objective_families", "K"]
        check(f"{id_}: has all required fields", all(k in spec for k in required))
        check(f"{id_}: data_ref starts with gen:instance_family", spec["data_ref"].startswith("gen:instance_family?"))
        margins = spec["margins"]
        check(f"{id_}: length(margins) == K", len(margins) == spec["K"])
        check(f"{id_}: targets keys == margins", set(spec["targets"].keys()) == set(margins))
        for m in margins:
            check(f"{id_}: targets[{m}] sums to 1", near(sum(spec["targets"][m].values()), 1.0))
        check(f"{id_}: bounds.min == 0", near(spec["bounds"]["min"], 0.0))

    print("== gen: resolution through the WU-2 loader (load_problem_spec) ==")
    ifam.install_gen_resolver()
    small_ids = [gid for gid, n in zip(res["ids"], n_vals) if n <= 10000]
    for id_ in small_ids[:3]:
        path = Path(res["out_dir"]) / f"{id_}.json"
        problem = problem_io.load_problem_spec(path)
        with open(path) as f:
            spec_raw = json.load(f)
        expected_n = int(re.search(r"[?&]n=(\d+)", spec_raw["data_ref"]).group(1))
        check(f"{id_}: loaded n matches spec n query param", len(problem["data"]) == expected_n)
        check(f"{id_}: K margin columns present", all(m in problem["data"].columns for m in problem["margins"]))
        check(f"{id_}: design_weights default to all-ones", bool((problem["design_weights"] == 1.0).all()))

print("== unrecognized gen: recipe still falls through to WU-2's not-implemented guard ==")
gen_err = None
try:
    problem_io._resolve_data_ref("gen:toy_recipe", None)
except NotImplementedError as e:
    gen_err = str(e)
check(
    "gen:toy_recipe (unregistered) still raises WU-3-scope error after install_gen_resolver()",
    gen_err is not None and "WU-3" in gen_err,
)

if len(sys.argv) >= 2:
    r_dump_path = Path(sys.argv[1])
    print(f"\n== R<->Py exact parity check (vs {r_dump_path}) ==")
    with open(r_dump_path) as f:
        r_dumps = json.load(f)
    for spec in r_dumps:
        n, k, c, s, rho, seed = spec["n"], spec["K"], spec["C"], spec["s"], spec["rho"], spec["seed"]
        py_mat, py_zp = ifam._generate_categories(n, k, c, s, rho, seed)
        # R's jsonlite serializes the n x K matrix column-major (as a flat
        # array of arrays-of-columns via toJSON on a matrix: rows of the
        # JSON array correspond to matrix rows since toJSON on a base-R
        # matrix emits it as a list of row-vectors after row/col traversal
        # -- verified below by shape + value comparison rather than assumed).
        r_mat = spec["mat"]
        label = f"n={n},K={k},C={c},s={s},rho={rho},seed={seed}"
        check(f"{label}: matrix shape matches R", len(r_mat) == n and len(r_mat[0]) == k)
        check(f"{label}: matrix values EXACT match R", r_mat == py_mat)
        check(
            f"{label}: target Zipf probs EXACT match R",
            all(near(a, b, tol=0.0) or abs(a - b) < 1e-15 for a, b in zip(spec["probs"], py_zp[0])),
        )

print(f"\nRESULT: {failures} failure(s)")
if failures > 0:
    sys.exit(1)
