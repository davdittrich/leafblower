#!/usr/bin/env python3
"""test_problem_io.py -- loader smoke test for benchmarks/study/common/problem_io.py

Exercises the data_ref origins reachable from Python:
    inline -> spec/toy_inline.json
    file:  -> spec/stepstone_unbounded.json, spec/stepstone_bounded.json
    gen:   -> asserted to raise the documented WU-3-not-implemented error
    pkg:   -> asserted to raise the documented "no Python resolver" error
              (spec/canonical_survey_apistrat.json is survey::apistrat, an
              R-only package -- single-arm home-turf by design, DESIGN.md §3)

Usage: python/.venv/bin/python benchmarks/study/common/test_problem_io.py [roundtrip_r_summary_path]
Set OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1 before `import numpy`
(CLAUDE.md determinism rule) -- this script does not import numpy until after
those are read from the environment by the interpreter's C extensions, but the
caller must still export them for parity with the R run.

When roundtrip_r_summary_path is given, cross-checks this run's independently
computed values for the cross-language-portable specs (toy_inline,
stepstone_unbounded) against the R loader's summary -- the R<->Py round-trip
identical-problem check (WU-2 DoD).
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from problem_io import _resolve_data_ref, load_problem_spec  # noqa: E402

SPEC_DIR = Path("benchmarks/study/spec")

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


print("== toy_inline (data_ref='inline') ==")
toy = load_problem_spec(SPEC_DIR / "toy_inline.json")
check("id", toy["id"] == "toy_inline")
check("n rows", len(toy["data"]) == 4)
check("K", toy["K"] == 1)
check("margins", toy["margins"] == ["grp"])
check("targets sum to 1", near(sum(toy["targets"]["grp"].values()), 1.0))
check(
    "targets values",
    near(toy["targets"]["grp"]["A"], 0.5) and near(toy["targets"]["grp"]["B"], 0.5),
)
check("design_weights (inline array)", list(toy["design_weights"]) == [1.0, 1.0, 2.0, 2.0])
check("bounds min", near(toy["bounds"]["min"], 0.0))
check("bounds max", near(toy["bounds"]["max"], 10.0))

print("== canonical_survey_apistrat (data_ref='pkg:survey::apistrat') ==")
pkg_err = None
try:
    load_problem_spec(SPEC_DIR / "canonical_survey_apistrat.json")
except NotImplementedError as e:
    pkg_err = str(e)
check(
    "pkg: raises 'no Python resolver' (R-only package, single-arm by design)",
    pkg_err is not None and "no Python resolver registered" in pkg_err,
)

print("== stepstone_unbounded (data_ref='file:...parquet') ==")
step_u = load_problem_spec(SPEC_DIR / "stepstone_unbounded.json")
check("id", step_u["id"] == "stepstone_unbounded")
check("n rows", len(step_u["data"]) == 1_582_732)
check("K", step_u["K"] == 9)
check("margins count", len(step_u["margins"]) == 9)
check("design_weights all ones", bool((step_u["design_weights"] == 1.0).all()))
check("design_weights length", len(step_u["design_weights"]) == len(step_u["data"]))
check("bounds max unbounded", math.isinf(step_u["bounds"]["max"]))
for nm in step_u["margins"]:
    check(f"targets '{nm}' sum to 1", near(sum(step_u["targets"][nm].values()), 1.0))

print("== stepstone_bounded (data_ref='file:...parquet', bounds.max=5) ==")
step_b = load_problem_spec(SPEC_DIR / "stepstone_bounded.json")
check("id", step_b["id"] == "stepstone_bounded")
check("n rows", len(step_b["data"]) == 1_582_732)
check("bounds max = 5", near(step_b["bounds"]["max"], 5.0))

print("== gen: origin (WU-3 not-yet-implemented guard) ==")
gen_err = None
try:
    _resolve_data_ref("gen:toy_recipe", None)
except NotImplementedError as e:
    gen_err = str(e)
check("gen: raises WU-3-scope error", gen_err is not None and "WU-3" in gen_err)

if len(sys.argv) >= 2:
    r_summary_path = Path(sys.argv[1])
    print(f"\n== R<->Py round-trip check (vs {r_summary_path}) ==")
    with open(r_summary_path) as f:
        r_summary = json.load(f)

    r_toy = r_summary["toy_inline"]
    check("toy_inline: id matches R", r_toy["id"] == toy["id"])
    check("toy_inline: n matches R", r_toy["n"] == len(toy["data"]))
    check("toy_inline: K matches R", r_toy["K"] == toy["K"])
    check("toy_inline: margins match R", sorted(r_toy["margins"]) == sorted(toy["margins"]))
    check(
        "toy_inline: targets match R (rtol=1e-6)",
        all(
            near(r_toy["targets"]["grp"][cat], toy["targets"]["grp"][cat], tol=1e-6)
            for cat in toy["targets"]["grp"]
        ),
    )
    check(
        "toy_inline: design_weights sum matches R",
        near(r_toy["design_weights_sum"], float(toy["design_weights"].sum()), tol=1e-6),
    )

    r_step = r_summary["stepstone_unbounded"]
    check("stepstone_unbounded: id matches R", r_step["id"] == step_u["id"])
    check("stepstone_unbounded: n (incidence dims) matches R", r_step["n"] == len(step_u["data"]))
    check("stepstone_unbounded: K matches R", r_step["K"] == step_u["K"])
    check(
        "stepstone_unbounded: margins match R",
        sorted(r_step["margins"]) == sorted(step_u["margins"]),
    )
    check(
        "stepstone_unbounded: design_weights (d_i) sum matches R",
        near(r_step["design_weights_sum"], float(step_u["design_weights"].sum()), tol=1e-6),
    )
    targets_match = True
    for nm in step_u["margins"]:
        r_t = r_step["targets"][nm]
        py_t = step_u["targets"][nm]
        if set(r_t.keys()) != set(py_t.keys()):
            targets_match = False
            break
        for cat in py_t:
            if not near(r_t[cat], py_t[cat], tol=1e-6):
                targets_match = False
                break
    check("stepstone_unbounded: targets (T_kj) match R (rtol=1e-6)", targets_match)

print(f"\nRESULT: {failures} failure(s)")
if failures > 0:
    sys.exit(1)
