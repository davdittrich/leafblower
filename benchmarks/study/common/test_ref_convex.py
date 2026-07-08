#!/usr/bin/env python3
"""benchmarks/study/common/test_ref_convex.py

Cross-check for benchmarks/study/common/ref_convex.py (ticket leafblower-2ouc.5,
WU-4) against an INDEPENDENT hand derivation on spec/toy_inline.json -- Python
twin of test_ref_convex.R, same golden literals.

Hand derivation (toy_inline.json: n=4, margin grp in {A,B}, design_weights=
[1,1,2,2], targets{A:2,B:2} normalized to {A:0.5,B:0.5}):
  N = sum(d) = 6; T_A = T_B = 0.5*6 = 3.
  Single margin, closed-form KL-raking scale per group = T_group / D_group:
    group A: D_A = d1+d2 = 2, scale = 3/2 = 1.5  -> w1=w2=1.5
    group B: D_B = d3+d4 = 4, scale = 3/4 = 0.75 -> w3=w4=2*0.75=1.5
  All four weights converge to exactly 1.5 (unique feasible point on this
  toy -- single 2-level margin pins the whole vector regardless of family).
  weight_kl = sum(w*log(w/d)) = 2*1.5*log(1.5) + 2*1.5*log(0.75) = 0.35334910696915045

Run: python/.venv/bin/python benchmarks/study/common/test_ref_convex.py
(single-thread BLAS forced before numpy import, matching ref_convex.py)
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import math  # noqa: E402
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from metrics import weight_kl  # noqa: E402
from problem_io import load_problem_spec  # noqa: E402
from ref_convex import REF_MAX_N, solve_ref, store_ref  # noqa: E402

fail_count = 0


def check(desc, got, want, tol=1e-8):
    ok = math.isclose(got, want, rel_tol=0, abs_tol=tol) or math.isclose(got, want, rel_tol=tol)
    print(f"[{'PASS' if ok else 'FAIL'}] {desc:<70} got={got!r} want={want!r}")
    global fail_count
    if not ok:
        fail_count += 1
    return ok


def check_true(desc, got):
    ok = bool(got)
    print(f"[{'PASS' if ok else 'FAIL'}] {desc:<70} got={got!r}")
    global fail_count
    if not ok:
        fail_count += 1
    return ok


repo_root = Path(__file__).resolve().parents[3]
spec_path = repo_root / "benchmarks" / "study" / "spec" / "toy_inline.json"
problem = load_problem_spec(spec_path)

W_HAND = np.array([1.5, 1.5, 1.5, 1.5])
WEIGHT_KL_HAND = 0.35334910696915045

print("== kl family (WU-1 golden toy cross-check) ==")
res_kl = solve_ref(problem, "kl")
check_true("kl mode == weight_vector", res_kl["mode"] == "weight_vector")
check("kl weights[0]", res_kl["weights"][0], 1.5)
check("kl weights[1]", res_kl["weights"][1], 1.5)
check("kl weights[2]", res_kl["weights"][2], 1.5)
check("kl weights[3]", res_kl["weights"][3], 1.5)
check_true("kl weights match hand-derived vector (max abs diff < 1e-8)",
           float(np.max(np.abs(res_kl["weights"] - W_HAND))) < 1e-8)
wkl = weight_kl(res_kl["weights"], problem["design_weights"], family="kl")["weight_kl"]
check("weight_kl(ref kl weights) matches hand-derived 0.3533489267", wkl, WEIGHT_KL_HAND)

print("\n== chi2 family (same toy: unique feasible point == same 1.5 vector) ==")
res_chi2 = solve_ref(problem, "chi2")
check_true("chi2 mode == weight_vector", res_chi2["mode"] == "weight_vector")
check_true("chi2 weights match hand-derived vector",
           float(np.max(np.abs(res_chi2["weights"] - W_HAND))) < 1e-8)

print("\n== logit family (same toy: bounds [0,10] finite, unique feasible point) ==")
res_logit = solve_ref(problem, "logit")
check_true("logit mode == weight_vector", res_logit["mode"] == "weight_vector")
check_true("logit weights match hand-derived vector",
           float(np.max(np.abs(res_logit["weights"] - W_HAND))) < 1e-8)

print("\n== minimax family (Blocker G: objective_value anchor, not weight-vector) ==")
res_mm = solve_ref(problem, "minimax")
check_true("minimax mode == objective_value (never weight_vector)",
           res_mm["mode"] == "objective_value")
check("minimax achieved margin_linf ~ 0 (single margin fully saturable)", res_mm["obj_val"], 0.0, tol=1e-8)

print("\n== REF_MAX_N scope guard (stepstone 1.58M has no anchor, DESIGN.md Sec.6) ==")
fake_big = dict(id="fake_stepstone", data=pd.DataFrame({"x": range(REF_MAX_N + 1)}))
guard_err = None
try:
    solve_ref(fake_big, "kl")
except ValueError as e:
    guard_err = str(e)
check_true("solve_ref refuses n > REF_MAX_N with an explicit message (not a silent hang/fake)",
           guard_err is not None and "REF_MAX_N" in guard_err)

print("\n== store_ref: pseudo-solver row persistence ==")
out_dir = str(repo_root / "benchmarks" / "study" / "results")
p_kl = store_ref(problem, "kl", res_kl, out_dir)
p_chi2 = store_ref(problem, "chi2", res_chi2, out_dir)
p_logit = store_ref(problem, "logit", res_logit, out_dir)
p_mm = store_ref(problem, "minimax", res_mm, out_dir)
check_true("kl weights parquet written", p_kl.exists())
check_true("chi2 weights parquet written", p_chi2.exists())
check_true("logit weights parquet written", p_logit.exists())
check_true("minimax objective-value JSON written (NOT a weights parquet)", p_mm.exists())
check_true("minimax anchor path is under ref_objective/, not weights/", "ref_objective" in str(p_mm))

reread = pd.read_parquet(p_kl)
check_true("re-read kl parquet round-trips the weight vector",
           float(np.max(np.abs(reread["weight"].to_numpy() - res_kl["weights"]))) < 1e-12)

print(f"\n{'GOLDEN PASS' if fail_count == 0 else 'GOLDEN FAIL'}: {fail_count} assertion(s) failed.")
if __name__ == "__main__":
    sys.exit(0 if fail_count == 0 else 1)
