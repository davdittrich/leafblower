#!/usr/bin/env bash
# WU-11 (leafblower-2ouc.12) SCORED CORE run orchestrator. Non-frozen execution
# artifact (like --reps/--jobs): lives OUTSIDE benchmarks/study/spec, disclosed
# via git + environment.json. NOT the rehearsal launcher (rehearsal/launch_v8.sh,
# unscored reps=1) -- this is the article's canonical data-producing run.
#
# ARCHITECTURE (v2, 2026-07-10): the driver now executes cells IN PARALLEL
# (--jobs N; parallel::mclapply / ThreadPoolExecutor) with per-cell checkpoint
# (result_<key>.json) + resume-skip. This replaces the v1 serial-per-arm design
# whose 3600s DNF ceiling blocked the queue an hour per pathological competitor
# (scipy_trust_constr/cvxpy_linf DNF at 1h on ~every cell -> 40-60h projected).
# Now DNF cells run CONCURRENTLY across cores and a killed run RESUMES.
#
# Protocol (DESIGN.md §4 CORE / §5): 1 thread; reps=10 nonheavy / 5 heavy(n=1.58M);
# warmups=2 discarded; min-total-duration 0.05s; cell-timeout 3600s (run_config
# cell_time_budget_s -> status='dnf', right-censored). Randomized cell order via
# --seed. WU-9 pre-run gate: --assert-runnable-tag hard-stops a dirty/drifted tree.
# Governor NOT root-settable here (powersave) -> driver WARN-not-fail, recorded in
# environment.json (disclosed limitation; WU-12 robust stats absorb the jitter).
#
# Memory gate (run_config.json memory_gate; 60 GB box):
#   - Nonheavy (n<=100k, small RSS): --jobs 20 (= default_concurrency). Arms run
#     SEQUENTIALLY (R then Py) so total concurrent workers stay <=20.
#   - Heavy non-sbw (n=1.58M): --jobs 3 (conservative). sbw's dense phantom is the
#     ONLY measured OOM; the other heavy solvers are sparse conic/iterative, but
#     run_config flags optweight_linf/cvxpy_linf RSS "unmeasured, do NOT assume
#     light" -> cap at 3 concurrent (<=~45 GB even at 15 GB/worker). Arms sequential.
#   - sbw @ if_n1580000_* (measured ~60 GB dense): --jobs 1, serial-ALONE (R arm).
#   - sbw @ stepstone_{unbounded,bounded}: BEYOND CAPACITY (~10 TB dense -> OOM
#     >56 GB deterministic; run_config beyond_capacity_cells). PRE-CLASSIFIED as a
#     status='error' row in the merge WITHOUT running (never spawned anywhere).
set -uo pipefail
cd /home/dd/Gemini/leafblower || exit 2

TAG=benchmark-runnable-freeze-v12
OUT=benchmarks/study/results
SHARDS="$OUT/_shards"
PY="python/.venv/bin/python benchmarks/study/python/run_arm.py"
R="Rscript benchmarks/study/R/run_arm.R"
SEED=20260710
COMMON="--threads 1 --warmups 2 --min-total-duration 0.05 --cell-timeout 3600 \
        --seed $SEED --assert-runnable-tag $TAG"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

HEAVY='if_n1580000_|stepstone'
NONHEAVY='canonical_|toy_inline|if_n1000_|if_n10000_|if_n100000_'

# Heavy stage sub-split (opus review): optweight_linf + cvxpy_linf are
# run_config suspected_heavy_unmeasured -- same dense QP/KKT class as sbw (60 GB),
# peak RSS NOT yet measured. They MUST run at concurrency 1 on n=1.58M (never
# co-scheduled), else OOM-kills the box OR OOM-kills a would-succeed solver and
# writes it as a false status='error' into the SCORED dataset. NONSBW_LIGHT =
# every solver EXCEPT {sbw, optweight_linf, cvxpy_linf} -> safe at jobs>1; the two
# QP suspects run in their own jobs=1 stage (HEAVY_QP).
NONSBW_LIGHT=$(python/.venv/bin/python -c "import json;ex={'sbw','optweight_linf','cvxpy_linf'};ks=[k for k in json.load(open('benchmarks/study/registry.json'))['solvers'] if k not in ex];print('^('+'|'.join(ks)+')\$')")
[ -n "$NONSBW_LIGHT" ] && [[ "$NONSBW_LIGHT" == '^('*')$' ]] || { echo "[run_all] FATAL: NONSBW_LIGHT empty/malformed ('$NONSBW_LIGHT')"; exit 3; }
case "|$NONSBW_LIGHT|" in *"|sbw|"*|*"|optweight_linf|"*|*"|cvxpy_linf|"*|*"(sbw|"*|*"(optweight_linf|"*|*"(cvxpy_linf|"*) echo "[run_all] FATAL: a heavy-QP/sbw solver leaked into NONSBW_LIGHT"; exit 3;; esac
HEAVY_QP='^(optweight_linf|cvxpy_linf)$'

mkdir -p "$SHARDS" "$OUT/weights"
# NOTE: no wipe of results/weights or _shards -- the per-cell result_<key>.json
# checkpoints make this run RESUMABLE (a re-invocation skips completed cells). To
# force a clean run, `rm -rf benchmarks/study/results/{_shards,weights} runs.parquet`
# manually first. results/ref_objective (WU-4 convex reference) is never touched.

echo "[run_all] START $(date -Iseconds)  seed=$SEED  tag=$TAG"

# ---------------------------------------------------------------- Stage A: NONHEAVY
echo "[run_all] Stage A NONHEAVY R (reps=10, jobs=20)  $(date -Iseconds)"
$R  $COMMON --reps 10 --jobs 20 --problems "$NONHEAVY" --out "$SHARDS/r_nonheavy" \
    > "$SHARDS/r_nonheavy.log" 2>&1
echo "[run_all]   r_nonheavy exit=$?  $(date -Iseconds)"
echo "[run_all] Stage A NONHEAVY Py (reps=10, jobs=20)  $(date -Iseconds)"
$PY $COMMON --reps 10 --jobs 20 --problems "$NONHEAVY" --out "$SHARDS/py_nonheavy" \
    > "$SHARDS/py_nonheavy.log" 2>&1
echo "[run_all]   py_nonheavy exit=$?  $(date -Iseconds)"

# ---------------------------------------------------------------- Stage B: HEAVY light (non-sbw, non-QP)
echo "[run_all] Stage B HEAVY light Py (reps=5, jobs=3)  $(date -Iseconds)"
$PY $COMMON --reps 5 --jobs 3 --problems "$HEAVY" --solvers "$NONSBW_LIGHT" \
    --out "$SHARDS/py_heavy" > "$SHARDS/py_heavy.log" 2>&1
echo "[run_all]   py_heavy exit=$?  $(date -Iseconds)"
echo "[run_all] Stage B HEAVY light R (reps=5, jobs=3)  $(date -Iseconds)"
$R  $COMMON --reps 5 --jobs 3 --problems "$HEAVY" --solvers "$NONSBW_LIGHT" \
    --out "$SHARDS/r_heavy" > "$SHARDS/r_heavy.log" 2>&1
echo "[run_all]   r_heavy exit=$?  $(date -Iseconds)"

# ---- Heavy QP-suspects (optweight_linf R / cvxpy_linf Py): jobs=1, serial (unmeasured
#      RSS ~ sbw dense class; run_config suspected_heavy_unmeasured). Never co-scheduled. ----
echo "[run_all] Stage C HEAVY QP-suspect Py cvxpy_linf (reps=5, jobs=1)  $(date -Iseconds)"
$PY $COMMON --reps 5 --jobs 1 --problems "$HEAVY" --solvers "$HEAVY_QP" \
    --out "$SHARDS/py_heavy_qp" > "$SHARDS/py_heavy_qp.log" 2>&1
echo "[run_all]   py_heavy_qp exit=$?  $(date -Iseconds)"
echo "[run_all] Stage C HEAVY QP-suspect R optweight_linf (reps=5, jobs=1)  $(date -Iseconds)"
$R  $COMMON --reps 5 --jobs 1 --problems "$HEAVY" --solvers "$HEAVY_QP" \
    --out "$SHARDS/r_heavy_qp" > "$SHARDS/r_heavy_qp.log" 2>&1
echo "[run_all]   r_heavy_qp exit=$?  $(date -Iseconds)"

# ---- sbw @ if_n1580000_* : measured-heavy dense, serial-ALONE (R arm only) ----
echo "[run_all] Stage B sbw @ if_n1580000 (reps=5, jobs=1, serial-alone)  $(date -Iseconds)"
$R  $COMMON --reps 5 --jobs 1 --problems 'if_n1580000_' --solvers '^sbw$' \
    --out "$SHARDS/r_sbw_if" > "$SHARDS/r_sbw_if.log" 2>&1
echo "[run_all]   r_sbw_if exit=$?  $(date -Iseconds)"

# ---------------------------------------------------------------- Merge + sbw@stepstone synth
echo "[run_all] merge + sbw@stepstone beyond-capacity rows  $(date -Iseconds)"
python/.venv/bin/python - "$SHARDS" "$OUT" <<'PY'
import sys, glob, os
import pandas as pd, numpy as np
shards, out = sys.argv[1], sys.argv[2]
paths = sorted(glob.glob(f"{shards}/*/runs.parquet"))
frames, found = [], set()
for p in paths:
    d = pd.read_parquet(p); sh = p.split("/")[-2]; d["shard"] = sh; frames.append(d); found.add(sh)
    print(f"  shard {sh:14s} rows={len(d)}")
if not frames:
    print("FATAL: no shard runs.parquet found"); sys.exit(1)
# A single-arm driver-FATAL crash (venv/import death before writing runs.parquet)
# would silently drop that arm's rows and still exit 0. Require every expected shard.
expected = {"py_nonheavy", "r_nonheavy", "py_heavy", "r_heavy",
            "py_heavy_qp", "r_heavy_qp", "r_sbw_if"}
missing = expected - found
if missing:
    print(f"FATAL: expected shard(s) missing (driver arm crashed?): {sorted(missing)}"); sys.exit(1)
m = pd.concat(frames, ignore_index=True)
cols = [c for c in m.columns if c != "shard"]

# sbw @ stepstone_{unbounded,bounded}: beyond-capacity dense phantom OOM (>56 GB;
# run_config beyond_capacity_cells). Pre-classified error row, NOT run. R arm only.
# NaN-sentinel weight parquet so weights_ref resolves (mirrors run_arm.R::.crash_row).
wdir = os.path.join(out, "weights"); os.makedirs(wdir, exist_ok=True)
synth = []
for prob in ("stepstone_unbounded", "stepstone_bounded"):
    ref = f"benchmarks/study/results/weights/sbw__{prob}__t1__na.parquet"
    pd.DataFrame({"weight": [np.nan]}).to_parquet(ref, index=False)
    row = {c: None for c in cols}
    row.update(dict(
        solver="sbw", problem=prob, thread=1, build="na", rep=0, weights_ref=ref,
        iterations=None, status="error", converged=False,
        error_message=("beyond memory capacity: sbw dense ATT phantom (1.58M x ~836 "
                       "real-name dummies ~10 TB) OOM-kills >56 GB on 60 GB box every "
                       "attempt (exit 137, deterministic); pre-classified without "
                       "running per run_config.json beyond_capacity_cells"),
        wall_time_s=float("nan"), peak_rss_bytes=float("nan"), trajectory_ref=None,
    ))
    synth.append(row); print(f"  synth  sbw@{prob} -> status=error (beyond_capacity)")
m = pd.concat([m[cols], pd.DataFrame(synth)[cols]], ignore_index=True)

runs_path = os.path.join(out, "runs.parquet")
m.to_parquet(runs_path, index=False)
cells = m[["solver", "problem", "thread", "build"]].drop_duplicates().shape[0]
print(f"MERGED rows={len(m)} cells={cells} -> {runs_path}")
print("status:", m["status"].value_counts().to_dict())
print(f"weight store: {len(glob.glob(os.path.join(wdir,'*.parquet')))} files")
PY
echo "[run_all] DONE $(date -Iseconds)"
