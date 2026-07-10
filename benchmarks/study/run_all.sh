#!/usr/bin/env bash
# WU-11 (leafblower-2ouc.12) SCORED CORE run orchestrator. Non-frozen execution
# artifact (like --reps/--threads): lives OUTSIDE benchmarks/study/spec, disclosed
# via git + environment.json. NOT the rehearsal launcher (rehearsal/launch_v8.sh,
# unscored reps=1) -- this is the article's canonical data-producing run.
#
# Protocol (DESIGN.md §4 CORE / §5):
#   - 1 thread (CORE); single-thread BLAS (OMP/OPENBLAS/MKL=1) for determinism.
#   - reps=10 on nonheavy (canonical + small/medium instance family);
#     reps=5 on heavy (n=1.58M: if_n1580000_* + stepstone_{unbounded,bounded}).
#   - warmups=2 discarded; min-total-duration 0.05s (sub-ms cells inner-batched);
#     cell-timeout 3600s (run_config.json cell_time_budget_s -> status='dnf').
#   - CPU pinning via --pin-core (taskset); randomized cell order via --seed
#     (rm_build_matrix rng_seed). Governor is NOT root-settable on this box
#     (powersave, recorded in environment.json as a disclosed limitation);
#     robust stats (median + 5/95 + bootstrap CI, WU-12) absorb governor jitter.
#   - WU-9 pre-run gate: --assert-runnable-tag hard-stops on a dirty/drifted tree.
#
# Memory gate (run_config.json memory_gate; 60 GB box):
#   - Nonheavy: R || Python concurrent (small RSS, pinned to distinct cores).
#   - Heavy (n=1.58M): STRICTLY serial, one arm at a time, one cell at a time
#     (the driver is serial within an arm; we never co-run the two arms). Max one
#     heavy worker alive => OOM-killer (if any) targets that worker only.
#   - sbw @ stepstone_{unbounded,bounded}: BEYOND CAPACITY (dense ATT phantom =
#     1.58M x 836 real-name dummies ~= 10 TB dense -> OOM >56 GB every attempt,
#     deterministic; run_config beyond_capacity_cells). PRE-CLASSIFIED as a
#     status='error' (exit 137 / OOM) row WITHOUT running, to avoid thrashing the
#     box. sbw @ if_n1580000_* (fewer categories, measured ~60 GB peak) IS run,
#     serial-alone -- it either completes or the driver records an honest crash row.
set -uo pipefail
cd /home/dd/Gemini/leafblower || exit 2

TAG=benchmark-runnable-freeze-v11
OUT=benchmarks/study/results
SHARDS="$OUT/_shards"
PY="python/.venv/bin/python benchmarks/study/python/run_arm.py"
R="Rscript benchmarks/study/R/run_arm.R"
SEED=20260710
# Uniform scored knobs. reps set per-tier below.
COMMON="--threads 1 --warmups 2 --min-total-duration 0.05 --cell-timeout 3600 \
        --seed $SEED --assert-runnable-tag $TAG"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

HEAVY='if_n1580000_|stepstone'
# Nonheavy = everything not matching HEAVY. run_arm has no --exclude-problems, so
# nonheavy is expressed as a positive alternation of the non-1.58M problem-id
# prefixes actually present in spec/ (canonical + toy + if_n{1000,10000,100000}).
NONHEAVY='canonical_|toy_inline|if_n1000_|if_n10000_|if_n100000_'

# Non-sbw solver alternation (^(a|b|...)$) for the heavy stage: run_arm --solvers
# is an inclusive regex, so we list every registry solver EXCEPT sbw. Cells whose
# arm lacks a given id simply don't match -- harmless. Anchored to avoid partial
# (substring) matches (e.g. 'oris' vs 'oris_soft').
NONSBW=$(python/.venv/bin/python -c "import json;ks=[k for k in json.load(open('benchmarks/study/registry.json'))['solvers'] if k!='sbw'];print('^('+'|'.join(ks)+')\$')")
# HARD guard: an empty NONSBW (python/registry failure) is NOT caught by set -u
# (var is set, just empty); --solvers "" => grepl("",solver)=TRUE for ALL solvers
# incl sbw => sbw@stepstone runs => >56 GB OOM (the exact catastrophe we pre-skip).
[ -n "$NONSBW" ] && [[ "$NONSBW" == '^('*')$' ]] || { echo "[run_all] FATAL: NONSBW alternation empty/malformed ('$NONSBW')"; exit 3; }
case "|$NONSBW|" in *"|sbw|"*|*"(sbw|"*|*"|sbw)"*) echo "[run_all] FATAL: sbw leaked into NONSBW"; exit 3;; esac

rm -rf "$SHARDS"; mkdir -p "$SHARDS"
# Clean scored weight store (gitignored working data); keep results/ref_objective
# (WU-4 convex reference, produced separately -- NOT a run output).
rm -f "$OUT"/runs.parquet "$OUT"/weights/*.parquet 2>/dev/null || true
mkdir -p "$OUT/weights"

echo "[run_all] START $(date -Iseconds)  seed=$SEED  tag=$TAG"

# ---------------------------------------------------------------- Stage A: NONHEAVY
# R || Python concurrent, reps=10, pinned to distinct cores.
echo "[run_all] Stage A NONHEAVY (reps=10, R||Py)  $(date -Iseconds)"
$PY $COMMON --reps 10 --pin-core 3 --problems "$NONHEAVY" --out "$SHARDS/py_nonheavy" \
    > "$SHARDS/py_nonheavy.log" 2>&1 &
PY_PID=$!
$R  $COMMON --reps 10 --pin-core 2 --problems "$NONHEAVY" --out "$SHARDS/r_nonheavy" \
    > "$SHARDS/r_nonheavy.log" 2>&1 &
R_PID=$!
wait $PY_PID; PY_A=$?
wait $R_PID;  R_A=$?
echo "[run_all] Stage A done  py_exit=$PY_A r_exit=$R_A  $(date -Iseconds)"

# ---------------------------------------------------------------- Stage B: HEAVY
# Strictly serial (one worker alive at a time). reps=5. sbw excluded (handled below).
# Order: Python arm fully, then R arm.  Pin all heavy to one core.
echo "[run_all] Stage B HEAVY non-sbw (reps=5, serial)  $(date -Iseconds)"
$PY $COMMON --reps 5 --pin-core 2 --problems "$HEAVY" --solvers "$NONSBW" \
    --out "$SHARDS/py_heavy" > "$SHARDS/py_heavy.log" 2>&1
echo "[run_all]   py heavy done $(date -Iseconds)"
$R  $COMMON --reps 5 --pin-core 2 --problems "$HEAVY" --solvers "$NONSBW" \
    --out "$SHARDS/r_heavy" > "$SHARDS/r_heavy.log" 2>&1
echo "[run_all]   r heavy done $(date -Iseconds)"

# ---- sbw @ if_n1580000_* : measured-heavy, RUN serial-alone (R arm only) ----
echo "[run_all] Stage B sbw @ if_n1580000 (serial-alone)  $(date -Iseconds)"
$R  $COMMON --reps 5 --pin-core 2 --problems 'if_n1580000_' --solvers '^sbw$' \
    --out "$SHARDS/r_sbw_if" > "$SHARDS/r_sbw_if.log" 2>&1
echo "[run_all] Stage B done  $(date -Iseconds)"

# ---------------------------------------------------------------- Merge + sbw@stepstone synth
echo "[run_all] merge + sbw@stepstone beyond-capacity rows  $(date -Iseconds)"
python/.venv/bin/python - "$SHARDS" "$OUT" <<'PY'
import sys, glob, os
import pandas as pd, numpy as np
shards, out = sys.argv[1], sys.argv[2]
paths = sorted(glob.glob(f"{shards}/*/runs.parquet"))
frames = []
found = set()
for p in paths:
    d = pd.read_parquet(p); sh = p.split("/")[-2]; d["shard"] = sh; frames.append(d); found.add(sh)
    print(f"  shard {sh:14s} rows={len(d)}")
if not frames:
    print("FATAL: no shard runs.parquet found"); sys.exit(1)
# A single-arm driver-FATAL crash (venv/import death before writing runs.parquet)
# would silently drop that whole arm's rows and still exit 0. Assert every expected
# shard is present so a lost arm is a hard failure, not an inflated "success".
expected = {"py_nonheavy", "r_nonheavy", "py_heavy", "r_heavy", "r_sbw_if"}
missing = expected - found
if missing:
    print(f"FATAL: expected shard(s) missing (driver arm crashed?): {sorted(missing)}")
    sys.exit(1)
m = pd.concat(frames, ignore_index=True)
cols = [c for c in m.columns if c != "shard"]

# sbw @ stepstone_{unbounded,bounded}: beyond-capacity (dense phantom OOM >56 GB;
# run_config beyond_capacity_cells). Pre-classified error row, NOT run. R arm only
# (sbw is an R competitor). Write the NaN-sentinel weight parquet so weights_ref
# resolves, mirroring run_arm.R::.crash_row.
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
    synth.append(row)
    print(f"  synth  sbw@{prob} -> status=error (beyond_capacity)")
m = pd.concat([m[cols], pd.DataFrame(synth)[cols]], ignore_index=True)

runs_path = os.path.join(out, "runs.parquet")
m.to_parquet(runs_path, index=False)
cells = m[["solver", "problem", "thread", "build"]].drop_duplicates().shape[0]
print(f"MERGED rows={len(m)} cells={cells} -> {runs_path}")
print("status:", m["status"].value_counts().to_dict())
nweights = len(glob.glob(os.path.join(wdir, "*.parquet")))
print(f"weight store: {nweights} files in {wdir}")
PY
echo "[run_all] DONE $(date -Iseconds)"
