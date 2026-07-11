#!/usr/bin/env bash
# STUDY-BRANCH-ONLY-DO-NOT-MERGE
# WU-12c (leafblower-2ouc.42) bounded trajectory capture for RQ3 convergence
# curves. leafblower solvers x representative problems x {R,Py}. Cold path
# (byte-identical weights, proven). Outputs CSV under benchmarks/study/results/.
set -euo pipefail
cd "$(dirname "$0")/../../.."          # repo root
STUDY=benchmarks/study
OUTDIR=$STUDY/results/trajectories
mkdir -p "$OUTDIR"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LBW_TRAJECTORY_AT="1,2,3,5,10,20,50,100,200,500,1000,2000"

METHODS=(raking sinkhorn greenkhorn newton_kl logit greg chebyshev oris oris_soft)
PROBLEMS=(stepstone_bounded stepstone_unbounded \
          if_n10000_K4_medium_moderate_well_loose_s0 \
          if_n10000_K9_high_extreme_ill_tight_s0 \
          if_n1000_K2_high_extreme_ill_tight_s0)
PY=python/.venv/bin/python

resolve_spec() {   # echo the spec path for a problem id
  local pid="$1"
  if [ -f "$STUDY/spec/$pid.json" ]; then echo "$STUDY/spec/$pid.json"
  elif [ -f "$STUDY/spec/instance_family/$pid.json" ]; then echo "$STUDY/spec/instance_family/$pid.json"
  else echo ""; fi
}

for pid in "${PROBLEMS[@]}"; do
  spec=$(resolve_spec "$pid")
  if [ -z "$spec" ]; then echo "SKIP $pid (no spec)"; continue; fi
  for m in "${METHODS[@]}"; do
    for lang in r py; do
      out="$OUTDIR/${m}__${pid}__${lang}.csv"
      export LBW_TRAJECTORY_OUT="$out"
      if [ "$lang" = "r" ]; then
        Rscript "$STUDY/analysis/capture_one.R" "$m" "$spec" >/dev/null 2>&1 \
          && echo "OK  $out ($(wc -l <"$out" 2>/dev/null || echo 0) ln)" \
          || echo "ERR $m/$pid/r"
      else
        "$PY" "$STUDY/analysis/capture_one.py" "$m" "$spec" >/dev/null 2>&1 \
          && echo "OK  $out ($(wc -l <"$out" 2>/dev/null || echo 0) ln)" \
          || echo "ERR $m/$pid/py"
      fi
    done
  done
done
echo "trajectory capture complete -> $OUTDIR"
