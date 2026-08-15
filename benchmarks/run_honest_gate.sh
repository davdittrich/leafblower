#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

# Single-thread BLAS envelope (CLAUDE.md Build & Test) plus the opt-in gate
# switch — exported once here, not repeated per step below, so there is
# exactly one place a future reader can weaken the determinism envelope and
# exactly one place to audit it.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export LBW_BENCH_GATE=1
# testthat::test_dir() does not set NOT_CRAN itself (only devtools::test()'s
# pkgload/test_check machinery does); every bench-gate test opens with
# skip_on_cran(), so without this every assertion below silently reports
# SKIP ("On CRAN") instead of actually asserting, regardless of
# LBW_BENCH_GATE. Exporting it here is what makes this an executed gate
# rather than a no-op.
export NOT_CRAN=true

echo '=== stepstone regression gate ===' && Rscript -e "testthat::test_dir('tests/testthat', filter='bench-gate', stop_on_failure=TRUE)"
echo '=== oris_soft vs survey::calibrate (medium_100k_5margins) ===' && Rscript benchmarks/oris_soft_vs_competitors.R
echo '=== greenkhorn/sinkhorn vs POT (k2_margin_pot_equiv) ===' && python/.venv/bin/python benchmarks/greenkhorn_sinkhorn_vs_pot.py

echo "Wrote benchmarks/results/oris_soft_vs_competitors.csv and benchmarks/results/oris_soft_vs_competitors_env.txt"
echo "Wrote benchmarks/results/greenkhorn_sinkhorn_vs_pot.csv and benchmarks/results/greenkhorn_sinkhorn_vs_pot_env.txt"
