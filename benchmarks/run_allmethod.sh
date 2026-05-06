#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
echo '=== R benchmark ===' && Rscript benchmarks/allmethod_bench.R
echo '=== Python benchmark ===' && python benchmarks/allmethod_bench.py
