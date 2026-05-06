#!/usr/bin/env bash
set -e
echo '=== R benchmark ===' && Rscript benchmarks/allmethod_bench.R
echo '=== Python benchmark ===' && python benchmarks/allmethod_bench.py
