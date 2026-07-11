# Study branch: `study/benchmark-instrumented-DO-NOT-MERGE`

Benchmark study only (epic **leafblower-2ouc**). **Never merges to `main`/`master`.**

## Why this branch exists
The article-grade benchmark needs two things forbidden on the shippable package:

1. **Trajectory instrumentation** — per-iteration convergence-metric emission
   from `harvest()` → C API → solver loops, so WU-12 can render RQ3 convergence
   curves. Production `harvest()` deliberately exposes no per-iteration trace.
2. **`metrics.py` native-divergence extension** — chi2 (Pearson `Σ(w−d)²/d`) and
   logit (binary-entropy) native objectives, so each family's work-precision is
   judged on its own divergence. Re-frozen as **benchmark-runnable-freeze-v14**
   (study-branch tag).

Base: `847aad0` (= benchmark-runnable-freeze-v13, master).

## Safeguards (local `.git/hooks`; repo has no remote)
- `.NEVER_MERGE_TO_MAIN` sentinel at repo root (this branch only).
- Every study-only source edit carries `STUDY-BRANCH-ONLY-DO-NOT-MERGE`.
- `pre-commit` + `pre-merge-commit` ABORT if either the sentinel or the marker
  would land on `main`/`master`.
- Do not delete the sentinel; do not `--no-verify`.

## What is safe to merge
Nothing from this branch. Downstream deliverables (results parquets, report,
gap analysis) are data/report artifacts under `benchmarks/study/` produced by
running this branch; they are consumed as files, not merged as code.
