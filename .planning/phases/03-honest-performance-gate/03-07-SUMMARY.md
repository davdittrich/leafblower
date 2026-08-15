---
phase: 03-honest-performance-gate
plan: 07
subsystem: benchmarks
tags: [honest-performance-gate, greenkhorn, sinkhorn, pot, python, gap-closure]
requirements: [US-003, KPI-04]
gap_closure: true
gap_ids: [G-03-1]
dependency-graph:
  requires:
    - benchmarks/oris_soft_vs_competitors.R (CSV schema, provenance-capture pattern)
    - python/leafblower/_harvest.py (harvest() Python API)
  provides:
    - benchmarks/greenkhorn_sinkhorn_vs_pot.py
    - benchmarks/results/greenkhorn_sinkhorn_vs_pot.csv
    - benchmarks/results/greenkhorn_sinkhorn_vs_pot_env.txt
  affects:
    - benchmarks/run_honest_gate.sh
tech-stack:
  added:
    - "POT 0.9.7.post1 (Python Optimal Transport, benchmark-only, uv-managed venv, not in pyproject.toml)"
  patterns:
    - "M=-log(prior), reg=1 kernel-equivalence construction to make POT's 2-marginal entropic-OT solve the identical KL(X||X_init) objective leafblower's greenkhorn/sinkhorn solve for K=2"
    - "implied per-observation weight T[i,j]/K_prior[i,j] to give POT's cell-only solution the same deff/n_eff/max_w/min_w columns as leafblower's per-observation output"
key-files:
  created:
    - benchmarks/greenkhorn_sinkhorn_vs_pot.py
    - benchmarks/results/greenkhorn_sinkhorn_vs_pot.csv
    - benchmarks/results/greenkhorn_sinkhorn_vs_pot_env.txt
  modified:
    - benchmarks/run_honest_gate.sh
decisions:
  - "Same skewed target (0.40/0.30/0.20/0.10, normalized) applied to both m1 and m2 margins, matching the plan's own example proportions."
  - "max_weight column left NA for pot_* rows (POT has literally no bounds parameter) rather than reporting the same N cap leafblower's arms use; the cap only has meaning on the leafblower side."
  - "POT rows' max_w/min_w/deff/n_eff derived from an implied per-observation weight T[i,j]/K_prior[i,j], reusing leafblower.design_effect/effective_sample_size on that vector, rather than leaving those 4 columns NA -- the K=2 IPF equivalence this plan demonstrates makes observations within a cell exchangeable, so this is the natural (not invented) extension, not an additional analytical degree of freedom."
metrics:
  duration: "~35m"
  completed: "2026-08-16"
status: complete
actuals:
  tokens: 32000
  tasks: 3
  commits: 3
---

# Phase 03 Plan 07: greenkhorn/sinkhorn vs POT (K=2 gap closure) Summary

Built a K=2-margin fixture and an `M=-log(prior)`, `reg=1` kernel-equivalence
construction that makes POT's `ot.bregman.greenkhorn`/`ot.sinkhorn` solve the
identical `KL(X||X_init)` objective leafblower's `greenkhorn`/`sinkhorn` solve
in the pure 2-marginal case, closing G-03-1's per-method coverage gap for the
two Python-reachable leafblower methods.

## What was built

`benchmarks/greenkhorn_sinkhorn_vs_pot.py` builds a `n=10000`, 4-level-per-margin
fixture with a skewed target (0.40/0.30/0.20/0.10, normalized) for both `m1`
and `m2`, computes the observed cross-tabulation `K_prior` (leafblower's
implicit `X_init` prior), and constructs `M = -log(K_prior)` so that POT's
kernel `exp(-M/reg=1)` equals `K_prior` exactly -- the same prior leafblower's
own solvers start from. Both leafblower's arms (`max_weight` set to `n`,
effectively unbounded) and POT's arms (which have no bounds mechanism at all)
run 4 total measurements: `leafblower_greenkhorn`, `leafblower_sinkhorn`,
`pot_greenkhorn`, `pot_sinkhorn`.

**Measured result confirms the equivalence claim empirically, not just
theoretically**: all four arms converge to the identical `max_w=2.424146...`,
`n_eff=7110.59...` solution, with `max_error` at or near machine precision
for three of the four arms (`leafblower_sinkhorn`: 1.7e-16, `pot_greenkhorn`:
6.3e-15, `pot_sinkhorn`: 5.6e-17; `leafblower_greenkhorn`: 5.2e-9, limited by
its own `improvement`-rule convergence tolerance, not the construction).

The BLAS single-thread determinism guard runs before any `numpy`/`ot`
import (`import ot` transitively imports `numpy`) and refuses to measure
(raising `RuntimeError` naming the offending var) unless
`OMP_NUM_THREADS`/`OPENBLAS_NUM_THREADS`/`MKL_NUM_THREADS` are all `"1"` --
verified as a negative control (fires when unset).

POT has no per-observation identity, only a solved 4x4 cell-mass table, so
`max_w`/`min_w`/`deff`/`n_eff` for the `pot_*` rows are derived from an
implied per-observation weight `T[i,j]/K_prior[i,j]` (observations sharing a
cell are exchangeable under this margin-only objective -- the same
cell-uniform-weight structure leafblower's own unbounded unit-mode
water-fill produces), reusing `leafblower.design_effect`/
`leafblower.effective_sample_size` on that vector rather than leaving those
columns NA. `max_weight` itself is left NA for `pot_*` rows since POT has no
such parameter at all.

The unbounded-comparison caveat ("POT has no bounds mechanism; both sides
run effectively unbounded ... NOT a test of leafblower's normal bounded,
K>2 workload") is embedded verbatim in every row's `note`, per the plan's
instruction and the threat register's T-03-07-02 mitigation.

`benchmarks/run_honest_gate.sh` gained one new step invoking this script,
reusing the script's already-exported single-thread BLAS envelope (no
re-export inside the new line).

## Verification

- `bash benchmarks/run_honest_gate.sh` runs all three steps (stepstone
  regression gate, R competitor script, this new Python script) and exits 0.
- `benchmarks/results/greenkhorn_sinkhorn_vs_pot.csv` has a header + 4 data
  rows, same 17-column schema as `oris_soft_vs_competitors.csv`
  (`input_class="k2_margin_pot_equiv"`, `n_margins=2`).
- `git diff --stat python/pyproject.toml` is empty -- POT stayed
  benchmark-only (`uv pip install POT` into the venv, never added as a
  shipped dependency).
- Python DoD (`cd python && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
  MKL_NUM_THREADS=1 .venv/bin/python -m pytest`): 161 passed, 0 failed.

## Deviations from Plan

### Auto-fixed Issues

None -- plan executed as written.

### Discrepancies reported (not silently reconciled)

**1. Python DoD test count: 161 passed, not the plan's cited "160 passed."**
- **Found during:** Task 3 verification.
- **Detail:** The plan's `<verification>` section states "Python DoD: 160
  passed, 0 failed." The actual run (both before and after this plan's
  changes; no test file was touched by any of this plan's 3 tasks) reports
  161 passed, 0 failed. This is pre-existing drift from a later, unrelated
  test addition in the repo's history since the plan was authored, not a
  regression introduced here. Reported per the project's own convention of
  surfacing measured-vs-cited discrepancies rather than silently forcing a
  match.
- **Files modified:** none (no test file touched).
- **Commit:** n/a (observation only).

**2. `oris_soft_vs_competitors.csv` regenerated as a side effect of running
`bash benchmarks/run_honest_gate.sh` (Task 3's own verify step), producing a
timing-noise-only diff (wall_s values only; every accuracy/weight field
byte-identical) in a file outside this plan's `files_modified` scope.**
- **Action taken:** Restored via `git checkout -- benchmarks/results/oris_soft_vs_competitors.csv`
  rather than committing an out-of-scope, noise-only change under this
  plan's history.

## Threat Flags

None -- the threat register's T-03-07-SC (package legitimacy), T-03-07-01
(information disclosure), and T-03-07-02 (repudiation/equivalence-claim
verifiability) were all pre-mitigated per the plan; no new surface was
introduced beyond what the threat model anticipated.

## Self-Check

- `benchmarks/greenkhorn_sinkhorn_vs_pot.py` — FOUND
- `benchmarks/results/greenkhorn_sinkhorn_vs_pot.csv` — FOUND
- `benchmarks/results/greenkhorn_sinkhorn_vs_pot_env.txt` — FOUND
- `benchmarks/run_honest_gate.sh` (modified) — FOUND
- Commit `6e9c0f1` (Task 1) — FOUND in `git log`
- Commit `7636bf0` (Task 2) — FOUND in `git log`
- Commit `78103f8` (Task 3) — FOUND in `git log`

## Self-Check: PASSED
