# Deferred Items — Phase 1

Out-of-scope discoveries logged during plan execution, per the executor's
scope-boundary rule (only auto-fix issues directly caused by the current
task's changes).

## From 01-01 (weight-vector parity relocation + extension)

- **`python/leafblower/test_trajectory_csv_smoke.py::test_trajectory_csv_smoke` fails
  pre-existing, unrelated to this plan.** `assert lines[0] == "iter,errRp"` fails because
  the actual CSV header is `iter,errRp,marginal_kl` — `src/oris_trajectory.cpp`'s writer
  gained a `marginal_kl` column (last touched in commit `f4b56bd`, before this plan's HEAD)
  and the test assertion was never updated to match. Confirmed pre-existing: the file this
  plan modified (`python/leafblower/test_parity_weights.py`) is unrelated, and
  `test_trajectory_csv_smoke.py` / `src/oris_trajectory.cpp` were untouched by any of this
  plan's three commits. Not fixed here — 01-01's file scope is
  `tests/test_parity_weights.py` / `python/leafblower/test_parity_weights.py` only, no
  `src/` change is permitted in this phase (test-layer only). Needs a beads ticket: either
  widen the assertion to `"iter,errRp,marginal_kl"` or `.startswith("iter,errRp")`.
