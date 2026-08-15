# Testing Patterns

**Analysis Date:** 2026-08-15

## Test Framework

**R - testthat v3:**
- Runner: `devtools::test()` or `R CMD INSTALL --preclean .` then `Rscript -e "devtools::test()"`
- Config: none required; tests discovered in `tests/testthat/test-*.R`
- Run commands:
  ```bash
  Rscript -e "devtools::test()"     # Run all tests
  Rscript -e "devtools::test(filter='raking')"  # Run specific file pattern
  ```

**Python - pytest:**
- Runner: `pytest` (installed via uv)
- Config: `python/conftest.py` sets single-threaded BLAS and removes local path
- Run commands:
  ```bash
  cd python && uv pip install -e . --reinstall-package leafblower
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python -m pytest
  python -m pytest -v                          # Verbose
  python -m pytest python/leafblower/test_solver_parity.py -v  # Specific file
  ```

**Assertion Library:**
- R: `testthat::expect_*` functions (expect_true, expect_equal, expect_match, expect_error)
- Python: `pytest` assertions + `np.allclose()` for numerical parity

## Test File Organization

**Location:**
- R: `tests/testthat/test-*.R` — one test file per feature/module
- Python: `python/leafblower/test_*.py` — one test file per feature/module
- Shared: test data fixtures in `tests/testthat/fixtures/` (JSON, .rds, R scripts)

**Naming Convention:**
- R: `test-<feature>.R` (e.g., `test-raking.R`, `test-oris.R`)
- Python: `test_<feature>.py` (e.g., `test_solver_parity.py`, `test_harvest_na_parity.py`)
- Edge case tests: `tests/testthat/_problems/test-<issue-id>.R`

**Structure - R:**
```r
test_that("descriptive test name", {
  # Setup
  set.seed(42)
  df <- data.frame(x = factor(...))
  tgt <- list(x = c(a=0.5, b=0.5))
  
  # Execute
  result <- leafblower::harvest(df, tgt, method = "raking", ...)
  
  # Assert
  expect_true(attr(result, "algorithm") == "raking")
  expect_equal(mean(result$weights), 1.0, tolerance = 1e-10)
})
```

**Structure - Python:**
```python
def test_harvest_returns_copy():
    """weights_out must be a copy, not a view into input."""
    # Setup
    n = 100
    weights = np.ones(n, dtype=np.float64)
    gids = [np.zeros(n, dtype=np.int32)]
    
    # Execute
    status, weights_out, res = calibrate(...)
    weights_out[0] = 9999.0
    
    # Assert
    assert weights[0] != 9999.0, "weights_out must be a copy"
```

## Test Structure

**Suite Organization - R:**
- One `test_that()` block per logical test case
- Descriptive name includes: test subject + expected behavior
- Setup/Execute/Assert pattern (Arrange/Act/Assert)
- Example: `test_that("raking respects max_weight=2 on tight bounds", { ... })`

**Suite Organization - Python:**
- One `def test_*():` function per logical test case
- Docstring explains R↔Python parity or test motivation
- Same Arrange/Act/Assert pattern
- Example from `test_python.py`:
  ```python
  def test_sor_omega_max_is_wired():
      """eb79.1: sor_omega_max must reach the C solver (regression guard)."""
  ```

**Setup:**
- R: `set.seed(42)` for reproducibility; build minimal data fixtures inline
- Python: embedded fixture constants (e.g., `_A`, `_B` arrays in `test_solver_parity.py`)
- Both: single-threaded BLAS enforced at test runner start

**Teardown:**
- R: testthat auto-cleans up; no explicit cleanup needed
- Python: pytest fixture cleanup via `conftest.py`

**Assertion Patterns - R:**
```r
expect_true(condition)
expect_equal(actual, expected, tolerance = 1e-10)
expect_match(string, pattern)
expect_lt(x, threshold)
expect_error(expr, pattern)
expect_silent(expr)
```

**Assertion Patterns - Python:**
```python
assert condition
assert np.allclose(actual, expected, rtol=1e-6, atol=0.0)
with pytest.raises(ValueError, match="pattern"):
    some_function()
np.testing.assert_allclose(actual, expected, rtol=1e-6)
```

## Mocking

**Framework - R:**
- testthat does not use mocks; uses real solvers and verification instead
- Fixtures validated against R subprocess calls (parity tests)

**Framework - Python:**
- pytest fixtures for test data
- No mock framework in use; real C++ solver calls
- Subprocess validation for cross-language parity via `subprocess.run()`

**Patterns:**
- Direct solver invocation (`calibrate()`, `harvest()`)
- Verify expected errors via exception matching
- Subprocess parity tests: run R script, parse JSON output, compare

**What to Mock:**
- Nothing; all tests call real C++ core

**What NOT to Mock:**
- Solver logic (calibrate, harvest)
- Parameter parsing/validation (mirror R in Python)
- Convergence checks

## Fixtures and Factories

**Test Data - R:**
Location: `tests/testthat/fixtures/`

Examples:
- `oris_baseline_snapshot.R`: shell script that generates `oris_fixed_omega_baseline.json`
- `stepstone_reference_run.R`: runs stepstone benchmark
- `oris_shipgate_fixture.R`: builds problem instances for ship-gate validation

Pattern:
```r
run_case <- function(df, tgt, extra = list()) {
  base_args <- c(list(data = df, target = tgt, method = "oris"), extra)
  res <- do.call(harvest, base_args)
  list(status = as.integer(info$status), iters = info$iterations, ...)
}
```

**Test Data - Python:**
Location: embedded constants in test files

Examples from `test_solver_parity.py`:
```python
_A = ["x", "x", "x", ...]  # 300-element category vector
_B = ["q", "p", "p", ...]  # 300-element category vector

_TARGETS = {
    "a": {"x": 1/3, "y": 1/3, "z": 1/3},
    "b": {"p": 0.5, "q": 0.5},
}

_CONV_PY = {"rule": "improvement", "tol": 0.001}
_CONV_R  = 'list(rule="improvement", tol=0.001)'
```

**Factory Pattern:**
None used; tests build minimal data inline:
```python
df = pd.DataFrame({"x": ["a","b","a","b"]})
tgts = {"x": {"a": 0.5, "b": 0.5}}
```

## Coverage

**Requirements:** 
- No explicit coverage target in `.coverage-thresholds.json` (gate uses behavioral validation)
- BEHAVIORAL gate: R tests pass (0 FAIL) + Python tests pass (0 FAIL) + stepstone benchmark shows no regression
- No pytest-cov or covr integration

**View Coverage:**
```bash
# R: none configured
# Python: none configured
# Validation gate:
R CMD INSTALL --preclean .
Rscript -e "devtools::test()"
cd python && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python -m pytest
LBW_BENCH_GATE=1 ./tools/stepstone-benchmark.R  # optional regression gate
```

## Test Types

**Unit Tests - R:**
- Located: `tests/testthat/test-*.R`
- Scope: individual solver functions, parameter parsing, edge cases
- Examples: `test-raking.R`, `test-oris.R`, `test-convergence-criteria.R`
- Pattern: call public API (`harvest()`), verify output and attributes

**Unit Tests - Python:**
- Located: `python/leafblower/test_*.py`
- Scope: individual API functions, Python-layer logic, parameter validation
- Examples: `test_python.py`, `test_design_effect.py`
- Pattern: call public API (`harvest()`, helper functions), verify outputs

**Integration Tests - R:**
- Located: `tests/testthat/test-*.R` (mixed with unit tests)
- Scope: multi-method comparison, bounds mode interaction, bounds + convergence
- Examples: `test-harvest-bounds-mode.R`, `test-oris-bounds-mode.R`

**Integration Tests - Python:**
- Located: `python/leafblower/test_*.py`
- Scope: end-to-end calibration, parity with R
- Examples: `test_solver_parity.py`, `test_harvest_na_parity.py`

**E2E / Parity Tests:**
- Framework: pytest with R subprocess backend
- Pattern: run identical problem in Python and R, compare weights
- Tolerance: `rtol=1e-6, atol=0.0` (numpy.allclose)
- Pre-check: both must converge (max_error < 0.01) before asserting parity
- Examples from `test_solver_parity.py`:
  ```python
  def test_newton_kl_parity():
      r = subprocess.run(["Rscript", "test_solver_parity.R"], capture_output=True)
      w_r = json.loads(r.stdout)
      w_py = harvest(df, _TARGETS, method="newton_kl", convergence=_CONV_PY)
      np.testing.assert_allclose(w_r, w_py, rtol=1e-6, atol=0.0)
  ```

**Problem-Specific / Regression Tests:**
- Located: `tests/testthat/_problems/test-<issue-id>.R`
- Scope: specific bugs fixed, edge cases discovered in production
- Examples: `test-ieppa-nonuniform-d-28.R`, `test-raking-89.R`
- Pattern: reproduce exact problem, verify it no longer occurs

## Common Patterns

**Async Testing:**
Not applicable (single-threaded, deterministic solvers)

**Error Testing - R:**
```r
test_that("min_weight > max_weight raises error", {
  expect_error(
    harvest(df, tgt, min_weight = 2.0, max_weight = 1.0),
    "min_weight must be <= max_weight"
  )
})
```

**Error Testing - Python:**
```python
def test_min_weight_badarg_python():
    with pytest.raises(Exception):
        harvest(df, tgts, min_weight=5.0, max_weight=5.0)
```

**Determinism / Reproducibility:**
- All randomness seeded: `set.seed(42)` in R, embedded constants in Python
- Single-threaded BLAS enforced:
  ```bash
  export OMP_NUM_THREADS=1
  export OPENBLAS_NUM_THREADS=1
  export MKL_NUM_THREADS=1
  ```
- Python `conftest.py` sets these before numpy import

**Convergence Verification:**
- Tests use explicit convergence specs to lock stopping point
- Example:
  ```python
  convergence = {"rule": "improvement", "tol": 0.001}
  ```
- Precheck: verify max_error < tolerance before asserting parity
- Stall detection: verify that constrained optima are detected, not false convergence

**Tolerance Strategy:**
- Weight parity: `rtol=1e-6, atol=0.0` (relative 1e-6)
- Convergence precheck: max_error < 0.01 (covers all solver types)
- Absolute error tolerance: depends on metric (1e-6 for max_err, 1e-10 for individual weight constraints)

---

*Testing analysis: 2026-08-15*
