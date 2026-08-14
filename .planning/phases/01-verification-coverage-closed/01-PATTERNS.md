# Phase 1: Verification Coverage Closed - Pattern Map

**Mapped:** 2026-08-15
**Files analyzed:** 5
**Analogs found:** 5 / 5

## File Classification

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|--------------------|------|-----------|-----------------|----------------|
| `python/leafblower/test_parity_weights.py` (relocated from `tests/test_parity_weights.py`, `leafblower-x7n8`) | test | file-I/O (subprocess) | `python/leafblower/test_solver_parity.py` | exact (same directory target, same subprocess-parity role) |
| `python/leafblower/test_parity_weights.py` — parametrize extension (+chebyshev, +greg, +oris_soft) | test | request-response (subprocess) | itself, `tests/parity/run_oris_soft_r.R` + `run_chebyshev_r.R` helper pattern | exact |
| `python/leafblower/test_solver_parity.py` — new `test_raking_*`, `test_sinkhorn_*` | test | request-response (subprocess) | `test_logit_default_rule_parity` / `_assert_parity` in same file | exact |
| `tests/testthat/test-bound-property.R` (new, KPI-02, D-01..D-05) | test | batch (50-dataset property sweep) | `tests/testthat/test-clamp-contract.R` + `tests/testthat/test-unit-bounds-status-consistency.R` | role-match (bound-assertion idiom + mixture-weight generator idiom) |
| `DESCRIPTION` (SC5, `+Config/testthat/edition: 3`) | config | — | n/a (DCF field addition) | n/a — no analog needed |

## Pattern Assignments

### `python/leafblower/test_parity_weights.py` (relocated + extended) (test, request-response/file-I/O)

**Analog:** `python/leafblower/test_solver_parity.py` (sibling file, same directory, same role)

**Single-thread BLAS env-guard pattern** — MUST be set before `numpy`/`leafblower` import (current file, lines 10-15, keep verbatim after move):
```python
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pandas as pd
import pytest
```

**REPO_ROOT path fix required on relocation (Pitfall 2 from RESEARCH.md)** — current (pre-move) computation, `tests/test_parity_weights.py:29-30`:
```python
REPO_ROOT = Path(__file__).resolve().parent.parent
R_HELPER  = REPO_ROOT / "tests" / "parity" / "run_parity_r.R"
```
After `git mv` to `python/leafblower/test_parity_weights.py`, `.parent.parent` resolves to `.../python`, not repo root — `R_HELPER`/`_ORIS_SOFT_R_HELPER`/`_CHEBYSHEV_R_HELPER` (lines 30, 103, 254) must be recomputed to still point at `<repo_root>/tests/parity/*.R`. One extra `.parent` (`Path(__file__).resolve().parent.parent.parent`) is the fix — verify with `R_HELPER.exists()` assertion or a collection-count check, not by trusting the arithmetic.

**Parametrize extension pattern (+chebyshev, +greg, +oris_soft) — SC1** — current list, `tests/test_parity_weights.py:73-75`:
```python
@pytest.mark.skipif(not RSCRIPT_AVAILABLE, reason="Rscript not found")
@pytest.mark.parametrize("method", [
    "greenkhorn", "logit", "raking", "oris", "sinkhorn", "newton_kl",
])
def test_weight_parity(method, tmp_path):
```
Extend the list to the full 9-solver `rk_algorithm_t` non-AUTO set. Note `chebyshev` and `oris_soft` already have DEDICATED parity tests further down this file (`test_oris_soft_default_tol_parity` line 147, `test_chebyshev_default_tol_parity` line 266) using bespoke correlated fixtures and R helper scripts (`run_oris_soft_r.R`, `run_chebyshev_r.R`) — decide whether adding them to the generic `test_weight_parity` matrix duplicates or complements those; `greg` has no existing R-helper wiring in `tests/parity/`, check `run_parity_r.R` supports method="greg" before adding.

**Skip-on-missing-Rscript guard (every R-subprocess test in repo uses this):**
```python
RSCRIPT_AVAILABLE = shutil.which("Rscript") is not None
...
if result.returncode == 2:
    pytest.skip(f"R package not available: {result.stderr.strip()}")
if result.returncode != 0:
    raise RuntimeError(f"Rscript failed:\n{result.stderr}")
```

**Tolerance ternary needing diagnosis-driven comment (SC3, D-06/D-07)** — `tests/test_parity_weights.py:93-96`:
```python
tol = 1e-6 if method == "logit" else 1e-10
assert diff < tol, (
    f"{method}: max|w_py - w_r| = {diff:.2e} (threshold {tol:.0e})"
)
```
Do not touch the number without first measuring per Pitfall 3 — annotate with the mechanism found (build `-O3` asymmetry vs. Newton/Armijo conditioning) or replace with a value tied to a filed ticket.

---

### `python/leafblower/test_solver_parity.py` — new raking/sinkhorn tests (test, request-response)

**Analog:** same file, `test_logit_default_rule_parity` (lines 212-242) as the template to mirror per RESEARCH.md's explicit recommendation.

**Core four-step protocol to reuse verbatim** (`python/leafblower/test_solver_parity.py:171-186`):
```python
def _assert_parity(method: str, conv_py: dict = _CONV_PY, conv_r: str = _CONV_R):
    w_py, ri_py = _run_py(method, conv_py)
    assert ri_py.get("max_error", 1.0) < _CONV_TOL, (
        f"Python {method} did not converge: max_error={ri_py.get('max_error')}"
    )
    r_out = _run_r(method, conv_r)
    assert r_out["max_error"] < _CONV_TOL, (
        f"R {method} did not converge: max_error={r_out['max_error']}"
    )
    w_r = np.array(r_out["weights"])
    assert len(w_py) == len(w_r), (
        f"{method}: length mismatch Python={len(w_py)} R={len(w_r)}"
    )
    assert np.allclose(w_r, w_py, rtol=1e-6, atol=0.0), (
        f"{method} R↔Python mismatch: max|Δw|={np.max(np.abs(w_r - w_py)):.3e}"
    )
```

**Simple per-method test wrapper pattern** (lines 245-247):
```python
def test_chebyshev_parity():
    """chebyshev: R↔Python weights match to rtol=1e-6."""
    _assert_parity("chebyshev")
```
Copy this shape for `test_raking_parity()` / `test_sinkhorn_parity()`.

**Default-rule lock pattern** (lines 231-242, when a per-method default resolution needs pinning — likely relevant since RESEARCH.md flags raking's Boyle-Dykstra→KL-Bregman geometry change):
```python
_, ri_py = _run_py("logit", conv_py={})
assert ri_py.get("convergence_used", {}).get("rule") == "improvement", (...)
r_out = _run_r("logit", conv_r="list()")
assert r_out.get("rule") == "improvement", (...)
_assert_parity("logit", conv_py={}, conv_r="list()")
```

**Shared fixture and _CONV_TOL rationale (lines 78-107)** — the embedded categorical fixture (`_A`, `_B`, `_TARGETS`) and the per-method `_CONV_TOL = 0.01` are already calibrated against measured max_error values for 5 solvers; adding raking/sinkhorn requires measuring and documenting their max_error on this SAME fixture in the same comment block, not assuming 0.01 covers them blind.

---

### `tests/testthat/test-bound-property.R` (new file, KPI-02) (test, batch/property)

**Analog 1 — bound-assertion idiom:** `tests/testthat/test-clamp-contract.R:51-58`
```r
for (m in c("raking", "sinkhorn", "greg", "oris", "chebyshev")) {
  w <- as.numeric(harvest(df, target, method = m, design_weights = dw,
                          max_weight = 2.0, min_weight = 0.3, bounds_mode = "unit",
                          max_iterations = 1000L, attach_weights = FALSE))
  expect_true(all(w <= 2.0 + 1e-9),  info = sprintf("%s: max=%.6f > 2.0", m, max(w)))
  expect_true(all(w >= 0.3 - 1e-9),  info = sprintf("%s: min=%.6f < 0.3", m, min(w)))
}
```
KPI-02 requires `1e-10` tolerance, not this file's `1e-9` — do not copy the tolerance constant, only the loop/assertion shape. Note `bounds_mode = "unit"` is mandatory (Pattern 2: cell-mode is not unconditionally bound-enforcing).

**Analog 2 — skewed/fat-tailed design-weight generator idiom:** `tests/testthat/test-clamp-contract.R:48-49`
```r
dw <- 2 ^ (rnorm(n))            # log-normal design weights: heavy skew per cell
dw[sample(n, 60L)] <- 50.0      # extreme outliers force cascading clamps
```
and the project's `rlnorm`-based mixture idiom already in-repo, `tests/testthat/test-unit-bounds-status-consistency.R:17-19`:
```r
d <- rep(1, n); d[x == "B"] <- rlnorm(sum(x == "B"), meanlog = 0, sdlog = 1.6)
```
D-03 requires lognormal bulk + heavy-tailed contaminant fraction (not Gaussian) — combine `rlnorm(n, meanlog, sdlog)` for the bulk with `rcauchy`/`rt(df=1)` (base `stats`, already `Imports: stats` per `DESCRIPTION:18`) replacing a small index subset, mirroring the `dw[sample(n, 60L)] <- 50.0` outlier-injection shape above rather than a second full-vector draw.

**Analog 3 — sparse-cells / skewed-marginal generator (D-04's second stratum):** `test-unit-bounds-status-consistency.R:16-21` (`.drift_fixture`):
```r
.drift_fixture <- function() {
  set.seed(101L); n <- 2000L
  x <- factor(sample(c("A", "B"), n, TRUE, prob = c(0.85, 0.15)))
  y <- factor(sample(c("P", "Q"), n, TRUE, prob = c(0.5, 0.5)))
  d <- rep(1, n); d[x == "B"] <- rlnorm(sum(x == "B"), meanlog = 0, sdlog = 1.6)
  list(df = data.frame(x = x, y = y), d = d,
       tgt = list(x = c(A = 0.55, B = 0.45), y = c(P = 0.5, Q = 0.5)))
}
```
This is the near-empty-cell / skewed-category-marginal shape (`prob = c(0.85, 0.15)`) for the "sparse-cells-only" stratum in D-04.

**Fixed-seed requirement (D-02):** every fixture-generating helper in this file MUST call `set.seed(<fixed literal>)` per dataset (see `.drift_fixture`'s `set.seed(101L)`, `.feasible_fixture`'s `set.seed(3L)`) — 50 distinct fixed literals, never a loop-derived or randomized seed.

**File header/naming convention** (no `context()` — new files should NOT add the deprecated call; `test-clamp-contract.R` and `test-unit-bounds-status-consistency.R` both already omit it, confirming file-based grouping is the current convention for new files even before the edition-3 flip).

---

### `DESCRIPTION` (SC5, config)

**Current relevant fields** (`DESCRIPTION:18-19`):
```
Imports: stats
Suggests: autumn (>= 0.2.0), testthat (>= 3.0.0), bench, lhs, DiceKriging, ggplot2, rprojroot,
```
Add `Config/testthat/edition: 3` as a new DCF field (no existing analog in this file — it is a bare field addition). Sequence per D-12: fix all deprecated constructs first (see Shared Patterns below), commit green under 2e, THEN add this field as its own commit.

## Shared Patterns

### Deprecated `context()` calls — must be removed/replaced before the 3e flip (D-12)
**Source:** `tests/testthat/test-oris.R:1`
```r
context("oris (faithful algBCD)")
```
**Apply to:** 11 files flagged in RESEARCH.md (`test-alm-config-grouping.R`, `test-oris-nonuniform-d.R`, `test-oris.R`, `test-oris-sraa-log-path.R`, `test-oris-sraa.R`, and 6 more — grep `^context(` across `tests/testthat/*.R` for the authoritative list at execution time). Edition 3 auto-derives context from filename; the fix is deleting the `context(...)` line, not replacing it.

### R-side single-thread BLAS / subprocess env propagation
**Source:** `tests/test_parity_weights.py:167-172` (Chebyshev/oris_soft R subprocess calls)
```python
_single_thread_env = {
    **os.environ,
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
}
result = subprocess.run([...], env=_single_thread_env, timeout=180)
```
**Apply to:** any new R-subprocess call added for raking/sinkhorn parity tests or the relocated file.

### Cross-language float comparison idiom
**Source:** used identically throughout `test_solver_parity.py` and `test_parity_weights.py`
```python
max_abs_diff = np.max(np.abs(w_r - w_py))
assert max_abs_diff < tol, f"... max|Δw|={max_abs_diff:.3e}"
```
**Apply to:** all new parity assertions — do not hand-roll a different epsilon-diff helper (RESEARCH.md "Don't Hand-Roll" table).

## No Analog Found

None — all 5 files/changes have a direct in-repo analog. This phase adds no new architectural role (no new controller/service/model), only test-layer extensions and one config field, consistent with CONTEXT.md's "Test-layer only" phase boundary.

## Metadata

**Analog search scope:** `tests/testthat/`, `python/leafblower/test_*.py`, `tests/test_parity_weights.py`, `tests/parity/*.R`, `DESCRIPTION`
**Files read this session:** `01-CONTEXT.md`, `01-RESEARCH.md`, `python/leafblower/test_solver_parity.py`, `tests/test_parity_weights.py`, `tests/testthat/test-clamp-contract.R` (via prior research context), `tests/testthat/test-unit-bounds-status-consistency.R`, `tests/testthat/test-oris.R` (header), `DESCRIPTION` (targeted grep)
**Pattern extraction date:** 2026-08-15
