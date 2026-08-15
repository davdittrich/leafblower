---
phase: 2
slug: one-engine-not-two
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-08-15
---

# Phase 2 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | R: testthat edition 3 (`DESCRIPTION:22`); Python: pytest |
| **Config file** | `DESCRIPTION` (R); `python/pyproject.toml` + `python/conftest.py` (Python) |
| **Quick run command** | R: `Rscript -e 'library(testthat); library(leafblower); test_dir("tests/testthat", filter="<name>")'`; Python: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -m pytest -k <name> -q` |
| **Full suite command** | `R CMD INSTALL --preclean . && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript -e 'library(testthat); library(leafblower); out <- test_dir("tests/testthat", stop_on_failure=TRUE)' && cd python && uv pip install -e . --reinstall-package leafblower && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -m pytest -q` (`.coverage-thresholds.json`'s `enforcement.command`, verbatim) |
| **Estimated runtime** | ~180 seconds |

---

## Sampling Rate

- **After every task commit:** targeted solver's own R testthat file(s) + `python/leafblower/test_solver_parity.py` / `test_parity_weights.py` for that method (per D-01 solver-by-solver migration order)
- **After every plan wave:** Full DoD gate (`.coverage-thresholds.json`'s `enforcement.command`)
- **Before `/gsd-verify-work`:** Full suite must be green, plus `LBW_BENCH_GATE=1` stepstone benchmark on any commit touching a TU boundary
- **Max feedback latency:** 180 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 02-01-01 | 01 | 0 | SC4 | — | N/A | regression-prevention | `python/leafblower/test_core_sources_sync.py` (new, per D-05) | ❌ Wave 0 | ⬜ pending |
| 02-01-02 | 01 | 0 | SC3 | — | N/A | verification | grep/reflection assertion that each solver `.cpp` calls `finalize_weights[_buf]` (new, per D-03) | ❌ Wave 0 | ⬜ pending |
| 02-0N-0M | TBD | TBD | SC1 (dispatch table + field-gap survival) | — | N/A | regression (existing) | `tests/testthat/test-newton-kl.R`, `test-newton-tsvd-projection.R`, `test-newton-kl-tsvd-ratio.R`, `test-cr-d5-auto-fallback-fields.R`, `python/leafblower/test_solver_parity.py`, `test_parity_weights.py` | ✅ | ⬜ pending |
| 02-0N-0M | TBD | TBD | SC2 (opt-level symmetry, documented if asymmetric) | — | N/A | build-config check | manual + DoD gate | ✅ (existing gate) | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*
*Planner fills exact Task ID / Plan / Wave columns once the solver migration order (D-01) is sequenced into waves.*

---

## Wave 0 Requirements

- [ ] `python/leafblower/test_core_sources_sync.py` (or R-side equivalent) — asserts `src/*.cpp` (minus `r_bridge.cpp`) == `python/CMakeLists.txt` `CORE_SOURCES`. Covers SC4.
- [ ] New SC3 verification test/script asserting all solvers (7, or 8 pending the `newton_kl` bounds-mode scope decision — see RESEARCH.md Pitfall 2) route through `calib_dispatch.hpp::finalize_weights`/`finalize_weights_buf`. Covers SC3.
- [ ] `leafblower-rywn`'s own Step 1 field inventory (a `bd comment`, not a test file) runs before any solver migration task starts.

---

## Manual-Only Verifications

*None — all phase behaviors have automated verification (existing R/Python parity suites + two new Wave 0 tests).*

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 180s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
