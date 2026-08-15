---
phase: 2
slug: one-engine-not-two
status: validated
nyquist_compliant: true
wave_0_complete: true
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
| 01-T1 | 01 | 1 | SC1 (field inventory, no code change) | T-02-01 | N/A | bd comment | `git status --porcelain src/` empty | ✅ | ⬜ pending |
| 01-T2 | 01 | 1 | SC1 (tracer: shared dispatch table + neutral struct) | T-02-01 | N/A | tdd/regression | `tests/testthat/test-sinkhorn*.R` + `python/leafblower/test_solver_parity.py -k sinkhorn` | ✅ | ⬜ pending |
| 01-T3 | 01 | 1 | SC1 (tracer rollback proof, D-02) | T-02-01 | N/A | manual + regression | `git revert` dry-run + full suite green | ✅ | ⬜ pending |
| 02-T1 | 02 | 1 | SC4 (build-list sync) | — | N/A | regression-prevention | `python/leafblower/test_core_sources_sync.py` (new, per D-05) | ❌ Wave 0 | ⬜ pending |
| 02-T2 | 02 | 1 | SC3 (water-fill single source; `newton_kl` gap ticketed, per D-03) | — | N/A | verification | grep/reflection assertion each solver `.cpp` calls `finalize_weights[_buf]` (new) | ❌ Wave 0 | ⬜ pending |
| 02-T3 | 02 | 1 | SC2 (opt-level symmetry, documented per CRAN `PKG_CXXFLAGS` constraint) | — | N/A | doc check | manual + DoD gate | ✅ (existing gate) | ⬜ pending |
| 03-T1..T3 | 03 | 2 | SC1 (greg/greenkhorn/logit migration) | T-02-01 | N/A | regression | per-solver R testthat + `test_solver_parity.py -k {greg,greenkhorn,logit}` | ✅ | ⬜ pending |
| 04-T1..T2 | 04 | 3 | SC1 (chebyshev/raking migration, warm-start divergence tabulated) | T-02-01 | N/A | regression | per-solver R testthat + `test_solver_parity.py -k {chebyshev,raking}` | ✅ | ⬜ pending |
| 05-T1..T2 | 05 | 4 | SC1 (oris/oris_soft migration; `sraa_demoted`, ALM/SRAA fields; `estimate_M_cell` dup measured) | T-02-01 | N/A | regression | R testthat + `test_solver_parity.py -k oris` + `test_parity_weights.py` | ✅ | ⬜ pending |
| 06-T1..T2 | 06 | 5 | SC1 (newton_kl migration; `n_projected_dims`, `lm_mu_final`, `stall_kind` survive) | T-02-01 | N/A | regression | `test-newton-kl.R`, `test-newton-tsvd-projection.R`, `test-newton-kl-tsvd-ratio.R`, `test-cr-d5-auto-fallback-fields.R`, `test_solver_parity.py -k newton_kl` | ✅ | ⬜ pending |
| 07-T1..T3 | 07 | 6 | SC1 (AUTO + chain collapse, `estimate_M_cell` consolidated, two `kAlgNames` tables reconciled) | T-02-01 | N/A | regression | full R testthat + full `python/leafblower/test_*.py` | ✅ | ⬜ pending |
| 08-T1..T3 | 08 | 7 | SC5 (full DoD gate + stepstone; empty/encoding probe rows closed; phase gate, `autonomous: false`) | T-02-01 | N/A | full-suite + benchmark | `.coverage-thresholds.json` `enforcement.command` + `LBW_BENCH_GATE=1` stepstone run | ✅ | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky — flips to ✅/⬜ per-task as `/gsd-execute-phase` runs each wave.*

---

## Wave 0 Requirements

- [x] `python/leafblower/test_core_sources_sync.py` (or R-side equivalent) — asserts `src/*.cpp` (minus `r_bridge.cpp`) == `python/CMakeLists.txt` `CORE_SOURCES`. Covers SC4. Scheduled as plan 02 task 1 (wave 1, parallel to the 01 tracer — not a literal pre-wave-1 gate, but has no dependency on the tracer and is safe to land first).
- [x] New SC3 verification test/script asserting all solvers route through `calib_dispatch.hpp::finalize_weights`/`finalize_weights_buf` (`newton_kl`'s gap explicitly ticketed, not silently absorbed, per RESEARCH.md Pitfall 2). Covers SC3. Scheduled as plan 02 task 2.
- [x] `leafblower-rywn`'s own Step 1 field inventory (a `bd comment`, not a test file) scheduled as plan 01 task 1, runs before the tracer task (01 task 2) which is the first solver migration.

---

## Manual-Only Verifications

*None — all phase behaviors have automated verification (existing R/Python parity suites + two new Wave 0 tests).*

---

## Validation Sign-Off

- [x] All tasks have `<automated>` verify or Wave 0 dependencies
- [x] Sampling continuity: no 3 consecutive tasks without automated verify
- [x] Wave 0 covers all MISSING references (test_core_sources_sync.py, SC3 verification test)
- [x] No watch-mode flags
- [x] Feedback latency < 180s
- [x] `nyquist_compliant: true` set in frontmatter

**Approval:** approved 2026-08-15 (gsd-plan-checker: 0 blockers, 8 plans verified against ROADMAP SC1-SC5 + US-004)
