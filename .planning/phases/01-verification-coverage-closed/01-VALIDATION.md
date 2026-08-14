---
phase: 1
slug: verification-coverage-closed
# status lifecycle: draft (seeded by plan-phase) → validated (set by validate-phase §6)
# audit-milestone §5.5 distinguishes NOT-VALIDATED (draft) from PARTIAL (validated + nyquist_compliant: false) (#2117)
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-08-15
---

# Phase 1 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | testthat >= 3.0.0 (R, installed 3.3.2; edition flip is SC5's deliverable) + pytest (Python, via `python/.venv/bin/python -m pytest`) |
| **Config file** | `DESCRIPTION` (R — `Config/testthat/edition`, currently absent); `python/pyproject.toml` (Python) |
| **Quick run command** | `Rscript -e 'testthat::test_file("tests/testthat/<new-file>.R")'` (R); `.venv/bin/python -m pytest python/leafblower/test_solver_parity.py -q` (Python) |
| **Full suite command** | `.coverage-thresholds.json`'s `enforcement.command` verbatim — `R CMD INSTALL --preclean .` + testthat + `uv pip install -e .` + `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -m pytest` |
| **Estimated runtime** | ~180 seconds (R build + full R/Python suites) |

---

## Sampling Rate

- **After every task commit:** Run the single new/changed file's quick command
- **After every plan wave:** Full `.coverage-thresholds.json` `enforcement.command`
- **Before `/gsd-verify-work`:** Full suite must be green; Python collected-test count must read **149** (141 pre-existing + 8 relocated from `leafblower-x7n8`), not 141
- **Max feedback latency:** 180 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 01-01 | 01 | 0 | leafblower-x7n8 | — / N/A | `git mv tests/test_parity_weights.py python/leafblower/test_parity_weights.py` with `REPO_ROOT` fix, 8 tests collected | integration (collection-count) | `cd python && .venv/bin/python -m pytest --collect-only -q \| tail -1` | ❌ Wave 0 | ⬜ pending |
| 01-02 | 01 | 1 | SC1 (9-solver parity matrix) | — / N/A | R↔Python weight-vector `rtol=1e-6` parity per method, `chebyshev`/`greg`/`oris_soft` added | integration (subprocess) | `.venv/bin/python -m pytest python/leafblower/test_parity_weights.py -q` | ❌ Wave 0 — extend | ⬜ pending |
| 01-03 | 01 | 1 | SC2 (raking/sinkhorn convergence-rule parity) | — / N/A | per-method default-rule resolution matches across bindings | integration (subprocess) | `.venv/bin/python -m pytest python/leafblower/test_solver_parity.py -q` | ❌ Wave 0 — new functions | ⬜ pending |
| 01-04 | 01 | 1 | SC3 (logit tolerance) | — / N/A | tolerance tightened to `1e-10` or comment names conditioning mechanism | integration + investigation | same file as 01-02, tolerance line only | ✅ existing | ⬜ pending |
| 01-05 | 01 | 1 | KPI-02 (weight-bound property test) | — / N/A | `min(w) >= min_weight` and `max(w) <= max_weight` within `1e-10`, 50 fixed-seed mixture-distributed datasets, `bounds_mode="unit"` | property-based | `Rscript -e 'testthat::test_file("tests/testthat/<new-file>.R")'` | ❌ Wave 0 — new file | ⬜ pending |
| 01-06 | 01 | 1 | SC5 (testthat edition) | — / N/A | full 94-file R suite green under 3e | full-suite regression | `.coverage-thresholds.json`'s R half | ✅ existing files, config-only | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `git mv tests/test_parity_weights.py python/leafblower/test_parity_weights.py` + `REPO_ROOT` path fix (`leafblower-x7n8`, must land first per D-09; silent-skip failure mode otherwise)
- [ ] New R file under `tests/testthat/` for the KPI-02 property test (name per planner — e.g. `test-bound-property.R`)
- [ ] No framework install needed — testthat and pytest both already present at required versions

---

## Manual-Only Verifications

*All phase behaviors have automated verification.*

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 180s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
