---
phase: 4
slug: truthful-surface
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-08-15
---

# Phase 4 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | testthat 3 (`Config/testthat/edition: 3` in DESCRIPTION) |
| **Config file** | `tests/testthat.R` (standard harness, not modified by this phase) |
| **Quick run command** | `Rscript -e "testthat::test_dir('tests/testthat', filter='harvest-rval')"` (basename filter, project convention) |
| **Full suite command** | `Rscript -e "devtools::test()"` |
| **Estimated runtime** | ~30-90s quick filter; full suite per project baseline |

---

## Sampling Rate

- **After every task commit:** relevant filtered `testthat::test_dir(..., filter=...)` run
  (`harvest-rval`, `logit`) for code-touching tasks; SC1/SC4 verification greps for doc-only tasks.
- **After every plan wave:** full `Rscript -e "devtools::test()"` — run `R CMD INSTALL --preclean .`
  first if `src/leafblower.h` was touched (SC2 forces a full rebuild of both R and Python
  extensions per CLAUDE.md's "Two build sites" note).
- **Before `/gsd-verify-work`:** Full DoD — R tests 0 FAIL + Python parity tests 0 FAIL
  (`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`). Stepstone benchmark
  (`LBW_BENCH_GATE=1`) not required — no task touches a hot path.
- **Max feedback latency:** ~120s (full suite).

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|--------------------|-------------|--------|
| 04-01-01 | 01 | 0 | SC3 | V5 | RVAL.4 test exists asserting `design_weights` error message | unit | `Rscript -e "testthat::test_dir('tests/testthat', filter='harvest-rval')"` | ❌ W0 | ⬜ pending |
| 04-01-02 | 01 | 1 | SC1 | — | docs/raking.md §8.2/§12 attribute correctly | manual | `grep -o 'L-BFGS-B' docs/raking.md \| wc -l` = 0; `grep -ci Gurobi docs/raking.md` = 0 | ✅ | ⬜ pending |
| 04-01-03 | 01 | 1 | SC2 | — | `rk_algorithm_t` slot 7 documented as removed grake slot | build | `git diff src/leafblower.h` comment-only; full suite 0 FAIL | ✅ | ⬜ pending |
| 04-01-04 | 01 | 1 | SC3 | V5 | `harvest(weights=)` errors naming `design_weights=`; `eb79.18` renamed | unit | `Rscript -e "testthat::test_dir('tests/testthat', filter='harvest-rval')"` + `filter='logit'` | ✅ (edit) | ⬜ pending |
| 04-01-05 | 01 | 1 | SC4 | — | README/NEWS/man/docs audit — zero true positives outside raking.md | manual | Re-run session's grake/lbfgsb/cp greps | N/A | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `tests/testthat/test-harvest-rval.R` — add `test_that("RVAL.4: bare weights= errors naming design_weights", ...)` per RVAL.1-3/META.2 numbering convention.

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|--------------------|
| docs/raking.md no longer misattributes outer entropic-proximal-point loop | SC1 | Prose fix, no code path to assert against | `grep -o 'L-BFGS-B' docs/raking.md \| wc -l` = 0; `grep -ci Gurobi docs/raking.md` = 0; visually confirm §8.2/§12 attribute to the paper's own name or are removed |
| README/NEWS/man/docs free of grake/lbfgsb/cp-as-method references | SC4 | Prose/doc audit, not a code path | Re-run researcher's greps for `grake`, `lbfgsb`, `\bcp\b` across `README*`, `NEWS.md`, `man/*.Rd`, `docs/*.md`; every hit must match the researcher's false-positive table (grake_norm, survey::grake, Gurobi/CPLEX comparison, withdrawn-spike cp path, removal-assertion test) |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references (RVAL.4)
- [ ] No watch-mode flags
- [ ] Feedback latency < 120s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
