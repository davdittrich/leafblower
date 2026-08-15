---
phase: 03-honest-performance-gate
verified: 2026-08-15T16:20:00Z
status: passed
score: 4/4 must-haves verified
behavior_unverified: 0
overrides_applied: 0
---

# Phase 3: Honest Performance Gate Verification Report

**Phase Goal:** The package's headline performance claim is a measured fact about a solver
that exists, expressed as a gate a user can run and a maintainer can regress against.
**Verified:** 2026-08-15
**Status:** passed
**Re-verification:** No — initial verification.

## Goal Achievement

### Observable Truths (Success Criteria from ROADMAP/PLAN)

| # | Truth | Status | Evidence |
|---|---|---|---|
| 1 | A user reading package documentation sees a large-scale performance figure measured on this codebase with a live solver, on a stated input class and machine — not an lbfgsb-inherited target | ✓ VERIFIED | `README.md` states `oris_soft` calibrates in 0.0427s/max_error 3.35e-05/n_eff 67,489.4 on the medium class, linking to `docs/performance.md`. That page's `large_stepstone_fulldata` section states `n=1,582,732`, `wall_s=3.5459`, machine "AMD Ryzen 9 9950X3D 16-Core Processor" and R/BLAS/LAPACK versions, all transcribed from `benchmarks/results/oris_soft_vs_competitors_env.txt` (verified byte-matching by direct read). No `lbfgsb`/`test-lbfgsb.R` reference survives as a live claim anywhere in `README.md` or `docs/performance.md` (grep negative). |
| 2 | The medium-scale target states ONE number, not the 1s/2s contradiction, and names the measuring artefact | ✓ VERIFIED | `tasks/prd-leafblower-core.md` lines 30, 181, 673 (the three `100K rows, 5 margins` sites, confirmed via `grep -n -E`) each now carry an in-place `**Superseded 2026-08-15**` marker pointing to `docs/performance.md`; none was deleted (PRD is a frozen historical doc). `README.md`/`docs/performance.md` publish exactly one medium-scale figure (wall_s=0.0427), naming `benchmarks/oris_soft_vs_competitors.R` and the `honest gate:` testthat block as the measuring artefacts. |
| 3 | The K=20 uniform-random / M_cell/n=1.0 class is documented as a known limit with the structural reason, not a silently failing promise | ✓ VERIFIED | `docs/performance.md` § "Known limit" states `m_cell/n = 1.0000` (measured, confirmed against the live CSV), explains zero cell-compression benefit, gives both `leafblower_oris_soft` (max_error=5.229e-03) and `leafblower_raking_accelerated` (max_error=1.516e-03) figures, states plainly the composite `<30s AND <1e-6` gate is confirmed structurally unachievable, and cross-links to both investigation docs without restating them. `.planning/REQUIREMENTS.md`'s US-003/KPI-04 entries restate this as a documented known limit, not an open blocker. No gate assertion exists for this class (`grep -c 'known_limit_k20_uniform' tests/testthat/test-bench-gate.R` = 0, confirmed). |
| 4 | A maintainer can re-run the benchmark that produced every published figure with one command; KPI rows name a live measuring artefact | ✓ VERIFIED | `CI=1 bash benchmarks/run_honest_gate.sh` executed live during this verification: exit 0, regenerated `benchmarks/results/oris_soft_vs_competitors.csv`/`_env.txt` with figures matching the published ones within run-to-run wall-time noise (medium wall_s 0.0427→0.0427, large 3.5459→3.5321, known-limit 7.3934→7.3687). `.planning/REQUIREMENTS.md` US-003/KPI-04 name `benchmarks/oris_soft_vs_competitors.R` and the `honest gate:` assertion in `tests/testthat/test-bench-gate.R` (both confirmed present and passing), replacing the void `test-lbfgsb.R`. |

**Score:** 4/4 truths verified (0 present, behavior-unverified)

### Required Artifacts

| Artifact | Expected | Status | Details |
|---|---|---|---|
| `benchmarks/oris_soft_vs_competitors.R` | Fresh measurement script, 3 input classes, 4+ arms | ✓ VERIFIED | Exists, substantive (single-thread-BLAS guard confirmed to exit non-zero when `OPENBLAS_NUM_THREADS` unset — ran live), no `autumn` mention, `tryCatch` x5, `requireNamespace` x7. |
| `benchmarks/run_honest_gate.sh` | One-command wrapper | ✓ VERIFIED | Executable, ran live to completion (exit 0), regenerates both result artefacts, exports all three thread vars + `LBW_BENCH_GATE` + `NOT_CRAN`. |
| `benchmarks/results/oris_soft_vs_competitors.csv` | 17-column, 3 input classes, 10 rows | ✓ VERIFIED | Present on disk, schema matches, content matches every number quoted in `docs/performance.md` and `README.md`. |
| `benchmarks/results/oris_soft_vs_competitors_env.txt` | Machine/version provenance | ✓ VERIFIED | Present, content matches `docs/performance.md`'s "Machine and versions" table verbatim. |
| `tests/testthat/test-bench-gate.R` | Opt-in honest gate + re-gated kk1204 block | ✓ VERIFIED | Ran live: default filtered run (`LBW_BENCH_GATE` unset) completes in 0.34s reporting 4 SKIP, 0 PASS (no heavy solve paid). Gated run (`LBW_BENCH_GATE=1 NOT_CRAN=true CI=1`) produces `honest gate: wall_s=0.0427 max_error=3.347e-05 max_w=3.0000 n_eff=67489.4` and `kk1204 gate: status=0 iters=10 best_error=-7.376e-14 time=1.6s`, 0 FAIL / 10 PASS. |
| `docs/performance.md` | Methodology page, all sections | ✓ VERIFIED | 215 lines, contains headline claim, results tables (every wall_s row carries max_error/max_w/n_eff), input classes, machine/versions, methodology, reproduction, known limit (cross-linked, not restated), competitors (linked to `docs/methods/oris.md`, no restated bibliography). No `uuid`/`/home/dd/stepstone` disclosure (grep negative). |
| `README.md` | One-line headline claim + link | ✓ VERIFIED | 12 lines, states `oris_soft` paired claim (wall_s/max_error/bound/n_eff), links to `docs/performance.md`, no install/badge/CRAN content, no `autumn` mention. |
| `.planning/REQUIREMENTS.md` | US-003/KPI-04 name live artefacts | ✓ VERIFIED | Both entries rewritten, name `benchmarks/oris_soft_vs_competitors.R` and the `honest gate:` assertion; Traceability rows show `Partial` (honestly, not rounded up — see Notable Honesty below). |
| `tasks/prd-leafblower-core.md` | Contradiction sites marked superseded | ✓ VERIFIED | All 3 sites (lines 30, 181, 673) carry in-place superseded markers; `test-lbfgsb.R` explicitly annotated as void rather than silently left standing. |

### Key Link Verification

| From | To | Via | Status | Details |
|---|---|---|---|---|
| `run_honest_gate.sh` | `oris_soft_vs_competitors.R` | invocation | ✓ WIRED | Confirmed by live execution — script step 2 runs it, produces fresh CSV. |
| `oris_soft_vs_competitors.R` | `test-bench-gate.R`'s `honest gate:` block | CSV read | ✓ WIRED | `csv_path` anchored on `testthat::test_path()`; live gated run reads the same CSV the script writes and asserts on it. |
| `harvest(bounds_mode="unit")` | `survey::calibrate(bounds=...)` | comparable problem | ✓ WIRED | CSV confirms both arms received `[0,3]` per-observation bound; `icarus`/`ReGenesees` arms use their own bound-honouring methods (`method='logit'`, `calfun='raking'`). |
| `docs/performance.md` | `docs/methods/oris.md`, investigation docs | citation links | ✓ WIRED | Confirmed present in the page body (grep + read); bibliography not duplicated. |
| `.planning/REQUIREMENTS.md` KPI-04 | live artefacts | naming | ✓ WIRED | `benchmarks/oris_soft_vs_competitors.R` and `test-bench-gate.R` named and confirmed to exist and run. |

### Behavioral Spot-Checks (executed live during this verification, not from SUMMARY claims)

| Behavior | Command | Result | Status |
|---|---|---|---|
| R build gate | `R CMD INSTALL --preclean .` | `* DONE (leafblower)` | ✓ PASS |
| R test suite (DoD) | `Rscript -e "devtools::test()"` (single-thread BLAS) | `[ FAIL 0 \| WARN 141 \| SKIP 13 \| PASS 1836 ]` — matches 03-04-SUMMARY's claimed baseline exactly | ✓ PASS |
| Python parity suite | `.venv/bin/python -m pytest` (single-thread BLAS) | `160 passed, 0 failed` | ✓ PASS |
| Default bench-gate filtered run | `testthat::test_dir(filter='bench-gate')`, `LBW_BENCH_GATE` unset | 0.34s, `[ FAIL 0 \| WARN 0 \| SKIP 4 \| PASS 0 ]` — no heavy solve paid by default | ✓ PASS |
| Gated bench-gate run | `CI=1 LBW_BENCH_GATE=1 NOT_CRAN=true testthat::test_dir(filter='bench-gate', stop_on_failure=TRUE)` | `honest gate:`/`kk1204 gate:` lines printed, `[ FAIL 0 \| WARN 0 \| SKIP 2 \| PASS 10 ]` | ✓ PASS |
| One-command reproduction | `CI=1 bash benchmarks/run_honest_gate.sh` | exit 0, regenerated CSV with figures matching published ones within noise | ✓ PASS |
| Determinism guard | `OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript benchmarks/oris_soft_vs_competitors.R` (OPENBLAS unset) | exit 1, `Error: refusing to measure: OPENBLAS_NUM_THREADS must be set to "1"...` | ✓ PASS |
| D-05 negative check | `grep -v '^\s*#' ... \| grep -ci autumn` on `benchmarks/oris_soft_vs_competitors.R` | 0 | ✓ PASS |
| DESCRIPTION/pyproject.toml scoping | `grep -i 'icarus\|regenesees' DESCRIPTION python/pyproject.toml` | no matches | ✓ PASS |
| Task commits present | `git log --oneline` for all 12 task commits across 03-01..03-04 | all present, matches every SUMMARY's claimed hash | ✓ PASS |

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|---|---|---|---|---|
| US-003 | 03-01..03-04 | 1M+ obs / 20+ margins < 30s | ✓ HONESTLY PARTIAL | Rewritten entry names live artefacts, states measured figures on 3 classes, and explicitly states no single fixture combines "1M+ rows AND 20+ margins" as the PRD literally demands — status kept `Partial`, not rounded up to Implemented. This is correct, not a gap: the PRD's literal AND-combination genuinely was never demonstrated by one measurement, and the entry says so plainly. |
| KPI-04 | 03-01..03-04 | Large-scale 1M rows/20 margins/<30s | ✓ HONESTLY PARTIAL | Same honest-partial treatment; names both live artefacts; large-scale (1.58M rows, single-digit margins) and the 20-margin known-limit shape are measured separately and reported as such, not conflated. |

### Anti-Patterns Found

None blocking. No `TBD`/`FIXME`/`XXX` markers found in phase-modified files. The one deliberately-left inconsistency — the kk1204 block's fixture uses a uniform target (converges trivially, `best_error≈-7e-14`) rather than the skewed target that reproduces the actual degenerate known limit — is explicitly flagged in-code (lines 78-96 of `test-bench-gate.R`) and in `03-04-SUMMARY.md`'s "Fixture Mismatch Finding" section as a reported-not-fixed, out-of-scope-for-this-task observation, not hidden. This is the honest disclosure pattern the phase goal demands, not a defect.

### Human Verification Required

None. All four success criteria are independently verifiable via file content and live command execution, which was done during this verification (not merely re-stated from SUMMARY.md).

### Notable Honesty Checks (per verification_focus)

1. **kk1204 fixture-mismatch note (03-04):** Confirmed flagged both in-comment (`test-bench-gate.R:78-96`) and in `03-04-SUMMARY.md`'s dedicated "Fixture Mismatch Finding" section — not hidden. The block is correctly relabelled as a "shape-matched regression floor," not a test of the documented degenerate case.
2. **REQUIREMENTS.md US-003/KPI-04 left "Partial"/"Open→Partial" (03-04):** Confirmed via direct read of `.planning/REQUIREMENTS.md` lines 102-160 — both entries explicitly state the PRD's literal "1M+ rows AND 20+ margins" combination was never demonstrated by one single measurement, and status markers (`- [ ]`, `Partial`) were kept rather than rounded up to Implemented/complete. `requirements.mark-complete` reportedly no-op'd (per 03-04-SUMMARY, not independently re-run here since it is a planning-tool no-op already evidenced by the unchanged status markers in the file).

### Gaps Summary

None. All four ROADMAP success criteria are independently verified against the live codebase (file content plus executed commands), not merely SUMMARY.md claims. The R build, R test suite (0 FAIL), and Python parity suite (0 failed) all pass under the mandated single-thread BLAS envelope. The one-command reproduction (`benchmarks/run_honest_gate.sh`) was executed live and confirmed to regenerate every published figure. The determinism guard was executed live and confirmed to refuse measurement under a partially-set thread-variable environment.

---

*Verified: 2026-08-15*
*Verifier: Claude (gsd-verifier)*
