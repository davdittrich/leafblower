---
phase: 03-honest-performance-gate
plan: 05
subsystem: testing
tags: [benchmark, R, survey-calibrate, oris, raking, newton_kl, greg, logit]

# Dependency graph
requires:
  - phase: 03-honest-performance-gate (plans 01-04)
    provides: benchmarks/oris_soft_vs_competitors.R harness, medium_100k_5margins fixture, arm_row()/margin_max_error()/normalize_to_n() helpers, run_honest_gate.sh wrapper
provides:
  - leafblower_oris, leafblower_raking, leafblower_newton_kl arms on medium_100k_5margins, reusing the fixture's existing survey_calibrate/icarus_calibration/ReGenesees_e_calibrate competitor rows
  - leafblower_greg, leafblower_logit arms paired with NEW distance-matched survey_calibrate_linear/survey_calibrate_logit competitor rows
affects: [03-06-PLAN.md (chebyshev/optweight), 03-07-PLAN.md (greenkhorn/sinkhorn/POT), 03-08-PLAN.md (docs/performance.md per-method table, consumes this CSV)]

actuals:
  tokens: 3592
  tasks: 3
  commits: 3

tech-stack:
  added: []
  patterns:
    - "lb_only_arm_row(method_name, arm_label): leafblower-only arm helper reusing an already-built fixture (df/tgt/n/K/nj/m_cell_medium) and an already-computed competitor set, for methods sharing an existing competitor's objective"
    - "survey_calfun_arm_row(calfun_name, arm_label): tryCatch-isolated survey::calibrate(calfun=) competitor helper, parameterized by distance function, for methods needing a NEW objective-matched competitor row"

key-files:
  created: []
  modified:
    - benchmarks/oris_soft_vs_competitors.R
    - benchmarks/results/oris_soft_vs_competitors.csv

key-decisions:
  - "oris/raking/newton_kl reuse the medium fixture's existing survey_calibrate/icarus_calibration/ReGenesees_e_calibrate rows (D-08 narrow-not-duplicate) rather than triggering three more competitor solves, since all three attack the identical K-margin box-bounded KL objective per docs/methods/oris.md, raking.md, newton_kl.md's own Practitioner-implementations tables"
  - "greg gets calfun='linear' (chi-square/linear objective match, docs/methods/greg.md), logit gets calfun='logit' (docs/methods/logit.md names survey::calibrate(calfun='logit') first) -- distinct new competitor rows survey_calibrate_linear/survey_calibrate_logit, never the existing raking-calfun arm, which would be an objective mismatch"
  - "greg's measured max_error=0.109 on this skewed/tight-bound fixture (with the package's own runtime warning naming K=5/max_weight=3 as cause) is left as-is -- an honest measured limitation, not a bug to hide behind a friendlier fixture"

requirements-completed: [US-003, KPI-04]

coverage:
  - id: D1
    description: "leafblower_oris/leafblower_raking/leafblower_newton_kl each have a fresh row on medium_100k_5margins, reusing the fixture's existing survey/icarus/ReGenesees competitor rows"
    requirement: "US-003"
    verification:
      - kind: other
        ref: "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript benchmarks/oris_soft_vs_competitors.R && grep -c '\"leafblower_oris\",\\|\"leafblower_raking\",\\|\"leafblower_newton_kl\",' benchmarks/results/oris_soft_vs_competitors.csv"
        status: pass
    human_judgment: false
  - id: D2
    description: "leafblower_greg/leafblower_logit each have a fresh row on medium_100k_5margins, paired with a NEW distance-matched survey_calibrate_linear/survey_calibrate_logit competitor row (never the raking-calfun arm)"
    requirement: "US-003"
    verification:
      - kind: other
        ref: "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript benchmarks/oris_soft_vs_competitors.R && grep -c '\"leafblower_greg\",\\|\"leafblower_logit\",\\|\"survey_calibrate_linear\",\\|\"survey_calibrate_logit\",' benchmarks/results/oris_soft_vs_competitors.csv"
        status: pass
    human_judgment: false
  - id: D3
    description: "CSV/env provenance regenerated fresh; DoD unbroken (R CMD INSTALL, R+Python test suites 0 FAIL, DESCRIPTION diff empty)"
    verification:
      - kind: other
        ref: "R CMD INSTALL --preclean . && devtools::test() (0 FAIL/141 WARN/13 SKIP/1837 PASS) && python -m pytest (161 passed) && git diff --stat DESCRIPTION"
        status: pass
    human_judgment: false

duration: 25min
completed: 2026-08-16
status: complete
---

# Phase 03 Plan 05: Per-Method Competitor Coverage (oris/raking/newton_kl/greg/logit) Summary

**medium_100k_5margins now carries 11 measured rows: 5 leafblower methods (oris_soft, oris, raking, newton_kl, greg, logit) each graded against a doc-grounded `survey::calibrate` competitor on the same wall_s/max_error/bound/n_eff metrics — closing the G-03-1 gap for the zero-new-dependency slice of leafblower's method surface.**

## Performance

- **Duration:** ~25 min
- **Started:** 2026-08-16T00:00:00Z (approx)
- **Completed:** 2026-08-16
- **Tasks:** 3
- **Files modified:** 2

## Accomplishments

- Added `lb_only_arm_row()` — reuses the medium fixture's already-measured `survey_calibrate`/`icarus_calibration`/`ReGenesees_e_calibrate` rows for `leafblower_oris`, `leafblower_raking`, `leafblower_newton_kl` (all three attack the identical K-margin box-bounded KL objective per `docs/methods/oris.md`, `raking.md`, `newton_kl.md`'s own Practitioner-implementations tables). Each call uses `convergence = list()` (per-method natural default, R/harvest.R:424), never the borrowed `marginal_kl` pin used by `oris_soft`.
- Added `survey_calfun_arm_row()` — a `tryCatch()`-isolated, `calfun`-parameterized `survey::calibrate` competitor row, paired with `leafblower_greg` (`calfun="linear"`, chi-square/linear-distance match) and `leafblower_logit` (`calfun="logit"`, the implementation `docs/methods/logit.md` names first). Two new competitor rows (`survey_calibrate_linear`, `survey_calibrate_logit`) added; the existing `survey_calibrate` (implicit `calfun="raking"`) row is untouched.
- Regenerated `benchmarks/results/oris_soft_vs_competitors.csv` and confirmed `_env.txt` byte-identical (same machine/package versions). `medium_100k_5margins` now has 11 rows (up from 4); `large_stepstone_fulldata`/`known_limit_k20_uniform` classes untouched.
- Confirmed the Definition of Done unbroken: `R CMD INSTALL --preclean .` succeeds; `devtools::test()` reports `0 FAIL/141 WARN/13 SKIP/1837 PASS` (matches the 03-04 baseline exactly); Python pytest `161 passed/0 failed`; `DESCRIPTION` diff empty (`survey` was already `Suggests:`, no new package).

## Task Commits

1. **Task 1: leafblower_oris / leafblower_raking / leafblower_newton_kl arms, reusing existing competitor rows** - `3d57fef` (feat)
2. **Task 2: leafblower_greg / leafblower_logit arms with distance-matched NEW competitor rows** - `f10c98c` (feat)
3. **Task 3: Regenerate CSV/env provenance, verify Definition of Done unbroken** - `eb25cfe` (docs)

## Files Created/Modified

- `benchmarks/oris_soft_vs_competitors.R` - added `lb_only_arm_row()` and `survey_calfun_arm_row()` helpers plus 7 new arm-row calls on the medium fixture (104 lines added)
- `benchmarks/results/oris_soft_vs_competitors.csv` - regenerated; `medium_100k_5margins` grew from 4 to 11 rows

## Measured Results (medium_100k_5margins, max_weight=3)

| arm | wall_s (median) | max_error | max_w | n_eff |
|---|---|---|---|---|
| leafblower_oris_soft | 0.042 | 3.35e-05 | 3.000 | 67489 |
| survey_calibrate (calfun=raking) | 0.499 | 1.49e-04 | 2.998 | 67434 |
| icarus_calibration | 0.426 | 6.20e-09 | 2.873 | 69075 |
| ReGenesees_e_calibrate | 0.905 | 1.10e-08 | 3.000 | 67478 |
| **leafblower_oris** (new) | 0.090 | 7.59e-11 | 3.000 | 67478 |
| **leafblower_raking** (new) | 0.025 | 5.00e-16 | 3.000 | 67523 |
| **leafblower_newton_kl** (new) | 0.051 | 4.22e-09 | 6.210 | 66189 |
| **leafblower_greg** (new) | 0.023 | 1.09e-01 | 1.064 | 99924 |
| **leafblower_logit** (new) | 0.025 | 2.19e-08 | 2.873 | 69075 |
| **survey_calibrate_linear** (new) | 0.229 | 5.55e-17 | 3.000 | 69704 |
| **survey_calibrate_logit** (new) | 0.292 | 1.17e-15 | 2.873 | 69075 |

Two genuine findings worth flagging for the docs plan (03-08) that consumes this CSV:
- `leafblower_newton_kl` measures `max_w=6.210` (above the `max_weight=3` bound) — matches its documented shipped contract (`leafblower-73d7`, cited in PROJECT.md): newton_kl *reports* bound violations (`RK_ERR_NOCONV`) rather than clamping, unlike the other 8 solvers.
- `leafblower_greg` measures `max_error=0.109` on this skewed/tight-bound fixture. The package itself emits a runtime warning at solve time: `"greg converged but max_err=0.1089 (>5%). greg may be unreliable for K=5 margins or tight bounds (max_weight=3)."` This is GREG's well-known linear-approximation degradation under tight bounds on strongly skewed targets — an honest, measured limitation of the method, not a bug. No fixture change was made to hide it (consistent with D-08/this phase's honesty mandate).

## Decisions Made

- Reuse (not recompute) the three existing competitor rows for oris/raking/newton_kl — see key-decisions in frontmatter.
- `calfun="linear"` for greg, `calfun="logit"` for logit — distance-matched per each method's own doc, never the raking-calfun arm.
- greg's high measured error is reported as-is, cross-referenced for the downstream docs plan rather than silently patched.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Task 1/2's literal `<verify>` grep pattern didn't account for CSV quoting**
- **Found during:** Task 1 verification
- **Issue:** The plan's literal verify command (`grep -c "leafblower_oris,\|leafblower_raking,\|leafblower_newton_kl,"`) returns 0 because the CSV quotes string fields (`"leafblower_oris",` not `leafblower_oris,`) — the pattern's comma never directly follows the bare method name.
- **Fix:** Verified with a quote-aware pattern (`grep -c '"leafblower_oris",\|...'`) instead, confirming the actual row presence the `<done>` criterion requires. No code change; verification-methodology fix only.
- **Verification:** `grep -c '"leafblower_oris",\|"leafblower_raking",\|"leafblower_newton_kl",' benchmarks/results/oris_soft_vs_competitors.csv` → 3; same pattern applied for Task 2's four new rows → 4.
- **Committed in:** n/a (verification-only, no file change)

**2. [Rule 3 - Blocking] `design`/`formula_list`/`population_list` for the new competitor rows rebuilt, not literally reused**
- **Found during:** Task 2 implementation
- **Issue:** The plan instructed reusing the medium fixture's already-built `design`/`formula_list`/`population_list` objects "from the earlier survey_calibrate block" — but those are local variables inside `run_input_class()`'s function scope, unreachable from the top-level script scope where Task 2's new arms are added.
- **Fix:** Rebuilt `design_medium`/`formula_list_medium`/`population_list_medium` at top-level scope, using the identical construction (same `margin_cols`, `tgt`, `n`, `df`) as the original block — this is cheap, deterministic object declaration, not a re-solve of the calibration itself, so it does not violate D-08's "measure fresh, don't duplicate the expensive part" intent.
- **Files modified:** `benchmarks/oris_soft_vs_competitors.R`
- **Verification:** `survey_calibrate_linear`/`survey_calibrate_logit` rows produced correct control totals matching `tgt` (`max_error` at machine precision — 5.55e-17 and 1.17e-15 respectively — confirms identical population totals to the original `survey_calibrate` arm's construction).
- **Committed in:** `f10c98c` (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (both Rule 3 — blocking issues preventing literal task completion, neither changing measured behavior)
**Impact on plan:** No scope creep. Deviation 1 is verification-methodology only. Deviation 2 reconstructs identical declarative objects from the same already-computed fixture inputs; it does not re-run any solver or recompute an expensive result.

## Issues Encountered

None beyond the two deviations above.

## User Setup Required

None - no external service configuration required. `survey` was already an installed `Suggests:` dependency before this plan (D-09); no new package added.

## Next Phase Readiness

- G-03-1's zero-new-dependency slice (oris, raking, newton_kl, greg, logit — all mappable to `survey::calibrate`) is closed. Remaining G-03-1 scope: chebyshev (optweight, 03-06-PLAN.md) and greenkhorn/sinkhorn (POT, cross-language, 03-07-PLAN.md).
- `benchmarks/results/oris_soft_vs_competitors.csv` now has the full per-method row set 03-08-PLAN.md's `docs/performance.md` competitors table needs for the oris-family/raking/greg/logit rows; chebyshev/greenkhorn/sinkhorn rows still pending from 03-06/03-07.
- No blockers.

---
*Phase: 03-honest-performance-gate*
*Completed: 2026-08-16*

## Self-Check: PASSED

All claimed files found on disk; all claimed commit hashes found in `git log --all`.
