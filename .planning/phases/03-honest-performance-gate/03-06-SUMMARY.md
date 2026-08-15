---
phase: 03-honest-performance-gate
plan: 06
subsystem: testing
tags: [benchmark, R, chebyshev, optweight, minimax, LP]

# Dependency graph
requires:
  - phase: 03-honest-performance-gate (plans 01-05)
    provides: benchmarks/oris_soft_vs_competitors.R harness, medium_100k_5margins fixture, arm_row()/margin_max_error()/normalize_to_n()/lb_only_arm_row() helpers, run_honest_gate.sh wrapper
provides:
  - leafblower_chebyshev arm on medium_100k_5margins, reusing lb_only_arm_row() (same convergence=list() convention as oris/raking/newton_kl)
  - optweight_linf competitor arm (NEW dependency: optweight, benchmark-scoped only per D-09), paired with an honest objective-mismatch caveat and ok=NA (no max.w to verify against)
affects: [03-08-PLAN.md (docs/performance.md per-method table, consumes this CSV)]

actuals:
  tokens: 5150
  tasks: 2
  commits: 2

tech-stack:
  added: []
  patterns:
    - "optweight_linf follows the same requireNamespace()-guarded, tryCatch()-isolated competitor-row structure as icarus_calibration/ReGenesees_e_calibrate, but sets ok=NA (not TRUE/FALSE) when the competitor package cannot express the bound being measured -- a new but minimal extension of the existing arm-row convention for honestly reporting unverifiable claims"

key-files:
  created: []
  modified:
    - benchmarks/oris_soft_vs_competitors.R
    - benchmarks/results/oris_soft_vs_competitors.csv
    - benchmarks/results/oris_soft_vs_competitors_env.txt

key-decisions:
  - "leafblower_chebyshev reuses lb_only_arm_row('chebyshev', ...) verbatim -- no new helper needed, since chebyshev's convergence=list() natural-default call shape is identical to oris/raking/newton_kl's existing arm"
  - "optweight_linf's targets vector is built as unlist(tgt[margin_cols], use.names=FALSE) -- ALL categories per margin, not 'non-reference category only' as the plan's literal text specified. Verified empirically this session against the installed optweight package (process_targets() run on a fixture matching this shape, plus its own Rd: 'a target value must be specified for each level of the factor... must add up to 1') and against optweight's internal model.covs expansion (m1_a,m1_b,m1_c,m1_d,m2_a,... -- every level kept as a dummy column, none dropped as a reference). This is a [Rule 3] blocking-issue auto-fix: the plan's literal instruction would have produced a target vector of the wrong length (K*(nj-1)=15 instead of K*nj=20), which optweight would reject."
  - "optweight_linf's ok is set to NA (not TRUE), per plan instruction, since optweight.svy() has no max.w argument (confirmed against its own formals: formula, data, tols, targets, s.weights, b.weights, norm, min.w, verbose, ...) and there is no bound to certify compliance against"
  - "optweight_linf's max_w is reported on the same normalize_to_n()-scale every other competitor row in this file already uses, for comparability -- consistent with the file's established convention (competitor weights normalized to sum(w)==n; leafblower's own arms never are), rather than reporting a raw un-normalized value the plan's language left ambiguous"

requirements-completed: [US-003, KPI-04]

coverage:
  - id: D1
    description: "leafblower_chebyshev has a fresh row on medium_100k_5margins, compared against optweight's norm='linf' minimax-flavored objective"
    requirement: "US-003"
    verification:
      - kind: other
        ref: "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript benchmarks/oris_soft_vs_competitors.R && grep -c '\"leafblower_chebyshev\"\\|\"optweight_linf\"' benchmarks/results/oris_soft_vs_competitors.csv"
        status: pass
    human_judgment: false
  - id: D2
    description: "the objective mismatch between optweight's weight-deviation minimax and chebyshev's margin-error minimax is stated explicitly in the CSV note, not silently treated as a like-for-like replication; optweight's absence of a max.w argument is recorded as an honest gap (ok=NA), not silently assumed bound-compliant"
    requirement: "US-003"
    verification:
      - kind: other
        ref: "grep '\"optweight_linf\"' benchmarks/results/oris_soft_vs_competitors.csv | grep -c 'UNVERIFIABLE on this arm' -- 1 match; ok field literal value is NA"
        status: pass
    human_judgment: false
  - id: D3
    description: "CSV/env provenance regenerated fresh; DESCRIPTION diff empty (optweight stays benchmark-scoped, D-09); DoD unbroken (R CMD INSTALL, R+Python test suites 0 FAIL)"
    verification:
      - kind: other
        ref: "R CMD INSTALL --preclean . && devtools::test() (0 FAIL/141 WARN/13 SKIP/1837 PASS) && python -m pytest (161 passed) && git diff --stat DESCRIPTION | wc -l -> 0"
        status: pass
    human_judgment: false

duration: 55min
completed: 2026-08-16
status: complete
---

# Phase 03 Plan 06: chebyshev vs. optweight (norm="linf") Summary

**leafblower_chebyshev now has a doc-grounded, honestly-caveated competitor row against `optweight::optweight.svy(norm="linf")` on `medium_100k_5margins` -- closing the last R-side method in G-03-1's per-method coverage gap, and surfacing a genuine finding: optweight's LP formulation takes ~660s (median, single-threaded) to solve a 100,000-unit problem chebyshev solves in 0.028s.**

## Performance

- **Duration:** ~55 min (dominated by the optweight LP solve itself: 2 bench::mark iterations x ~660s each)
- **Started:** 2026-08-16 (approx)
- **Completed:** 2026-08-16
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments

- Installed `optweight` (v2.0.1) into the ambient R library, matching how `icarus`/`ReGenesees` were made available in 03-02 -- an ambient install, not a script-embedded auto-install.
- Added `leafblower_chebyshev` arm via the existing `lb_only_arm_row("chebyshev", "leafblower_chebyshev")` helper -- no new code needed, since chebyshev's `convergence=list()` call shape already matches the pattern established for oris/raking/newton_kl in 03-05.
- Added `optweight_linf` competitor arm: `optweight::optweight.svy(formula, data=df, targets=<all-category proportion vector>, tols=0, norm="linf", min.w=1e-8)`, `requireNamespace()`-guarded and `tryCatch()`-isolated per the file's established competitor-row structure. `ok` is `NA` (not `TRUE`), since `optweight.svy()` has no `max.w` argument and there is nothing to certify bound compliance against.
- The CSV `note` field records, verbatim, the objective mismatch specified in the plan: optweight's minimax is over *weight deviation from a base weight* (margin balance enforced via `tols=0` as a *constraint*), while chebyshev minimizes *max margin error directly* as the LP objective (weights unconstrained in the objective) -- a related but non-identical minimax formulation. The note also states optweight's `max.w`-argument absence explicitly.
- Regenerated `benchmarks/results/oris_soft_vs_competitors.csv` (13 rows on `medium_100k_5margins`, up from 11) and `_env.txt` (added `optweight: 2.0.1` to the provenance vector, same capture pattern as `survey`/`icarus`/`ReGenesees`).
- Confirmed the Definition of Done unbroken: `R CMD INSTALL --preclean .` succeeds; `devtools::test()` reports `0 FAIL/141 WARN/13 SKIP/1837 PASS` (matches the 03-05 baseline exactly); Python `pytest` `161 passed/0 failed`; `DESCRIPTION` diff empty and contains zero `optweight` mentions (benchmark-scoped only, per D-09).

## Task Commits

1. **Task 1: leafblower_chebyshev + optweight_linf arms** - `985a49e` (feat)
2. **Task 2: Regenerate CSV/env, confirm optweight stays benchmark-only, DoD check** - `a24e64b` (docs)

## Files Created/Modified

- `benchmarks/oris_soft_vs_competitors.R` - added the `leafblower_chebyshev` call and the `optweight_linf` competitor block (88 lines); extended `competitor_pkgs` with `"optweight"`
- `benchmarks/results/oris_soft_vs_competitors.csv` - regenerated; `medium_100k_5margins` grew from 11 to 13 rows
- `benchmarks/results/oris_soft_vs_competitors_env.txt` - regenerated; now records `optweight: 2.0.1`

## Measured Results (medium_100k_5margins, max_weight=3, new rows this plan)

| arm | wall_s (median) | max_error | max_w | ok | n_eff |
|---|---|---|---|---|---|
| leafblower_chebyshev | 0.028 | 3.65e-10 | 3.000 | TRUE | 67478.2 |
| optweight_linf | 660.388 | 2.01e-12 | 1.811 | **NA** | 60307.8 |

**Genuine finding worth flagging for the docs plan (03-08) that consumes this CSV:** `optweight_linf`'s LP solve (HiGHS backend) took ~660 seconds (median of 2 runs) to converge at n=100,000 -- roughly 24,000x slower than `leafblower_chebyshev`'s 0.028s on the identical fixture. This is not a fixture artifact: a smaller-scale sanity check (n=2,000, same K/nj/skew) confirmed `optweight.svy()` runs in ~1-2s, so the LP's cost is genuinely scaling badly with the per-observation weight-variable count as n grows into the hundreds of thousands. `optweight`'s own accuracy (`max_error=2.01e-12`) is excellent -- the finding is purely about runtime scaling of its LP reformulation at this problem size, an honest measured limitation reported as-is (consistent with D-08/this phase's honesty mandate), not a defect in leafblower's own solver.

## Decisions Made

See `key-decisions` in frontmatter: reuse of `lb_only_arm_row()` for chebyshev; the empirically-corrected (all-category, not non-reference-only) target vector construction; `ok=NA` for the unverifiable max.w bound; `max_w` reported on the file's established `normalize_to_n()` scale for cross-arm comparability.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Plan's literal "non-reference category" target-vector spec would have produced the wrong vector length**
- **Found during:** Task 1 implementation
- **Issue:** The plan's action text specified `targets = <margin-proportion vector, one entry per (margin, non-reference category), in the SAME order process_targets() would produce>`. Empirically running `optweight::process_targets()` against a fixture matching this file's medium-fixture shape (installed package this session, v2.0.1) showed the package requires a target value for **every** level of every factor margin (confirmed both by the package's own Rd: "a target value must be specified for each level of the factor, and these values must add up to 1", and by the internal `model.covs` column expansion: `m1_a, m1_b, m1_c, m1_d, m2_a, ...` -- no reference level is dropped). Passing a non-reference-only vector (`K*(nj-1)=15` entries instead of the required `K*nj=20`) would have errored or silently mis-mapped values to the wrong dummy columns.
- **Fix:** Built the target vector as `unlist(tgt[margin_cols], use.names = FALSE)` (all categories, all margins, in `margin_cols` order matching `df`'s factor level order) and verified it against a live `process_targets()` call on a matching fixture shape before writing the production call.
- **Files modified:** `benchmarks/oris_soft_vs_competitors.R`
- **Verification:** `optweight_linf` row's `max_error=2.006e-12` (machine-precision-adjacent) confirms the target vector correctly matched the achieved population totals; a mis-ordered or wrong-length vector would have failed the LP or produced a large `max_error`.
- **Committed in:** `985a49e` (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (Rule 3 -- blocking issue preventing literal task completion; the fix is a correction to match the actual, empirically-verified package contract, not a scope change)
**Impact on plan:** No scope creep. The corrected target-vector construction is what the plan's own stated intent ("in the SAME order `process_targets()` would produce") requires once checked against the real installed package -- the plan's prose description of "non-reference category" was simply inaccurate for this package's actual API.

## Issues Encountered

`optweight_linf`'s LP solve took ~11 minutes per `bench::mark` iteration at n=100,000 (22 minutes total for `iterations=2`), making this plan's total wall-clock time dominated by that single measurement. This was not a bug or a stuck process -- confirmed via direct process monitoring (steady ~98% CPU throughout) and cross-checked against a fast small-n sanity run. The slow runtime is itself the genuine finding documented above for 03-08.

## User Setup Required

None - no external service configuration required. `optweight` was installed into the ambient R library this session (D-09 benchmark-scoped install, same treatment as `icarus`/`ReGenesees`); no new package added to `DESCRIPTION`.

## Next Phase Readiness

- G-03-1's R-reachable method coverage is now fully closed: `oris`, `oris_soft`, `raking`, `newton_kl`, `greg`, `logit`, and `chebyshev` all have a measured, doc-grounded competitor row on `medium_100k_5margins`.
- Remaining G-03-1 scope: `greenkhorn`/`sinkhorn` vs. Python `POT` (cross-language, 03-07-PLAN.md).
- `benchmarks/results/oris_soft_vs_competitors.csv` now carries the `chebyshev`/`optweight_linf` rows 03-08-PLAN.md's `docs/performance.md` per-method table needs, including the honest objective-mismatch caveat and the optweight LP-scaling finding.
- No blockers.

---
*Phase: 03-honest-performance-gate*
*Completed: 2026-08-16*

## Self-Check: PASSED

All claimed files found on disk; all claimed commit hashes found in `git log --all`.
