---
phase: 03-honest-performance-gate
plan: 08
subsystem: docs
tags: [docs, performance, benchmarks, competitors, gap-closure]

# Dependency graph
requires:
  - phase: 03-honest-performance-gate (plans 05-07)
    provides: "benchmarks/results/oris_soft_vs_competitors.csv (13 rows, incl. leafblower_chebyshev/optweight_linf) and benchmarks/results/greenkhorn_sinkhorn_vs_pot.csv (4 rows, k2_margin_pot_equiv class)"
provides:
  - "docs/performance.md's Competitors section rewritten into a 9-method table (oris, oris_soft, raking, sinkhorn, greenkhorn, chebyshev, greg, logit, newton_kl), each row citing its docs/methods/*.md practitioner-implementations source"
  - "docs/performance.md Results section gains two subsections transcribing every new CSV row verbatim: medium_100k_5margins per-method competitor coverage, and k2_margin_pot_equiv greenkhorn/sinkhorn vs POT"
affects: []

actuals:
  tokens: 2900
  tasks: 2
  commits: 1

tech-stack:
  added: []
  patterns:
    - "Competitors table columns (Method | Closest competitor | Source | Objective match | Caveats) transcribe CSV note fields verbatim rather than paraphrasing, so a reader can diff the published sentence against the CSV note to confirm nothing was softened in transit (T-03-08-02 mitigation)"

key-files:
  created: []
  modified:
    - docs/performance.md

key-decisions:
  - "greg/logit rows point to their own survey::calibrate(calfun=...) competitor rather than re-citing the oris-family trio, since GREG and logit are distance-matched arms (chi-square / logit distance), not raking-calfun matches -- keeps the Objective match column honest per method, not blanket-copied from the headline claim"
  - "sinkhorn/greenkhorn's Objective match column states 'equivalent at K=2 only', not a blanket 'same' -- the K>2 case has no POT-side counterpart at all, so overstating equivalence here would misrepresent the actual coverage"

requirements-completed: [US-003, KPI-04]

coverage:
  - id: D1
    description: "docs/performance.md's Competitors section covers all 9 leafblower methods, each against its own doc-named closest competitor, with the source docs/methods/*.md citation named per row"
    requirement: "US-003"
    verification:
      - kind: other
        ref: "for m in oris oris_soft raking sinkhorn greenkhorn chebyshev greg logit newton_kl; do grep -qi \"$m\" docs/performance.md || echo MISSING; done -- zero MISSING lines"
        status: pass
    human_judgment: false
  - id: D2
    description: "every new figure in the two new Results subsections is transcribed verbatim from the regenerated CSVs, never estimated -- spot-checked 3 numbers (leafblower_chebyshev row, optweight_linf row, pot_greenkhorn row) against the actual CSV cells"
    requirement: "US-003"
    verification:
      - kind: other
        ref: "manual diff of docs/performance.md's chebyshev/optweight_linf/pot_greenkhorn table cells against benchmarks/results/oris_soft_vs_competitors.csv and benchmarks/results/greenkhorn_sinkhorn_vs_pot.csv -- exact match on all 3"
        status: pass
    human_judgment: false
  - id: D3
    description: "both honest scope caveats (optweight's objective mismatch + no max.w; POT's K=2-only scope + no bounds) are stated on the published page, not left only in CSV note fields"
    requirement: "KPI-04"
    verification:
      - kind: other
        ref: "grep -ci optweight docs/performance.md -> 3; grep -ci 'POT\\|Python Optimal Transport' docs/performance.md -> 11"
        status: pass
    human_judgment: false
  - id: D4
    description: "full project DoD re-verified green across all four gap-closure plans (05-08) combined: R CMD INSTALL, R testthat (0 FAIL), Python pytest (0 failed, single-thread BLAS), run_honest_gate.sh exit 0, DESCRIPTION/pyproject.toml diff empty"
    requirement: "KPI-04"
    verification:
      - kind: other
        ref: "R CMD INSTALL --preclean . (DONE); devtools::test() -> 0 FAIL/141 WARN/13 SKIP/1837 PASS; cd python && uv pip install -e . --reinstall-package leafblower && pytest -> 161 passed/0 failed; bash benchmarks/run_honest_gate.sh -> all 3 steps clean, no errors; git diff --stat DESCRIPTION python/pyproject.toml -> empty"
        status: pass
    human_judgment: false

duration: ~70min
completed: 2026-08-16
status: complete
---

# Phase 03 Plan 08: docs/performance.md 9-method competitor coverage Summary

**Rewrote `docs/performance.md`'s 3-package `oris_soft`-only Competitors section into a full 9-method table (oris, oris_soft, raking, sinkhorn, greenkhorn, chebyshev, greg, logit, newton_kl), each citing its own `docs/methods/*.md` practitioner-implementations source, and published every new figure from plans 03-05/06/07's regenerated CSVs verbatim -- closing G-03-4 and, with it, G-03-1 end to end.**

## Performance

- **Duration:** ~70 min (dominated by re-running `benchmarks/run_honest_gate.sh`'s `optweight_linf` LP solve as part of Task 2's DoD re-verification: ~651s median wall time, matching 03-06-SUMMARY.md's measurement)
- **Started:** 2026-08-16 (approx)
- **Completed:** 2026-08-16
- **Tasks:** 2
- **Files modified:** 1

## Accomplishments

- Added `### medium_100k_5margins -- per-method competitor coverage` under `## Results`: a 9-row table transcribing `leafblower_oris`/`raking`/`newton_kl`/`greg`/`logit`/`chebyshev` alongside their competitor rows (`survey_calibrate_linear`, `survey_calibrate_logit`, `optweight_linf`), each row's `wall_s`/`max_error`/`max_w`/`min_w`/`n_eff` transcribed verbatim from `benchmarks/results/oris_soft_vs_competitors.csv`.
- Added `### k2_margin_pot_equiv -- greenkhorn/sinkhorn vs POT (K=2-margin special case)` under `## Results`: transcribes all 4 rows from `benchmarks/results/greenkhorn_sinkhorn_vs_pot.csv`, with the K=2-only scope, the unbounded-comparison caveat, and the `docs/methods/oris.md` Axis-1 citation stated in prose above the table.
- Rewrote `## Competitors` into a `Method | Closest competitor | Source | Objective match | Caveats` table with one row per each of the 9 leafblower methods, citing the exact `docs/methods/<name>.md` "Practitioner implementations & use cases" passage each choice came from (`oris.md`, `newton_kl.md`, `greg.md`, `logit.md`, `chebyshev.md`, `sinkhorn.md`, `greenkhorn.md`). Both honest scope caveats are transcribed verbatim in the table, not paraphrased: `optweight`'s objective-mismatch statement + no-`max.w`-argument statement (from the `optweight_linf` CSV row's note), and POT's no-bounds-mechanism + K=2-only-equivalence statement (from the `k2_margin_pot_equiv` CSV rows' notes).
- Re-verified the full project Definition of Done across all four gap-closure plans (05-08) combined: `R CMD INSTALL --preclean .` succeeds; `devtools::test()` reports `0 FAIL/141 WARN/13 SKIP/1837 PASS` (matches 03-06/03-07's own baseline); Python `pytest` under single-thread BLAS reports `161 passed/0 failed`; `bash benchmarks/run_honest_gate.sh` ran all 3 steps (stepstone regression gate, R competitor script including the `optweight_linf` arm, the new Python POT script) with zero errors; `git diff --stat DESCRIPTION python/pyproject.toml` is empty (`optweight`/POT stayed benchmark-only across plans 05-08, per D-09).
- Spot-checked 3 transcribed numbers against their CSV source cells before marking the plan done (T-03-08-01 mitigation): `leafblower_chebyshev` row (`wall_s=0.0278`, `max_error=3.647e-10`, `max_w=3.0000`, `min_w=0.1298`, `n_eff=67478.2`), `optweight_linf` row (`wall_s=660.3880`, `max_error=2.006e-12`, `max_w=1.8113`, `min_w=0.1887`, `n_eff=60307.8`), `pot_greenkhorn` row (`wall_s=0.000279`, `max_error=6.328e-15`, `max_w=2.4241`, `min_w=0.1669`, `n_eff=7110.6`) -- all matched exactly, no transcription slip found.

## Task Commits

1. **Task 1: Rewrite the Competitors section into a full 9-method table; extend Results** - `4270041` (docs)
2. **Task 2: Final Definition-of-Done verification across all four gap-closure plans** - verification-only, no code/docs changed, no commit (files: none per plan)

## Files Created/Modified

- `docs/performance.md` - added two new `## Results` subsections (per-method competitor coverage on `medium_100k_5margins`; greenkhorn/sinkhorn vs POT on `k2_margin_pot_equiv`) and rewrote `## Competitors` from a 3-package `oris_soft`-only section into a 9-method table

## Decisions Made

See `key-decisions` in frontmatter: greg/logit competitor rows point to their own distance-matched `survey::calibrate` arm rather than the raking-family trio; sinkhorn/greenkhorn's Objective match column states "equivalent at K=2 only" rather than a blanket "same", since the K>2 case has no POT-side measurement at all.

## Deviations from Plan

None -- plan executed exactly as written. Task 2's `run_honest_gate.sh` re-run regenerated `benchmarks/results/oris_soft_vs_competitors.csv` and `benchmarks/results/greenkhorn_sinkhorn_vs_pot.csv` as an unavoidable side effect of the verify step (same as 03-07's Task 3); the diff was `wall_s`-only timing noise (every accuracy/weight/n_eff field byte-identical), so both files were restored via `git checkout --` rather than committed as an out-of-scope change under this plan's `files_modified` scope (`docs/performance.md` only).

## Issues Encountered

Task 2's `run_honest_gate.sh` re-run was initially raced against a still-in-progress `R CMD INSTALL --preclean .` in a parallel shell, so the first attempt failed with `Error in library(leafblower) : there is no package called 'leafblower'`. Not a deviation from the plan's own scope -- re-ran after confirming the install completed (`EXIT: 0`), and the second run passed clean end to end (~70 min wall time, dominated by `optweight_linf`'s two `bench::mark` LP-solve iterations at ~651s median, consistent with 03-06-SUMMARY.md's own measurement of the same arm).

## User Setup Required

None -- no external service configuration required.

## Next Phase Readiness

- G-03-4 closed: `docs/performance.md`'s Competitors section covers all 9 leafblower methods against their own doc-named closest competitor, every figure traceable to a regenerated CSV, both honest scope caveats stated on the page itself.
- G-03-1 closed end to end: measured in plans 03-05/06/07, documented here in 03-08.
- This is the final plan in the 03-honest-performance-gate gap-closure wave (plans 05-08). Full DoD confirmed green across the combined batch. No blockers.

---
*Phase: 03-honest-performance-gate*
*Completed: 2026-08-16*

## Self-Check: PASSED
