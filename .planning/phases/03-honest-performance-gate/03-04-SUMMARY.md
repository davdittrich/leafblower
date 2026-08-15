---
phase: 03-honest-performance-gate
plan: 04
subsystem: testing
tags: [r, benchmarking, documentation, requirements, honest-performance-gate]

# Dependency graph
requires:
  - phase: 03-honest-performance-gate
    provides: "03-01: tracer measurement + benchmarks/run_honest_gate.sh + honest-gate opt-in block. 03-02: full measured table (medium/large/known-limit classes, full D-07 competitor set). 03-03: developer-selected paired ceiling (wall_s<=0.5s, n_eff>=60000) and docs/performance.md, the single linked methodology page."
provides:
  - "README.md: the package's one-line headline claim, linked to docs/performance.md, with no autumn mention (D-05) and no Phase 5 distribution scope."
  - "tasks/prd-leafblower-core.md: all three '100K rows, 5 margins' contradiction sites marked superseded in place, naming docs/performance.md; the void test-lbfgsb.R reference explicitly annotated as void rather than left standing."
  - ".planning/REQUIREMENTS.md: US-003 and KPI-04 rewritten to state what was measured, on which classes, naming benchmarks/oris_soft_vs_competitors.R and the honest gate: assertion in tests/testthat/test-bench-gate.R; Traceability and Superseded tables updated."
  - "tests/testthat/test-bench-gate.R: kk1204 block re-gated onto the file's own LBW_BENCH_GATE convention (no more CI-inverted default execution), relabelled as a shape-matched regression floor with a cross-reference to docs/performance.md's known-limit section."
affects: []

actuals:
  tokens: 3480
  tasks: 3
  commits: 3

tech-stack:
  added: []
  patterns:
    - "A PRD statement that contradicts a later measurement is superseded in place (marker on the same line, never a deletion) so the historical promise stays legible while a live reader lands on the current figure — matches the project's existing 'frozen historical document' convention for tasks/prd-leafblower-core.md."
    - "A test fixture that shares SIZE parameters (n/K/seed/max_weight) with a documented degenerate case but not its TARGET SKEW does not exercise that degenerate case — the two must not be conflated in a test's own description, even when the plan that requested the change assumed they were the same fixture."

key-files:
  created:
    - README.md
  modified:
    - .planning/REQUIREMENTS.md
    - tasks/prd-leafblower-core.md
    - tests/testthat/test-bench-gate.R

key-decisions:
  - "README.md's orienting sentence does not name autumn anywhere (D-05, user-emphatic, verified twice this phase in 03-CONTEXT.md) — rephrased around DESCRIPTION's own Title/Description wording instead of the drop-in-replacement framing that would have required naming it."
  - "The kk1204 block's fixture uses a uniform 1/5 target (test-bench-gate.R's own literal parameters), not the skewed 0.3/0.175x4 target that 03-02 used to reproduce the documented degenerate case. This phase does NOT change the block's target to match — that's a fixture-content change outside Task 3's explicit scope (skip-guard convention + relabelling only) — but the mismatch is now stated in-comment so a future reader is not misled into thinking this block tests the known limit it cites."
  - "requirements.mark-complete was run for US-003 and KPI-04 and correctly no-oped (both reported not_found, zero bytes written) — their Traceability Status stays 'Partial' by design, since the PRD's literal 'AND' target (1M+ rows AND 20+ margins together) was never demonstrated by one measurement and the verbose=1 routing-reason clause was not re-exercised by Phase 3."

patterns-established:
  - "A supersede marker sits on the SAME line as the statement it retires (not on a following line), so it survives regardless of later line-number drift elsewhere in the document and is trivially greppable/verifiable."

requirements-completed: []  # US-003/KPI-04 remain 'Partial' by design — see key-decisions; requirements.mark-complete confirmed no-op (not_found, 0 bytes written).

coverage:
  - id: D1
    description: "README.md created: package name, one orienting sentence (agreeing with DESCRIPTION, no autumn mention), the paired headline claim (wall_s/max_error/bound/n_eff) for oris_soft/medium_100k_5margins, and a link to docs/performance.md. At most 25 lines; no Phase 5 distribution scope (install/CRAN/badges)."
    requirement: "US-003"
    verification:
      - kind: unit
        ref: "Rscript -e stopifnot chain over README.md (<=25 lines, docs/performance.md link present, oris_soft present, no install/badge/CRAN strings) -> 'OK readme'"
        status: pass
      - kind: unit
        ref: "grep -ci 'install.packages|pip install|badge|shields.io|CRAN' README.md -> 0; grep -ci 'autumn' README.md -> 0"
        status: pass
    human_judgment: false
  - id: D2
    description: "tasks/prd-leafblower-core.md's three '100K rows, 5 margins' contradiction sites (Sec.1 Success Criteria, US-006 L-BFGS-B AC, Sec.11 KPI row) each marked superseded in place, naming docs/performance.md; the Sec.11 row's test-lbfgsb.R reference explicitly annotated as void (file does not exist, solver withdrawn) rather than deleted."
    requirement: "US-003"
    verification:
      - kind: unit
        ref: "Rscript -e stopifnot chain: grep('100K rows, 5 margins') matches >=2 lines, each carries 'Superseded' within its own line, docs/performance.md referenced -> 'OK prd markers'"
        status: pass
      - kind: unit
        ref: "grep -c 'test-lbfgsb' tasks/prd-leafblower-core.md: 6 before edit, 6 after edit (unchanged -- lines marked, none deleted)"
        status: pass
    human_judgment: false
  - id: D3
    description: ".planning/REQUIREMENTS.md's US-003 and KPI-04 entries rewritten to state what Phase 3 actually measured (three input classes, live oris_soft solver), name both live artefacts (benchmarks/oris_soft_vs_competitors.R, the honest gate: assertion in tests/testthat/test-bench-gate.R) in place of the void test-lbfgsb.R, and restate the K=20 structural-unachievability finding as a documented known limit rather than an open blocker. Traceability rows and the Superseded table updated."
    requirement: "KPI-04"
    verification:
      - kind: unit
        ref: "Rscript -e stopifnot chain: REQUIREMENTS.md names oris_soft_vs_competitors.R, LBW_BENCH_GATE, docs/performance.md -> 'OK bookkeeping'"
        status: pass
      - kind: unit
        ref: "grep -c 'oris_soft_vs_competitors.R' -> 3, grep -c 'test-bench-gate.R' -> 3, grep -c 'LBW_BENCH_GATE' -> 2, all >= 1 as required"
        status: pass
      - kind: unit
        ref: "test -f benchmarks/oris_soft_vs_competitors.R (exit 0) && CI=1 LBW_BENCH_GATE=1 NOT_CRAN=true Rscript -e testthat::test_dir(filter='bench-gate', stop_on_failure=TRUE) exits 0"
        status: pass
      - kind: unit
        ref: "git status --porcelain .beads/ shows only pre-existing working-tree noise (.beads/plans/active-plan.md, modified before this session started) -- leafblower-ylsy cited only, not modified"
        status: pass
    human_judgment: false
  - id: D4
    description: "tests/testthat/test-bench-gate.R's kk1204 block re-gated: skip_if(Sys.getenv('CI') != '') replaced with the file's own skip_if(Sys.getenv('LBW_BENCH_GATE') == ''), test_that() description and a new cross-referencing comment rewritten to say what the block actually measures (a shape-matched regression floor, not the documented degenerate known limit -- the block's uniform target converges trivially, unlike the skewed target that reproduces the degeneracy)."
    requirement: "US-003"
    verification:
      - kind: unit
        ref: "Pre-edit observation (CI unset, NOT_CRAN=true): status=0 iters=10 best_error=-7.376e-14 elapsed=1.6s, PASS 3/3 assertions"
        status: pass
      - kind: unit
        ref: "grep -c 'LBW_BENCH_GATE' tests/testthat/test-bench-gate.R -> 4; grep -c 'Sys.getenv(\"CI\")' -> 0"
        status: pass
      - kind: unit
        ref: "Default filtered run (LBW_BENCH_GATE unset): Rscript -e devtools::test(filter='bench-gate') -> [ FAIL 0 | WARN 0 | SKIP 4 | PASS 0 ] in 0.539s (kk1204 block now SKIP, not a paid 500k-row solve)"
        status: pass
      - kind: unit
        ref: "Gated run: CI=1 LBW_BENCH_GATE=1 NOT_CRAN=true testthat::test_dir(filter='bench-gate', stop_on_failure=TRUE) -> [ FAIL 0 | WARN 0 | SKIP 2 | PASS 10 ], kk1204 gate: status=0 iters=10 best_error=-7.376e-14 time=1.6s printed; regate verify script -> 'OK regate'"
        status: pass
    human_judgment: false
  - id: D5
    description: "Project Definition of Done unbroken after all three tasks: R CMD INSTALL --preclean . succeeds, devtools::test() reports 0 FAIL (WARN 141, SKIP 13 [+1 vs 03-03's baseline 12, from the kk1204 block now skipping by default], PASS 1836 [-3, the three kk1204 assertions no longer running under a plain devtools::test()]), Python parity 160 passed/0 failed, and benchmarks/run_honest_gate.sh still exits 0 and now also runs the re-gated kk1204 block."
    requirement: "US-003"
    verification:
      - kind: unit
        ref: "R CMD INSTALL --preclean . -> * DONE (leafblower)"
        status: pass
      - kind: unit
        ref: "OMP/OPENBLAS/MKL_NUM_THREADS=1 Rscript -e devtools::test() -> [ FAIL 0 | WARN 141 | SKIP 13 | PASS 1836 ]"
        status: pass
      - kind: unit
        ref: "OMP/OPENBLAS/MKL_NUM_THREADS=1 .venv/bin/python -m pytest (python/) -> 160 passed, 0 failed"
        status: pass
      - kind: unit
        ref: "CI=1 bash benchmarks/run_honest_gate.sh -> exit 0; regenerated CSV/env.txt and tests/testthat/fixtures/oris_shipgate_eval_results.json (unrelated test side-effect) both restored via git checkout -- to keep the diff scoped to this plan's four files"
        status: pass
    human_judgment: true
    rationale: "The WARN/SKIP/PASS count shift (SKIP +1, PASS -3) is the intended behavioral change from re-gating kk1204 -- confirming it reads as 'the gate correctly stopped running by default' rather than 'something silently broke' benefits from a human's read of the before/after numbers, not just the FAIL=0 automated check."

duration: ~30min
completed: 2026-08-15
status: complete
---

# Phase 3 Plan 4: Close the Loop — README, Retired Contradictions, Live Bookkeeping, Re-gated kk1204 Summary

**`README.md` now carries the package's one-line paired performance claim linked to `docs/performance.md`; all three `100K rows, 5 margins` contradiction sites in the frozen PRD are marked superseded in place; `.planning/REQUIREMENTS.md`'s US-003 and KPI-04 entries name the live measuring artefacts (`benchmarks/oris_soft_vs_competitors.R`, the `honest gate:` assertion) in place of the void `test-lbfgsb.R`; and the kk1204 block in `tests/testthat/test-bench-gate.R` is re-gated onto the file's own `LBW_BENCH_GATE` convention so an ordinary `devtools::test()` no longer pays for a 500,000-row solve.**

## Performance

- **Duration:** ~30 min
- **Tasks:** 3 of 3 complete
- **Files modified:** 4 (`README.md` created; `.planning/REQUIREMENTS.md`, `tasks/prd-leafblower-core.md`, `tests/testthat/test-bench-gate.R` modified)

## Accomplishments

- **Task 1 (README.md):** Created at 12 lines — package name, one orienting sentence (agreeing with `DESCRIPTION`'s Title/Description, no `autumn` mention per D-05), the paired headline claim (`wall_s=0.0427`, `max_error=3.35e-05`, `[0,3]` bound, `n_eff=67,489.4`), and a link to `docs/performance.md`. No Phase 5 distribution scope (no install instructions, badges, or CRAN references) leaked in.
- **Task 2 (retire contradictions, name live artefacts):** All three `100K rows, 5 margins` sites in `tasks/prd-leafblower-core.md` (§1 Success Criteria line 30, US-006 L-BFGS-B AC line 181, §11 KPI row line 673) marked superseded in place, each naming `docs/performance.md`; none deleted (`test-lbfgsb` line count unchanged at 6 before/after). `.planning/REQUIREMENTS.md`'s US-003 and KPI-04 entries rewritten to state what was actually measured (three input classes on a live `oris_soft` solver), name both live artefacts, and restate the K=20 structural-unachievability finding as a documented known limit rather than an open blocker. Traceability rows updated (KPI-04: Open → Partial); three new Superseded-table rows added.
- **Task 3 (re-gate kk1204):** Pre-edit state observed and recorded (status=0, iterations=10, `best_error=-7.376e-14`, elapsed=1.6s, PASS on all three assertions — the block was NOT failing before this plan). Skip guard flipped from the backwards `skip_if(Sys.getenv("CI") != "")` to the file's own `skip_if(Sys.getenv("LBW_BENCH_GATE") == "")` idiom; `test_that()` description and a new cross-referencing comment rewritten to state plainly that this fixture shares size parameters with `docs/performance.md`'s `known_limit_k20_uniform` class but NOT its target skew, so it does not exercise the documented degenerate case — it is a shape-matched regression floor. Existing thresholds (`best_error<=1e-3`, `elapsed<=30`) kept unchanged, reconfirmed against the pre-edit observation with large headroom.

## US-003 — Quoted Before/After (required by output spec)

**Before:**
> - [ ] **US-003**: 1M+ observations across 20+ margins calibrate in under 30 seconds,
>       single-threaded, so a census microsimulation researcher can iterate on synthetic
>       population models. `verbose = 1` prints the selected algorithm and routing reason.
>       *Status: **Partial**.* Routing observability shipped. The performance targets have
>       **never been verified against a live solver** and are internally contradictory:
>       §1 says medium-scale 100K rows / 5 margins < 1 s, §11 says < 2 s for the same shape,
>       and both were written against the removed `lbfgsb` whose measuring artefact
>       (`test-lbfgsb.R`) is void. Investigation `docs/investigations/2026-04-23-kk1204-convergence.md`
>       (commit `3effd3a`) found the composite gate "<30 s AND <1e-6" **structurally
>       unachievable** on K=20 uniform-random input (M_cell/n = 1.0 → zero compression
>       benefit; extrapolated ~1.76M iterations). → `leafblower-kk1.20.4` demands a REFRAME.

**After:**
> - [ ] **US-003**: 1M+ observations across 20+ margins calibrate in under 30 seconds,
>       single-threaded, so a census microsimulation researcher can iterate on synthetic
>       population models. `verbose = 1` prints the selected algorithm and routing reason.
>       *Status: **Partial**.* Routing observability shipped (not re-exercised by Phase 3 —
>       that clause's evidence predates this phase and is carried forward as-is). Phase 3
>       measured `oris_soft` against a live solver on three input classes, transcribed in
>       `docs/performance.md`: `medium_100k_5margins` (100K rows/5 margins, wall_s=0.0427,
>       max_error=3.35e-05), `large_stepstone_fulldata` (1,582,732 rows/9 margins, a real
>       salary-survey fixture, wall_s=3.5459 — clears the <30s large-scale budget on a real
>       shape exceeding the PRD's 1M-row target), and `known_limit_k20_uniform` (500,000
>       rows/20 margins, `m_cell/n = 1.0000`). No single measured fixture combines "1M+ rows
>       AND 20+ margins" as the PRD literally states — the 20-margin shape is measured
>       separately and found **structurally unable to clear an accuracy floor**:
>       `max_error = 5.229e-03` at wall_s=7.3934 under a bounded 500-iteration budget,
>       `m_cell/n = 1.0000` meaning every one of the 500,000 observations lands in its own
>       cell, so ORIS's cell-compression advantage yields nothing at this shape (confirmed
>       by `docs/investigations/2026-04-23-kk1204-convergence.md`, commit `3effd3a`, and
>       `leafblower-ylsy`'s closed research — cited only, not reopened, per D-03). This is
>       now a documented known limit (`docs/performance.md` § Known limit), not an open
>       blocker: the "<30 s AND <1e-6" composite gate on K=20 uniform-random input is
>       confirmed structurally unachievable with any solver currently implemented, and
>       `leafblower-kk1.20.4`'s REFRAME resolved on **dropping it as the headline basis**
>       (D-01) — not on any of the ticket's three originally-worded options. Live measuring
>       artefacts: `benchmarks/oris_soft_vs_competitors.R` (the measurement) and the
>       `honest gate:` assertion in `tests/testthat/test-bench-gate.R` (the regression gate,
>       `LBW_BENCH_GATE=1`), replacing the void `test-lbfgsb.R`.

Status marker unchanged (`- [ ]` / Partial): the `verbose = 1` clause was not re-exercised by this phase, and no single measurement demonstrates the PRD's literal "1M+ rows AND 20+ margins" combination together — both are stated honestly rather than rounded up to Implemented.

## KPI-04 — Quoted Before/After (required by output spec)

**Before:**
> - [ ] **KPI-04**: Large-scale — 1M rows, 20 margins, `max_weight = 3`, < 30 s. **Blocked on
>       US-003 reframe.** Its named measuring artefact was the Phase-2 gate against `lbfgsb`.

**After:**
> - [ ] **KPI-04**: Large-scale — 1M rows, 20 margins, `max_weight = 3`, < 30 s.
>       *Status: **Partial**.* Measured on the real 1,582,732-row/9-margin
>       `large_stepstone_fulldata` class (`max_weight = 3`): `oris_soft` calibrates in
>       wall_s=3.5459, comfortably inside the 30s budget, on a fixture larger than the PRD's
>       literal 1M-row target. The 20-margin corner of the original target
>       (`known_limit_k20_uniform`, 500,000 rows/20 margins) is measured separately and
>       documented as a known limit (`m_cell/n = 1.0000`, `max_error = 5.229e-03` at
>       wall_s=7.3934 — does not clear an accuracy floor at this shape; see
>       `docs/performance.md` § Known limit). Live artefacts: `benchmarks/oris_soft_vs_competitors.R`
>       (measurement) and the `honest gate:` assertion in `tests/testthat/test-bench-gate.R`,
>       run via `LBW_BENCH_GATE=1 CI=1 NOT_CRAN=true Rscript -e "testthat::test_dir('tests/testthat',
>       filter='bench-gate', stop_on_failure=TRUE)"`. No other performance-adjacent row in
>       this KPI list was found still naming a non-existent measuring artefact (checked
>       KPI-01 → `test-harvest.R`, KPI-02 → property test in `test-harvest.R`, KPI-03 →
>       `check_convergence` in `src/calib_dispatch.hpp:204`, KPI-05 → `R CMD check --as-cran`,
>       KPI-06 → no CI pipeline exists at all, which is KPI-06's own recorded open status, not
>       a stale artefact reference — all name something that exists or is honestly marked
>       absent, not a void file).

## PRD Supersede Markers (as written, required by output spec)

1. **§1 Success Criteria (line 30):**
   `| Performance — medium | 100K rows, 5 margins (3–5 cats each) < 1 s — **Superseded 2026-08-15**: written against the withdrawn L-BFGS-B solver and contradicted by §11's `< 2 s` row below; the live measured figure is in [docs/performance.md](../docs/performance.md). |`

2. **US-006 L-BFGS-B acceptance checkbox (line 181):**
   `- [ ] `rk_calibrate(..., algorithm=RK_ALG_LBFGSB)` converges on: 100K rows, 5 margins within 1 s — **Superseded 2026-08-15**: L-BFGS-B was withdrawn entirely (`RK_ALG_LBFGSB` is a permanently reserved enum hole, `src/leafblower.h:44`); this acceptance criterion has no live solver to measure. The live medium-scale figure, on `oris_soft`, is in [docs/performance.md](../docs/performance.md).`

3. **§11 KPI table row (line 673):**
   `| L-BFGS-B convergence | 100K rows, 5 margins < 2 s, `max_error < 1e-6` | `test-lbfgsb.R` Phase 1 gate — **Superseded 2026-08-15**: `test-lbfgsb.R` does not exist and the solver it measured was withdrawn (`RK_ALG_LBFGSB` is a permanently reserved enum hole). Live artefacts: `benchmarks/oris_soft_vs_competitors.R` and the `honest gate:` assertion in `tests/testthat/test-bench-gate.R`; see [docs/performance.md](../docs/performance.md). |`

## README.md — Number-to-CSV-Cell Mapping (required by output spec)

| README.md value | `docs/performance.md` source |
|---|---|
| `0.0427s` wall time | Headline claim / Results table, `medium_100k_5margins`, `leafblower_oris_soft` row: `wall_s = 0.0427` |
| `max margin error 3.35e-05` | Same row: `max_error = 3.347e-05` (displayed rounded as `3.35e-05`, matching `docs/performance.md`'s own headline rounding) |
| `[0,3]` weight bound | Same row: `max_w = 3.0000`, and the fixture's stated `[0, 3]` per-observation bound in the Headline claim paragraph |
| `67,489.4` effective observations | Same row: `n_eff = 67489.4` |
| `oris_soft` (solver name) | Headline claim: "`oris_soft` calibrates in 0.0427s..." |

Every numeric value in `README.md` is transcribed verbatim from `docs/performance.md`'s own already-verified transcription of the measured CSV — no new number was typed.

## kk1204 Block — Pre-edit Observation (required by output spec)

Run in isolation before any edit (`CI` unset, `NOT_CRAN=true` — required for `skip_on_cran()` to not short-circuit the run before reaching the guard under test):

```
kk1204 gate: status=0 iters=10 best_error=-7.376e-14 time=1.6s
[ FAIL 0 | WARN 0 | SKIP 3 | PASS 3 ]
```

**Verdict: PASS on all three assertions** (`status == 0L`, `best_error <= 1e-3`, `elapsed <= 30`). The block was NOT failing before this plan — no threshold was relaxed to hide a pre-existing failure.

**No threshold changed.** The existing `best_error <= 1e-3` and `elapsed <= 30` values are kept, reconfirmed against this plan's own pre-edit observation (both cleared by several orders of magnitude: `best_error` at `-7.376e-14` vs. a `1e-3` floor, `elapsed` at `1.6s` vs. a `30s` ceiling) rather than against 03-02's `known_limit_k20_uniform` measurement — that measurement (`max_error = 5.229e-03`, which would NOT clear this block's own `1e-3` threshold) used a different, skewed target and is not the same fixture. See the next section.

## Fixture Mismatch Finding (reported, not fixed — per Task 3's scope)

This plan's `read_first` step assumed the kk1204 block's fixture parameters were identical to 03-02's `known_limit_k20_uniform` measurement. On inspection they share `n=500000`/`K=20`/`seed=1204`/`max_weight=3`, but **not** the per-margin target: this in-repo block uses a literal uniform `1/5` target (`setNames(rep(1/cats, cats), ...)`), while 03-02 used the original 2026-04-23 kk1204 investigation's skewed `0.3/0.175x4` target — a decision 03-02-SUMMARY.md itself already recorded ("uniform targets on this exact data converge trivially ... which cannot back the 'known limit is unachievable' claim"). Under the uniform target, this block converges essentially exactly (`best_error ≈ -7e-14`), so it does **not** exercise the degenerate, near-infeasible behavior that `docs/performance.md`'s Known-limit section and the kk1204 investigation document.

Per Task 3's explicit action scope (skip-guard convention + relabelling only — no fixture-content change requested), this plan does **not** change the block's target to match the skewed one. Instead, the block's own comment now states this mismatch plainly, cross-referencing `docs/performance.md` and the investigation doc, so a future reader is not misled into thinking a passing kk1204-shape regression floor also confirms the documented known limit. If a future phase wants a live in-repo test of the actual degenerate case, it would need to adopt 03-02's skewed target — that is out of this plan's scope and is not filed as a ticket here (a documentation-only observation, not a defect).

## Task Commits

1. **Task 1: README.md — one headline line, one link, nothing else** - `5f5dc66` (docs)
2. **Task 2: Retire the contradictory numbers and name live artefacts** - `44e63a4` (docs)
3. **Task 3: Bring the kk1204 assertion under the one gate convention** - `4395dc4` (test)

## Files Created/Modified

- `README.md` (new) - one-line paired headline claim linked to `docs/performance.md`; no `autumn` mention; no Phase 5 scope.
- `tasks/prd-leafblower-core.md` - three `100K rows, 5 margins` sites marked superseded in place, naming `docs/performance.md`; no lines deleted.
- `.planning/REQUIREMENTS.md` - US-003 and KPI-04 entries rewritten; Traceability row for KPI-04 updated (Open → Partial); three new Superseded-table rows added.
- `tests/testthat/test-bench-gate.R` - kk1204 block re-gated onto `LBW_BENCH_GATE`; description and a new cross-referencing comment rewritten; thresholds unchanged.

## Decisions Made

See `key-decisions` in frontmatter, and the Fixture Mismatch Finding section above.

## Deviations from Plan

None — plan executed exactly as written for all three tasks. The fixture-mismatch discovery in Task 3 (see above) is a finding surfaced by the plan's own "establish the block's current state by observation" instruction, not a deviation from it; no code beyond what Task 3's action text specifies was changed as a result.

## Issues Encountered

**`skip_on_cran()` requires `NOT_CRAN=true`, not just `CI` unset, to actually run a guarded block (self-caught, no plan/code impact).** The first attempt to observe the kk1204 block's pre-edit state (with only `CI` left unset) reported all four `bench-gate` blocks as `SKIP: On CRAN`, because `testthat::test_dir()` — unlike `devtools::test()`'s own `pkgload`/`test_check` machinery — does not set `NOT_CRAN` itself, and every block in this file opens with `skip_on_cran()`. This matches a decision already carried forward from 03-02/03-03 (`benchmarks/run_honest_gate.sh` exports `NOT_CRAN=true` explicitly for the same reason). Re-ran with `NOT_CRAN=true` added and got the real pre-edit observation reported above. No incorrect data reached the summary; caught before any threshold reasoning was done.

**`benchmarks/run_honest_gate.sh` and `devtools::test()` regenerate two files as side effects, both restored via `git checkout --` to keep the diff scoped:** `benchmarks/results/oris_soft_vs_competitors.csv` (run-to-run wall-time noise, same pattern as 03-03) and, newly observed this plan, `tests/testthat/fixtures/oris_shipgate_eval_results.json` (a ship-gate evaluation fixture some test in the full suite regenerates with a fresh `generated_at` timestamp and flipped `ship_decision`/`omega_mode_default` values — unrelated to any of this plan's three tasks; last touched by an unrelated prior commit `3fa915e`). Both restored via `git checkout --` before committing.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Phase 3 (Honest Performance Gate) is complete: all 4 plans (03-01..03-04) executed. The headline claim is measured, gated, published on `README.md` and `docs/performance.md`, and every performance-adjacent bookkeeping row (`US-003`, `KPI-04`, and the checked sibling KPI rows) names a live artefact. The kk1204 assertion runs under the same `LBW_BENCH_GATE` convention as every other heavy gate in the file.
- `leafblower-kk1.20.4`'s REFRAME resolved on D-01 (dropped as headline basis) — confirmed not reopened; `leafblower-ylsy` cited only, confirmed not reopened or commented on (`git status --porcelain .beads/` shows only pre-existing, session-start working-tree noise on `active-plan.md`).
- Carry-forward for any future phase touching `tests/testthat/test-bench-gate.R`'s kk1204 block: its fixture target is uniform, not the skewed target that reproduces the documented degenerate case (see Fixture Mismatch Finding above) — a future change wanting a live regression test of the actual known limit needs to adopt 03-02's skewed target, which this plan intentionally left unchanged.
- Phase 4 (Truthful Surface) and Phase 5 (CRAN + PyPI Release) are next in the roadmap; neither depends on anything left open by this plan.

---
*Phase: 03-honest-performance-gate*
*Completed: 2026-08-15*

## Self-Check: PASSED

All 3 task commits (`5f5dc66`, `44e63a4`, `4395dc4`) confirmed present via `git log --oneline --all`. All 4 claimed files confirmed present on disk: `README.md`, `.planning/REQUIREMENTS.md`, `tasks/prd-leafblower-core.md`, `tests/testthat/test-bench-gate.R`.
