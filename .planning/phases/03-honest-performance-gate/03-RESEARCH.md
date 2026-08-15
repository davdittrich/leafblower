# Phase 3: Honest Performance Gate - Research

**Researched:** 2026-08-15
**Domain:** In-repo documentation/benchmark truth-telling (no new external tech stack — this
phase restates an existing performance claim as a measured, regressable gate)
**Confidence:** HIGH (all load-bearing claims verified by reading source/docs/beads this
session; the one open item is whether the existing `test-bench-gate.R` kk1204 assertion
currently passes, which needs an executor-time run, not research-time inference)

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

**REFRAME scope — kk1204 split-off**
- D-01: The `<30s AND <1e-6` gate on kk1204 (K=20 uniform-random, `M_cell/n=1.0`) is
  confirmed structurally unachievable and is dropped as the basis for any headline claim.
  kk1204 is NOT part of Phase 3's restated performance number in any form.
- D-02: kk1204 is not "ongoing algorithmic research" to track going forward —
  `leafblower-ylsy` (the ticket meant to carry that research) is already CLOSED
  (2026-05-03) with a firm conclusion: kk1204 is a near-infeasible **degenerate fixture**
  (DEFF 8000-14000, n_eff 71-118 across ALL solvers — ieppa, newton_kl, auto), not an
  algorithm gap. See `docs/investigations/2026-05-02-ylsy-cp-ipm-spike-result.md`.
  Reversibility: reversible — documentation framing choice, not code.
- D-03: No beads work needed on `leafblower-ylsy` — just cite its close reason and the
  2026-05-02 report as the SC3 "known limit" explanation. Do not reopen or comment on it.

**Headline claim — solver and metric**
- D-04: `oris_soft` is the primary headline solver, NOT `autumn` and NOT raking+SQUAREM
  alone. raking+SQUAREM stays cited as the fallback path for zero-compression-benefit
  inputs (K=20-uniform-random-style), not as the flagship number.
- D-05: User correction — do not mention `autumn` in any comparative/relevant position
  anywhere in this phase's output. Autumn was never released to CRAN and is not a real
  competitor; using it as the baseline ("Nx faster than autumn") is dishonest framing this
  phase exists specifically to eliminate.
- D-06: Headline metric (speed/wall-time vs. design-effect/effective-observations) is
  Claude's discretion — pick whichever number is the strongest, most honest differentiator
  once the fresh benchmark actually runs.

**Competitor selection — doc-grounded, not ad hoc**
- D-07: Competitors are NOT picked freely from `benchmarks/study/`. Each solver's
  `docs/methods/*.md` file already names its direct competitors in both R and Python. For
  `oris_soft`, pull from `docs/methods/oris.md`'s existing "How leafblower deviates"
  comparison table (bounds-handling row: `survey::calibrate`, `icarus`, `ReGenesees`).
- D-08: The fresh benchmark run reuses `benchmarks/` infrastructure conventions but is a
  NEW, narrow run — not a read of existing `benchmarks/study/report/tables/*.csv` output.
- D-09: New R/Python competitor packages needed for the fresh run are added as
  **benchmark-only dependencies**, scoped under `benchmarks/` — NOT added to the package's
  own `DESCRIPTION` or `pyproject.toml`. Reversibility: reversible.

**Claim mechanics**
- D-10: The restated performance number is a **hard gate** (pass/fail, regressable),
  matching the existing stepstone `LBW_BENCH_GATE=1` pattern.
- D-11: The stepstone gate stays **opt-in** (`LBW_BENCH_GATE=1`); Phase 3 does not change
  that contract.
- D-12: SC4's "one command" requirement is satisfied by **one wrapper command** that runs
  both the existing stepstone gate and the new oris_soft-vs-competitors comparison.
- D-13: Claim location: **README headline (one line) + linked docs/ page**. Not everything
  inlined into README.

**kk1204 documentation (SC3)**
- D-14: SC3's "known limit" documentation stays in
  `docs/investigations/2026-04-23-kk1204-convergence.md` (already exists), cross-linked
  from the `leafblower-ylsy` close reason and the 2026-05-02 follow-up report. No new
  user-facing mention in `docs/methods/*.md` is needed.

### Claude's Discretion
- Which metric (speed vs. design-effect/ESS) headlines the restated claim (D-06).
- Exact wrapper-command shape/location for the combined benchmark entry point (D-12).
- Exact README wording and linked-doc structure, within the "one line + linked page"
  shape (D-13).

### Deferred Ideas (OUT OF SCOPE)
- Full comparative benchmark study against every package in `benchmarks/study/` — belongs
  to the `leafblower-2ouc` epic (article-length deliverable), not Phase 3. Phase 3's fresh
  benchmark is a narrow subset (1-2 doc-named competitors for the headlining solver).
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| US-003 | 1M+ obs / 20+ margins calibrate <30s single-threaded; `verbose=1` prints routing. Currently Partial: internally contradictory medium-scale target, gate never verified against a live solver. | See "The 1s/2s Contradiction" and "The kk1204 Gate" sections below — both source locations verified, both quoted verbatim. |
| KPI-04 | Large-scale — 1M rows, 20 margins, `max_weight=3`, <30s. Blocked on US-003 reframe; named measuring artefact was the Phase-2 gate against removed `lbfgsb`. | See "KPI Table Location" and "Existing Gate Infrastructure" — the row to rewrite is identified by file:line, and a live replacement artefact (`test-bench-gate.R`) already exists and is inventoried. |
</phase_requirements>

## Summary

This phase requires no new library research — it is a documentation-and-benchmark
truth-telling exercise entirely within this repo. Every fact the planner needs is a
grep-and-read away, and this research session read every source: the PRD's contradictory
§1/§11 numbers, the kk1204 investigation and its closure, the doc-named competitor table,
the existing (but under-used) `test-bench-gate.R` gate, and the state of `README.md`
(**it does not exist — git history confirms it has never existed in this repo**).

Three findings materially affect planning beyond what CONTEXT.md already locked:

1. **There is no README.md today.** D-13 says "README headline + linked docs page," but
   the plan must **create** README.md from scratch, not edit an existing one. This is a
   bigger task than a one-line edit — CRAN packages conventionally need a README (badges,
   install instructions, quick example) and Phase 5 (Distribution/CRAN) will also touch it.
   Coordinate scope: Phase 3 only needs the performance headline; do not let this phase's
   plan balloon into "write the whole package README" (that risk is real given no baseline
   exists).

2. **A live kk1204 gate already exists in code** (`tests/testthat/test-bench-gate.R:28-54`),
   but it does not match the ticket's stated `<30s AND <1e-6` framing: it runs `n=500,000`
   (not 1,000,000), asserts `best_error <= 1e-3` (not `1e-6`), and is skipped by
   `Sys.getenv("CI") != ""` (not the documented `LBW_BENCH_GATE` pattern). Per D-01,
   kk1204 is dropped as the headline basis anyway — but this test still exists, still runs
   in non-CI environments, and its relationship to the "known limit" documentation (SC3)
   needs an explicit planning decision: keep it as a regression floor on the *current*
   `oris`/`ieppa` best-effort number, relabel its intent in a comment, or leave it
   untouched. Do not assume it currently passes — this needs an executor-time run.

3. **All three doc-named R competitors are already installed and two are already
   `Suggests:`-declared** (`survey` and, implicitly usable, the base R environment has
   `icarus` and `ReGenesees` installed too, per this session's `requireNamespace()` check).
   `DESCRIPTION` already lists `survey` under `Suggests`. This simplifies D-09
   considerably: R-side, no NEW packages need to be added under `benchmarks/`-scoped deps
   at all — `icarus` and `ReGenesees` merely need to be added to `Suggests` (or
   benchmark-scoped, per D-09's letter) since they are not yet declared anywhere in
   `DESCRIPTION`, even though they resolve locally.

**Primary recommendation:** Treat this as a **decision-and-documentation phase with one
small new benchmark script**, not a solver-engineering phase. The bulk of the work is: (a)
run a narrow, fresh `oris_soft` vs. `survey::calibrate()`/`icarus::calibration()`/
`ReGenesees::e.calibrate()` benchmark on a realistic large-scale fixture reusing
`benchmarks/` conventions (pattern: `benchmarks/newton_kl_bench.R`'s `bench::mark(...,
iterations=2, check=FALSE)` loop, NOT a from-scratch harness); (b) create `README.md` with
one headline line; (c) create a linked `docs/` methodology page; (d) rewrite the KPI-04 row
in `.planning/REQUIREMENTS.md` (and, if the PRD source-of-truth pattern is followed,
`tasks/prd-leafblower-core.md` — see "KPI Table Location") to name the new live artefact;
(e) wire a single wrapper command that runs both the stepstone `LBW_BENCH_GATE=1` gate and
the new comparison script.

## Architectural Responsibility Map

The web-app tier taxonomy (Browser/SSR/API/CDN/DB) does not apply to this domain — leafblower
is a compiled C++ core with R/Python bindings and a benchmark/documentation layer. Mapped to
this project's actual tiers:

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Solver performance (the number itself) | C++ core (`src/oris.cpp`, `oris_soft` = `RK_ALG_ORIS_SOFT`) | — | The number is a property of the compiled solver; nothing in this phase changes solver code. |
| Benchmark measurement | R benchmark script under `benchmarks/` | testthat gate (`tests/testthat/test-bench-gate.R`) | Measurement lives in `benchmarks/` (ad hoc, human-run); the pass/fail *gate* lives in `tests/testthat/` per the existing `LBW_BENCH_GATE=1` pattern (D-10). |
| Regression protection | `tests/testthat/test-bench-gate.R` (opt-in, `LBW_BENCH_GATE=1`) | — | D-11: stays opt-in; Phase 3 does not fold it into the default DoD gate. |
| Headline claim surface | `README.md` (new file) | linked `docs/` page | D-13: one line in README, full methodology in `docs/`. |
| KPI bookkeeping | `.planning/REQUIREMENTS.md` §"Project-level acceptance" | `tasks/prd-leafblower-core.md` §11 (historical source, itself superseded per REQUIREMENTS.md) | SC4 requires the KPI table's performance rows to each name a live measuring artefact — REQUIREMENTS.md is the currently-maintained copy; the PRD is frozen/historical. |
| Known-limit documentation | `docs/investigations/2026-04-23-kk1204-convergence.md` (existing) | `leafblower-ylsy` beads close-reason (citation only) | D-14: no new file, cross-link only. |

## Standard Stack

No new external libraries are introduced by this phase. Everything needed is already a
declared or resolvable dependency:

### Core (already present)
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| `bench` | (R, already `Suggests:` in `DESCRIPTION`) | `bench::mark()` — the project's established micro-benchmark timer | Already used in `benchmarks/newton_kl_bench.R:29` and the study harness (`benchmarks/study/spec/contract.md:150`); matches CLAUDE.md's determinism-protocol vocabulary. |
| `survey` | (R, already `Suggests:` in `DESCRIPTION`, installed — `requireNamespace("survey")` = TRUE this session) [VERIFIED: DESCRIPTION:16, R session `requireNamespace` check] | Doc-named competitor #1 (`survey::calibrate()`) | Named in `docs/methods/oris.md:218,261` "How leafblower deviates" table. |
| `icarus` | CRAN, installed locally (`requireNamespace("icarus")` = TRUE this session) [VERIFIED: R session check] but **NOT yet in `DESCRIPTION` `Suggests:`** [VERIFIED: DESCRIPTION full-file read] | Doc-named competitor #2 (`icarus::calibration()`) | Named in `docs/methods/oris.md:221,261`. |
| `ReGenesees` | CRAN, installed locally (`requireNamespace("ReGenesees")` = TRUE this session) [VERIFIED: R session check] but **NOT yet in `DESCRIPTION` `Suggests:`** [VERIFIED: DESCRIPTION full-file read] | Doc-named competitor #3 (`ReGenesees::e.calibrate()`) | Named in `docs/methods/oris.md:222,261`. |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| `arrow` (R) | already used (`benchmarks/regression_timer.R:4`) | Parquet I/O for benchmark fixtures | Reuse the stepstone-fulldata data-generation pattern (`benchmarks/stepstone_benchmark.R`) if a fresh large fixture is needed, rather than inventing a new data format. |
| `jsonlite` (R) | already used | Target-margin JSON I/O | Same reuse rationale. |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `bench::mark()` timing | `microbenchmark` package | Not used anywhere in this repo; `bench` is already the established convention (D-08's "reuse `benchmarks/` conventions" makes this a non-choice — `bench` wins by existing precedent, not by fresh evaluation). |
| Doc-named R competitors (`survey`/`icarus`/`ReGenesees`) | Ad hoc pick from `benchmarks/study/` roster (autumn, ipfn, sbw, gecal, nonprobsvy, etc.) | Explicitly rejected by user correction D-07/D-08 — do not invent a comparison basis outside the doc-named set. |

**Installation:** None required for R — all three competitor packages already resolve via
`requireNamespace()` in this environment. If `DESCRIPTION` `Suggests:` is extended (to keep
`R CMD check` honest about what the benchmark script needs), add:
```
Suggests: ..., icarus, ReGenesees
```
Per D-09, if the planner instead treats these as strictly benchmark-scoped (not
package-level `Suggests:`), a `benchmarks/DESCRIPTION`-style scoping file or an in-script
`requireNamespace()` guard with an informative skip message is the alternative — pick
whichever the plan's benchmark script structure calls for; both are consistent with D-09's
"benchmark-only, not in the shipped package's own manifest" intent, since `Suggests:` in the
package's own `DESCRIPTION` is arguably NOT benchmark-scoped enough to satisfy D-09's letter.
**Recommendation: do NOT add to package `DESCRIPTION`; guard with `requireNamespace()` +
`skip_if_not_installed()`/informative message inside the `benchmarks/` script**, matching
D-09 exactly (`survey` is a pre-existing exception already in `Suggests:` for unrelated
reasons — do not use it as precedent to add the other two).

**Version verification:** Not applicable — no new package installs; both `icarus` and
`ReGenesees` resolved at their currently-installed local versions this session (exact
version strings not captured; the planner's benchmark script should print
`packageVersion()` for all three competitors into the published docs page for
reproducibility, since that is exactly the kind of number a stale claim would silently
drift on).

## Package Legitimacy Audit

**Not applicable in the SLOP/SUS/OK sense** — no new packages are being installed from a
registry search or WebSearch discovery. All three competitor packages (`survey`, `icarus`,
`ReGenesees`) are:
1. Named explicitly in this repo's own `docs/methods/oris.md` (a first-party, already-cited
   source with DOI/journal citations for each — `[zardetto2015regenesees]`,
   `[rebecq2017icarus]`), not discovered via search this session.
2. Already resolvable via `requireNamespace()` in the current R environment (verified by
   running the check, not by registry lookup alone).
3. `survey` is already a declared `Suggests:` dependency in `DESCRIPTION` (`DESCRIPTION:16`,
   read this session).

| Package | Registry | Age | Downloads | Source Repo | Verdict | Disposition |
|---------|----------|-----|-----------|-------------|---------|-------------|
| `survey` | CRAN | Long-established (Lumley, cited 2010 paper) | High (foundational R survey package) | CRAN canonical | OK | Approved — already `Suggests:` |
| `icarus` | CRAN | Established (Rebecq, 2017 CRAN release cited) | — | CRAN canonical | OK | Approved — add to benchmark scope, not package `Suggests:` per D-09 |
| `ReGenesees` | CRAN | Established (Zardetto 2015, Italian national statistics office) | — | CRAN canonical | OK | Approved — add to benchmark scope, not package `Suggests:` per D-09 |

**Extended 2026-08-15 (gap closure, G-03-1/G-03-4):** per-method competitor coverage
(UAT gap closure) requires two more benchmark-only packages, each named explicitly in this
repo's own `docs/methods/*.md` (not discovered via ad hoc search this session) and verified
directly against the registry (CRAN PDF manual fetch / PyPI JSON), not by registry-listing
alone:

| Package | Registry | Age | Downloads | Source Repo | Verdict | Disposition |
|---------|----------|-----|-----------|-------------|---------|-------------|
| `optweight` | CRAN | Established (v2.0.1, published 2026-03-19; author Noah Greifer, ORCID 0000-0003-3067-7154 — also maintains `WeightIt`/`cobalt`, well-known R causal-inference ecosystem) | — | CRAN canonical (`cran.r-project.org/package=optweight`) | OK | Approved — add to benchmark scope only (chebyshev's minimax-flavored competitor per `docs/methods/chebyshev.md` §Practitioner implementations), not package `Suggests:`, per D-09. Verified via direct CRAN PDF manual fetch this session: `norm = "linf"` argument confirmed real (`f(w,b,s) = max_i|w_i-b_i|`); no `max.w` argument exists (only `min.w` floor) — bound-compliance is unverifiable on this competitor, documented as a caveat, not silently assumed compliant. |
| `POT` (import `ot`) | PyPI | Established (v0.9.7.post1; JMLR-published, `requires_python >= 3.7`) | High (standard Python OT library) | PyPI canonical (`pypi.org/project/POT`) | OK | Approved — add to benchmark scope only (greenkhorn/sinkhorn's competitor per `docs/methods/greenkhorn.md` and `docs/methods/sinkhorn.md` §Practitioner implementations, both already citing `[flamary2021pot]` in-repo), not added to `pyproject.toml`, per D-09's Python analogue. Verified via direct PyPI JSON metadata fetch this session. `ot.bregman.greenkhorn(a, b, M, reg, ...)` confirmed to take a cost matrix `M` + scalar `reg`, no arbitrary-kernel or bounds argument — the K=2-margin fixture and `M = -log(prior)`, `reg = 1` construction (recovering kernel = prior exactly) is required to make this a faithful comparison; documented as scope, not silently generalized to K>2. |

**Packages removed due to `[SLOP]` verdict:** none.
**Packages flagged as suspicious `[SUS]`:** none.

## Architecture Patterns

### System Architecture Diagram

```
                         ┌─────────────────────────────┐
                         │  benchmarks/ (new script)    │
                         │  oris_soft_vs_competitors.R  │
                         └──────────────┬───────────────┘
                                         │ single-thread BLAS env
                                         │ (OMP/OPENBLAS/MKL_NUM_THREADS=1)
                                         ▼
        ┌────────────────────────────────────────────────────┐
        │  bench::mark(run = harvest(..., method="oris_soft"),│
        │              iterations = N, check = FALSE)         │
        │  ── vs (same fixture, sequential loop, not the      │
        │     "interleaved before/after" case — this is an    │
        │     absolute cross-solver comparison, not a         │
        │     before/after regression of the SAME code) ──    │
        │  survey::calibrate() / icarus::calibration() /      │
        │  ReGenesees::e.calibrate()                          │
        └───────────────────────┬──────────────────────────────┘
                                 │ writes wall_s, max_err, converged
                                 ▼
                 ┌───────────────────────────────┐
                 │ published figure: headline #   │
                 └───────────┬─────────┬──────────┘
                              │         │
                 ┌────────────▼─┐   ┌───▼─────────────────────┐
                 │ README.md    │   │ docs/<new page>.md      │
                 │ (new file —  │   │ full methodology,       │
                 │  1-line hdl) │   │ machine spec, versions   │
                 └───────────────┘   └──────────────────────────┘

        ┌──────────────────────────────────────────────────────┐
        │ tests/testthat/test-bench-gate.R                      │
        │  — EXISTING kk1204 assertion (n=500k, best_error<=1e-3,│
        │    skip_if(CI!="")) — pre-dates this phase; D-01 drops│
        │    kk1204 as headline basis but does not mandate      │
        │    deleting this test. Planner must decide its fate.  │
        │  — new/updated assertion for the oris_soft headline   │
        │    number, gated LBW_BENCH_GATE=1 per D-10/D-11        │
        └──────────────────────────────────────────────────────┘
                                 │
                                 ▼
                 ┌────────────────────────────────────┐
                 │ benchmarks/run_honest_gate.sh (new) │
                 │ ONE command (D-12): runs stepstone   │
                 │ LBW_BENCH_GATE=1 test AND the new    │
                 │ oris_soft-vs-competitors script      │
                 └────────────────────────────────────┘
```

### Recommended Project Structure
```
benchmarks/
├── oris_soft_vs_competitors.R   # new: the narrow, doc-grounded fresh benchmark (D-08)
├── run_honest_gate.sh           # new: SC4's one-command wrapper (D-12), mirrors run_allmethod.sh's shape
docs/
├── methods/oris.md              # existing — cited as competitor-selection source, not edited unless metric requires it
└── performance.md               # new (exact name/location is Claude's discretion, D-13) — full methodology
README.md                        # new — does not currently exist; one headline line + link
tests/testthat/
└── test-bench-gate.R            # existing — extend or add adjacent test for the new gate (D-10)
.planning/REQUIREMENTS.md        # KPI-04 row rewritten to name the new live artefact (SC4)
```

### Pattern 1: Cross-solver absolute-timing loop (NOT interleaved before/after)
**What:** Loop over a fixed list of methods/competitors, call `bench::mark(iterations=2,
check=FALSE, memory=FALSE, filter_gc=FALSE)` once per method against the same fixture,
collect median wall-time + reported error into a `data.frame`.
**When to use:** This is the correct pattern for "solver A vs. solver B on the same input,"
which is what SC1/SC2/D-06 need. It is distinct from CLAUDE.md's "interleaved before/after
in one `bench::mark()` call" rule, which governs a *different* scenario — measuring the
effect of a code change on the SAME algorithm (a regression check), where sequential
before/after runs are confounded by machine-state drift between the two runs. Comparing two
different algorithms on the same fixture has no such confound to interleave away; the
existing codebase's own pattern (`benchmarks/newton_kl_bench.R:28-33`) already does this
non-interleaved, one-`bench::mark()`-call-per-method loop, and it is not a CLAUDE.md
violation.
**Example:**
```r
# Source: benchmarks/newton_kl_bench.R:28-33 (verified read this session)
for (m in methods) {
  res <- bench::mark(
    run = harvest(df, tgt, method = m, max_weight = 3,
                  max_iterations = 50, accelerate = (m == "oris")),
    iterations = 2, check = FALSE, memory = FALSE, filter_gc = FALSE
  )
  ...
}
```
**Caveat:** if the plan ALSO includes a "does this change regress the number we published
last time" check (i.e. `LBW_BENCH_GATE=1` regression protection, D-10/D-11), THAT check is
the one CLAUDE.md's interleaved-before/after rule governs — it compares this run's number
against a stored baseline, and if the comparison is done live (not baseline-file lookup),
it must interleave, not run sequentially before/after in separate processes.

### Pattern 2: Single-thread BLAS determinism envelope
**What:** `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1` set together, R-side
before the benchmark script runs (shell env, not in-script) — matches
`benchmarks/stepstone_all_methods.R:3`'s comment ("OMP_NUM_THREADS=1 set by caller") and
CLAUDE.md's own DoD command line.
**When to use:** Every timed run this phase produces, without exception — CLAUDE.md:90
states this explicitly and calls out that skipping it causes R↔Python and benchmark drift.
**Example:**
```bash
# Source: CLAUDE.md build/test section (verified read this session)
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript benchmarks/oris_soft_vs_competitors.R
```

### Anti-Patterns to Avoid
- **Building a parallel harness:** D-08 and the ROADMAP notes explicitly forbid a new
  benchmark framework — reuse `benchmarks/`'s existing script style (a plain `Rscript`,
  `bench::mark()`, `data.frame` output), not `benchmarks/study/`'s registry/adapter
  machinery (that belongs to `leafblower-2ouc`).
- **Citing `autumn` anywhere comparative:** D-05 is explicit and was a repeated user
  correction — do not resurrect the old "Nx faster than autumn" framing in any new copy,
  even implicitly (e.g., a benchmark script variable named `autumn_baseline`).
- **Picking competitors from `benchmarks/study/`'s registry:** D-07/D-08 — the registry
  contains ~17 packages (`sbw`, `gecal`, `nonprobsvy`, `ott-jax`, etc., per
  `leafblower-2ouc`'s ticket body) scoped for the separate article-length study; none of
  them are the doc-named `oris_soft` competitors.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Wall-time measurement | Custom `proc.time()` diffing with manual GC warmup | `bench::mark()` | Already the project's convention (`benchmarks/newton_kl_bench.R`, `benchmarks/study/spec/contract.md:150`); handles GC/warmup/robust-statistics concerns leafblower's own ad hoc timers (`regression_timer.R`) approximate manually and imperfectly. |
| Competitor calibration calls | Reimplementing `survey`/`icarus`/`ReGenesees`'s calibration APIs | Call the packages directly | These are the doc-named, cited, peer-reviewed reference implementations — the entire point of SC1 is comparing against a REAL competitor, not a reimplementation. |
| Regression gate wiring | A new env-var / CI convention | `LBW_BENCH_GATE=1` (existing, D-10/D-11) | One opt-in-heavy-gate convention already exists and is documented in CLAUDE.md and `tests/testthat/test-bench-gate.R`; a second convention would fragment the DoD story. |

**Key insight:** every mechanical piece this phase needs (timer, competitor APIs, gate
convention, single-thread-BLAS discipline) already exists in this repo. The actual work is
choosing what to measure and being honest about the answer — not writing new
infrastructure.

## Common Pitfalls

### Pitfall 1: Restating the same 1s/2s contradiction with a new pair of numbers
**What goes wrong:** A plan "fixes" US-003/SC2 by writing a new headline number without
retiring BOTH of the old ones. PRD §1 ("Performance — medium") and §11 (KPI table) state
DIFFERENT numbers for the identical shape (100K rows, 5 margins):
```
tasks/prd-leafblower-core.md:30:  | Performance — medium | 100K rows, 5 margins (3–5 cats each) < 1 s |
tasks/prd-leafblower-core.md:673: | L-BFGS-B convergence | 100K rows, 5 margins < 2 s, `max_error < 1e-6` | `test-lbfgsb.R` Phase 1 gate |
```
[VERIFIED: tasks/prd-leafblower-core.md:30,673 — read this session, quoted verbatim]. Both
are measured against `lbfgsb`, which is WITHDRAWN (`RK_ALG_LBFGSB=2` permanently reserved
hole, `src/leafblower.h:44`, per REQUIREMENTS.md's Superseded table) and `test-lbfgsb.R`
does not exist as a live gate.
**Why it happens:** The PRD is frozen (Draft v3, 2026-04-18) and REQUIREMENTS.md correctly
marks it historical, but nothing forces a plan to touch BOTH sites — a plan that only edits
README.md/docs/ could leave the PRD's own internal contradiction unresolved as a stale,
still-readable artifact.
**How to avoid:** SC2 says "states ONE number" — the plan must either (a) add a note to
`tasks/prd-leafblower-core.md` marking §1/§11 superseded (consistent with how
REQUIREMENTS.md already treats other PRD sections, e.g. the "§7 `-O3`" and "§6 enum" rows
in its Superseded table), or (b) rely on REQUIREMENTS.md/KPI table being the sole
maintained source going forward and say so explicitly in the new docs page. Either is
consistent with D-13 as long as a reader lands on exactly one number.
**Warning signs:** grep for "100K" or "< 1 s"/"< 2 s" after the plan lands — if both PRD
lines are still there unmarked, SC2 is not actually met.

### Pitfall 2: Treating kk1204's existing code gate as already retired
**What goes wrong:** Assuming D-01 ("kk1204 dropped as headline basis") means the kk1204
test in code is gone or irrelevant. It is not — it is live, currently executed outside CI:
```r
# tests/testthat/test-bench-gate.R:28-54 (verified read this session, quoted verbatim)
test_that("kk1204 gate: n=500k K=20 converges in <30s with best_error<1e-3", {
  skip_on_cran()
  skip_if(Sys.getenv("CI") != "")
  ...
  expect_lte(elapsed, 30, label = "elapsed: speed gate <30s on n=500k K=20")
  expect_lte(r$best_error, 1e-3, label = "best_error: quality gate below 1e-3")
})
```
Note this ALREADY differs from both the PRD's `<30s AND <1e-6` framing (it asserts `1e-3`,
not `1e-6`) and from the ticket's stated fixture (`n=500,000`/`K=20`, not the investigation
doc's `n=1,000,000`/`K=20`). This is a THIRD variant of the "kk1204 gate," distinct from
both the PRD number and the investigation's finding — the plan needs to reconcile or
retire it, not just leave it silently diverging further from whatever gets published.
**Why it happens:** the test was evidently adjusted at some point after the 2026-04-23
investigation (which used n=1M and found `<1e-6` unreachable) without a matching update to
the PRD or beads ticket language, and without matching the investigation's own recommended
reframe options ("K ≤ 10" or "<1e-4" — this test picked neither; it kept K=20 but relaxed
to `1e-3` and shrank `n`).
**How to avoid:** Have the plan explicitly decide this test's fate (keep as an internal
regression floor with a corrected label; delete since kk1204 is out of headline scope per
D-01; or fold into the SC3 "known limit" documentation as the live artefact that
demonstrates the ceiling). Do not assume it's dead code — it currently runs whenever
`Sys.getenv("CI") == ""`, i.e., in ordinary local dev sessions.
**Warning signs:** running `Rscript -e "testthat::test_dir('tests/testthat', filter='bench-gate')"` locally and getting either a silent 30s+ hang or a failure the plan didn't anticipate.

### Pitfall 3: Conflating the "known limit" input class with the headline benchmark's input class
**What goes wrong:** Using K=20/uniform-random/M_cell/n=1.0 (the input class SC3 says is a
documented known limit) as the SAME fixture for the SC1/SC2 headline number. The whole
point of D-01/D-04 is that the headline claim is measured on a DIFFERENT, non-degenerate
input class where `oris_soft` actually performs well (per the ylsy close reason, DEFF is
only pathological at that specific K=20-uniform fixture — not universally).
**Why it happens:** kk1204 is the most extensively-benchmarked fixture in the repo
(`benchmarks/results/newton_kl_kk1204.csv`, `benchmarks/results/tsvd_kk1204_K20.csv`,
multiple investigation docs) — it's the path of least resistance to reach for, but it's
exactly the wrong fixture for a positive claim.
**How to avoid:** Use `benchmarks/stepstone_bench_data.parquet`-style fixtures (9 margins,
835 categories, n=200K, per `benchmarks/stepstone_benchmark.R:1-15` — realistic
survey-shaped data, not K=20 uniform-random) or a purpose-built large-scale fixture with
non-degenerate `M_cell/n`, and state that input class explicitly in the docs page per SC1's
"on a stated input class" requirement.
**Warning signs:** the published headline number matching any of the kk1204 CSVs in
`benchmarks/results/`.

## Code Examples

### Existing single-thread BLAS + `bench::mark()` cross-solver pattern (reuse as-is)
```r
# Source: benchmarks/newton_kl_bench.R (verified read this session)
for (m in methods) {
  res <- bench::mark(
    run = harvest(df, tgt, method = m, max_weight = 3,
                  max_iterations = 50, accelerate = (m == "oris")),
    iterations = 2, check = FALSE, memory = FALSE, filter_gc = FALSE
  )
  r <- harvest(df, tgt, method = m, max_weight = 3,
               max_iterations = 50, accelerate = (m == "oris"))
  R <- attr(r, "result")
  cat(sprintf("  %-12s wall=%6.2fs status=%d max_err=%.3e iters=%d\n",
              m, as.numeric(res$median), R$status, R$max_error, R$iterations))
}
```

### Existing wrapper-command shape to mirror for SC4's "one command"
```bash
# Source: benchmarks/run_allmethod.sh (verified read this session, full file)
#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
echo '=== R benchmark ===' && Rscript benchmarks/allmethod_bench.R
echo '=== Python benchmark ===' && python benchmarks/allmethod_bench.py
```

### Existing opt-in gate pattern (D-10/D-11 — do not invent a new one)
```r
# Source: tests/testthat/test-bench-gate.R:1-17 (verified read this session, quoted verbatim)
test_that("stepstone-fulldata AB config meets merge floor (errRp) + Pearson agreement", {
  skip_on_cran()
  skip_if(Sys.getenv("LBW_BENCH_GATE") == "")
  rpt_path <- "benchmarks/stepstone_fulldata_homotopy_report.rds"
  skip_if(!file.exists(rpt_path))
  ...
})
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|---------------|--------|
| `iEPPA` solver name / L-BFGS-B as one of two algorithms | `oris`/`oris_soft` (renamed), 8 solvers total, `lbfgsb` withdrawn | NEWS.md "development" entry: "BREAKING: solver method 'ieppa'/'ieppa_soft' renamed to 'oris'/'oris_soft'"; REQUIREMENTS.md US-006 marked WITHDRAWN | Any performance claim tied to `lbfgsb` or citing "iEPPA" is stale terminology as well as stale numbers. |
| PRD §1/§11 targets vs. removed `lbfgsb` | No live measuring artefact for the medium-scale claim at all | Ongoing (this phase's SC2) | The plan must supply a NEW artefact, not just fix the number — `test-lbfgsb.R` is void. |

**Deprecated/outdated:**
- `test-lbfgsb.R`: does not exist; was the PRD §11 KPI table's named measuring artefact for
  the medium-scale target. [VERIFIED: REQUIREMENTS.md:105-108, read this session — quoted:
  "both were written against the removed `lbfgsb` whose measuring artefact
  (`test-lbfgsb.R`) is void."]
- `autumn` as a comparison baseline anywhere in new copy — D-05, user-mandated removal.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|----------------|
| A1 | `icarus` and `ReGenesees` should be added as benchmark-scoped dependencies rather than package `Suggests:` to satisfy D-09's intent, even though `survey` (already `Suggests:`) sets a precedent for the opposite | Standard Stack / Installation | Low — this is a phase-3 implementation-detail recommendation, not a locked decision; the planner/user can choose either path without breaking SC1-SC4. Flag for confirmation if the plan diverges from `survey`'s existing precedent. |
| A2 | The stepstone-style fixture (`benchmarks/stepstone_benchmark.R`, 9 margins/n=200K) is a suitable non-degenerate input class for the headline number, as opposed to a purpose-built 1M-row/20-margin fixture that still avoids `M_cell/n≈1.0` | Common Pitfalls / Pattern 3 | Medium — if the chosen fixture doesn't scale to something recognizably "large-scale" (SC1 says "large-scale performance figure"), the headline claim may look under-powered versus the PRD's original 1M-row ambition. The planner should size the fixture deliberately, citing `M_cell/n` for whatever is chosen. |
| A3 | `docs/methods/oris.md`'s bounds-handling competitor row (written for base ORIS) is the intended grounding source for `oris_soft` specifically, since the file covers both `RK_ALG_ORIS` and `RK_ALG_ORIS_SOFT` under one document (`docs/methods/oris.md:3` — "Enum: `RK_ALG_ORIS = 1` (+ variant `RK_ALG_ORIS_SOFT = 8`)") | Standard Stack, Architecture Patterns | Low — CONTEXT.md D-07 explicitly names this exact table for oris_soft, so this is confirming (not inventing) the mapping; risk only arises if a future doc split separates oris_soft into its own file before this phase lands. |

**If this table is empty:** N/A — see entries above; none of them are load-bearing facts
(all load-bearing facts in this document carry `[VERIFIED: ...]` tags with file:line and
verbatim quotes). These three are phase-3-specific implementation judgment calls flagged
for the planner/user, not unverified factual claims.

## Open Questions

1. **Does `tests/testthat/test-bench-gate.R`'s existing kk1204 test (lines 28-54) currently
   pass or fail?**
   - What we know: it asserts `best_error <= 1e-3` and `elapsed <= 30` at `n=500,000, K=20`,
     gated by `skip_if(Sys.getenv("CI") != "")` — i.e., it runs in ordinary local
     development. The 2026-04-23 investigation found `errRp≈1.15e-3` after 500 iterations
     at `n=1,000,000` (double this test's `n`) using default routing — at half the row
     count this test's `oris` (method="oris" explicitly, not AUTO) may or may not clear
     `1e-3` within its `max_iterations=500` budget.
   - What's unclear: this session did not execute the test (a 500k-row, 500-iteration
     `oris` solve is a real wall-clock cost not worth spending research budget on when the
     planner/executor will run it anyway as part of normal DoD verification).
   - Recommendation: the plan's Wave 0 (or an early task) should run this specific test in
     isolation (`Rscript -e "testthat::test_dir('tests/testthat', filter='bench-gate')"`,
     ensuring `CI` is unset) and record the actual result before deciding whether to keep,
     relabel, or retire it per Common Pitfall 2.

2. **Where exactly should the new `docs/` methodology page live, and does it need a
   References-style citation for `survey`/`icarus`/`ReGenesees` matching
   `docs/methods/oris.md`'s existing bibliography style?**
   - What we know: `docs/methods/oris.md` already has full academic citations for all three
     competitors (`[zardetto2015regenesees]`, `[rebecq2017icarus]`, `[lumley2010survey]`
     via the `docs/methods/references.bib` shared bibliography, `docs/methods/oris.md:216-226`).
   - What's unclear: whether the new performance-claim doc page should link back to
     `docs/methods/oris.md` for citations (avoiding duplication) or restate them.
   - Recommendation: link, don't duplicate — `docs/methods/oris.md` is the canonical
     citation source per D-07; the new page's job is methodology/numbers, not bibliography.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|-------------|-----------|---------|----------|
| R (`Rscript`) | Benchmark script, gate test | ✓ | — (project-standard toolchain, unchanged by this phase) | — |
| `bench` (R pkg) | `bench::mark()` timing | ✓ [VERIFIED: `Suggests:` in DESCRIPTION, read this session] | — | — |
| `survey` (R pkg) | Competitor #1 | ✓ [VERIFIED: `Suggests:` in DESCRIPTION + `requireNamespace()` TRUE this session] | — | — |
| `icarus` (R pkg) | Competitor #2 | ✓ [VERIFIED: `requireNamespace()` TRUE this session] — not yet in `DESCRIPTION` | — | — |
| `ReGenesees` (R pkg) | Competitor #3 | ✓ [VERIFIED: `requireNamespace()` TRUE this session] — not yet in `DESCRIPTION` | — | — |
| `arrow`, `jsonlite` (R pkgs) | Fixture I/O if a new large fixture is built | ✓ (already used by `benchmarks/regression_timer.R`, `stepstone_benchmark.R`) | — | — |

**Missing dependencies with no fallback:** none.
**Missing dependencies with fallback:** none — everything needed resolves locally today.

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | testthat edition 3 [VERIFIED: `DESCRIPTION:22` — "Config/testthat/edition: 3", read this session] |
| Config file | `DESCRIPTION` (testthat config lives here, no separate `tests/testthat.R`-level override found) |
| Quick run command | `Rscript -e "devtools::test()"` [VERIFIED: CLAUDE.md build/test section, read this session] |
| Full suite command | `R CMD INSTALL --preclean . && Rscript -e "devtools::test()"` (behavioral DoD, per CLAUDE.md) |
| Opt-in perf gate | `LBW_BENCH_GATE=1 Rscript -e "devtools::test()"` (or `test_dir` filtered to `bench-gate`) |

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|---------------------|--------------|
| US-003 / SC1 | Large-scale figure measured on a live solver, stated input class | integration (opt-in perf) | `LBW_BENCH_GATE=1 Rscript -e "testthat::test_dir('tests/testthat', filter='bench-gate')"` | ✅ file exists (`test-bench-gate.R`); ❌ new assertion for the `oris_soft` headline number — Wave 0 |
| US-003 / SC2 | Medium-scale target states ONE number with a named live artefact | integration (opt-in perf) or unit, depending on chosen scale | new — Wave 0 | ❌ Wave 0 |
| KPI-04 / SC3 | kk1204 known-limit documented | documentation-only (no automated test — a doc/reference check) | manual — cross-link verification | N/A — doc task, not a test |
| KPI-04 / SC4 | One command reproduces every published figure; KPI table names live artefacts | integration (shell wrapper) | new `benchmarks/run_honest_gate.sh` (or similar), plus `.planning/REQUIREMENTS.md` KPI-04 row edit | ❌ Wave 0 |

### Sampling Rate
- **Per task commit:** default DoD (`R CMD INSTALL --preclean .` + testthat, 0 FAIL) — the
  perf gate itself is NOT part of per-commit sampling, per D-11 (stays opt-in/heavy).
- **Per wave merge:** run `LBW_BENCH_GATE=1` once to confirm the new gate assertion is
  wired correctly (not necessarily to enforce a numeric threshold at every merge).
- **Phase gate:** full DoD green + `LBW_BENCH_GATE=1` run producing the actual published
  figure, captured into the docs page, before `/gsd-verify-work`.

### Wave 0 Gaps
- [ ] New assertion in `tests/testthat/test-bench-gate.R` (or a new adjacent file) for the
  `oris_soft` headline number, gated `LBW_BENCH_GATE=1` — covers SC1/SC2.
- [ ] `benchmarks/oris_soft_vs_competitors.R` — the actual measurement script — covers
  SC1/SC4.
- [ ] `benchmarks/run_honest_gate.sh` (or equivalent single entry point) — covers SC4.
- [ ] Resolve Open Question 1 (does the existing kk1204 test pass today) before deciding
  its fate.

## Security Domain

No new attack surface. This phase adds a benchmark script (developer-run, not
user-invoked at runtime), a README, and a docs page — no new user input parsing, no new
network calls, no new deserialization paths. `security_enforcement` is not explicitly
disabled in `.planning/config.json` (the file does not exist in this repo, so per the
governing instruction absence = enabled), but the ASVS categories below are effectively
N/A for this phase's actual deliverables.

### Applicable ASVS Categories

| ASVS Category | Applies | Standard Control |
|----------------|---------|-------------------|
| V2 Authentication | No | No auth surface touched. |
| V3 Session Management | No | N/A. |
| V4 Access Control | No | N/A. |
| V5 Input Validation | No | Benchmark script reads fixed, repo-local fixtures; no external/user input parsed. |
| V6 Cryptography | No | N/A. |

### Known Threat Patterns for this stack
None applicable — a benchmark script invoked manually by a maintainer, reading local
parquet/JSON fixtures already present in `benchmarks/`, is not a new threat surface.

## Sources

### Primary (HIGH confidence — read directly this session)
- `.planning/phases/03-honest-performance-gate/03-CONTEXT.md` — full file, locked decisions
- `.planning/REQUIREMENTS.md` — US-003, KPI-04, Superseded table, Traceability
- `.planning/STATE.md` — Blockers/Concerns, carried decisions
- `.planning/ROADMAP.md` — Phase 3 section (lines 151-182)
- `bd show leafblower-kk1.20.4` — REFRAME options, investigation reference
- `bd show leafblower-kk1.20` — Phase 2 parent, close reason
- `bd show leafblower-ylsy` — CLOSED, verbatim close reason (DEFF/n_eff quote)
- `bd show leafblower-2ouc` — epic scope, confirms 47/49 children complete, out of Phase 3 scope
- `tasks/prd-leafblower-core.md:30,673` — the §1/§11 1s/2s contradiction, quoted verbatim
- `docs/investigations/2026-04-23-kk1204-convergence.md` — full file, kk1204 rate/infeasibility findings
- `docs/investigations/2026-05-02-ylsy-cp-ipm-spike-result.md` — full file, ylsy closure evidence
- `docs/methods/oris.md` — competitor table (lines 216-226, 253-267), oris_soft description (lines 54-56)
- `tests/testthat/test-bench-gate.R` — full file, existing gate patterns (LBW_BENCH_GATE, CI-skip variant)
- `benchmarks/newton_kl_bench.R` — cross-solver `bench::mark()` loop pattern
- `benchmarks/regression_timer.R`, `benchmarks/stepstone_benchmark.R`, `benchmarks/run_allmethod.sh` — existing conventions
- `benchmarks/study/spec/contract.md` — `wall_time_s` measurement definition
- `DESCRIPTION` — full file, `Suggests:` list, testthat edition
- `NEWS.md` — solver rename history
- `CLAUDE.md` — determinism protocol, DoD commands, LBW_BENCH_GATE contract
- R session: `requireNamespace("survey"|"icarus"|"ReGenesees")` — all TRUE, run this session
- `git log --all --full-history -- README.md` — empty; confirms README.md has never existed

### Secondary (MEDIUM confidence)
None — no external web sources were needed; this phase's research surface is entirely
inside the repo.

### Tertiary (LOW confidence)
None.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — no new packages; all three competitors confirmed installed via
  direct `requireNamespace()` check, not inference.
- Architecture: HIGH — every pattern cited is an existing, read file in this repo.
- Pitfalls: HIGH — the 1s/2s contradiction and the kk1204 test's actual current assertions
  were both read directly, not recalled or assumed.

**Research date:** 2026-08-15
**Valid until:** No external dependency drift risk (nothing here tracks upstream releases);
re-check only if `tests/testthat/test-bench-gate.R` or `docs/methods/oris.md` change before
planning executes.
