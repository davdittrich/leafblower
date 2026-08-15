# Phase 3: Honest Performance Gate - Context

**Gathered:** 2026-08-15
**Status:** Ready for planning

<domain>
## Phase Boundary

Restate leafblower's headline performance claim as a measured fact about a solver that
exists, expressed as a gate a user can run and a maintainer can regress against. This
phase does NOT include the full comparative benchmark study against all packages in
`benchmarks/study/` (that's the separate `leafblower-2ouc` epic) — it needs a narrow,
targeted benchmark scoped to exactly the numbers SC1/SC2 require.

</domain>

<decisions>
## Implementation Decisions

### REFRAME scope — kk1204 split-off
- **D-01:** The `<30s AND <1e-6>` gate on kk1204 (K=20 uniform-random, M_cell/n=1.0) is
  confirmed structurally unachievable and is dropped as the basis for any headline claim.
  kk1204 is NOT part of Phase 3's restated performance number in any form.
- **D-02:** **Correction mid-discussion:** kk1204 is not "ongoing algorithmic research" to
  track going forward — `leafblower-ylsy` (the ticket meant to carry that research) is
  already CLOSED (2026-05-03) with a firm conclusion: kk1204 is a near-infeasible
  **degenerate fixture** (DEFF 8000-14000, n_eff 71-118 across ALL solvers — ieppa,
  newton_kl, auto), not an algorithm gap. See
  `docs/investigations/2026-05-02-ylsy-cp-ipm-spike-result.md`. — **Reversibility:**
  reversible — this is a documentation framing choice, not code; if new evidence emerges
  the ticket can be reopened.
- **D-03:** No beads work needed on `leafblower-ylsy` — just cite its close reason and the
  2026-05-02 report as the SC3 "known limit" explanation. Do not reopen or comment on it.

### Headline claim — solver and metric
- **D-04:** `oris_soft` is the primary headline solver, NOT `autumn` and NOT
  raking+SQUAREM alone. raking+SQUAREM stays cited as the fallback path for
  zero-compression-benefit inputs (K=20-uniform-random-style), not as the flagship number.
- **D-05:** **User correction — do not mention `autumn` in any comparative/relevant
  position anywhere in this phase's output.** Autumn was never released to CRAN and is
  not a real competitor; using it as the baseline ("Nx faster than autumn") is dishonest
  framing that this phase exists specifically to eliminate.
- **D-06:** Headline metric (speed/wall-time vs. design-effect/effective-observations) is
  Claude's discretion — pick whichever number is the strongest, most honest differentiator
  once the fresh benchmark actually runs. Both framings are legitimate; don't force one.

### Competitor selection — doc-grounded, not ad hoc
- **D-07:** **User correction — competitors are NOT picked freely from `benchmarks/study/`.**
  Each solver's `docs/methods/*.md` file already names its direct competitors in both R and
  Python. For `oris_soft`, pull from `docs/methods/oris.md`'s existing "How leafblower
  deviates" comparison table (bounds-handling row: `survey::calibrate`, `icarus`,
  `ReGenesees` — bounds folded into every iteration via bisection/logit, vs. oris_soft's
  ALM/ADMM-inside-the-loop bounds). Research/planning must locate and use the doc-named
  competitor set for whichever solver ends up headlining, not invent one.
- **D-08:** The fresh benchmark run reuses `benchmarks/` infrastructure conventions
  (per PROJECT.md's explicit "reuse `benchmarks/`, do not build a parallel harness"
  instruction) but is a NEW, narrow run — not a read of existing
  `benchmarks/study/report/tables/*.csv` output, since that existing study wasn't scoped
  to oris_soft specifically and mixes in packages beyond the doc-named competitor set.
- **D-09:** New R/Python competitor packages (`survey`, `icarus`, `ReGenesees`, etc.)
  needed for the fresh run are added as **benchmark-only dependencies**, scoped under
  `benchmarks/` — NOT added to the package's own `DESCRIPTION` or `pyproject.toml`.
  — **Reversibility:** reversible — benchmark-scoped deps don't affect the shipped
  package's dependency footprint.

### Claim mechanics
- **D-10:** The restated performance number is a **hard gate** (pass/fail, regressable),
  matching the existing stepstone `LBW_BENCH_GATE=1` pattern — not a reported measurement
  with no threshold.
- **D-11:** The stepstone gate stays **opt-in** (`LBW_BENCH_GATE=1`), matching the existing
  CLAUDE.md Definition-of-Done split (behavioral gate is default; stepstone regression is
  opt-in/heavy/local-only). Phase 3 does not change that contract.
- **D-12:** SC4's "one command" requirement is satisfied by **one wrapper command** that
  runs both the existing stepstone gate and the new oris_soft-vs-competitors comparison —
  not two separately-documented commands.
- **D-13:** Claim location: **README headline (one line) + linked docs/ page** with full
  methodology, machine spec, input class, and competitor detail. Not everything inlined
  into README.

### kk1204 documentation (SC3)
- **D-14:** SC3's "known limit" documentation stays in
  `docs/investigations/2026-04-23-kk1204-convergence.md` (already exists), cross-linked
  from the `leafblower-ylsy` close reason and the 2026-05-02 follow-up report. No new
  user-facing mention in `docs/methods/*.md` is needed — kk1204 is not part of the shipped
  claim, so it doesn't need to preempt a user's expectation there.

### Claude's Discretion
- Which metric (speed vs. design-effect/ESS) headlines the restated claim — pick the
  stronger, more honest number once the fresh benchmark runs (D-06).
- Exact wrapper-command shape/location for the combined benchmark entry point (D-12).
- Exact README wording and linked-doc structure, within the "one line + linked page"
  shape (D-13).

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### REFRAME decision and kk1204
- `leafblower-kk1.20.4` — the REFRAME beads ticket with the three original options (this
  phase resolved on Option 3 only, per D-01/D-02)
- `docs/investigations/2026-04-23-kk1204-convergence.md` — original kk1204 slow-rate
  finding, commit `3effd3a`
- `leafblower-ylsy` (CLOSED) — kk1204 degeneracy conclusion; cite the close reason verbatim
- `docs/investigations/2026-05-02-ylsy-cp-ipm-spike-result.md` — the investigation that
  closed `leafblower-ylsy`, concluding kk1204 is a degenerate fixture not an algorithm gap

### Competitor grounding and comparative infra
- `docs/methods/oris.md` §"How leafblower deviates" — the doc-named competitor comparison
  table (`survey::calibrate`, `icarus`, `ReGenesees` for bounds handling) that MUST ground
  the fresh benchmark's competitor selection; do not pick competitors from elsewhere
  without checking the relevant `docs/methods/*.md` first
- `docs/methods/00-overview.md` and the other `docs/methods/*.md` files — check for the
  doc-named competitor set for whichever solver ends up headlining, if not oris_soft
- `benchmarks/study/report/tables/quality.csv`,
  `benchmarks/study/report/tables/rss_by_language_package.csv` — existing comparative
  study output; informative for context but NOT the source of Phase 3's fresh benchmark
  numbers (per D-08)
- `study/benchmark-instrumented-DO-NOT-MERGE` (git branch) — where the existing
  comparative study infrastructure lives; reuse its conventions, do not duplicate
- `leafblower-2ouc` — the separate epic that owns the FULL comparative study; Phase 3's
  fresh run must stay narrow and not encroach on this epic's scope

### Existing gate and roadmap context
- `.planning/ROADMAP.md` §"Phase 3: Honest Performance Gate" — SC1-SC4, requirements
  (US-003, KPI-04), beads references
- `CLAUDE.md` — Definition of Done split (behavioral gate default, `LBW_BENCH_GATE=1`
  opt-in); single-thread BLAS determinism protocol; "no cancellations" numeric philosophy
- `benchmarks/` (existing stepstone fixtures, `regression_timer.R`, `parity_bench.R`) —
  reuse per PROJECT.md's explicit instruction, do not build a parallel harness

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `benchmarks/` — existing stepstone fixtures and timer/parity harness; the new
  oris_soft-vs-competitors run should follow the same conventions (single-thread BLAS,
  interleaved before/after `bench::mark()` calls per CLAUDE.md's determinism protocol)
- `benchmarks/study/` — has a working pattern for invoking real R/Python competitor
  packages (`benchmarks/study/python/leafblower_adapter.py` and friends); worth reading
  even though its aggregate output isn't the source of Phase 3's numbers

### Established Patterns
- `LBW_BENCH_GATE=1` env-var gate pattern (stepstone) — the new combined wrapper command
  should extend this pattern, not replace it
- CLAUDE.md's benchmark determinism protocol: single-thread BLAS, interleaved
  before/after comparison in one `bench::mark()` call, never sequential measurement

### Integration Points
- README.md — gets the one-line headline claim
- A new linked docs/ page (exact location is Claude's discretion, D-13) — gets full
  methodology/machine/competitor detail
- `docs/methods/oris.md` — already has the competitor comparison table to pull from; no
  edit needed unless the headline metric changes what's cited

</code_context>

<specifics>
## Specific Ideas

- The user was explicit and emphatic that autumn must not appear anywhere as a comparison
  baseline — it was never CRAN-released and using it as "the" competitor is exactly the
  dishonest framing this phase is meant to fix.
- The user framed oris_soft's actual innovation precisely: "high speed, bounded weights,
  small design effect -> more effective observations" — this phrasing should inform the
  headline claim's actual content, not just its metric choice.
- The user's repeated correction pattern (twice) was: don't invent a comparison basis —
  find what's already documented per-method and use that. This is a general instruction
  for how downstream research should approach any solver-vs-competitor claim in this
  phase, not just the oris_soft case.

</specifics>

<deferred>
## Deferred Ideas

- Full comparative benchmark study against every package in `benchmarks/study/` —
  belongs to the `leafblower-2ouc` epic (an article-length deliverable), not Phase 3.
  Phase 3's fresh benchmark is a narrow subset (1-2 doc-named competitors for the
  headlining solver), not a duplicate of that epic's scope.

### Reviewed Todos (not folded)
None — discussion stayed within phase scope.

</deferred>

---

*Phase: 3-Honest Performance Gate*
*Context gathered: 2026-08-15*
