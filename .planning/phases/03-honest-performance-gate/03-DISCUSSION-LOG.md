# Phase 3: Honest Performance Gate - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-08-15
**Phase:** 3-Honest Performance Gate
**Areas discussed:** REFRAME decision, Restated claim wording, kk1204 documentation, Benchmark command / maintainer workflow

---

## REFRAME decision

| Option | Description | Selected |
|--------|-------------|----------|
| Option 1: raking+SQUAREM | Accept raking+SQUAREM as the fast path (stepstone 1.46e-3, 1.6s) | |
| Option 2: gate on stepstone | Redefine the gate to measure on stepstone instead of kk1204 | |
| Option 3: track kk1204 separately | kk1204 splits off entirely as algorithmic research (leafblower-ylsy) | ✓ |

**User's choice:** Option 3 only.
**Notes:** First round of answers was self-contradictory (picked Option 3, then picked
"strictly Option 1 only" on the combine question). A direct clarifying question resolved
it to Option 3 only. Mid-discussion in the "kk1204 documentation" area, a beads lookup
revealed `leafblower-ylsy` is already CLOSED with a conclusion that kk1204 is a degenerate
fixture, not an algorithm gap — this meant "track separately" needs no new work, just
citation. See below.

---

## Restated claim wording

| Question | Options presented | Selected |
|---|---|---|
| Benchmark basis | stepstone fixture / new fixture / you decide | stepstone fixture |
| Gate type | hard gate / reported measurement | hard gate |
| Doc location (1st ask) | README+doc / README only / you decide | superseded by correction below |
| Stepstone gate scope | keep opt-in / promote to default | keep opt-in |
| Headline solver | oris_soft only / both, oris_soft primary / you decide | both, oris_soft primary |
| Competitor source | reuse study output / fresh narrow run / you decide | fresh narrow run |
| Competitor selection method | (free response) | doc-grounded (see correction) |
| Doc-grounded confirmation | yes / not quite | yes, generalized to R+Python |
| Doc location (re-asked) | README+doc / README only | README + linked doc |
| Headline metric | speed / DEFF-effective-obs / you decide | you decide (most dramatic number) |

**User's choice:** oris_soft is the primary headline solver (not autumn, not
raking+SQUAREM alone); competitors must come from each solver's `docs/methods/*.md`
doc-named comparison table, not picked freely from `benchmarks/study/`; a fresh narrow
benchmark run is needed (existing study output isn't scoped correctly); claim lives as a
README one-liner + linked doc page.

**Notes / corrections:**
1. **Correction 1:** When asked where the claim should live (README vs. linked doc), the
   user instead corrected the entire premise: autumn must never appear in a comparative
   position (never CRAN-released, not a real competitor). Real competitors are the ones
   in `benchmarks/study/`. oris_soft — not raking+SQUAREM — was flagged as the actual
   innovation to headline (speed, bounded weights, small design effect → more effective
   observations).
2. **Correction 2:** When asked whether to reuse existing study output or run fresh, and
   told the answer was "fresh," a follow-up on HOW to pick fresh competitors got another
   correction: don't pick freely from `benchmarks/study/` — each method's own doc already
   names its direct competitors (oris.md's comparison table: `survey::calibrate`,
   `icarus`, `ReGenesees`). Confirmed and generalized: every method doc names competitors
   in both R and Python.
3. The original "doc location" question never got a clean answer on the first ask (the
   slot was consumed by Correction 1) — re-asked afterward and resolved to README + linked
   doc.

---

## kk1204 documentation

| Question | Options presented | Selected |
|---|---|---|
| Where SC3's finding lives | investigations-only / investigations + method-doc note / you decide | investigations-only, linked from ylsy |
| ylsy ticket work | update/create as part of Phase 3 / out of scope | update/create — see correction below |

**User's choice:** Document in `docs/investigations/` only, cite ylsy's close reason.
**Notes:** The user's initial answer ("update/create leafblower-ylsy as part of Phase 3")
was given before a beads lookup surfaced that `leafblower-ylsy` is already CLOSED
(2026-05-03) with a firm conclusion: kk1204 is a near-infeasible degenerate fixture
(DEFF 8000-14000 across all solvers), not an algorithm gap — investigation report already
exists (`docs/investigations/2026-05-02-ylsy-cp-ipm-spike-result.md`). Surfaced this to the
user directly; they resolved it to "just cite it, no ticket work needed."

---

## Benchmark command / maintainer workflow

| Question | Options presented | Selected |
|---|---|---|
| One command or two | one wrapper command / two separate commands | one wrapper command |
| Competitor package deps | benchmark-only deps / you decide | benchmark-only deps, scoped under benchmarks/ |

**User's choice:** One wrapper command covering stepstone + the new oris_soft comparison;
new R/Python competitor packages added as benchmark-scoped dependencies only (not the
shipped package's own DESCRIPTION/pyproject).
**Notes:** None beyond the selections above.

---

## Claude's Discretion

- Which metric (speed vs. design-effect/effective-observations) headlines the restated
  claim — pick the stronger, more honest number once the fresh benchmark actually runs.
- Exact shape/location of the combined wrapper command.
- Exact README wording and linked-doc structure, within the "one line + linked page" shape.

## Deferred Ideas

- Full comparative benchmark study against every package in `benchmarks/study/` —
  belongs to the `leafblower-2ouc` epic, not Phase 3.
