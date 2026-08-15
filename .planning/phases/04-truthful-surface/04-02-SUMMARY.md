---
phase: 04-truthful-surface
plan: 02
subsystem: docs
tags: [documentation, truthfulness, enum-annotation]
status: complete
requirements: []
dependency-graph:
  requires: ["04-01"]
  provides: ["SC1", "SC2", "SC4"]
  affects: ["docs/raking.md", "src/leafblower.h", "CLAUDE.md"]
tech-stack:
  added: []
  patterns:
    - "Exact-substring boundary deletion for prose-only docs with no markdown headings (docs/raking.md is 10 physical lines / 43KB — section markers are inline prose, edited via byte-offset-anchored substring find, not line numbers)"
    - "Enum-slot removal annotation convention: `/* N = removed (was RK_ALG_X) */` comment inserted at the gap, matching the pre-existing slot-2 pattern"
key-files:
  created: []
  modified:
    - docs/raking.md
    - src/leafblower.h
    - CLAUDE.md
decisions:
  - "Deleted the false ORIS/L-BFGS-B capability claims outright rather than replacing with a correct description — docs/methods/oris.md:127 is already the single authoritative description (D-01); duplicating it into raking.md would create a second source of truth to drift"
  - "Left the 9 non-scoped L-BFGS-B mentions in raking.md's §6.2/§11 untouched — general classical-numerical-technique discussion, not an ORIS/leafblower capability claim, explicitly out of D-01's narrow §8.2/§12 scope"
metrics:
  duration: "~15min"
  completed: 2026-08-15
actuals:
  tokens: 9500
  tasks: 2
  commits: 2
---

# Phase 04 Plan 02: Truthful Surface Cleanup (SC1/SC2/SC4) Summary

Deleted docs/raking.md's ORIS/L-BFGS-B misattribution passages and annotated the
rk_algorithm_t enum's slot-7 gap, closing the phase's three lower-risk documentation
success criteria with zero code-behavior changes.

## What Was Built

**Task 1 (SC1):** Deleted two passages from `docs/raking.md` by exact-substring boundary
match — the file has zero markdown headings (10 physical lines / ~43KB of inline prose),
so line-number edits were not viable; located and deleted via a Python script anchored
on the exact section-title/closing-sentence substrings from the plan (each confirmed
`count==1` before deletion to guard against accidental duplicate matches):

- §8.2 ("The ORIS Solver (renamed from iEPPA)") — falsely attributed the paper's headline
  contribution (the unimplemented outer inexact entropic-proximal-point loop) to ORIS, and
  claimed ORIS "routinely and conclusively outperforms ... Gurobi".
- §12 ("Synthesis and Final Conclusions") — called L-BFGS-B "the definitive
  state-of-the-art" and recommended offloading L-BFGS-B/ORIS numerical routines, both
  false for this package (L-BFGS-B was removed; ORIS's inner solver is not L-BFGS-B).

No replacement prose was added — the deletion is the fix, per D-01 (the ground-truth
description already lives at `docs/methods/oris.md:127`). §8.1 and §9–§11 are byte-for-byte
unchanged; the file now reads directly from "...strictly respecting the capacity
constraints." into "9. Real-World Architectures...".

**Task 2 (SC2 + SC4):** Inserted `/* 7 = removed (was RK_ALG_GRAKE) */` in
`src/leafblower.h`'s `rk_algorithm_t` enum, between `RK_ALG_GREG = 6` and
`RK_ALG_ORIS_SOFT = 8`, matching the existing slot-2 comment convention exactly. No enum
value changed. Extended `CLAUDE.md`'s existing slot-2 sentence to name both reserved slots
in one line.

Re-ran the SC4 grake/lbfgsb/cp audit sweep across the seven in-scope globs (`README.md`,
`NEWS.md`, `man/*.Rd`, `docs/methods/*.md`, `docs/raking.md`, `docs/performance.md`;
`README.Rmd` does not exist in this repo). All hits matched one of the six pre-classified
false-positive classes from 04-RESEARCH.md's sweep: `grake_norm` (live convergence-metric
field name), `survey::grake` (comparative reference to the R `survey` package's own
function), `L-BFGS-B`/`LBFGSB` in `docs/methods/00-overview.md` (correctly documents
removed status) and in `docs/raking.md`'s untouched §6.2/§11 (9 mentions, exactly matching
the plan's predicted count, confirming §8.2/§12 deletion did not disturb them), the `cp`
hit in `docs/performance.md:176` (the withdrawn `cp`-ipm research-spike report). No new
true-positive references found; no additional edits required.

## Verification

- `docs/raking.md`: all five acceptance greps (`renamed from iEPPA`, `Gurobi`,
  `definitive state-of-the-art`, `8.2 The ORIS Solver`, `12. Synthesis and Final
  Conclusions`) return 0; `strictly respecting the capacity constraints` and `external
  hardware accelerators` both still return 1 (§8.1/§11 untouched).
- `src/leafblower.h`: `git diff` shows exactly one added comment line, zero enum value
  changes; `grep -c "RK_ALG_GRAKE"` returns 1 (comment only).
- `CLAUDE.md`: both slot 2 and slot 7 named in one sentence.
- `R CMD INSTALL --preclean .`: succeeded (leafblower.h edit forces a full rebuild, per
  CLAUDE.md's "Two build sites" note).
- `Rscript -e "devtools::test()"`: **0 FAIL / 1837 PASS / 141 WARN / 13 SKIP**.
- Python parity suite not re-run per-task (no `.py`/`.cpp` file changed this plan, per the
  plan's own acceptance criteria); deferred to the phase-level DoD gate before
  `/gsd:verify-work`.

## Deviations from Plan

None — plan executed exactly as written. Both tasks matched their `<action>` specs
precisely; the SC4 re-audit found zero new true-positive hits, so no additional file
edits beyond the two named in the plan were needed.

## Known Stubs

None.

## Threat Flags

None — this plan's threat register (T-04-03, T-04-04) both dispositioned in the plan
itself (mitigate via deletion; accept, no security surface). No new surface introduced.

## Self-Check: PASSED

- FOUND: docs/raking.md (modified, verified via grep)
- FOUND: src/leafblower.h (modified, verified via git diff)
- FOUND: CLAUDE.md (modified, verified via grep)
- FOUND: commit 3143594 (docs(raking): delete ORIS/L-BFGS-B misattribution)
- FOUND: commit 2a06f1f (docs(leafblower.h,CLAUDE.md): annotate rk_algorithm_t slot 7)
