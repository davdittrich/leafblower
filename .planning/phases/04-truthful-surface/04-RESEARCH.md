# Phase 4: Truthful Surface - Research

**Researched:** 2026-08-15
**Domain:** Documentation/API truthfulness audit — no new architecture, no new dependencies
**Confidence:** HIGH (every finding below is file:line-verified this session; no package research,
no external ecosystem research applies to this phase)

## Summary

This is a four-item defect-fix phase, not a build phase. All four success criteria were traced
to exact source locations this session. Three of the four (docs/raking.md misattribution, the
enum slot-7 hole, the `weights=` silent-swallow) are narrow, mechanical, single-file-or-two-file
edits with pre-existing style templates to mirror (`docs/methods/oris.md` for the correct ORIS
description; `src/leafblower.h`'s slot-2 comment for slot 7; `R/harvest.R`'s own `design_weights`
warning idiom for the new `stop()` message). The fourth (grake/lbfgsb/cp audit) turned out to
have **zero true positives** beyond the two already-ticketed items — the audit surfaced several
near-miss strings that must NOT be touched (`grake_norm`, `survey::grake` in a comparative-methods
doc, a `cp` file-path reference to a withdrawn research spike) and confirmed `README.md` and
`NEWS.md` are already clean.

The single most consequential finding, not named in the phase notes, is that **the `weights=`
guard (D-04) will break one existing test on implementation**: `tests/testthat/test-logit.R:195`
(`eb79.18` "HETEROGENEOUS design weights") passes `weights = base_w` to `harvest()` today, wrapped
in `suppressWarnings()` — which currently hides the RVAL.2 generic-unknown-argument warning and
silently drops the vector, so the test does NOT exercise heterogeneous design weights at all
despite its own name. Converting the guard from warn to `stop()` will make this test **fail**
outright unless it is fixed in the same commit. This must be a task, not an afterthought — see
Pitfall 1.

**Primary recommendation:** Four independent, narrow-scope tasks — (1) delete the docs/raking.md
§8.2/§12 passage per D-01 (delete, not rewrite); (2) add the slot-7 comment to `src/leafblower.h`
and CLAUDE.md; (3) add a `stop()` guard for bare `weights=` in `R/harvest.R`'s existing RVAL.2
block AND fix `test-logit.R:195` to use `design_weights=` in the same task; (4) close
`leafblower-05ha`/`leafblower-x2iq` and file the `weights=` ticket — no README/NEWS.md changes
needed (already clean), no other grake/lbfgsb/cp true positives found.

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| docs/raking.md content correction | Documentation | — | Pure prose edit, no code path |
| `rk_algorithm_t` slot-7 annotation | C header (shared R+Python ABI) | Documentation (CLAUDE.md) | Header is the single source both language bridges compile against |
| `harvest(weights=)` guard | R API surface (`R/harvest.R`) | Test suite (`tests/testthat/test-logit.R`) | R-level argument validation; guard and its one affected caller must land together |
| grake/lbfgsb/cp doc audit | Documentation | — | Prose-only; confirmed zero code-path risk |

## User Constraints

This phase has a CONTEXT.md; the planner MUST honor the decisions below verbatim.

<user_constraints>
### Locked Decisions

- **D-01:** Delete the docs/raking.md §8.2/§12 passage entirely rather than rewrite it under
  the paper's own name. Do not leave an ORIS-labelled description of the unimplemented outer
  entropic-proximal-point loop anywhere in the file. Ticket: `leafblower-05ha`.
- **D-02:** File a follow-up ticket (separate from this phase) noting the rest of
  `docs/raking.md` (§1-7, §9-11) lacks citations and includes off-topic content. Do NOT act on
  this in Phase 4 — flagged for a maintainer look only.
- **D-03:** The README/NEWS.md/man/docs audit (SC4) covers `README.md`, `NEWS.md`, `man/`,
  `docs/methods/`, `docs/raking.md`, `docs/performance.md`, and any vignettes.
  `docs/superpowers/specs/` is explicitly OUT of scope (historical planning records).
- **D-04:** `harvest(..., weights = w)` must raise a hard `stop()` naming `design_weights=`, not
  a warning. Reversibility: costly — a `checkpoint:decision` is required before the implementing
  task, per phase 4's own truthfulness goal. The researcher/planner must confirm no in-repo
  test or benchmark currently passes `weights=` expecting the old silent-fallback behavior
  before this ships. **Research finding: one test does — see Pitfall 1. It must be fixed in the
  same task, not merely flagged.**

### Claude's Discretion

- Exact wording of the `stop()` message for `weights=` (beyond naming `design_weights=`
  explicitly) — planner/executor's call. See Code Examples for a template matching the
  package's existing idiom.
- How to document the slot-7 hole in `rk_algorithm_t` — match slot 2's existing comment style
  (`/* N = removed (was RK_ALG_X) */`). Research confirmed the reason: slot 7 was
  `RK_ALG_GRAKE`, dropped in commit `9a67891` ("remove(grake): drop RK_ALG_GRAKE=7 enum value —
  no release, no ABI constraint").

### Deferred Ideas (OUT OF SCOPE)

- **docs/raking.md broader quality review** (D-02) — file as a new ticket; not Phase 4.
</user_constraints>

## Phase Requirements

<phase_requirements>
No PRD-derived requirement IDs apply — this phase is defect-driven (roadmap states
"Requirements: (none — defect-driven)"). The four Success Criteria in the phase description
function as the requirement set; each is mapped to research support below.

| SC | Description | Research Support |
|----|-------------|------------------|
| SC1 | docs/raking.md §8.2/§12 misattribution removed | Exact passage text extracted verbatim below (Section "docs/raking.md misattribution"); `docs/methods/oris.md:127` gives the ground-truth correction language already in the repo |
| SC2 | `rk_algorithm_t` slot 7 documented like slot 2 | `src/leafblower.h:41-53` read in full; slot 2's exact comment text confirmed; slot 7's gap location confirmed |
| SC3 | `harvest(weights=)` raises informative error | `R/harvest.R:266-307` (signature + RVAL.2 block) and `:510-515` (design_weights warning idiom) read in full; one breaking test identified |
| SC4 | grake/lbfgsb/cp audit of README/NEWS/man/docs | Full grep sweep executed; results and false-positive exclusions tabulated below |
</phase_requirements>

## Project Constraints (from CLAUDE.md)

- R build gate: `R CMD INSTALL --preclean .` (NOT `devtools::install`). R tests:
  `Rscript -e "devtools::test()"`.
- Definition of Done: R tests pass + Python parity tests pass + stepstone benchmark shows no
  regression (stepstone/`LBW_BENCH_GATE=1` gate is opt-in/heavy — not required for a docs+guard
  phase unless the executor judges the header touch warrants it; the slot-7 header edit forces
  a full rebuild of both R and Python extension per `leafblower-x2iq`'s own ticket body).
- Commit with explicit pathspec, NEVER `git add -A` (a bd/graphify hook re-stages
  `.beads/issues.jsonl`).
- No git remote — local-only. Complete = committed locally + gates green.
- `CLAUDE.md` itself states: "Algorithm slot 2 is reserved (LBFGSB removed). Do not reuse in
  `rk_algorithm_t` enum." — this line needs updating to also name slot 7, per the phase's own
  canonical_refs and D-04's discretion note.
- Beads-only task tracking — no TodoWrite/TaskCreate/markdown TODO lists.
- `bd` warns `.beads has permissions 0755 (recommended: 0700)` on every invocation in this repo
  — pre-existing, unrelated to this phase, not an action item here.

## Findings by Success Criterion

### SC1 — docs/raking.md §8.2/§12 misattribution

**File shape** [VERIFIED: docs/raking.md, `wc -l -c` this session]: the file is 10 physical
lines / 43,826 bytes with zero markdown headings — section numbers like "8.2" and "12" are
inline prose markers, not `#`-headers. Line-based edits will be coarse; the executor must
locate passages by quoted text, not line number (confirmed independently of the ticket body).

**False claim block ("8.2 The ORIS Solver")** [VERIFIED: docs/raking.md, extracted verbatim
this session via `grep -o`]:

> "8.2 The ORIS Solver (renamed from iEPPA)To solve high-dimensional CCOT problems with
> guaranteed stability, contemporary research points to the inexact Entropic Proximal Point
> Algorithm — now renamed ORIS. ... It achieves this by wrapping a proximal point framework
> around an inner solver. The algorithm executes the following sequence:Entropic Proximal
> Step: ... Dual Block Coordinate Descent Subsolver: ... Inexact Stopping Criteria: ...
> Extensive benchmarking on massive capacity-constrained Multi-Marginal Optimal Transport
> problems demonstrates that ORIS routinely and conclusively outperforms both stabilized
> Dykstra's algorithms and premium commercial linear programming solvers like Gurobi."

Counts confirmed this session: `L-BFGS-B` occurs **12** times in the file (not 4 as the ticket
body states — the ticket's count of 4 undercounts; re-verify before the executor's before/after
grep, do not trust the ticket's cached number), `Gurobi` occurs 1 time, `renamed ORIS` occurs 1
time.

**False claim block ("12. Synthesis and Final Conclusions")** [VERIFIED: docs/raking.md,
extracted verbatim this session]:

> "12. Synthesis and Final Conclusions ... For standard, high-dimensional survey tasks
> requiring strict weight clamping, the L-BFGS-B algorithm stands as the definitive
> state-of-the-art. ... ORIS (renamed from iEPPA), armed with a dual block coordinate descent
> subsolver, delivers unmatched computational resilience. ... offloading the L-BFGS-B or ORIS
> numerical routines into statically typed, compiled languages like Rust or C++ ..."

**Ground-truth correction, already live in the repo** [VERIFIED: docs/methods/oris.md:127, read
in full this session]:

> "**Not the iEPPA paper.** The solver was originally (mis)named "iEPPA" after Chu, Liang, Toh &
> Yang, *"An efficient implementable inexact entropic proximal point algorithm…"* (arXiv:2011.14312)
> [chu2022ieppa]. That paper's **headline contribution — an outer inexact entropic-proximal-point
> loop with a re-centered Bregman term to solve a transport-cost LP (`min⟨C,X⟩`) without `ε → 0`
> blow-up — is NOT implemented here.** Survey calibration has no transport cost (`C = 0`), at
> which the outer proximal loop is mathematically inert ... So citing arXiv:2011.14312 as the
> basis over-claims a contribution the code does not contain."

Per D-01, docs/raking.md must NOT duplicate this passage — delete the false §8.2/§12 blocks
outright; `docs/methods/oris.md` remains the sole authoritative description. Rewrite the "12.
Synthesis" passage only enough to stop recommending L-BFGS-B and stop repeating the ORIS
outer-loop claim (per ticket `leafblower-05ha` step 5) — the synthesis section's other content
(CPU/hardware discussion) is out of scope per D-02.

**Rename history for context** [VERIFIED: `git log --oneline`, this session]: `77d0614
refactor(oris)!: rename ieppa C++ core to oris (enum values frozen)` — confirms enum value 1
(`RK_ALG_ORIS`) is unchanged across the rename; only the string/doc name moved.

### SC2 — `rk_algorithm_t` slot 7

**Current enum, verbatim** [VERIFIED: src/leafblower.h:41-53, read in full this session]:

```c
typedef enum {
    RK_ALG_AUTO   = 0,
    RK_ALG_ORIS  = 1,
    /* 2 = removed (was RK_ALG_LBFGSB) */
    RK_ALG_RAKING    = 3,
    RK_ALG_SINKHORN  = 4,
    RK_ALG_CHEBYSHEV = 5,
    RK_ALG_GREG      = 6,
    RK_ALG_ORIS_SOFT = 8,   /* oris + ADMM soft capacity enforcement */
    RK_ALG_GREENKHORN = 9,   /* greedy coordinate-descent IPF (autumn::harvest style) */
    RK_ALG_LOGIT      = 10,  /* Deville-Sarndal 1992 logit Newton calibration (autumn::calibrate style) */
    RK_ALG_NEWTON_KL  = 11   /* Newton-KL smooth dual calibration (zero-compression regime) */
} rk_algorithm_t;
```

Gap confirmed: nothing between `RK_ALG_GREG = 6,` (line 48) and `RK_ALG_ORIS_SOFT = 8,` (line
49) — no comment marking the value-7 hole, unlike slot 2's inline `/* 2 = removed (was
RK_ALG_LBFGSB) */` at line 44.

**Reason for the hole, from git history** [VERIFIED: `git show 9a67891 --stat`, this session]:

```
commit 9a67891410251de3da1d3b5b1ad30492633680fe
    remove(grake): drop RK_ALG_GRAKE=7 enum value — no release, no ABI constraint

    Pre-release codebase; no ABI compatibility obligation. Remove the
    deprecated sentinel entirely rather than leaving dead enum values.

 .beads/issues.jsonl | 44 ++++++++++++++++++++++----------------------
 src/leafblower.h    |  1 -
 2 files changed, 22 insertions(+), 23 deletions(-)
```

Slot 7 was `RK_ALG_GRAKE = 7`. The fix is a one-line addition mirroring slot 2's exact wording:
`/* 7 = removed (was RK_ALG_GRAKE) */` inserted between the `RK_ALG_GREG` and `RK_ALG_ORIS_SOFT`
lines. `docs/methods/00-overview.md:81` [VERIFIED, read this session] already contains a
matching comment for slot 2 in a doc-code excerpt (`// 2 = removed (was LBFGSB)`) — no separate
doc file needs a parallel slot-7 line beyond CLAUDE.md, since 00-overview.md's excerpt is a
direct quote of the header and will pick up the fix if/when re-synced (out of scope to force
that resync in this phase; not a stated success criterion).

**CLAUDE.md update needed** [VERIFIED: CLAUDE.md's own text, this session's file read]:
current line reads "Algorithm slot 2 is reserved (LBFGSB removed). Do not reuse in
`rk_algorithm_t` enum." — must be extended to also name slot 7 (GRAKE), per D-04's discretion
note and the phase's own canonical_refs.

**Build-invalidation note** (from ticket `leafblower-x2iq`, cross-checked against CLAUDE.md's
own "Two build sites" footgun): editing `src/leafblower.h` forces a full rebuild of BOTH R
(`R CMD INSTALL --preclean .`) and Python (`cd python && uv pip install -e . --reinstall-package
leafblower`) extensions even though only a comment changes — the header is included by both
build sites. Both test suites must be re-run under single-thread BLAS per the project's
parity-determinism rule.

### SC3 — `harvest(weights=)` guard

**Signature and RVAL.2 block, verbatim** [VERIFIED: R/harvest.R:266-307, read in full this
session]:

```r
harvest <- function(
  data,
  target,
  min_weight       = 0,
  max_weight       = 5,
  capacity_penalty = NULL,
  alm_penalty      = NULL,
  method           = "oris",
  verbose          = 0,
  max_iterations   = 500,
  start_weights    = NULL,
  attach_weights   = TRUE,
  weight_column    = "weights",
  convergence      = list(),
  sor              = NULL,
  bounds_mode      = "cell",
  homotopy_levels       = 1L,
  homotopy_start_factor = 1.0,
  homotopy_end_factor   = 1.0,
  homotopy_budget_p     = 0.5,
  scheduler             = c("round_robin", "greedy"),
  eta_schedule          = c("fixed", "tang_dynamic"),
  eta_start             = 1.0,
  eta_end               = 1.0,
  eta_schedule_power    = 0.5,
  accelerate       = FALSE,
  add_na_proportion = FALSE,
  auto_collapse    = FALSE,
  collapse_vars    = NULL,
  target_map       = NULL,
  design_weights   = NULL,
  newton_tsvd_ratio = 1e-8,
  ridge_lambda = 0.0,
  ...
) {
  # RVAL.2: warn on unknown ... args (typos / removed params)
  dots <- list(...)
  if (length(dots) > 0L)
    warning("harvest: unknown argument(s) ignored: ",
            paste(names(dots), collapse = ", "), call. = FALSE)
```

There is no `weights=` formal — a bare `weights=w` call falls into `...` (`dots`) and triggers
only the generic RVAL.2 warning naming it among "unknown argument(s) ignored". This is the
exact silent-fallback D-04 describes: the call proceeds, `w` is discarded, `harvest()` returns
a plausible unweighted result.

**Existing, in-repo idiom to model the new `stop()` on** [VERIFIED: R/harvest.R:510-515, read
this session]:

```r
  # design_weights: used as start_weights when supplied (normalized to mean=1 by normalize_start_weights)
  if (!is.null(design_weights)) {
    if (!is.null(start_weights))
      warning("leafblower: both design_weights and start_weights supplied; design_weights ignored")
    else
      start_weights <- design_weights
  }
```

**Package-wide error idiom** [VERIFIED: DESCRIPTION `Imports: stats` only — no `cli` dependency
— and R/harvest.R:380-508, multiple `stop()` calls read this session]: every existing `stop()`
in `R/harvest.R` uses **base R** `stop("leafblower: <message>", call. = FALSE)` (e.g. line 508:
`stop("leafblower: 'data' must be a non-empty data.frame", call. = FALSE)`; lines 447, 451, 455
similarly). The `cli` package is NOT an installed/imported dependency of this package — do not
introduce `cli::cli_abort()`; use base `stop()` with the `"leafblower: "` prefix and
`call. = FALSE` to match every sibling error in this file.

**Insertion point**: the RVAL.2 block (lines 303-307) is the natural, single insertion point —
intercept `"weights" %in% names(dots)` there (before or as part of the generic warning),
`stop()` with the specific message, and let the generic RVAL.2 warning continue to cover any
other unrecognized `...` names. Do not add a parallel/duplicate mechanism.

### SC4 — grake/lbfgsb/cp audit

**Full-repo grep executed this session** across `README.md`, `README.Rmd`, `NEWS.md`,
`man/*.Rd`, `docs/methods/*.md`, `docs/raking.md`, `docs/performance.md`, and (checked, none
exist) `vignettes/`.

**True positives requiring no further action beyond SC1/SC2 (already covered above):**
- `docs/raking.md` — the SC1 passage (L-BFGS-B ×12, Gurobi ×1) — covered by D-01.

**False positives — DO NOT TOUCH, confirmed by reading source:**

| Hit | File:Line | Why it's NOT a violation |
|-----|-----------|---------------------------|
| `grake_norm` | `man/harvest.Rd:271` ("`grake_norm`: survey-grake normalized residual max_k\|misfit\|/(1+\|pop\|)") | Live convergence metric field, confirmed `[VERIFIED: src/greg.cpp:153]` `res.base.grake_norm = m.grake_norm;` and `[VERIFIED: src/logit_calib.cpp:558]` `res.base.grake_norm = abs_grake;` — a different object from the removed `RK_ALG_GRAKE` solver, correct as documented. |
| `grake` (survey pkg internal fn) | `docs/methods/logit.md:112,167,175,179` | Describes the R **`survey`** package's own internal `grake` function as a comparative reference implementation ("The Newton–Raphson solver is exposed as the internal `grake` function" — INSEE/Lumley `survey::calibrate`), not a claim that `leafblower` offers a `grake` method. Legitimate comparative-methods documentation; `docs/methods/logit.md` is explicitly in D-03 scope but this text is accurate as written. |
| `L-BFGS-B` / `LBFGSB` | `docs/methods/00-overview.md:3,81` | Correctly documents the removed status: "enum slot `2` is removed (was LBFGSB)" and a direct header-comment quote `// 2 = removed (was LBFGSB)`. Matches src/leafblower.h. No fix needed. |
| `Gurobi` / `CPLEX` | `docs/methods/chebyshev.md:115,119` | Comparative discussion of commercial LP solvers' algorithmic properties (Barrier method, simplex), not a claim leafblower includes or outperforms them. Legitimate. |
| `cp` (file-path reference) | `docs/performance.md:176` | References `docs/investigations/2026-05-02-ylsy-cp-ipm-spike-result.md`, a withdrawn research spike ("tried two additional solver families (Chambolle-Pock primal-dual and an interior-point Newton method) ... neither did [break the ceiling]"). Correctly reports `cp` was tried and rejected, consistent with REQUIREMENTS.md's own "Epic-K `cp` (slot 12): Landed `00a3f10`, reverted `3fac1d6`, same day. Never live." Not a capability claim. |
| `method="grake"` in a test | `tests/testthat/test-calibration-solvers.R:1-11` (test `T1a`) | Test asserts `expect_error(harvest(..., method="grake", ...), regexp = "must be exactly one of")` — this is a regression test confirming `grake` is REMOVED, not a doc claiming it's available. Test files are code, not the SC4 audit surface (`README`/`NEWS.md`/`man/`/`docs/`) per D-03; no action needed. |

**Clean files confirmed by direct grep this session (zero hits for grake/lbfgsb/cp-as-method):**
`README.md`, `README.Rmd`, `NEWS.md`. No vignettes directory exists in this repo. No `man/*.Rd`
file references `grake`/`lbfgsb` except the `grake_norm` false positive above.

**Conclusion:** SC4's audit surface has exactly ONE true-positive location (docs/raking.md,
already covered by SC1/D-01) and zero others. No README.md, NEWS.md, or man/ edits are needed
for SC4 beyond what SC1 already fixes.

## Common Pitfalls

### Pitfall 1: The `weights=` guard breaks a live test on implementation (not hypothetical)

**What goes wrong:** Implementing D-04's `stop()` guard, without also fixing the one caller that
relies on the current silent-fallback behavior, causes `R CMD check`/`devtools::test()` to fail
immediately.

**Exact location** [VERIFIED: tests/testthat/test-logit.R:182-202, read in full this session]:

```r
test_that("eb79.18: consistent collinear with HETEROGENEOUS design weights still reaches feasibility", {
  library(leafblower)
  set.seed(23)
  n <- 200
  base <- factor(c(rep("A", n/2), rep("B", n/2)))
  base2 <- factor(rep(c("P", "Q", "P", "Q"), length.out = n))
  data <- data.frame(v1 = base, v2 = base, v3 = base2, v4 = base2)
  target <- list(v1 = c(A = 0.5, B = 0.5), v2 = c(A = 0.5, B = 0.5),
                 v3 = c(P = 0.5, Q = 0.5), v4 = c(P = 0.5, Q = 0.5))
  base_w <- runif(n, 0.5, 2.0)
  w <- suppressWarnings(harvest(data, target, method = "logit", weights = base_w,
               min_weight = 0.05, max_weight = 20, max_iterations = 500L,
               attach_weights = FALSE))
  r <- attr(w, "result")
  expect_equal(r$status, 0L, label = "heterogeneous consistent collinear: converges (eb79.18)")
  expect_lt(r$max_error, 1e-4)
})
```

**Why it happens:** `weights = base_w` is not `design_weights=`. Today this falls into `...`,
triggers the RVAL.2 warning, and `suppressWarnings()` silently hides that warning — so the test
name's claim ("HETEROGENEOUS design weights") is currently **false**: the test actually runs
with uniform/default weighting, not `base_w`. A repo-wide grep this session
(`grep -rn -E '(^|[^_.a-zA-Z])weights[[:space:]]*=[[:space:]]*[a-zA-Z]' tests/testthat/*.R
R/*.R vignettes/*.R* benchmarks/*.R`, false-positive-filtered) found this is the **only**
occurrence anywhere in the repo of a bare `weights=` passed toward `harvest()`; two
`benchmarks/stepstone_fulldata_homotopy.R` hits are unrelated `list(..., weights = w)` return
values, not `harvest()` call arguments.

**How to avoid:** Fix `test-logit.R:195` to `design_weights = base_w` in the SAME task/commit
that adds the `stop()` guard — not a follow-up. After the rename, re-check whether
`suppressWarnings()` is still needed: the adjacent test in the same file
(`test-logit.R:164-179`, the INCONSISTENT-margins case) also wraps its `harvest()` call in
`suppressWarnings()` with no `weights=` argument at all, suggesting `logit`'s collinear/redundant-
margin path may independently emit a legitimate warning unrelated to the `weights=` typo — do
not blindly strip `suppressWarnings()`; confirm empirically (run the test after the rename and
observe whether any warning still fires) before deciding whether to keep or drop it.

**Warning signs:** If the planner scopes SC3 as "add one `if` block to harvest.R" without a
paired test-file task, `R CMD INSTALL --preclean .` + `Rscript -e "devtools::test()"` will show
a new FAIL in `test-logit.R` that did not exist before the change — this is the tell that the
paired fix was skipped.

### Pitfall 2: docs/raking.md is headingless — line-number edits will silently corrupt the file

**What goes wrong:** The file is 10 physical lines holding 43,826 bytes of running prose with
zero `#` headers (confirmed `wc -l -c` and `grep -c '^#'` this session). A diff/patch tool or an
executor assuming standard markdown section boundaries will either fail to match or, worse,
match a substring inside an unrelated sentence.

**Why it happens:** The doc was authored (per ticket `leafblower-05ha`'s own investigation) as a
single unstructured blob per numbered "section" with no line breaks between them.

**How to avoid:** Locate and replace by exact quoted string (as done in this research — the two
passages are extracted verbatim above), not by line/byte offset. Re-run the ticket's own
verification greps (`grep -o 'L-BFGS-B' docs/raking.md | wc -l`, `grep -c -i Gurobi
docs/raking.md`) before and after editing and paste both counts into the completion evidence —
this session found the actual "before" L-BFGS-B count is **12**, not the ticket body's cached
**4**; trust a fresh count, not the stale ticket number.

**Warning signs:** `git diff docs/raking.md` showing near-total file rewrite (expected, given
the format) vs. a targeted few-line diff (would indicate the section boundaries were
mis-located).

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|--------------|-----|
| Informative "unknown arg, did you mean X" error | A generic fuzzy-match/`did-you-mean` argument suggester | A specific, hard-coded `if ("weights" %in% names(dots)) stop(...)` check inside the existing RVAL.2 block | This is one known typo for one known correct name — a general suggestion engine is speculative scope the phase does not need; `weights` is the only argument name confirmed (by repo-wide grep) to actually trip this path today. |
| Roxygen/man regeneration | Hand-editing `man/harvest.Rd` directly | `roxygen2::roxygenize()` (project convention, per CLAUDE.md: "`man/` is roxygen2-generated: edit the roxygen block, regenerate") | Only applies if any SC3/SC4 fix touches an `@param`/`@details` roxygen block — current findings show none of the four SCs require a roxygen edit (the `design_weights=` doc block already exists correctly at `R/harvest.R:173`), but flag for the executor in case the `stop()` message needs a `@param` addendum. |

**Key insight:** All four fixes are edits to existing, correctly-patterned code/prose elsewhere
in the same files — no new abstraction, helper, or dependency is justified anywhere in this
phase.

## Code Examples

### `weights=` guard, matching the package's existing idiom

```r
# Insert inside the existing RVAL.2 block, R/harvest.R:303-307 — before or alongside the
# generic warning, using the same "leafblower: " prefix + call. = FALSE convention as every
# other stop() in this file (e.g. line 508).
dots <- list(...)
if ("weights" %in% names(dots))
  stop("leafblower: unrecognized argument 'weights' — did you mean 'design_weights'? ",
       "harvest() takes per-observation design weights via design_weights=, not weights=.",
       call. = FALSE)
if (length(dots) > 0L)
  warning("harvest: unknown argument(s) ignored: ",
          paste(names(dots), collapse = ", "), call. = FALSE)
```

Exact wording beyond naming `design_weights=` is Claude's discretion per CONTEXT.md.

### `rk_algorithm_t` slot-7 annotation

```c
// src/leafblower.h — insert between line 48 (RK_ALG_GREG = 6,) and line 49 (RK_ALG_ORIS_SOFT = 8,)
    RK_ALG_GREG      = 6,
    /* 7 = removed (was RK_ALG_GRAKE) */
    RK_ALG_ORIS_SOFT = 8,   /* oris + ADMM soft capacity enforcement */
```

Mirrors the existing slot-2 form exactly (`/* 2 = removed (was RK_ALG_LBFGSB) */`, line 44).

### NEWS.md breaking-change entry style (for the `weights=` stop(), if the planner elects to log it)

The project's own `## Breaking changes` heading pattern already exists in NEWS.md (used 4 times
this file), e.g.:

```
* BREAKING: solver method "ieppa"/"ieppa_soft" renamed to "oris"/"oris_soft" (ORIS — Over-Relaxed
  Iterative Scaling). No alias; update calls. Algorithm, numeric output, and enum values unchanged.
```

D-04 flags this as a "costly"/breaking reversibility item requiring a `checkpoint:decision` —
if the user confirms proceeding, a matching NEWS.md bullet under `## Breaking changes` (or a
new top dev-version entry, matching the file's existing top-of-file dev-version convention) is
the established idiom; this was NOT independently requested by the phase's success criteria, so
treat it as optional/discretionary unless the checkpoint decision calls for it.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| — | none | — | All claims in this research were verified this session via direct file reads, `git log`/`git show`, or repo-wide `grep` — no `[ASSUMED]`-tagged claims. This phase requires no package/ecosystem research (no new dependencies), so the package-legitimacy and version-verification protocols do not apply. |

**This table is empty:** All claims in this research were verified or cited — no user
confirmation needed beyond the pre-existing D-04 `checkpoint:decision` already specified in
CONTEXT.md (which this research resolves with a concrete finding: yes, one test relies on the
old behavior, and the fix is specified above).

## Open Questions

1. **Should `suppressWarnings()` be removed from `test-logit.R:195` once `weights=` is renamed
   to `design_weights=`?**
   - What we know: the sibling test in the same file (INCONSISTENT-margins case, lines 164-179)
     also uses `suppressWarnings()` with no `weights=` argument at all, suggesting `logit` may
     emit an independent, legitimate warning on certain collinear/heterogeneous inputs.
   - What's unclear: whether the HETEROGENEOUS-case test specifically triggers that same
     independent warning, or whether `suppressWarnings()` there exists ONLY to hide the
     `weights=` typo's RVAL.2 warning.
   - Recommendation: keep `suppressWarnings()` unless the executor empirically confirms (by
     running the test with `design_weights=` and no `suppressWarnings()`) that zero warnings
     fire — do not decide this by inspection alone.

2. **Does the D-04 `checkpoint:decision` (user must confirm before implementing the breaking
   change) block Task ordering within the phase plan?**
   - What we know: CONTEXT.md explicitly requires a `checkpoint:decision` before the
     implementing task, and this research now supplies the missing evidence that checkpoint was
     waiting on (one affected test, with a concrete fix).
   - What's unclear: whether the plan should sequence the checkpoint before or after presenting
     this research's finding to the user.
   - Recommendation: the planner should surface Pitfall 1's finding directly in the
     `checkpoint:decision` prompt so the user is deciding with full information (not deferring
     research to checkpoint time).

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | testthat 3 — `[VERIFIED: DESCRIPTION]` `Config/testthat/edition: 3`, `Suggests: testthat (>= 3.0.0)` |
| Config file | `tests/testthat.R` (standard testthat 3 harness; not modified by this phase) |
| Quick run command | `Rscript -e "testthat::test_dir('tests/testthat', filter='harvest-rval')"` (matches basenames only per project convention) |
| Full suite command | `Rscript -e "devtools::test()"` |

### Phase Requirements → Test Map
| SC | Behavior | Test Type | Automated Command | File Exists? |
|----|----------|-----------|-------------------|-------------|
| SC1 | docs/raking.md contains no L-BFGS-B/Gurobi/outer-loop claim | manual-only (prose, no code path) | `grep -o 'L-BFGS-B' docs/raking.md \| wc -l` must be 0; `grep -ci Gurobi docs/raking.md` must be 0 — the ticket's own verification greps, re-run post-edit | N/A — grep-based check, not a testthat file |
| SC2 | `rk_algorithm_t` slot 7 annotated; enum values unchanged | unit (build-level) | `git diff src/leafblower.h` shows only comment lines added (manual diff review); full `Rscript -e "devtools::test()"` + Python pytest must stay 0 FAIL (regression net — no dedicated test asserts a comment's presence) | ✅ (regression coverage exists; no NEW file needed) |
| SC3 | `harvest(weights=w)` raises an error naming `design_weights=` | unit | `Rscript -e "testthat::test_dir('tests/testthat', filter='harvest-rval')"` after adding a new `test_that("RVAL.4: bare weights= errors naming design_weights", ...)` block | ❌ Wave 0 — see gap below |
| SC3 (regression) | `test-logit.R:195`'s `eb79.18` HETEROGENEOUS case still passes after `weights=` → `design_weights=` rename | unit | `Rscript -e "testthat::test_dir('tests/testthat', filter='logit')"` | ✅ file exists, needs the one-line arg rename (Pitfall 1) |
| SC4 | No surviving grake/lbfgsb/cp-as-method reference in README/NEWS/man/docs | manual-only (prose audit) | Re-run this session's greps (see SC4 findings table) post-edit; expect the same zero-true-positive result outside docs/raking.md | N/A — grep-based check |

### Sampling Rate
- **Per task commit:** the relevant filtered `testthat::test_dir(..., filter=...)` run for the
  file(s) touched (`harvest-rval`, `logit`), plus the SC1/SC4 verification greps for doc-only
  tasks.
- **Per wave merge:** full `Rscript -e "devtools::test()"` (R CMD INSTALL --preclean . first if
  `src/leafblower.h` was touched) — the header edit (SC2) forces a full rebuild of both R and
  Python extensions per CLAUDE.md's "Two build sites" note.
- **Phase gate:** Full DoD — R tests 0 FAIL + Python parity tests 0 FAIL (single-thread BLAS,
  `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`) — before `/gsd:verify-work`.
  Stepstone benchmark (`LBW_BENCH_GATE=1`) is opt-in/heavy and not required by this phase's own
  Definition of Done unless the executor judges the header touch (SC2) warrants re-confirming no
  perf regression; none of SC1/SC3/SC4 touch any hot path.

### Wave 0 Gaps
- [ ] A new `test_that("RVAL.4: ...")` block in `tests/testthat/test-harvest-rval.R`, asserting
      `expect_error(harvest(df, tgt, weights = ...), regexp = "design_weights")` — covers SC3.
      No dedicated test for the bare-`weights=` case exists today; the file's own RVAL.1/RVAL.2/
      RVAL.3/META.2 numbering convention is the natural home for RVAL.4.
- [ ] `tests/testthat/test-logit.R:195` — the existing `eb79.18` HETEROGENEOUS test must be
      edited (`weights=` → `design_weights=`) in the SAME task as the RVAL.4 addition, or the
      full suite goes red (Pitfall 1). Not a new file — an edit to an existing one.
- No other gaps: SC1/SC2/SC4 are prose/comment-only changes covered by manual grep verification
  plus the existing full-suite regression net (already green today), not new test files.

## Security Domain

### Applicable ASVS Categories

| ASVS Category | Applies | Standard Control |
|---------------|---------|-------------------|
| V2 Authentication | no | N/A — no auth surface in this package |
| V3 Session Management | no | N/A |
| V4 Access Control | no | N/A |
| V5 Input Validation | yes | SC3 IS a V5 finding: an unrecognized/misspelled argument (`weights=`) is currently accepted silently and produces a plausible-but-wrong result instead of being rejected. The fix is exactly V5.1-class validation — reject unrecognized input at the trust boundary (the `harvest()` public API) rather than silently coercing/ignoring it. Standard control: explicit allow-list check against `names(dots)` before the generic warning, `stop()` on the known-typo case — implemented with base R `stop()`, consistent with every other input-validation `stop()` already in `R/harvest.R` (lines 447-508). No new validation library needed (`Don't Hand-Roll` table). |
| V6 Cryptography | no | N/A — no cryptographic material in this package |

### Known Threat Patterns for this stack

| Pattern | STRIDE | Standard Mitigation |
|---------|--------|-----------------------|
| Silent argument-name typo produces a plausible-but-wrong statistical result (the SC3 defect itself) | Tampering (of the effective computation, via an unnoticed input-validation gap — not a malicious-actor threat, but the same STRIDE category as an unvalidated input silently altering program behavior) | Explicit reject-on-unrecognized-name at the API boundary (SC3's `stop()` fix); this is a correctness/trust defect, not a network-facing security vulnerability — this package has no network, auth, or untrusted-multi-tenant surface, so no other STRIDE category applies. |
| Documentation over-claiming an unimplemented capability (SC1) | Spoofing (of capability — a reader trusts a false claim about what the software does) | Delete/correct the false claim (D-01); not a code-level mitigation, a documentation-truthfulness fix, but classified here because a downstream automated doc-consumer (e.g. an LLM agent reading `docs/raking.md` to decide how to call this package) is exactly the "attacker" this ASVS-adjacent framing protects against — a caller cannot be tricked into invoking a nonexistent capability. |

## Sources

### Primary (HIGH confidence — direct file reads / git history this session)
- `docs/raking.md` — full-text grep extraction of §8.2 and §12 passages, byte/line count.
- `docs/methods/oris.md` — read in full (lines 1-60), ground-truth ORIS description.
- `src/leafblower.h` — lines 30-59 read in full, `rk_algorithm_t` enum verbatim.
- `R/harvest.R` — lines 1-40, 260-334, 505-519 read in full.
- `tests/testthat/test-logit.R` — lines 160-220 read in full.
- `tests/testthat/test-calibration-solvers.R` — lines 1-20 read in full.
- `src/greg.cpp:148-156`, `src/logit_calib.cpp:553-562` — read in full, `grake_norm` field
  assignments.
- `man/harvest.Rd` — lines 265-276 read in full.
- `docs/methods/logit.md` — lines 105-119 read in full.
- `docs/methods/00-overview.md`, `docs/performance.md:168-181`, `docs/methods/chebyshev.md`
  — grep-located, contextually read.
- `git log --oneline --all`, `git show 9a67891 --stat` — slot-7/grake removal history.
- `bd show leafblower-05ha`, `bd show leafblower-x2iq` — full ticket bodies.
- `DESCRIPTION` — Imports/Suggests read, confirms no `cli` dependency.
- `NEWS.md` — read in full (138 lines), confirms zero grake/lbfgsb hits and documents the
  `## Breaking changes` style convention.

### Secondary / Tertiary
- None used — this phase required no external package research, no WebSearch, no Context7
  lookups. All findings are in-repo.

## Metadata

**Confidence breakdown:**
- SC1 (docs/raking.md): HIGH — exact passages extracted verbatim, ground-truth comparison text
  read in full.
- SC2 (enum slot 7): HIGH — header read in full, git history for the removal commit confirmed.
- SC3 (weights= guard): HIGH — signature, existing warning idiom, and the one affected test all
  read in full; repo-wide grep confirms no other caller is affected.
- SC4 (grake/lbfgsb/cp audit): HIGH — exhaustive grep across all in-scope files, every hit
  individually read and classified true/false positive.

**Research date:** 2026-08-15
**Valid until:** No expiry driver — this is a point-in-time audit of the current working tree,
not a claim about an external ecosystem. Re-run the greps if `docs/raking.md`,
`src/leafblower.h`, `R/harvest.R`, or the audited doc surfaces change before this phase is
planned/executed.
