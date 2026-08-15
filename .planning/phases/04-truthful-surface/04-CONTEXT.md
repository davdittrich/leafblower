# Phase 4: Truthful Surface - Context

**Gathered:** 2026-08-15
**Status:** Ready for planning

<domain>
## Phase Boundary

Defect-driven phase (no user stories) fixing four specific truth/silence violations in the
shipped surface: a misattributed algorithm description in `docs/raking.md`, an undocumented
hole in the `rk_algorithm_t` enum, a silently-swallowed `weights=` typo in `harvest()`, and
stale references to removed methods (`grake`, `lbfgsb`, `cp`) in shipped user-facing
documentation. Nothing new is built; existing false or silent surfaces are corrected.

</domain>

<decisions>
## Implementation Decisions

### docs/raking.md §8.2/§12 misattribution
- **D-01:** Delete the §8.2/§12 passage entirely rather than rewrite it under the paper's own
  name (Chu-Liang-Toh-Yang 2022 iEPPA, arXiv:2011.14312). Do not leave an ORIS-labelled
  description of the unimplemented outer entropic-proximal-point loop anywhere in the file.
  Ticket: `leafblower-05ha`.
- **D-02:** File a follow-up ticket (separate from this phase) noting that the rest of
  `docs/raking.md` (§1-7, §9-11) has zero bracketed citations — unlike `docs/methods/oris.md`'s
  59 real citations — and includes content unrelated to this codebase (ROAM medical-imaging
  architecture, a generic Rust-vs-Python microbenchmark, Statistics Canada GES tangents), and
  §12's synthesis calls L-BFGS-B "the definitive state-of-the-art" despite it being removed
  from the package. Do not act on this in Phase 4 — out of scope, flagged for a maintainer
  look only. — **Reversibility:** reversible — filing a ticket has no downstream lock-in.

### grake/lbfgsb/cp audit scope
- **D-03:** The README/NEWS.md/man/docs audit (SC4) covers shipped, user-facing surfaces only:
  `README.md`, `NEWS.md`, `man/`, `docs/methods/`, `docs/raking.md`, `docs/performance.md`,
  and any vignettes. `docs/superpowers/specs/` (internal design/planning documents) is
  explicitly OUT of scope — those documents are historical planning records expected to
  reference retired or future method names (`cp`, `lbfgsb`) and are not something a package
  user or CRAN reviewer reads.

### weights= guard
- **D-04:** `harvest(..., weights = w)` must raise a hard `stop()` naming `design_weights=`,
  not a warning. Currently the argument falls into `...` and triggers the generic RVAL.2
  warning ("unknown argument(s) ignored: weights") — non-fatal, easy to miss, and the call
  proceeds with a plausible-but-wrong unweighted result. — **Reversibility:** costly —
  changing an established R package's argument-handling from warn to error is a breaking
  change for any caller currently relying on (or silently tripping) the old fallback; a
  CRAN-facing package should flag this in NEWS.md as a behavior change. Needs a
  `checkpoint:decision` before the task implementing it, per phase 4's own truthfulness goal —
  the researcher/planner should confirm no in-repo test or benchmark currently passes
  `weights=` expecting the old silent-fallback behavior before this ships.

### Claude's Discretion
- Exact wording of the `stop()` message for `weights=` (beyond naming `design_weights=`
  explicitly) — planner/executor's call.
- How to document the slot-7 hole in `rk_algorithm_t` — match slot 2's existing comment style
  (`/* N = removed (was RK_ALG_X) */`) once research determines why slot 7 was never assigned
  (git history — no prior CLAUDE.md or ROADMAP note records the reason, unlike slot 2 which is
  explicitly called out in CLAUDE.md's Architecture section).

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### docs/raking.md misattribution
- `docs/raking.md` §8.2, §12 — the passage to delete (D-01).
- `docs/methods/oris.md` — states the opposite, correct claim: the paper's "headline
  contribution — an outer inexact entropic-proximal-point loop … is NOT implemented here",
  "mathematically inert" at `C = 0`. This is the authoritative description; do not duplicate
  it into `docs/raking.md`.
- `.planning/INGEST-CONFLICTS.md` lines 24-51 — the full prior conflict analysis (WARNING
  entry) that already traced this misattribution, confirmed the paper is
  Chu-Liang-Toh-Yang 2022 (arXiv:2011.14312), and recorded the resolution precedence
  (SPEC > DOC; DOC-vs-DOC broken by the rename spec).
- `docs/superpowers/specs/2026-05-30-oris-rename-design.md` §1, §8.10 — the original rename
  rationale and scope list that named `docs/raking.md` as a LIVE doc requiring the update.
- `tasks/prd-leafblower-core.md` § US-005 — where the "paper-faithful iEPPA" framing
  originates; superseded on this point, not to be treated as current truth.

### rk_algorithm_t slot 7
- `src/leafblower.h` lines 40-53 — the `rk_algorithm_t` enum. Slot 2 already has an inline
  comment (`/* 2 = removed (was RK_ALG_LBFGSB) */`); slot 7 has NO comment and is simply
  absent (enum jumps GREG=6 → ORIS_SOFT=8). Match slot 2's comment style once the reason for
  slot 7 is researched.
- `CLAUDE.md` — "Algorithm slot 2 is reserved (LBFGSB removed). Do not reuse in
  `rk_algorithm_t` enum." Slot 7 is NOT currently documented here either — the planner should
  consider whether this project convention file also needs the slot-7 note added, matching
  the existing slot-2 line.

### weights= guard
- `R/harvest.R` lines 267-310 — the `harvest()` signature (no `weights=` formal; falls into
  `...`) and the RVAL.2 generic-warning block right after the signature (`dots <- list(...)`
  / `warning("harvest: unknown argument(s) ignored: ...")`).
- `R/harvest.R` line 173, 298, 510-515 — the actual `design_weights` parameter and its
  existing warning for the *different* case of `design_weights` + `start_weights` both
  supplied — a working example of an informative, specific warning to model the new
  `stop()` message's tone on.

### grake/lbfgsb/cp audit
- `leafblower-x2iq` (P2, slot-7 ticket) and `leafblower-05ha` (P1, docs/raking.md ticket) —
  existing beads tickets this phase closes.
- Ticket needed: `weights=` guard (per ROADMAP notes — "file a ticket for the weights= guard").
- `.planning/REQUIREMENTS.md` lines 195, 271 — the two existing tracking rows for this phase's
  scope.

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `R/harvest.R`'s existing RVAL.2 unknown-arg warning mechanism (`dots <- list(...)`) is the
  natural insertion point for a specific `weights=` check — intercept before the generic
  warning fires, not a parallel mechanism.
- `src/leafblower.h`'s slot-2 inline-comment convention (`/* N = removed (was RK_ALG_X) */`)
  is the pattern to replicate for slot 7 — no new documentation convention needed.

### Established Patterns
- `docs/methods/oris.md`'s citation-rigorous style (bracketed `[key1234]` citations resolved
  against a real bibliography) is the standard the rest of `docs/methods/` follows;
  `docs/raking.md` does not currently meet it outside this phase's narrow fix, per D-02.

### Integration Points
- None — this phase touches documentation and one R-level argument guard; no new code paths
  or cross-module wiring.

</code_context>

<specifics>
## Specific Ideas

No specific requirements beyond the four ROADMAP success criteria and the decisions above —
open to standard approaches for exact wording/formatting.

</specifics>

<deferred>
## Deferred Ideas

- **docs/raking.md broader quality review** (D-02) — the non-§8.2/§12 sections lack citations
  and include off-topic content. File as a new ticket; a future phase or ad-hoc doc pass, not
  Phase 4.

### Reviewed Todos (not folded)
None — `todo.match-phase 4` returned zero matches.

</deferred>

---

*Phase: 4-Truthful Surface*
*Context gathered: 2026-08-15*
