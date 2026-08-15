# Decisions (ADR-derived)

**No ADR-class documents in the ingest set.** All 44 classified docs are `SPEC` (32),
`DOC` (11) or `PRD` (1); none carries `type: ADR` and none carries `locked: true`. This
file is therefore empty of decision entries by construction — per the extraction
contract, absent input is marked absent, not synthesised from lower-precedence material.

Design decisions in this corpus live inside SPECs and are extracted to `constraints.md`
with their source attribution. Where a later SPEC reverses an earlier one on the same
subject, the supersession is recorded in the constraint entry's `content` block and
summarised in `SYNTHESIS.md` § Supersession chains.

Two documents function as decision records without being ADRs:

- `docs/superpowers/specs/2026-08-14-removed-solver-slots-supersession.md` (Status:
  Accepted) fixes the dispositions of `grake` (slot 7), `lbfgsb` (slot 2) and `cp`
  (slot 12) against named commits. It is not `locked`, but it is the newest document in
  the set and explicitly declares itself authoritative over every earlier spec on those
  three subjects; it is applied as such. Its content is carried in `constraints.md` under
  "Live algorithm enum and permanently reserved slots (authoritative)" and "Removed solver
  dispositions (grake, lbfgsb, cp)".
- `tasks/prd-leafblower-core.md` § 12 Open Questions records five resolved decisions
  (OQ-1 … OQ-5) in ADR-like form. Three still stand — fixed `epsilon = 0.05` per paper
  with adaptive outer tolerances (OQ-3), `log_fn` callback for verbose output rather than
  `message()`/`print()` (OQ-4), and Windows CRAN binary resolved as a non-goal with
  `Makevars.win` deferred (OQ-5). Two (OQ-1 rake/nr mapping, OQ-2 logit degenerate case)
  resolved into L-BFGS-B content that has since been withdrawn. These are carried in
  `requirements.md`, not here, because their source is a PRD.

**Consequence for downstream consumers:** there is no locked-decision layer to enforce, so
LOCKED-vs-LOCKED contradiction is impossible in this set. Precedence is resolved by
(a) SPEC > PRD > DOC, then (b) document date (later supersedes earlier on the same
subject), as directed by the ingest brief. The PRD is the OLDEST document in the corpus,
so precedence and date agree everywhere they both apply.
