# Synthesis

Entry point for downstream consumers. Produced by `gsd-doc-synthesizer` from 44
classified documents. Mode: `new` (no pre-existing PROJECT.md / REQUIREMENTS.md /
ROADMAP.md / STATE.md). Third run: it supersedes the previous synthesis, which was
AWAITING USER on four warnings.

**STATUS: AWAITING USER — 0 blockers, 1 warning in `../INGEST-CONFLICTS.md`. The
remaining warning is a documentation defect that is already ticketed
(`leafblower-05ha`), not a contradiction that blocks planning. Safe to route once
acknowledged.**

What changed since run 2: the missing PRD was located at `tasks/prd-leafblower-core.md`
(outside the `docs/` convention that discovery walked), so `requirements.md` now carries
ten entries instead of being empty by construction. Three of run 2's four warnings are
closed — the PRD is present, and the `omega_mode_id` wiring and `kNCatsTotalMax` value
were verified against live source by the orchestrator between runs.

---

## Doc counts by type

- SPEC — 32 (`docs/superpowers/specs/`, dated 2026-04-18 … 2026-08-14)
- DOC — 11 (`docs/methods/` x 9, `docs/raking.md`, `docs/raking_bounds.md`)
- PRD — 1 (`tasks/prd-leafblower-core.md`, Draft v3, 2026-04-18)
- ADR — 0
- UNKNOWN — 0

All 44 classifications are `confidence: high` with `manifest_override: true`; none is
`locked`; none carries a per-doc `precedence` override.

## Decisions

- Locked decisions: **0** — there are no ADR-class documents in the set, so there is no
  locked-decision layer to enforce and LOCKED-vs-LOCKED contradiction is structurally
  impossible.
- `decisions.md` is empty of decision entries by construction and states why.
- Two documents function as decision records without being ADRs:
  `2026-08-14-removed-solver-slots-supersession.md` (Status: Accepted, commit-traced,
  self-declared authoritative over all earlier specs on `grake`/`lbfgsb`/`cp`) and the
  PRD's § 12 Open Questions table (OQ-1 … OQ-5). Applied on date + explicit
  self-declaration, not on a locked flag.
- All other design decisions live inside SPECs and are carried in `constraints.md` with
  source attribution and supersession notes.

## Requirements

**10 entries** in `requirements.md`, all from the single PRD. Entries keep the PRD's own
`US-xxx` identifier inside the `REQ-` slug, because the SPEC corpus cites those
identifiers verbatim:

- `REQ-us001-autumn-drop-in` — the seven exported functions and the honour / silently-
  ignore / deprecate parameter matrix (FR-29 … FR-35)
- `REQ-us002-min-weight-lower-bound` — the package's headline addition over `autumn`
- `REQ-us003-large-scale-performance` — 1M rows, 20 margins, < 30 s single-threaded
- `REQ-us004-c-api-contract` — FR-1 … FR-10: validation rules, error codes, memory and
  thread-safety contract, the Python GIL/trampoline contract
- `REQ-us005-oris-capacity-constrained-solver` — cell compression, log-space factors, the
  every-inner-step bounds invariant (partially superseded; see below)
- `REQ-us005b-classical-raking` — the "§ US-005b" the SPEC corpus cites
- `REQ-us008-python-pandas-interface` — FR-36 … FR-40, including the copy-never-view
  contract
- `REQ-us009-diagnostic-functions` — Kish (1992), Henry & Valliant (2015)
- `REQ-us010-cran-pypi-distribution`
- `REQ-kpi-success-metrics` — the § 11 KPI table

**The PRD is heavily superseded and must be read with its supersession notes.** It is
Draft v3 dated 2026-04-18 — the OLDEST document in the corpus — and describes a
TWO-algorithm package; eight solvers ship. Precedence (SPEC > PRD) and document date
agree, so nothing here is a conflict needing resolution. Eleven superseded/withdrawn
items are listed at the foot of `requirements.md`; the load-bearing ones:

- **US-006 and FR-20 … FR-28 (L-BFGS-B): WITHDRAWN entirely.** The PRD's § 6 enum
  (`RK_ALG_AUTO=0, RK_ALG_IEPPA=1, RK_ALG_LBFGSB=2`) is the ORIGIN of the slot-2 hole.
- **US-007 auto-routing: SUPERSEDED** — routing now spans eight solvers.
- **US-001 `method="rake"`/`"nr"` -> `"lbfgsb"`: SUPERSEDED** — the synonyms still need a
  destination; the PRD no longer supplies one.
- **§ 5 Non-Goals: four of nine have since shipped or been worked** (SQUAREM, bounded-IPF
  water-filling, `auto_collapse`, `add_na_proportion`) — superseded scope, NOT violated
  requirements. CPU-only and single-threaded remain in force.
- **§ 7 `PKG_CXXFLAGS ... -O3`: SUPERSEDED** — the R build sets no `-O` level by design.
- **US-005 paper-faithfulness: repudiated.** The outer entropic-proximal loop (FR-11,
  FR-12, FR-14) is not implemented and is inert at `C = 0`; that is the stated reason for
  the `ieppa` -> `oris` rename. The rest of US-005 is live and corroborated by
  `docs/methods/oris.md`.
- `iEPPA` throughout the PRD = the solver now named `oris`. Terminology, not a conflict.

Per-spec acceptance criteria (A1…A9, AC1…AC11, T1…T8 blocks) with their numeric gates
remain preserved separately as `nfr`-typed constraint entries. Every acceptance criterion
attached to a removed solver — A4 of `2026-04-25-calibration-solvers-design.md`, the
T1…T8 gate set of the Epic-K CP spec, the PRD's whole US-006 — does NOT carry forward.

## Constraints

**87 entries** in `constraints.md` (unchanged this run — SPEC-derived), by type:

- `protocol` — 50 (algorithm update rules, guard envelopes, projection/fallback ordering,
  state-reset checklists, test mechanics)
- `api-contract` — 23 (enum values, R/Python/C parameter surfaces, result-field contracts,
  status-code semantics, AUTO routing, removed-solver dispositions)
- `nfr` — 8 (merge gates, per-iteration cost targets, determinism/measurement protocol,
  invariant tolerances, build-list synchronisation)
- `schema` — 6 (`CellTable`, `SRAAState`, ABI tripwires, enum freeze, live enum + reserved
  slots)

Entries are grouped by subsystem inside the file for navigation only; every entry is
self-contained and carries its own `source:`. Three entries are re-titled "(WITHDRAWN
PROPOSAL — never landed)" — the Epic-K Chambolle-Pock trio. Four carry an in-place
supersession note from the 2026-08-14 record.

## Context topics

**19 topics** in `context.md` (unchanged this run — DOC-derived), covering the
eight-solver inventory and enum map, the shared input/exit contract, the
bounds-enforcement taxonomy, published selection guidance, the project's own "no
universal convergence proof" caveat, per-solver identity notes for all eight solvers, the
three-schemes/two-axes IPF taxonomy, why the ORIS theta2 estimator does not port to
sinkhorn or greenkhorn, literature failure modes, and how leafblower deliberately
deviates from standard implementations.

## Conflicts

- **0 blockers.**
- **1 warning** — `docs/raking.md` §8.2/§12 attribute the unimplemented outer
  entropic-proximal-point loop to ORIS, i.e. the exact claim the rename exists to
  repudiate. Ticketed as `leafblower-05ha`.
- **25 info** — seven prior blockers/warnings now resolved (grake, lbfgsb, cp,
  `RK_ALG_RAKING = 3`, the located PRD, `omega_mode_id`, `kNCatsTotalMax`); six PRD
  supersession records; the no-locked-layer note; the clean cycle-detection result; five
  supersession chains resolved by date; two explicit non-conflicts; ABI tripwire drift;
  and the fully superseded 2026-04-18 spec.

Full detail with per-claim source references: `../INGEST-CONFLICTS.md`.

## Resolved dispositions (authoritative — do not re-derive)

Source for the first four: `2026-08-14-removed-solver-slots-supersession.md`. The last
two were verified against live source by the orchestrator.

- **`grake` / slot 7** — removed pre-release, commit `9a67891` (2026-04-30). Slot
  permanently reserved. `src/grake.{cpp,hpp}` do not exist. All four specs treating it as
  live are historical. The `grake_norm` METRIC is live and unrelated
  (`src/greg.cpp:153`, `src/logit_calib.cpp:558`).
- **`lbfgsb` / slot 2** — removed, annotated at `src/leafblower.h:44`; doc drift purged in
  `7fa211c`. Four specs AND the PRD's US-006/FR-20 … FR-28 target it; all historical.
  Open follow-up: `alm_lambda`/`alm_mu` were reserved for its ALM and may now be dead.
- **`cp` / slot 12** — landed `00a3f10`, reverted `3fac1d6`, same day (2026-05-03). Never
  live. Epic-K is a withdrawn proposal, NOT pending work, unless deliberately revived.
- **Live enum** — AUTO 0, ORIS 1, RAKING 3, SINKHORN 4, CHEBYSHEV 5, GREG 6, ORIS_SOFT 8,
  GREENKHORN 9, LOGIT 10, NEWTON_KL 11. Eight solvers. Values frozen (`77d0614`). Slots 2
  and 7 stay holes; any positional `alg_names` table must carry those holes.
- **`omega_mode_id`** — FIXED. All three default sites agree at 2 (`src/types.hpp:74`,
  `src/c_api.cpp:213`, `R/harvest.R:1077`), and Python forwards it
  (`_harvest.py:593`, `_bindings.cpp:140-142`). No R<->Python divergence.
- **`kNCatsTotalMax`** — **2048** (`src/calib_validate.hpp:10`). The spec § 12 row saying
  8192 is the losing statement.

## Supersession chains (resolved, later wins)

1. **AUTO routing** — PRD § US-007 (complexity > 500K / max_weight < 3 / min_weight > 0,
   iEPPA vs L-BFGS-B) -> 04-20 benchmarked L-BFGS-B contour -> 04-23 unconditional iEPPA
   -> **05-01 three-way rule** (K >= 5, `M_cell/n >= 0.9`, `target_skew` vs 5.0).
2. **Convergence configuration** — PRD's four base codes -> 04-24 `pct` default +
   `criterion` enum -> **04-25 orthogonal `metric` + `rule`**, default
   `max_err + improvement + 0.001`; then per-solver metric defaults (`kl` for ORIS and
   sinkhorn) and **04-28 status codes 4/5 with a weight-KL stall**.
3. **Acceleration** — APVA/joint Anderson -> Tang 2024 (shelved) -> Halpern -> SQUAREM ->
   SQUAREM obs-level-alpha fix -> **SRAA-m** (+ global safeguard + adaptive-sort/
   outer-revert) -> **SRAA-m in log-factor space for ORIS**. (SQUAREM is the mechanism the
   PRD § 5 disclaimed as a non-goal.)
4. **ORIS optimal omega** — global spectral theta2 (NO-GO) -> free-subspace marginal
   residual (NO-SHIP) -> **free-coordinate iterate-change theta2 with block-root cadence
   recovery** (confirmed wired, mode 2).
5. **Raking projection geometry** — PRD § US-005b additive Boyle-Dykstra -> Euclidean/
   multiplicative hybrid critique -> **Bregman (multiplicative) Dykstra at cell level**.
6. **Solver roster** — PRD's two algorithms -> grake / lbfgsb / cp specified -> **all
   three removed or never landed; eight-solver roster with reserved slots 2, 7
   (2026-08-14)**.

Superseded entries are retained and marked rather than deleted: several carry measured
negative results (uncapped damping stress, cell-level CBB alpha, fixed-sort
`f_eval_sraa`, marginal-residual theta2, ORIS warm-start of Newton-KL) that a future plan
must not re-discover the hard way.

## Notes for the roadmapper

- **Read `requirements.md` supersession notes before scheduling anything from the PRD.**
  It is the requirement layer the SPECs cite, but it is four months older than the code
  and specifies a package that no longer exists in that shape. Roughly a third of it is
  withdrawn.
- The roster is eight. Do NOT schedule work for `grake`, `lbfgsb` or `cp`, and do NOT
  reuse slots 2, 7 or 12.
- **Two genuine open items surface from the PRD, both small:** the autumn `"rake"`/`"nr"`
  synonyms need a live destination now that `"lbfgsb"` is gone (US-001, FR-34
  `choosemethod`); and the medium-scale performance target is stated as both < 1 s (§ 1)
  and < 2 s (§ 11) against a removed solver, so it must be re-benchmarked before being
  used as a gate.
- One cheap source-of-truth check remains: whether `alm_lambda`/`alm_mu` still have a
  consumer now that lbfgsb is gone.
- One documentation defect is actionable as written and already ticketed
  (`leafblower-05ha`): `docs/raking.md` §8.2/§12.
- The corpus is dense with *pre-registered ship gates* (numeric, fixture-pinned,
  TIE = NO-SHIP). Where work is proposed, the gate is usually already written down and
  should be carried into the plan rather than re-invented. The PRD adds project-level
  gates on top: `R CMD check --as-cran` clean, the 50-random-dataset bound property test,
  and the Python 3.9–3.13 wheel matrix.
- Several constraints are explicitly guarded against "helpful" correction — the chebyshev
  Mehrotra cross-term `-Delta_s_aff * Delta_y_aff`, the ORIS ALM Newton step
  `X~(1-lambda+mu*z)/(1+rho)`, the lambda-capture ordering in `raking.cpp`, the
  normalize->bounds finalization order, and the SRAA best-iterate metric (re-introduced
  twice). These are recorded in `constraints.md` with their rationale.

## Files

- `requirements.md` — 10 PRD-derived requirements + 11 superseded/withdrawn records
- `decisions.md` — empty of entries by construction (0 ADRs), with rationale
- `constraints.md` — 87 SPEC-derived constraints, grouped, each with `source:`
- `context.md` — 19 DOC-derived topics, each with `source:`
- `../INGEST-CONFLICTS.md` — 0 blockers, 1 warning, 25 info
- `classifications/` — the 44 per-doc classification inputs (unmodified)
