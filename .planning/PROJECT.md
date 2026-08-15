# leafblower

## What This Is

A survey-weight calibration package: ONE C++17 numerical core (`src/`) exposed to BOTH R
(`src/r_bridge.cpp` + `R/`) and Python (pybind11 / scikit-build-core, `python/`). It is a
drop-in replacement for the R package `autumn`, adding a `min_weight` lower bound that
neither `autumn` nor `anesrake` exposes. Eight solvers ship: `oris`, `oris_soft`, `raking`,
`sinkhorn`, `chebyshev`, `greg`, `greenkhorn`, `logit`, `newton_kl`. Users are survey
methodologists, census microsimulation researchers and R/Python survey analysts.

**This is a brownfield project.** The package is largely built (v0.1.0, ~11k LOC of C++,
94+ R test files, 20 Python test files, 1478 closed beads tickets). This planning layer
covers what remains, not what shipped.

## Core Value

Calibrated weights that are numerically correct, bound-respecting, and **identical from R
and from Python** — because both call the same compiled core.

## Requirements

### Validated

Shipped and confirmed by the codebase map (`.planning/codebase/`) and the test suite:

- ✓ **US-001** autumn drop-in — 7 exported functions, honour/ignore/deprecate parameter
  matrix, `rake`/`nrake`/`nr` synonyms routed to live solvers (`R/harvest.R:975-987`)
- ✓ **US-002** `min_weight` lower bound — the headline addition over `autumn`
- ✓ **US-004** stable C ABI — `leafblower.h`, validation, status codes, memory and
  reentrancy contract (one clause outstanding, see Active)
- ✓ **US-005** ORIS capacity-constrained solver — cell compression, log-space factors,
  every-inner-step bounds invariant
- ✓ **US-005b** classical raking — now on multiplicative KL-Bregman Dykstra geometry
- ✓ **US-008** Python/pandas interface — copy-never-view contract, pybind11 module
  (wheel/CI clauses outstanding, see Active)
- ✓ **US-009** diagnostics — Kish (1992) and Henry & Valliant (2015)

### Active

- [ ] **US-003**: Large-scale performance target re-benchmarked against a *live* solver and
      restated as an achievable gate (the current statement was written against the removed
      `lbfgsb` and is contradictory: <1 s in §1 vs <2 s in §11)
- [ ] **US-004 residual**: both bridges reach the solvers through ONE dispatch path — today
      `src/r_bridge.cpp:654-899` dispatches on a method string and never calls
      `rk_calibrate` (`leafblower-rywn`, P0)
- [ ] **US-008 residual**: an installable, self-contained Python wheel on Python 3.9–3.13
- [ ] **US-010**: CRAN + PyPI distribution — the largest genuinely-unfinished requirement
- [ ] **KPI**: every row of the §11 success-metrics table has a *live* measuring artefact
      (KPI-02 validated in Phase 1 — see below)
- [ ] Known defects recorded in `.planning/codebase/CONCERNS.md` and tracked in beads

### Out of Scope

- **`grake` (slot 7), `lbfgsb` (slot 2), Epic-K `cp` (slot 12)** — removed or never landed
  per `docs/superpowers/specs/2026-08-14-removed-solver-slots-supersession.md`. Slots are
  permanently reserved. PRD US-006 and FR-20…FR-28 are WITHDRAWN with `lbfgsb`.
- **PRD US-007 two-way auto-routing** — superseded by the eight-solver three-way rule.
- **Multi-threading / GPU** — CPU-only, single-threaded by design; determinism for R↔Python
  parity requires single-thread BLAS, so internal parallelism would break the contract.
- **Variance estimation after calibration** — out of scope by design (`docs/methods/raking.md`).
- **`leafblower-2ouc`** (benchmark study vs competing R/Python calibration software, an
  article) — genuine open work, but a different deliverable with its own beads epic. Its
  `benchmarks/` infrastructure should be REUSED by Phase 3, not duplicated.
- **Line/branch coverage instrumentation** — deliberate: the quality gate is behavioural.
  Do NOT add `--cov-fail-under`; no `covr`/`pytest-cov` is wired and it would only cover
  the thin Python layer.

## Context

- **Repository:** local-only, NO git remote, branch `master`. Work is complete when
  committed locally and the gates are green; there is nothing to push to.
- **Task tracking:** `bd` (beads) is authoritative. 1495 tickets, 1478 closed, 13 open.
  These roadmap phases are a planning layer ABOVE beads, not a replacement — phases
  reference existing ticket IDs and never invent a parallel numbering scheme.
- **Document corpus:** 44 ingested docs (32 SPEC, 11 DOC, 1 PRD, **0 ADR**). The PRD
  (`tasks/prd-leafblower-core.md`, Draft v3, 2026-04-18) is the OLDEST document and
  describes a two-algorithm package; roughly a third of it is superseded or withdrawn. Read
  `.planning/intel/requirements.md` § Superseded / withdrawn before scheduling anything
  from it.
- **Pre-registered ship gates:** the SPEC corpus is dense with numeric, fixture-pinned,
  TIE = NO-SHIP gates. Where work is proposed, the gate is usually already written down in
  `.planning/intel/constraints.md` and should be carried into the plan, not re-invented.
- **Recorded negative results** (do not re-discover the hard way): uncapped damping stress,
  cell-level CBB alpha, fixed-sort `f_eval_sraa`, marginal-residual θ₂, ORIS warm-start of
  Newton-KL, global spectral θ₂.
- **No CI exists.** No `.github/workflows`, no Dockerfile. Gating is local and manual.

## Constraints

- **Philosophy**: SOTA, absolute statistical correctness, numeric stability, NO
  cancellations (`p*(1-p)`, never `p - p*p`), efficiency and speed before simplicity —
  from `CLAUDE.md`, non-negotiable.
- **Definition of Done**: `R CMD INSTALL --preclean .` succeeds + R testthat 0 FAIL +
  Python pytest 0 FAIL under single-thread BLAS (`OMP_NUM_THREADS` /
  `OPENBLAS_NUM_THREADS` / `MKL_NUM_THREADS` = 1, in Python set *before* `import numpy`) +
  stepstone benchmark no regression. R↔Python parity at `rtol=1e-6` is standing.
- **Tech stack**: C++17 core, R (testthat), Python ≥3.9 (pytest, uv-managed venv — never
  bare `pip`, never bare `python`/`pytest`; a stale `~/.local` shadow `.so` gets imported).
- **No LTO**: `-flto` is absent from `configure` and `Makevars`. Only COLD (once-per-solve)
  code may move to a new TU; hot per-iteration loops must stay with their callers.
- **No `-O` level in the R build**: CRAN's `tools:::.check_make_vars` rejects `-O` in
  `PKG_CXXFLAGS`. Only `python/CMakeLists.txt` sets `-O3`.
- **Two build source lists**: R auto-globs `src/*.cpp` (`PKG_SOURCES` is decorative); the
  Python build does NOT glob — a new `.cpp` MUST be added to `CORE_SOURCES` in
  `python/CMakeLists.txt` or the pybind11 link fails.
- **ABI frozen**: `static_assert` tripwires — `sizeof(rk_params_t) == 264`,
  `sizeof(rk_result_t) == 536`. Enum values frozen; slots 2 and 7 permanently reserved.
- **Guarded formulas — do NOT "fix"**: the chebyshev Mehrotra corrector's linear
  `y·Δs_aff` term, the ORIS ALM Newton step `X̃(1−λ+μz)/(1+ρ)`, the lambda-capture ordering
  in `raking.cpp`, the normalize→bounds finalization order, and the SRAA best-iterate
  metric (re-introduced twice).
- **Adding a solver is an 8-step checklist** (`CLAUDE.md`); adding a result field means
  editing both dispatch tables and both bridges.

## Key Decisions

**There is no locked-decision layer.** The ingest corpus contains **0 ADR-class documents**
and no document is marked `locked`, so this table is empty by construction rather than
back-filled from lower-precedence SPEC material. Design decisions live inside SPECs and are
carried, with source attribution and supersession notes, in
`.planning/intel/constraints.md` (87 entries).

Two documents function as decision records without being ADRs, and are applied on date +
explicit self-declaration rather than on a locked flag:

| Document | What it fixes |
|----------|---------------|
| `docs/superpowers/specs/2026-08-14-removed-solver-slots-supersession.md` (Status: Accepted, commit-traced) | Dispositions of `grake` (slot 7), `lbfgsb` (slot 2), `cp` (slot 12); the live eight-solver enum |
| `tasks/prd-leafblower-core.md` § 12 Open Questions | OQ-3 fixed `epsilon = 0.05`, OQ-4 `log_fn` callback for verbose output, OQ-5 Windows CRAN binary as a non-goal — all three still stand |

New decisions taken during roadmap execution get appended here.

**Phase 1 (2026-08-15):** KPI-02 (weight-bound property test) validated for 8 of the 9
non-`AUTO` solvers; `newton_kl` is a documented exception — it *reports* bound violations
(`RK_ERR_NOCONV` + `n_bounds_violated > 0`) rather than clamping them, a shipped contract
(`leafblower-73d7`) that conflicts with KPI-02's literal wording. Decision: pin the
reporting-contract assertion rather than force a `src/` fix inside a test-only phase; open
issue tracked on `leafblower-og7d.5`. R↔Python weight-vector parity now covers all 9
non-`AUTO` solvers (`chebyshev`, `greg`, `oris_soft` added); `logit`'s parity tolerance
tightened from `1e-6` to the uniform `1e-10` (measured divergence: 5.3e-15). Package now
opts into testthat edition 3.

---
*Last updated: 2026-08-15 after doc ingest + codebase map (new-project-from-ingest)*
