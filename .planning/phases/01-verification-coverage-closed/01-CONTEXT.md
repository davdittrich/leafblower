# Phase 1: Verification Coverage Closed - Context

**Gathered:** 2026-08-15
**Status:** Ready for planning

<domain>
## Phase Boundary

Make the Definition-of-Done gate actually detect the two failure classes it claims to cover
— R↔Python numerical divergence, and weight-bound violations — before Phase 2 touches the
code that could cause them.

**Test-layer only.** No solver, binding, or `src/` change lands in this phase. If a real
numerical divergence is found, it is measured, documented, and ticketed — not fixed here.
The one exception is `DESCRIPTION` (the testthat edition field) and test-file relocation.

</domain>

<decisions>
## Implementation Decisions

### Bound property test (SC4)

- **D-01:** The 50-dataset bound property test lives in **R testthat only**, not Python and
  not both. Bounds are enforced in the shared C++ core (`lbw::finalize_weights`), so one
  binding exercises the same code path; `tests/testthat/` already holds six targeted bound
  tests for it to sit beside. No new dependency.
- **D-02:** **Fixed seed, 50 fixed datasets.** The same 50 every run. Determinism is a hard
  project constraint (the same reason single-thread BLAS is mandatory) and a bound failure
  must be replayable without hunting a seed. Accepted cost: it is a fixed probe, not a
  search.
- **D-03:** **Design weights are NOT normally distributed.** Generate them as a
  **mixture — lognormal bulk plus a small heavy-tailed contaminant fraction.** Gaussian
  weights never push the calibrated result against the clamps, so a bound assertion on
  normal data proves nothing. The mixture models the real case: mostly well-behaved weights
  with a few extreme influentials, which is precisely what `min_weight`/`max_weight` exist
  to control.
- **D-04:** The 50 datasets are **stratified into three groups**: skewed-weights-only,
  sparse-cells-only (skewed category marginals, near-empty cells), and both. Each violation
  stays attributable to one axis while the interaction case still gets probed. Sparse cells
  are the documented `RK_ERR_INFEAS` trigger and a second independent failure axis.
- **D-05:** **Bounds are asserted unconditionally** — on every returned weight vector,
  regardless of status. No convergence precheck gates the assertion. PRD US-005 AC and FR-19
  require the invariant "at every inner iteration, not only at convergence", and docs state
  weights reflect the last iterate on `RK_ERR_NOCONV`. This is the strictest reading and
  catches a clamp that leaks only on the failure path. — **Reversibility:** reversible —
  the assertion can be softened to a precheck later without touching production code.

### logit parity tolerance (SC3)

- **D-06:** **Diagnose the conditioning before touching the number.** Do not tighten
  `tol = 1e-6 if method == "logit" else 1e-10` blind. Establish WHY logit would differ first
  — `src/logit_calib.cpp` runs a Newton solve with an Armijo line search, and FMA
  contraction at two sites is a known ~5e-8 scale source — then set the tolerance to match
  the identified mechanism and comment it.
- **D-07:** **If a real divergence above 1e-10 is found, file a ticket and keep this phase
  test-only.** Record the measured delta and suspected mechanism as its own beads ticket;
  set the tolerance to a documented value referencing it. Do NOT pull a solver/binding fix
  into Phase 1 — the suspected cause may itself be the dual dispatch path that Phase 2
  (`leafblower-rywn`) unifies.

### Parity harness (SC1, SC2)

- **D-08:** **`tests/test_parity_weights.py` moves to `python/leafblower/`** so a single
  pytest invocation collects everything and the gate command stays unchanged. Rejected:
  widening the gate command to sweep two roots. `python/pyproject.toml:34`
  (`wheel.exclude = ["leafblower/test_*.py"]`) already keeps such files out of the shipped
  wheel, so the relocation costs nothing on the artifact axis. — **Reversibility:** costly —
  the file's `REPO_ROOT`-relative path to the `tests/parity/*.R` helpers must be recomputed;
  getting it wrong makes the tests skip rather than fail.
- **D-09:** The gate-collection gap is **its own ticket — `leafblower-x7n8` (P0)**, filed
  during this discussion, not folded into the phase plan. It is a distinct defect class from
  the four coverage gaps in `CONCERNS.md`: existing, passing tests excluded from the blocking
  gate. It is P0 and sequences FIRST — adding `oris_soft` to a matrix the gate does not
  collect is a no-op until it lands.

### testthat edition (SC5)

- **D-10:** **Adopt edition 3.** Add `Config/testthat/edition: 3` to `DESCRIPTION`, making
  CLAUDE.md's existing "testthat v3" claim true. Rejected: staying on 2e and correcting the
  docs. waldo's strict comparison is the semantics a package whose stated philosophy is
  "absolute statistical correctness" should be running. — **Reversibility:** reversible —
  removing the DESCRIPTION field reverts the suite to 2e semantics.
- **D-11:** **Every newly-failing waldo assertion is a finding, not noise.** Investigate
  each; a test that only passed under `all.equal`'s tolerance was hiding something. Fix the
  code or document why the looseness is correct. No blanket "3e artifact" triage bucket —
  that judgment is where a real bug gets mislabeled.
- **D-12:** **Fix first, then enable.** Migrate deprecated constructs while still on 2e and
  green, then flip the edition last. Every commit stays green (satisfying "complete =
  committed locally + gates green"), and the flip commit isolates the true waldo fallout
  from mechanical deprecation churn.

### Claude's Discretion

- Exact lognormal sigma, contaminant fraction, and contaminant tail parameters (D-03) —
  derive from first principles so the mixture actually reaches the clamps at the default
  `max_weight = 5`; do not tune values to make tests pass.
- Dataset counts per stratum in D-04's three-way split.
- Which specific R RNG constructs implement the fixed-seed generator.

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Gate definition and test protocol
- `.coverage-thresholds.json` §`enforcement.command` — the blocking Definition-of-Done
  command, verbatim. The authority on what the gate actually runs. `"blocking": true`.
- `python/leafblower/test_solver_parity.py:1-25` — the four-step parity protocol
  (explicit convergence spec both sides → subprocess Rscript → **precheck both converged** →
  `np.allclose`). New parity tests follow it. Also documents why
  `test_logit_default_rule_parity` deliberately passes NO explicit rule.
- `tests/test_parity_weights.py` — the weight-vector parity file being relocated. Line 73
  is the parametrize list; line 93 is the logit tolerance ternary.
- `tests/parity/` — R-side helper scripts the parity tests shell out to
  (`run_oris_soft_r.R`, `run_chebyshev_r.R`). The relocation must keep resolving these.
- `conftest.py` (repo root) — strips `python/` from `sys.path` so the installed wheel is
  imported instead of the source tree (which lacks `_leafblower.so`). Any relocation must
  preserve this import behaviour.

### Coverage gaps this phase closes
- `.planning/codebase/CONCERNS.md` §"Test Coverage Gaps" — the four documented gaps:
  `oris_soft` outside the parity matrix, `raking`/`sinkhorn` absent from
  `test_solver_parity.py`, the unexplained logit tolerance, opt-in bench gate.
- `.planning/codebase/TESTING.md` — suite layout, run commands, single-thread BLAS
  requirement.
- `.planning/ROADMAP.md` §"Phase 1" — the five success criteria this context serves.

### Requirements
- `.planning/REQUIREMENTS.md` — KPI-02 (the bound property requirement, 50 random datasets,
  1e-10) is this phase's mapped requirement.
- `tasks/prd-leafblower-core.md` §US-005 AC, §FR-19 — the source of D-05: weights satisfy
  `[min_weight, max_weight]` at every inner iteration, not only at convergence. NOTE: this
  PRD is dated 2026-04-18 and is heavily superseded on algorithms — see below.
- `docs/superpowers/specs/2026-08-14-removed-solver-slots-supersession.md` — `grake`,
  `lbfgsb`, and Epic-K `cp` are removed/withdrawn. No test in this phase may target them.
  The PRD's US-006 and FR-20…FR-28 are withdrawn.

### Project constraints
- `CLAUDE.md` — build gate (`R CMD INSTALL --preclean .`, NOT devtools), uv-managed Python
  (no bare pip), single-thread BLAS trio, no-cancellation rule, explicit-pathspec commits.

### testthat 3e
- https://testthat.r-lib.org/articles/third-edition.html — waldo semantics, deprecations.

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `python/leafblower/test_solver_parity.py` — the four-step protocol and its embedded
  fixture-as-constants pattern (no filesystem side-effects). Reuse the protocol for the new
  `raking`/`sinkhorn` convergence-rule checks rather than inventing a second pattern.
- `tests/parity/run_*_r.R` — existing R-side subprocess helpers. Adding a solver to the
  matrix follows the same shape.
- `tests/testthat/` — six existing targeted bound tests give the assertion style for D-01.

### Established Patterns
- Parity tests fix the stopping point with an explicit convergence spec on BOTH sides
  (`rule="improvement"`, `tol=0.001`), so the same spec yields byte-identical weights from
  the shared core. Only the deliberate default-rule test omits it.
- `LBW_BENCH_GATE=1` is the repo's precedent for an env-gated heavier check — the shape to
  copy if any opt-in sweep is ever added.
- Test files live INSIDE `python/leafblower/` precisely so the editable install discovers
  them; `wheel.exclude` keeps them out of the artifact.

### Integration Points
- `.coverage-thresholds.json` `enforcement.command` is the single blocking gate — anything
  not collected by it is not enforced, regardless of whether it passes.
- `DESCRIPTION` `Config/testthat/edition` governs the semantics of all 94 test files at once.

### Measured facts (verified 2026-08-15, direct collection runs)
- `cd python && pytest --collect-only -q` → **141 tests**, all `leafblower/`-prefixed.
- `pytest tests/ --collect-only -q` → **8 tests**, none of which the gate runs.
- 94 files in `tests/testthat/`. `context()` appears in **13**; `expect_that`, `expect_is`,
  and `expect_equivalent` appear in **zero**. The 3e fallout is therefore mostly
  `expect_equal` under waldo, which cannot be grepped for — it surfaces only at run time.
- `hypothesis` is not a declared Python dependency.

</code_context>

<specifics>
## Specific Ideas

- "Data must not be generated by a normal distribution. Data needs to be skewed,
  fat-tailed." — the originating constraint behind D-03/D-04. Applies to BOTH axes: the
  design-weight distribution and the category marginals.

</specifics>

<deferred>
## Deferred Ideas

- **Adding `oris_soft` to the parametrized weight-parity matrix** stays in this phase (SC1),
  but is explicitly NOT part of `leafblower-x7n8` — that ticket is a pure relocation. Order:
  relocate first, extend second.
- **Making the stepstone bench gate non-opt-in** (`CONCERNS.md` flags it Medium) — belongs
  with Phase 2, which is the phase that moves TU boundaries under the no-LTO constraint.
- **Line/branch coverage instrumentation** — deliberate project non-goal, recorded in
  `CONCERNS.md` so it is not mistaken for an oversight. Not this phase, not any phase
  without an explicit decision to reverse it.
- **Fixing any real logit divergence found under D-06/D-07** — ticketed, executed after
  Phase 2 at the earliest.

</deferred>

---

*Phase: 1-Verification Coverage Closed*
*Context gathered: 2026-08-15*
