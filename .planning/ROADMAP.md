# Roadmap: leafblower

## Overview

leafblower is **already built**. Eight solvers ship over one C++17 core, exposed to R and
Python; 1478 beads tickets are closed. This roadmap therefore plans no features — it takes
the package from "works on the author's machine" to "correct, honest, and installable by a
stranger from CRAN and PyPI".

The journey is five phases in dependency order. First close the verification gaps so the
gate actually detects the divergence it claims to detect (Phase 1) — this is the safety net
for everything after it, and the project's standing TDD discipline requires the test before
the change. Then collapse the two hand-synced dispatch tables and the two divergent build
configurations into one engine, which is what makes the R and Python results the *same*
result rather than two results that happen to agree (Phase 2). Then replace the package's
one unverified headline claim — the performance target, written against a solver that no
longer exists — with a measured, achievable gate (Phase 3). Then make every documented
claim true (Phase 4). Then publish to r-universe and PyPI while keeping the CRAN tarball
proven check-clean but unsubmitted (Phase 5). Then write the vignettes and performance
comparisons a CRAN reviewer and a prospective user both expect (Phase 6). Then submit to
CRAN (Phase 7).

**Scope discipline.** Nothing here re-plans shipped work. Nothing here targets `grake`
(slot 7), `lbfgsb` (slot 2) or Epic-K `cp` (slot 12) — all removed or withdrawn, slots
permanently reserved. `bd` (beads) remains the authoritative task tracker; phases reference
existing ticket IDs and never introduce a competing numbering scheme.

## Tracked outside this roadmap

- **`leafblower-2ouc`** — EPIC: benchmark study vs competing R/Python calibration software
  (an article), with children `.15` (reproducibility bundle) and `.49` (logit `if_*`
  coverage). Genuine open work, different deliverable, own epic structure. Phase 3 must
  **reuse** its `benchmarks/` infrastructure rather than build a parallel harness.

## Phases

**Phase Numbering:**

- Integer phases (1, 2, 3): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

Decimal phases appear between their surrounding integers in numeric order.

- [x] **Phase 1: Verification Coverage Closed** - The gate detects R↔Python divergence and bound violations it currently cannot see (completed 2026-08-15)
- [x] **Phase 2: One Engine, Not Two** - R and Python reach the solvers through a single dispatch path, built the same way (completed 2026-08-15)
- [x] **Phase 3: Honest Performance Gate** - The performance claim is measured against a live solver and restated as an achievable gate (completed 2026-08-15)
- [x] **Phase 4: Truthful Surface** - Every documented claim matches shipped behaviour; every footgun errors instead of silently misleading (completed 2026-08-15)
- [ ] **Phase 5: r-universe + PyPI Release** - A stranger installs leafblower with `install.packages(repos = r-universe)` and `pip install`; the CRAN tarball stays proven check-clean but unsubmitted (SC1-5 complete 2026-08-15; SC6-7 added, not yet done)
- [ ] **Phase 6: Vignettes + Performance Comparisons** - Usage vignettes and a documented performance comparison against competing calibration software, the material a CRAN reviewer and a prospective user both expect
- [ ] **Phase 7: CRAN Submission** - The proven-clean tarball from Phase 5 is actually submitted to CRAN and accepted

## Phase Details

### Phase 1: Verification Coverage Closed

**Goal**: The Definition-of-Done gate actually detects the two failure classes it claims to
cover — R↔Python numerical divergence, and weight-bound violations — before any code that
could cause them is touched.
**Depends on**: Nothing (first phase)
**Requirements**: KPI-02
**Success Criteria** (what must be TRUE):

  1. Every one of the eight shipped solvers, `oris_soft` included, is compared R-vs-Python
     on the weight vector at `rtol=1e-6` in the parametrized parity matrix — a divergence
     in any solver fails the suite.

  2. `raking` and `sinkhorn` are covered by the convergence-rule / `max_error` parity
     checks, so a per-method default-rule divergence between bindings fails the suite
     (today only `logit` has this).

  3. The `logit` parity tolerance is either tightened to `1e-10` like every other method, or
     it carries a comment naming the conditioning mechanism that justifies the three-orders-
     of-magnitude relaxation — no unexplained tolerance survives.

  4. A property-based test asserts `min(w) >= min_weight` and `max(w) <= max_weight` within
     1e-10 across 50 randomly generated datasets (KPI-02's own wording).

  5. `Rscript -e 'testthat::test_dir("tests/testthat")'` runs under the testthat edition the
     project documents, with the `DESCRIPTION` field and `CLAUDE.md` in agreement.
**Beads**: `leafblower-og7d` (epic), `leafblower-x7n8` (P0, gate-collection gap, sequences
first). Filed at plan time from `.planning/codebase/CONCERNS.md` § Test Coverage Gaps:
`leafblower-og7d.1` (SC1 parity matrix), `.2` (SC2 raking/sinkhorn convergence-rule), `.3`
(SC3 logit tolerance), `.4` (KPI-02 property test). Also filed: `leafblower-og7d.5` — a
plan-time finding that `newton_kl` returns weights outside `[min_weight, max_weight]` in
`bounds_mode="unit"`, in tension with KPI-02 and scheduled no earlier than Phase 2.
**Notes for planning**: Moving to testthat 3e changes `expect_equal` to waldo semantics and
deprecates `expect_equivalent` across ~94 test files — expect fallout, and treat any newly
failing assertion as a finding, not noise. New parity tests follow the four-step protocol
documented at `python/leafblower/test_solver_parity.py:7-11`, including the
**precheck-both-sides-converged** step, so a non-converged run fails loudly instead of
silently skipping the comparison.

Measured at plan time (2026-08-15), superseding two assumptions in the criteria above:
SC1's "eight shipped solvers" is **nine** — the live `rk_algorithm_t` has nine non-AUTO
values and the set missing from the parity matrix is three methods (`chebyshev`, `greg`,
`oris_soft`), not one. SC5's edition-3 fallout is **zero failures**: the suite reports
`FAILED=0 ERROR=0 PASSED=1025` under both editions, and the only delta is five deprecation
warnings from the 11 `context()` files, so the "expect fallout across ~94 test files" note
overstates the risk. `expect_equivalent`, `expect_that` and `expect_is` appear in no file.
**Plans:** 4/4 plans complete

Plans:

- [x] 01-01-PLAN.md — Relocate the parity file into the blocking gate, extend the matrix to all nine solvers, retire the unexplained logit tolerance (SC1, SC3, `leafblower-x7n8`)
- [x] 01-02-PLAN.md — Convergence-rule and `max_error` parity for raking and sinkhorn (SC2)
- [x] 01-03-PLAN.md — Property-based weight-bound test over 50 fixed stratified datasets (KPI-02, SC4)
- [x] 01-04-PLAN.md — testthat edition 3, fix-first-then-flip (SC5, `leafblower-og7d`)

### Phase 2: One Engine, Not Two

**Goal**: An R user and a Python user get the same number because they ran the same code,
not because two independently maintained code paths happen to agree today.
**Depends on**: Phase 1
**Requirements**: US-004 (residual)
**Success Criteria** (what must be TRUE):

  1. Adding a solver, a result field, or a routing rule requires editing ONE dispatch site;
     the R bridge reaches the solvers through the same path as the C ABI rather than through
     its own method-string chain.

  2. R and Python binaries of the same source are built at the same optimization level, or
     the asymmetry is documented at the parity assertion as a deliberate, bounded decision.

  3. The per-cell unit-mode water-fill exists once; a fix to `bounds_mode="unit"` cannot
     land in one solver family and miss the other.

  4. A build-list divergence between `src/*.cpp` and `python/CMakeLists.txt:CORE_SOURCES`
     fails the test gate, instead of surfacing as an undefined-symbol link error after the
     R tests are already green.

  5. The full DoD gate is green after the change and the stepstone benchmark shows no
     regression.
**Beads**: `leafblower-rywn` (P0), `leafblower-qzto` (P1); file tickets for the water-fill
duplication and the `CORE_SOURCES` check.
**Notes for planning**: This is the highest-risk phase in the roadmap — it rewires the R
path for all eight solvers. Hard constraints: **no LTO**, so no hot per-iteration code may
move across a TU boundary; the ABI `static_assert` tripwires (264 / 536 bytes) must be
re-measured if any struct changes; enum slots 2 and 7 stay holes and any positional
`alg_names` table must carry them; the RAII-safe `Rf_error` unwinding convention
(`src/r_bridge.cpp:270, 642-648`) must be preserved. Run `LBW_BENCH_GATE=1` on every commit
that touches a TU boundary. The mirrored line ranges are already annotated in the source
(`src/c_api.cpp:458`, `:473`, `:505`) — use them.
**Plans:** 8/8 plans complete

Plans:

- [x] 02-01-PLAN.md — Tracer: shared dispatch table + neutral result struct in `calib_dispatch.hpp`, wired for `sinkhorn` on both bridges (SC1, `leafblower-rywn`)
- [x] 02-02-PLAN.md — Guards: `CORE_SOURCES` build-list sync test, shared finalize-helper verification test, optimization-level asymmetry documented (SC4, SC3, SC2, `leafblower-qzto`)
- [x] 02-03-PLAN.md — Migrate `greg`, `greenkhorn`, `logit` (no R-only fields) onto the shared table (SC1)
- [x] 02-04-PLAN.md — Migrate `chebyshev` (warm-start moved into the shared table) and `raking` (`sraa_demoted`) (SC1)
- [x] 02-05-PLAN.md — Migrate `oris` and `oris_soft` (SRAA / SOR / ALM diagnostics, capacity auto-resolution) (SC1)
- [x] 02-06-PLAN.md — Migrate `newton_kl` (`n_projected_dims`, `lm_mu_final`) (SC1)
- [x] 02-07-PLAN.md — Consolidate AUTO routing, unify the enum-to-name table, collapse the R bridge's method-string chain (SC1)
- [x] 02-08-PLAN.md — Phase gate: single-dispatch-site guard test, full DoD gate, stepstone no-regression, human verification (SC1, SC5, `leafblower-rywn` closed)

### Phase 3: Honest Performance Gate

**Goal**: The package's headline performance claim is a measured fact about a solver that
exists, expressed as a gate a user can run and a maintainer can regress against.
**Depends on**: Phase 2
**Requirements**: US-003, KPI-04
**Success Criteria** (what must be TRUE):

  1. A user reading the package documentation sees a large-scale performance figure that was
     measured on this codebase with a live solver, on a stated input class and machine — not
     a target inherited from the removed `lbfgsb`.

  2. The medium-scale target states ONE number, not the current 1 s / 2 s contradiction,
     and names the artefact that measures it.

  3. The input class on which the composite "<30 s AND <1e-6" gate is structurally
     unachievable (K=20 uniform-random, `M_cell/n = 1.0`, zero compression benefit) is
     documented as a known limit with the reason, rather than left as a silently failing
     promise.

  4. A maintainer can re-run the benchmark that produced every published figure with one
     command, and the KPI table's performance rows each name a live measuring artefact.
**Beads**: `leafblower-kk1.20.4` (the REFRAME decision, with three options already on the
ticket), `leafblower-kk1.20` context.
**Notes for planning**: The reframe is a **decision first, work second** — the ticket
already records the three options (accept raking+SQUAREM as the fast path; define the gate
on stepstone instead of kk1204; track kk1204 separately as algorithmic research). Reuse
`benchmarks/` and the existing stepstone fixtures; do NOT build a parallel harness — the
`leafblower-2ouc` epic owns the comparative-study infrastructure. Determinism protocol is
mandatory: single-thread BLAS, interleaved before/after comparison in one `bench::mark()`
call, never sequential measurement.
**Plans:** 8/8 plans complete

Plans:

- [x] 03-05-PLAN.md
- [x] 03-06-PLAN.md
- [x] 03-07-PLAN.md
- [x] 03-08-PLAN.md

- [x] 03-01-PLAN.md — Tracer: end-to-end medium-scale measurement (`oris_soft` vs `survey::calibrate`, accuracy and bound-compliance alongside every wall time), the D-12 one-command wrapper, and the opt-in `LBW_BENCH_GATE` assertion reading it (US-003, KPI-04)
- [x] 03-02-PLAN.md — Expansion: the 1.58M-row stepstone-fulldata large-scale class (SC1), the K=20 uniform-random known-limit class with `M_cell/n` measured (SC3), and the remaining doc-named competitors `icarus`/`ReGenesees` benchmark-scoped per D-09 (US-003, KPI-04)
- [x] 03-03-PLAN.md — Decide and publish: blocking `checkpoint:decision` on headline metric and gate ceiling (D-06/D-10), the completed hard gate, and `docs/performance.md` (US-003, KPI-04)
- [x] 03-04-PLAN.md — Close the loop: `README.md` one-line claim (D-13), the `< 1 s` / `< 2 s` contradiction retired in place at all three PRD sites (SC2), US-003 and KPI-04 rewritten to name live artefacts (SC4), and the kk1204 assertion re-gated onto `LBW_BENCH_GATE` (US-003, KPI-04)

### Phase 4: Truthful Surface

**Goal**: Nothing documented is untrue, and nothing the API silently swallows stays silent —
a reader can trust the docs and a user cannot get a plausible wrong answer by typo.
**Depends on**: Phase 2
**Requirements**: (none — defect-driven)
**Success Criteria** (what must be TRUE):

  1. No live document attributes the unimplemented outer entropic-proximal-point loop, or
     any other unshipped capability, to ORIS or to a removed solver — the claim the
     `ieppa` → `oris` rename exists to repudiate is gone from `docs/raking.md` §8.2/§12 and
     from anywhere else an audit finds it.

  2. A developer reading `rk_algorithm_t` sees why slot 7 is a hole, as they already do for
     slot 2, and cannot reuse either value by accident.

  3. `harvest(..., weights = w)` raises an informative error naming `design_weights=`
     instead of silently ignoring the argument and returning a plausible unweighted result.

  4. An audit of `README`, `NEWS.md`, `man/` and `docs/` finds no surviving reference to
     `grake`, `lbfgsb` or `cp` as available methods.
**Beads**: `leafblower-05ha` (P1), `leafblower-x2iq` (P2), `leafblower-lj8x` (P1, `weights=`
guard, filed at plan time).
**Notes for planning**: `grep grake` returns the live `grake_norm` **convergence metric**
(`src/greg.cpp:153`, `src/logit_calib.cpp:558`) — a different object, correct as-is, and
must NOT be removed. `docs/raking.md` is a legitimate research report on bounded raking;
only its §8.2/§12 misattribution is wrong, so rewrite the passage under the paper's own
name or delete it — do not delete the document. `man/` is roxygen2-generated: edit the
roxygen block, regenerate, and clean any stray `man/dot-*.Rd` before committing.
**Plans:** 2/2 plans complete

Plans:

- [x] 04-01-PLAN.md — Tracer: `checkpoint:decision` on the D-04 breaking-change confirmation,
  then the `harvest(weights=)` `stop()` guard + RVAL.4 test + `eb79.18` rename (SC3,
  `leafblower-lj8x`)

- [x] 04-02-PLAN.md — Expansion: delete `docs/raking.md` §8.2/§12 misattribution (SC1), annotate
  `rk_algorithm_t` slot 7 in `src/leafblower.h`/`CLAUDE.md` (SC2), re-audit grake/lbfgsb/cp
  (SC4, `leafblower-05ha`, `leafblower-x2iq`)

### Phase 5: r-universe + PyPI Release

**Goal**: A survey analyst who has never seen this repository installs leafblower from
r-universe in R and from PyPI in Python, and gets a working, self-contained package. The
CRAN tarball is proven check-clean so Phase 7 can submit it without redoing this work, but
submission itself is out of scope here.
**Depends on**: Phase 1, Phase 2, Phase 3, Phase 4
**Requirements**: US-010, US-008 (residual), KPI-05, KPI-06
**Success Criteria** (what must be TRUE):

  1. `R CMD check --as-cran` on the built tarball reports 0 errors and 0 warnings, with at
     most the new-submission note, and `cran-comments.md` explains any remaining note. This
     proves CRAN-readiness for Phase 7; it does not submit anything.

  2. The source tarball contains no development artifact — no `.cpp` snapshot copies, no
     nested tarball, no patch scripts, no log, no findings documents — and `git ls-files`
     shows no generated build output tracked in the repository.

  3. `pip wheel python/` produces a wheel that passes `twine check`, and
     `python -c "import leafblower; leafblower.harvest"` succeeds after installing it into a
     clean environment with no separately installed C library.

  4. The wheel imports and calibrates on Python 3.9 through 3.13.
  5. The R `DESCRIPTION` version and `python/pyproject.toml` version cannot silently drift —
     a mismatch is caught by a check, not by a reader.

  6. The package is registered and building on r-universe (a live `<owner>.r-universe.dev`
     feed), and `install.packages("leafblower", repos = "https://<owner>.r-universe.dev")`
     succeeds in a clean R session with no local source checkout.

  7. The sdist and wheels are actually uploaded to PyPI (a real `twine upload`, not just a
     local `twine check`), and `pip install leafblower` succeeds from a clean environment
     against the public PyPI index.
**Beads**: `leafblower-l6h0` (P1), `leafblower-dns3` (P2), `leafblower-kk1.24` (epic),
`leafblower-kk1.24.3` (T19 final gate).
**Notes for planning**: `.Rbuildignore` already excludes `cran-comments.md` — the file must
be created, and the exclusion is correct (CRAN wants it outside the tarball). Two build
constraints bite here: `PKG_CXXFLAGS` must carry **no** `-O` flag or
`tools:::.check_make_vars` rejects the package, and any new diagnostic print must be
`#ifndef LBW_NO_R`-guarded because CRAN forbids R packages writing to stderr. `LAPACK` is a
hard `find_package(... REQUIRED)` in the Python build with no fallback — a wheel cannot be
built without it, which constrains what "self-contained" can mean. `tests/testthat/fixtures`
is `.Rbuildignore`d, so fixture-backed tests must keep their `skip_if(!file.exists(...))`
guards or the CRAN run errors. There is no CI in this repository, so the 3.9–3.13 matrix is
either a new CI pipeline or a documented manual matrix — decide before planning.

Measured at plan-gate time (05-05, 2026-08-15), superseding this phase's original CI-only
limitation: SC1 and SC3/SC4 were closed with **real, executed CI evidence**, not local
structural verification. `.github/workflows/r-check.yml` ran on GitHub Actions
(`https://github.com/davdittrich/leafblower/actions/runs/31908234869`): 0 errors, 0
warnings, 2 NOTEs, both explained in `cran-comments.md` — SC1 met. `.github/workflows/
python-wheels.yml` ran on GitHub Actions
(`https://github.com/davdittrich/leafblower/actions/runs/31908234870`): wheels build, pass
`twine check`, and import + calibrate cleanly on `ubuntu-latest` + `macos-14` (arm64) across
Python 3.9–3.13 — SC3/SC4 met for those two platforms. `macos-13` (x86_64 macOS) was dropped
from the matrix after the runner failed to schedule across 25+ minutes on every CI run on
this account (a GitHub-side Intel-macOS-runner phase-out, not a package defect); x86_64
macOS wheel coverage remains unproven by CI, a residual honestly recorded rather than
silently folded into a blanket "macOS: done" claim. SC2 and SC5 were already closed locally
in 05-01/05-02 and re-confirmed unchanged.
**Plans:** 7 plans (SC1-5 closed by 05-01..05-05; SC6-7 planned in 05-06/05-07, not yet
executed)

Plans:

- [x] 05-01-PLAN.md — Tracer: git-hygiene strays removed, .Rbuildignore extended, local `R CMD check --as-cran` clean, cran-comments.md (SC1, SC2, D-01/D-04/D-05/D-06)
- [x] 05-02-PLAN.md — Version-sync test guarding DESCRIPTION vs. pyproject.toml drift (SC5, D-03)
- [x] 05-03-PLAN.md — R CI: .github/workflows/r-check.yml operationalizing the proven local check + hygiene guard (SC1, SC2, D-02)
- [x] 05-04-PLAN.md — Python wheel CI: cibuildwheel matrix config + local build/repair/import proof across Python 3.9-3.13 (SC3, SC4, D-02)
- [x] 05-05-PLAN.md — Phase gate: full local DoD gate, combined re-verification, human sign-off on scope closure
- [ ] 05-06-PLAN.md — r-universe registration (packages.json entry) + Rd2pdf pre-flight + clean-R-session install proof from the live feed (SC6, D-07/D-08, `leafblower-bl7g`)
- [ ] 05-07-PLAN.md — PyPI Trusted Publishing: `pypi` environment + `build-sdist`/`publish-to-pypi` jobs, human pending-publisher checkpoint, `v0.1.0` tag-publish + clean-env install proof (SC7, D-09/D-10/D-11, `leafblower-ej1n`)

### Phase 6: Vignettes + Performance Comparisons

**Goal**: A user landing on the package's documentation site finds worked usage vignettes
and a documented, reproducible performance comparison against competing R/Python
calibration software — the material both a CRAN reviewer and a prospective adopter expect
before trusting a new package.
**Depends on**: Phase 5
**Requirements**: TBD at planning time
**Success Criteria** (what must be TRUE):

  1. At least one R vignette (`vignettes/*.Rmd`) walks a realistic calibration workflow
     end-to-end and builds cleanly under `R CMD build`/`R CMD check`.
  2. A documented performance comparison against at least the competing packages already
     named in `docs/performance.md` (`survey::calibrate`, `icarus`, `ReGenesees`) exists,
     reusing the `leafblower-2ouc` benchmark infrastructure per this roadmap's scope
     discipline rather than a parallel harness.
  3. The comparison is reproducible with one documented command, consistent with the
     determinism protocol established in Phase 3 (single-thread BLAS, interleaved
     before/after measurement).
**Beads**: TBD at planning time — reuse `leafblower-2ouc` epic children where applicable.
**Notes for planning**: This phase does not re-plan `leafblower-2ouc`; it consumes its
`benchmarks/` infrastructure and existing stepstone fixtures the way Phase 3 did.

### Phase 7: CRAN Submission

**Goal**: The tarball proven check-clean in Phase 5 is actually submitted to CRAN and
accepted, with the vignettes and performance documentation from Phase 6 included.
**Depends on**: Phase 5, Phase 6
**Requirements**: TBD at planning time
**Success Criteria** (what must be TRUE):

  1. `R CMD check --as-cran` is re-run against the tarball including Phase 6's vignettes and
     passes with 0 errors/0 warnings (re-verifies Phase 5 SC1 is still true after Phase 6's
     additions).
  2. The package is submitted via the CRAN submission form / `devtools::release()`, and
     `cran-comments.md` reflects the actual submission (not just local-check readiness).
  3. CRAN accepts the submission, or every reviewer-requested change is tracked to
     resolution and re-submitted.
**Beads**: TBD at planning time.
**Notes for planning**: Do not re-run Phase 5's CI/tarball-hygiene work from scratch —
re-verify it still holds after Phase 6 lands, then submit.

## Progress

**Execution Order:**
Phases execute in numeric order: 1 → 2 → 3 → 4 → 5 → 6 → 7

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Verification Coverage Closed | 4/4 | Complete    | 2026-08-15 |
| 2. One Engine, Not Two | 8/8 | Complete    | 2026-08-15 |
| 3. Honest Performance Gate | 8/8 | Complete    | 2026-08-16 |
| 4. Truthful Surface | 2/2 | Complete    | 2026-08-15 |
| 5. r-universe + PyPI Release | 5/7 | SC1-5 done, SC6-7 pending | — |
| 6. Vignettes + Performance Comparisons | 0/? | Not planned | — |
| 7. CRAN Submission | 0/? | Not planned | — |

---
*Roadmap created: 2026-08-15 from `.planning/intel/` + `.planning/codebase/`*
*Granularity: standard (no `config.json` present) · Phase IDs: sequential*
