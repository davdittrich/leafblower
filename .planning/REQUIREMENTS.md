# Requirements: leafblower

**Defined:** 2026-08-15 (derived from `.planning/intel/requirements.md`, which extracted them
from the single PRD `tasks/prd-leafblower-core.md`, Draft v3, 2026-04-18)
**Core Value:** Calibrated weights that are numerically correct, bound-respecting, and
identical from R and from Python.

## How to read this file

**This is a brownfield project.** Requirements are recorded with their **current status
against the shipped code**, not as work to be scheduled. Only requirements with a residual
gap map to a roadmap phase.

IDs keep the PRD's own `US-xxx` identifiers because the 32-document SPEC corpus cites them
verbatim (`2026-04-23-ieppa-faithful-design.md` §3.1/§4.3/§10 cite "PRD § US-001 AC",
"§ US-005", "§ US-005b"). No new numbering scheme is introduced. Defects are referenced by
their existing **beads** ticket ID — `bd` is the authoritative task tracker.

Status values: `Implemented` · `Partial` · `Superseded` · `Withdrawn`

## v1 Requirements

### Compatibility surface

- [x] **US-001**: An R survey analyst replaces `library(autumn)` with `library(leafblower)`
      and calls `harvest(data, target)` identically, with no other script change.
      *Status: **Implemented**.* All 7 autumn functions exported (`NAMESPACE`); the
      honour / silently-ignore / deprecate parameter matrix is in `R/harvest.R`; FR-30…FR-35
      shipped. The `rake`/`nrake`/`nr` synonyms — which the PRD mapped to the now-removed
      `lbfgsb` — have live destinations: `rake`/`nrake` → `raking`, `nr` → `newton_kl`,
      both with a named warning (`R/harvest.R:975-987`, verified). The `auto_collapse` /
      `add_na_proportion` acceptance criterion ("raises an informative error") is stale —
      both have since shipped, which is superseded scope, not a violated requirement.

- [x] **US-009**: `diagnose_weights()`, `design_effect()` (Kish 1992 and Henry & Valliant
      2015), `effective_sample_size()`, `get_current_miss()`, `weighted_pct()` behave as in
      `autumn`; R implementations are pure R.
      *Status: **Implemented**.* `R/diagnose_weights.R`, `R/design_effect.R`,
      `R/current_miss.R`, `R/weighted_pct.R` + Python equivalents.

### Bounds and numerics

- [x] **US-002**: A survey methodologist specifies `min_weight` so calibrated weights never
      fall below a floor. The package's headline addition over `autumn`.
      *Status: **Implemented**.* `rk_params_t.min_weight`, R and Python surfaces,
      `min_weight >= max_weight` → `RK_ERR_BADARG`. The AC "`min_weight > 0` forces iEPPA"
      is superseded — routing now spans eight solvers.

- [x] **US-005**: Capacity-constrained solver with cell-compressed `O(M_cell·K)` inner
      cost: sort-based cell table, log-space Sinkhorn factors with the 700 clip, capacity
      clamp `X[c] = clamp(X̃[c]·W[c], L_c, U_c)`, `errRp < tol_abs` convergence, the
      `S_j < 1e-15·W` empty-category test, NA/OOV skip, and the every-inner-step bounds
      **invariant**.
      *Status: **Implemented** as `oris` (`src/oris.cpp`, enum value 1 unchanged).*
      FR-11 / FR-12 / FR-14 (the outer entropic-proximal-point loop) are **historical**: not
      implemented, mathematically inert at `C = 0`, and the stated reason for the
      `ieppa` → `oris` rename. Citing arXiv:2011.14312 over-claims.

- [x] **US-005b**: `method="raking"` provides classical cyclic IPF with box and hyperplane
      projection as an alternative to ORIS; `RK_ALG_RAKING` in the enum; pre-rev2 tests act
      as a regression guard.
      *Status: **Implemented**.* Projection geometry moved from additive Boyle-Dykstra to
      multiplicative KL-Bregman Dykstra at cell level
      (`2026-04-27-raking-bregman-dykstra-design.md`); existence and routing unaffected.

### Engine contract

- [x] **US-004**: A stable C API (`leafblower.h`) lets any language with C FFI call the
      engine without duplicating algorithm code — **both the R `.Call()` bridge and the
      Python pybind11 module route through the same shared dispatch table.**
      *Status: **Implemented** (leafblower-rywn, closed; Phase 2 plans 01-08 — SC1 now defended by `test_single_dispatch_site.py`).* FR-1…FR-10 shipped: C99-valid header, the full FR-4 validation
      set checked before any weight is modified, `snprintf`-only messages, the `log_fn` /
      Python GIL trampoline contract, caller-owned memory, reentrancy, and the
      `max_error` reporting formula. Base status codes are EXTENDED (not replaced) by
      `RK_ERR_BUDGET=4` / `RK_ERR_STALL=5`.
      **The one-symbol clause is satisfied via a shared dispatch table, not a literal shared
      `rk_calibrate()` call** — a deliberate architectural decision recorded during Phase 2
      research: `rk_result_t` (the C-ABI struct `rk_calibrate()` fills) is missing 4 fields
      R's `harvest()` result exposes and existing tests assert on (`n_projected_dims`,
      `lm_mu_final`, `sraa_demoted`, `convergence_stall_kind`), so a literal shared-symbol
      call would silently drop them. `lbw::dispatch_solver()` and `lbw::route_auto()`
      (`src/calib_dispatch.hpp`) are the single {enum → solver → result} table and single
      AUTO-routing decision both `src/r_bridge.cpp`'s `C_rk_calibrate` and
      `src/c_api.cpp`'s `rk_calibrate()` call — "same path" now means same dispatch table +
      same neutral result-extraction helper, each caller marshaling it into its own output
      shape (SEXP list vs. ABI-frozen `rk_result_t`). `r_bridge.cpp` has zero per-method
      `strcmp` branching left (plan 02-07).

- [x] **US-008**: A Python survey analyst calls `leafblower.harvest(df, targets)` with a
      pandas DataFrame over the same compiled core as R.
      *Status: **Implemented** (Phase 5, 05-04/05-05).* FR-36…FR-40 shipped, including the
      copy-never-view contract, the `float64` + C-contiguity enforcement, and the
      `convergence` dict semantics. `.github/workflows/python-wheels.yml` ran for real on
      GitHub Actions (`https://github.com/davdittrich/leafblower/actions/runs/31908234870`):
      wheels build, pass `twine check`, and `import leafblower; leafblower.harvest` succeeds
      after installing into a clean environment on `ubuntu-latest` (manylinux) and `macos-14`
      (arm64), across Python 3.9–3.13. **Residual, stated honestly:** the `macos-13` (x86_64
      macOS) runner never scheduled across 25+ minutes on every CI run on this account (a
      GitHub-side Intel-macOS-runner phase-out, not a package defect) and was dropped from
      the matrix — x86_64 macOS wheel coverage is unproven by CI, though the underlying
      build/repair mechanics were proven once locally without Docker in 05-04. US-008's PRD
      text does not name a specific macOS architecture, so this residual does not block
      closing US-008, but it is the one platform this phase cannot claim to have proven.
      *Note:* FR-38's stated mechanism (CMake EXCLUDES `r_bridge.cpp`) is the inverse of the
      live one (CMake does not glob; `CORE_SOURCES` is an explicit include list). Same
      outcome, different mechanism.

### Performance

- [ ] **US-003**: 1M+ observations across 20+ margins calibrate in under 30 seconds,
      single-threaded, so a census microsimulation researcher can iterate on synthetic
      population models. `verbose = 1` prints the selected algorithm and routing reason.
      *Status: **Partial**.* Routing observability shipped (not re-exercised by Phase 3 —
      that clause's evidence predates this phase and is carried forward as-is). Phase 3
      measured `oris_soft` against a live solver on three input classes, transcribed in
      `docs/performance.md`: `medium_100k_5margins` (100K rows/5 margins, wall_s=0.0427,
      max_error=3.35e-05), `large_stepstone_fulldata` (1,582,732 rows/9 margins, a real
      salary-survey fixture, wall_s=3.5459 — clears the <30s large-scale budget on a real
      shape exceeding the PRD's 1M-row target), and `known_limit_k20_uniform` (500,000
      rows/20 margins, `m_cell/n = 1.0000`). No single measured fixture combines "1M+ rows
      AND 20+ margins" as the PRD literally states — the 20-margin shape is measured
      separately and found **structurally unable to clear an accuracy floor**:
      `max_error = 5.229e-03` at wall_s=7.3934 under a bounded 500-iteration budget,
      `m_cell/n = 1.0000` meaning every one of the 500,000 observations lands in its own
      cell, so ORIS's cell-compression advantage yields nothing at this shape (confirmed
      by `docs/investigations/2026-04-23-kk1204-convergence.md`, commit `3effd3a`, and
      `leafblower-ylsy`'s closed research — cited only, not reopened, per D-03). This is
      now a documented known limit (`docs/performance.md` § Known limit), not an open
      blocker: the "<30 s AND <1e-6" composite gate on K=20 uniform-random input is
      confirmed structurally unachievable with any solver currently implemented, and
      `leafblower-kk1.20.4`'s REFRAME resolved on **dropping it as the headline basis**
      (D-01) — not on any of the ticket's three originally-worded options. Live measuring
      artefacts: `benchmarks/oris_soft_vs_competitors.R` (the measurement) and the
      `honest gate:` assertion in `tests/testthat/test-bench-gate.R` (the regression gate,
      `LBW_BENCH_GATE=1`), replacing the void `test-lbfgsb.R`.

### Distribution

- [x] **US-010**: leafblower is distributable on CRAN and PyPI so users install it with
      `install.packages()` and `pip install`.
      *Status: **Implemented** (Phase 5, 05-01 through 05-05).*
      Shipped: `configure` with C++17→C++14 fallback, no vendored dependency > 5 MB,
      `pyproject.toml` on scikit-build-core, a hygiene-clean git tree (11 tracked strays +
      120-file `leafblower.Rcheck/` removed, `leafblower-l6h0` closed), `cran-comments.md`
      (exists, documents the real CI result), a wheel that builds and passes `twine check`
      on CI (`leafblower-kk1.24.3` closed in substance), and a version-sync regression test
      (`python/leafblower/test_version_sync.py`) catching `DESCRIPTION`/`pyproject.toml`
      drift automatically rather than by hand. `R CMD check --as-cran` ran for real on
      GitHub Actions (`https://github.com/davdittrich/leafblower/actions/runs/31908234869`):
      0 errors, 0 warnings, 2 NOTEs, both explained in `cran-comments.md`. **Residual, stated
      honestly:** x86_64-macOS (`macos-13`) wheel coverage is unproven by CI (see US-008);
      this repository is a GitHub + r-universe intermediate release rather than a literal
      CRAN web-form submission (D-05, `cran-comments.md`'s own § Submission type), so
      `install.packages()` from the canonical CRAN mirror specifically is not yet live —
      the check that gates a real CRAN submission (`R CMD check --as-cran`, 0 errors/0
      warnings) is proven, the submission act itself is a separate, later step outside this
      phase's scope.
      The § 7 `PKG_CXXFLAGS ... -O3` statement is superseded — the R build sets no `-O`
      level by design; only `python/CMakeLists.txt` sets `-O3`.

### Project-level acceptance (PRD § 11 KPI table)

- [x] **KPI-01**: R API compat — all 7 autumn functions present. *Measured by `R CMD check`
      + `test-harvest.R`.* **Implemented.**
- [ ] **KPI-02**: Weight bound enforcement — `max(w) <= max_weight` and
      `min(w) >= min_weight` within 1e-10, **over 50 random datasets**, by a property-based
      test in `test-harvest.R`. **Not located.** Bounds are covered by ~6 targeted tests
      (`test-clamp-contract.R`, `test-harvest-bounds-mode.R`, `test-oris-bounds-mode.R`,
      `test-cr-d16-nbounds.R`, `test-unit-bounds-status-consistency.R`,
      `test-newton-bounds-write-guard.R`); no 50-dataset property test was found under
      `tests/testthat/`.
- [x] **KPI-03**: Convergence — `max_k max_j |Σw·1[g_k=j]/Σw − τ_j^(k)| < 1e-6` at reported
      convergence. **Implemented** (`check_convergence`, `src/calib_dispatch.hpp:204`).
- [ ] **KPI-04**: Large-scale — 1M rows, 20 margins, `max_weight = 3`, < 30 s.
      *Status: **Partial**.* Measured on the real 1,582,732-row/9-margin
      `large_stepstone_fulldata` class (`max_weight = 3`): `oris_soft` calibrates in
      wall_s=3.5459, comfortably inside the 30s budget, on a fixture larger than the PRD's
      literal 1M-row target. The 20-margin corner of the original target
      (`known_limit_k20_uniform`, 500,000 rows/20 margins) is measured separately and
      documented as a known limit (`m_cell/n = 1.0000`, `max_error = 5.229e-03` at
      wall_s=7.3934 — does not clear an accuracy floor at this shape; see
      `docs/performance.md` § Known limit). Live artefacts: `benchmarks/oris_soft_vs_competitors.R`
      (measurement) and the `honest gate:` assertion in `tests/testthat/test-bench-gate.R`,
      run via `LBW_BENCH_GATE=1 CI=1 NOT_CRAN=true Rscript -e "testthat::test_dir('tests/testthat',
      filter='bench-gate', stop_on_failure=TRUE)"`. No other performance-adjacent row in
      this KPI list was found still naming a non-existent measuring artefact (checked
      KPI-01 → `test-harvest.R`, KPI-02 → property test in `test-harvest.R`, KPI-03 →
      `check_convergence` in `src/calib_dispatch.hpp:204`, KPI-05 → `R CMD check --as-cran`,
      KPI-06 → no CI pipeline exists at all, which is KPI-06's own recorded open status, not
      a stale artefact reference — all name something that exists or is honestly marked
      absent, not a void file).
- [x] **KPI-05**: CRAN check — 0 errors, 0 warnings, via `R CMD check --as-cran`.
      **Implemented** (Phase 5, 05-05). Real GitHub Actions run:
      `https://github.com/davdittrich/leafblower/actions/runs/31908234869` — 0 errors, 0
      warnings, 2 NOTEs (`-mavx2` compilation flags, HTML-manual `tidy` absent), both
      explained in `cran-comments.md`.
- [x] **KPI-06**: Python wheel — installs on Linux/macOS, Python 3.9–3.13, via a CI matrix.
      **Implemented** (Phase 5, 05-04/05-05). Real GitHub Actions run:
      `https://github.com/davdittrich/leafblower/actions/runs/31908234870` — wheels build,
      pass `twine check`, import + calibrate cleanly on `ubuntu-latest` (manylinux) and
      `macos-14` (arm64), Python 3.9–3.13. **Residual:** `macos-13` (x86_64 macOS) dropped
      from the matrix — the runner never scheduled across 25+ minutes on every CI run on
      this account (GitHub-side Intel-macOS-runner phase-out); x86_64 macOS coverage is
      unproven by CI, though the wheel-build/repair mechanics were proven locally without
      Docker in 05-04. KPI-06's literal text ("Linux/macOS") does not name an architecture,
      so this residual does not block closing it, but is recorded rather than silently
      generalized away.
- **KPI-07** (scoping statement, not work): performance benchmarks are CI **artifacts, not
      gates** (§ 8); only the phase gates are binding.

## Known Defects (beads-tracked, not PRD-derived)

Sourced from `.planning/codebase/CONCERNS.md` and the live `bd` queue. **These keep their
beads IDs.** The roadmap schedules them; it does not renumber them.

| Ticket | P | Defect |
|--------|---|--------|
| `leafblower-rywn` | P0 | Two hand-synced solver dispatch tables — `src/r_bridge.cpp:654-899` (string) never calls `rk_calibrate`; `src/c_api.cpp:414+` (enum) |
| `leafblower-qzto` | P1 | `-O3` in `python/CMakeLists.txt:99` vs no `-O` R-side — parity compares differently-optimized binaries |
| `leafblower-og7d` | P1 | testthat runs edition 2; no `Config/testthat/edition: 3` in `DESCRIPTION` despite the "testthat v3" claim in `CLAUDE.md` (verified absent) |
| `leafblower-l6h0` | P1 | `.Rbuildignore` gaps ship dev artifacts in the source tarball (relates to US-010) |
| `leafblower-05ha` | P1 | `docs/raking.md` §8.2/§12 attribute the *unimplemented* outer entropic-proximal loop to ORIS — the exact claim the rename exists to repudiate |
| `leafblower-x2iq` | P2 | Enum slot 7 (ex-`RK_ALG_GRAKE`) is an undocumented hole; slot 2 has a comment, slot 7 does not |
| `leafblower-dns3` | P2 | Generated build artifacts tracked in git |
| `leafblower-kk1.20.4` | P1 | Phase-2 performance gate needs a REFRAME (see US-003) |
| `leafblower-kk1.24` / `.24.3` | P2 | Phase 4: Distribution — CRAN submission package + PyPI wheel artifacts |

Un-ticketed defects from `CONCERNS.md` picked up by the roadmap (file a beads ticket at
plan time — do not invent an ID here):

- Duplicated unit-mode water-fill: `src/calib_dispatch.hpp:359-383` vs
  `src/oris_finalize.cpp:19`, the second declaring itself a "mirror" of the first.
- No automated check that `src/*.cpp` and `python/CMakeLists.txt:CORE_SOURCES` agree.
- `oris_soft` is the one shipped solver absent from the R↔Python weight-parity matrix
  (`tests/test_parity_weights.py:73`) — **High** priority per the concerns audit.
- `raking` and `sinkhorn` absent from `python/leafblower/test_solver_parity.py` (the
  convergence-rule / `max_error` parity checks).
- Unexplained parity tolerance asymmetry: `tol = 1e-6 if method == "logit" else 1e-10`
  (`tests/test_parity_weights.py:93`), three orders of magnitude with no justifying comment.
- `harvest()` silently ignores `weights=` (the real argument is `design_weights=`),
  producing a plausible unweighted result.
- No version-sync check between `DESCRIPTION` and `python/pyproject.toml`.

## Superseded / Withdrawn — do NOT plan work for these

| Item | Disposition |
|------|-------------|
| **US-006 and FR-20…FR-28 (L-BFGS-B)** | **WITHDRAWN entirely.** Solver removed; `RK_ALG_LBFGSB = 2` is a permanently reserved hole (`src/leafblower.h:44`). Goes with it: Deville-Särndal link selection, the dual `φ(λ)` formulation, the L-BFGS 2-loop recursion and `lbfgs_m`, the Wolfe line search, the FR-28 final clamp. The Deville-Särndal logit *distance* survives in the separate `logit` solver (`RK_ALG_LOGIT = 10`) — a different solver, not a survival of US-006. |
| **PRD § 6 enum** `AUTO=0, IEPPA=1, LBFGSB=2` | **SUPERSEDED.** This PRD is the ORIGIN of the slot-2 hole. Live enum: AUTO 0, ORIS 1, RAKING 3, SINKHORN 4, CHEBYSHEV 5, GREG 6, ORIS_SOFT 8, GREENKHORN 9, LOGIT 10, NEWTON_KL 11. Values frozen (`77d0614`). Slots 2 and 7 stay holes; any positional `alg_names` table must carry them. |
| **US-007 two-way auto-routing** | **SUPERSEDED** by the three-way rule (`M_cell/n >= 0.9`, `K >= 5`, `target_skew` vs 5.0), implemented at both dispatch sites. FR-7's "never `RK_ALG_AUTO`" invariant survives; its enumeration of possible values does not. |
| **`grake` (slot 7)** | Removed pre-release, commit `9a67891`. `src/grake.{cpp,hpp}` do not exist. The live `grake_norm` **metric** is unrelated and unaffected. No grake acceptance criterion (A4 included) carries forward. |
| **Epic-K `cp` (slot 12)** | Landed `00a3f10`, reverted `3fac1d6`, same day. Never live. A **withdrawn proposal**, not pending work. Its T1…T8 gates do not carry forward. |
| **US-001 `rake`/`nr` → `lbfgsb`** | **SUPERSEDED**, and already resolved in code with different destinations (see US-001). |
| **US-005 paper-faithfulness (FR-11/12/14)** | **Repudiated.** See US-005. |
| **US-005b additive Boyle-Dykstra geometry** | **SUPERSEDED** by multiplicative KL-Bregman Dykstra at cell level. |
| **PRD § 5 non-goals: SQUAREM, bounded-IPF water-filling, `auto_collapse`, `add_na_proportion`** | **All four have since shipped or been worked.** Superseded *scope*, not violated requirements. CPU-only and single-threaded remain in force. |
| **PRD § 7 `PKG_CXXFLAGS ... -O3`** | **SUPERSEDED** — the R build sets no `-O` level by design. `PKG_SOURCES` in `Makevars.in` is likewise decorative. |
| **PRD § 9 Risks / § 10 Phased Rollout** | Phase 1 is L-BFGS-B and is void; Phases 2–4 describe delivered work. The TDD RED-phase discipline they encode outlives the phase plan and is the project's standing methodology. |
| **`iEPPA` throughout the PRD** | Terminology for the solver now named `oris`. Enum value 1 unchanged. Not a conflict. |
| **PRD § 1 Success Criteria, "Performance — medium ... < 1 s"** | **SUPERSEDED** (Phase 3, 2026-08-15) — written against the withdrawn L-BFGS-B solver and contradicted §11's own `< 2 s` row for the same shape. Live measured figure: `docs/performance.md` (`oris_soft`/`medium_100k_5margins`, wall_s=0.0427). |
| **PRD § US-006 AC, "converges on: 100K rows, 5 margins within 1 s"** | **SUPERSEDED** (Phase 3, 2026-08-15) — goes with US-006's whole-section withdrawal (L-BFGS-B never implemented, slot 2 permanently reserved); no live solver to measure this acceptance criterion. |
| **PRD § 11 KPI, "L-BFGS-B convergence ... `test-lbfgsb.R` Phase 1 gate"** | **SUPERSEDED** (Phase 3, 2026-08-15) — `test-lbfgsb.R` does not exist. Live artefacts: `benchmarks/oris_soft_vs_competitors.R` and the `honest gate:` assertion in `tests/testthat/test-bench-gate.R`. |

## Traceability

Shipped requirements carry no phase — the work predates this planning layer.

| Requirement | Phase | Status |
|-------------|-------|--------|
| US-001 | — (shipped pre-GSD) | Implemented |
| US-002 | — (shipped pre-GSD) | Implemented |
| US-003 | Phase 3 | Partial |
| US-004 | Phase 2 (residual only) | Implemented |
| US-005 | — (shipped pre-GSD) | Implemented |
| US-005b | — (shipped pre-GSD) | Implemented |
| US-008 | Phase 5 (residual only) | Implemented |
| US-009 | — (shipped pre-GSD) | Implemented |
| US-010 | Phase 5 | Implemented |
| KPI-01 | — (shipped pre-GSD) | Implemented |
| KPI-02 | Phase 1 | Open |
| KPI-03 | — (shipped pre-GSD) | Implemented |
| KPI-04 | Phase 3 | Partial |
| KPI-05 | Phase 5 | Implemented |
| KPI-06 | Phase 5 | Implemented |
| US-006, FR-20…FR-28 | — | **Withdrawn** |
| US-007 | — | **Superseded** |

Defect traceability:

| Ticket / defect | Phase |
|-----------------|-------|
| `oris_soft` parity gap, `raking`/`sinkhorn` rule parity, logit tolerance asymmetry | Phase 1 |
| `leafblower-og7d` (testthat edition) | Phase 1 |
| `leafblower-rywn` (P0, dual dispatch) | Phase 2 |
| `leafblower-qzto` (`-O` asymmetry) | Phase 2 |
| Duplicated water-fill, `CORE_SOURCES` sync check | Phase 2 |
| `leafblower-kk1.20.4` (gate reframe) | Phase 3 |
| `leafblower-05ha` (docs/raking.md), `leafblower-x2iq` (slot 7), `weights=` guard | Phase 4 |
| `leafblower-l6h0`, `leafblower-dns3`, `leafblower-kk1.24`, `leafblower-kk1.24.3`, version-sync check | Phase 5 |

**Coverage:**
- v1 requirements: **16** (10 PRD entries, of which the KPI entry is expanded into its 6
  measurable §11 rows)
- Fully satisfied by shipped code, no phase needed: **9**
- With residual scope, mapped to a phase: **7**
- Unmapped: **0** ✓
- Withdrawn / superseded, explicitly not planned: **US-006 + FR-20…FR-28, US-007**, plus the
  12 items in the table above

---
*Requirements defined: 2026-08-15*
*Source: `.planning/intel/requirements.md` · statuses cross-checked against
`.planning/codebase/` and live source*
