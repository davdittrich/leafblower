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
      *Status: **Implemented** (leafblower-rywn, Phase 2 plans 01-07).* FR-1…FR-10 shipped: C99-valid header, the full FR-4 validation
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

- [ ] **US-008**: A Python survey analyst calls `leafblower.harvest(df, targets)` with a
      pandas DataFrame over the same compiled core as R.
      *Status: **Partial**.* FR-36…FR-40 shipped, including the copy-never-view contract,
      the `float64` + C-contiguity enforcement, and the `convergence` dict semantics.
      **Outstanding:** `pip install leafblower` installing a self-contained wheel, and the
      Python 3.9–3.13 test matrix — no wheel artefact has ever been built
      (`python/dist/` and `python/wheelhouse/` do not exist) and no CI exists.
      *Note:* FR-38's stated mechanism (CMake EXCLUDES `r_bridge.cpp`) is the inverse of the
      live one (CMake does not glob; `CORE_SOURCES` is an explicit include list). Same
      outcome, different mechanism.

### Performance

- [ ] **US-003**: 1M+ observations across 20+ margins calibrate in under 30 seconds,
      single-threaded, so a census microsimulation researcher can iterate on synthetic
      population models. `verbose = 1` prints the selected algorithm and routing reason.
      *Status: **Partial**.* Routing observability shipped. The performance targets have
      **never been verified against a live solver** and are internally contradictory:
      §1 says medium-scale 100K rows / 5 margins < 1 s, §11 says < 2 s for the same shape,
      and both were written against the removed `lbfgsb` whose measuring artefact
      (`test-lbfgsb.R`) is void. Investigation `docs/investigations/2026-04-23-kk1204-convergence.md`
      (commit `3effd3a`) found the composite gate "<30 s AND <1e-6" **structurally
      unachievable** on K=20 uniform-random input (M_cell/n = 1.0 → zero compression
      benefit; extrapolated ~1.76M iterations). → `leafblower-kk1.20.4` demands a REFRAME.

### Distribution

- [ ] **US-010**: leafblower is distributable on CRAN and PyPI so users install it with
      `install.packages()` and `pip install`.
      *Status: **Partial** — the largest genuinely-unfinished requirement.*
      Shipped: `configure` with C++17→C++14 fallback, no vendored dependency > 5 MB,
      `pyproject.toml` on scikit-build-core.
      Outstanding: (a) `.Rbuildignore` does not exclude the tracked repository-root strays
      `cell_table_92c4f45.{cpp,hpp}`, `ieppa_92c4f45.cpp`, `patch_raking.py`,
      `patch_wolfe.py`, `test_output.log`, `leafblower_0.1.0.tar.gz`, `REVIEW_FINDINGS.md`,
      `code-review-findings.md` — `R CMD build` ships all of them (`leafblower-l6h0`);
      (b) `cran-comments.md` **does not exist** (it is already listed in `.Rbuildignore`);
      (c) no PyPI wheel has ever been built (`leafblower-kk1.24.3`);
      (d) `DESCRIPTION` and `python/pyproject.toml` versions (both `0.1.0`) are synced by
      hand with no check.
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
- [ ] **KPI-04**: Large-scale — 1M rows, 20 margins, `max_weight = 3`, < 30 s. **Blocked on
      US-003 reframe.** Its named measuring artefact was the Phase-2 gate against `lbfgsb`.
- [ ] **KPI-05**: CRAN check — 0 errors, 0 warnings, via `R CMD check --as-cran`. **Open.**
- [ ] **KPI-06**: Python wheel — installs on Linux/macOS, Python 3.9–3.13, via a CI matrix.
      **Open** — no CI pipeline exists in the repository at all.
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
| US-008 | Phase 5 (residual only) | Partial |
| US-009 | — (shipped pre-GSD) | Implemented |
| US-010 | Phase 5 | Partial |
| KPI-01 | — (shipped pre-GSD) | Implemented |
| KPI-02 | Phase 1 | Open |
| KPI-03 | — (shipped pre-GSD) | Implemented |
| KPI-04 | Phase 3 | Open |
| KPI-05 | Phase 5 | Open |
| KPI-06 | Phase 5 | Open |
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
