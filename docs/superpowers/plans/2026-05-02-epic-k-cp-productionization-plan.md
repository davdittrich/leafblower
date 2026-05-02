# Epic-K stepstone-CP Productionization — Implementation Plan

**Date:** 2026-05-02 (initial); 2026-05-03 (plan rev 1)
**Beads epic:** `leafblower-pcs9`
**Beads tasks:** `leafblower-pcs9.1` through `.7`
**Spec:** `docs/superpowers/specs/2026-05-02-epic-k-cp-productionization-design.md` rev 4 (design-review-gate iter 3 APPROVED 4/5 with override; CTO 6 textual blockers folded in rev 4)
**Predecessor:** Epic-J (`leafblower-y2ls`) FAIL at master `484e1e2`; spike CP at master `4e89769`. Investigation: `docs/investigations/2026-05-02-ylsy-cp-ipm-spike-result.md` Sec 6 (side-finding).

## Mechanism

**Target pattern:** opt-in `harvest(method="cp", ...)` solver wrapping Chambolle-Pock primal-dual (PDHG) — Algorithm 1 (vanilla, O(1/k)) and Algorithm 2 (accelerated, O(1/k²)) per Chambolle-Pock 2011 (J. Math. Imaging Vision 40:120–145). Cell-compressed by default; obs-level fallback at `M_cell/n > 0.9` or `bounds_mode="unit"`.

**Library dependencies:** RcppEigen (header-only, vendored at Epic-J WU-1), LAPACK-via-Makevars (existing). `RhpcBLASctl` already installed (Epic-J).

**ABI:** single-line `RK_ALG_CP = 12` enum addition to `src/leafblower.h`. NO `rk_params_t` field changes; `cp_safety_factor=1.05` hardcoded; `accelerate` reuses existing `st.accelerate` field. `EXPECTED_RK_PARAMS_BYTES` unchanged.

**Audit strategy:**
- Spec-compliance reviewer (Opus, per WU): each WU's DoD checklist verified against spec rev 4 sections.
- Code-quality reviewer (Opus, per WU): memory safety (std::vector / Eigen / no raw new/delete), R/C boundary (PROTECT/UNPROTECT mirror newton_kl pattern), numerical NaN propagation (status_code=2 detection).
- Mechanical CI gate (`tools/check_research_isolation.R`): forbidden-symbol list updated K-1 step 1 to remove cp_solve_R + cp_calibrate (now legitimately in src/leafblower.so) while keeping ipm_solve_R + ipm_calibrate (Epic-J FAIL artefacts).
- K-7 plan-review-gate convention (3 adversarial reviewers Feasibility/Completeness/Scope & Alignment) on cumulative Epic-K diff before close.

## Forbidden

- **No edits to** `src/Makevars.in`, `src/types.hpp` (CalibState/CalibResult), `src/cell_table.{hpp,cpp}` beyond what spec rev 4 explicitly permits.
- **No AUTO routing change** (Epic-K.2 deferred — fixture-class sweep required).
- **No Mehrotra predictor-corrector PDHG**, **no warm-start**, **no SIMD**, **no Python bindings**, **no vignette docs** in this epic.
- **No raw `new`/`delete`** in any C++ file — `std::vector` / Eigen only.
- **No Lambert-W in CP prox** (direct Newton with overflow-asymptote guard per Epic-J spike).
- **No `--no-verify` on commits** — pre-commit hook runs the isolation gate; fix root cause if hook rejects.
- **No bundling of WUs** — one bd ticket per atomic deliverable; K-1's 10 sub-steps are atomic-as-revertible-commit (sub-steps coupled by pre-commit isolation gate R17).
- **No staging** of `.beads/issues.jsonl`, `src/*.o`, `src/*.so`, `.tldr/`, `.wolf/` in any commit.
- **No closing of `leafblower-ylsy`** by K-7 directly — orchestrator stewards the long-standing kk1204 ticket separately.
- **No skipping plan-review-gate** (this skill's exit point).

## Audit (Spy/Mock Strategy)

- **Sanity recovery test** is K-2 inline T1 (K=3 fixture, status=0, max_err < 1e-6 = machine precision). Catches K-1 refactor errors.
- **Form discrimination** is K-2 inline T5 (cp ≢ greg by >1% rel diff; all cp weights > 0). Catches accidental KL-vs-chi2 mixing.
- **Algorithm 1 production parity** is K-2 inline T6 (stepstone Alg 1 max_err < 1e-4). Catches obs-level CP regression.
- **Algorithm 2 fallback** is K-3 inline T7 (max_weight=Inf → fell_back_to_pdhg=TRUE + algorithm_used="pdhg"). Catches γ=0 dispatch bug.
- **Cell-mode equivalence** is K-4 inline T3 (cell ≡ obs to 1e-10 weight diff on bounds_mode="cell"). Catches cell aggregation errors.
- **Headline win** is K-4 inline T2 (stepstone CP max_err ≤ 0.7× ieppa+sraa baseline). Locks Epic-J spike 2.2× win with safety margin.
- **Codified test suite** is K-5 (`tests/testthat/test-cp.R` T1-T7, runs under `devtools::test()`). Persistent CI regression detection.
- **Pre-commit isolation gate** is K-1 step 1 (`tools/check_research_isolation.R` updated forbidden list). Mechanical guard against accidental ipm productionization.
- **Adversarial review** is K-7 (3 fresh Opus reviewers on cumulative Epic-K diff). Catches cross-WU integration issues per-WU reviews missed.

## Work Units (one bd ticket each — see Epic `leafblower-pcs9`)

| WU | Bead | Title | Hard deps | Model | Wall budget |
|---|---|---|---|---|---|
| K-1 | `leafblower-pcs9.1` | Move + adapt cp_calib to src/ + isolation gate update + R wiring | — | Gemini | ~2.5h |
| K-2 | `leafblower-pcs9.2` | Validate Algorithm 1 obs-level production parity | K-1 | Gemini | ~2h |
| K-3 | `leafblower-pcs9.3` | Implement Algorithm 2 accelerated PDHG + fallbacks | K-2 | Opus | ~3h |
| K-4 | `leafblower-pcs9.4` | Cell-compressed CP with bounds_mode dispatch | K-2 | Opus | ~4h |
| K-5 | `leafblower-pcs9.5` | Test suite tests/testthat/test-cp.R T1-T7 + full regression | K-3, K-4 | Haiku | ~1h |
| K-6 | `leafblower-pcs9.6` | NEWS.md additive + harvest.Rd regen + harvest.R docstring | K-5 | Haiku | ~30min |
| K-7 | `leafblower-pcs9.7` | Code-review-gate (3 adversarial reviewers) + cleanup commit | K-6 | Opus | ~1h |

**Dependency graph (acyclic):**
```
K-1 ──► K-2 ──► K-3 ──► K-5 ──► K-6 ──► K-7
            │           ▲
            └─► K-4 ────┘
```

K-3 + K-4 both depend on K-2; K-5 depends on both. K-3 and K-4 may run in parallel (different model assignments) since both depend only on K-2 and edit different src/cp_calib.cpp regions.

**Total wall:** ~13–14h sequential; ~10–11h with K-3 ∥ K-4 parallel.

**Subagent routing**: Gemini for mechanical port + bench (K-1, K-2). Opus for algorithm correctness + cell compression + final review (K-3, K-4, K-7). Haiku for tests + docs (K-5, K-6). Per memory rule: avoid Sonnet; use Gemini for Tier-2 delegated work, Opus for Tier-3 design/review, Haiku for Tier-1 mechanical.

## Decision Rule (verbatim from spec rev 4 Sec 1)

Epic-K closes PASS iff:
- All 7 tests in `tests/testthat/test-cp.R` PASS via `devtools::test()`.
- `R CMD INSTALL --preclean .` clean.
- `harvest.Rd` + `NEWS.md` updated.
- AUTO routing untouched (existing `test-algo-selection.R` regression PASS).
- K-7 adversarial review (3 reviewers) APPROVED.

Any FAIL → halt, fix, re-test (no `--no-verify` bypass).

## Reversibility / FAIL artefact policy

Each WU = one revertible commit. K-1..K-7 are 7 sequential commits; revert any single SHA on bug. `research/cp_calib.{hpp,cpp}` stays as fossil with header comment "MOVED TO src/" (per K-1 sub-step 2). On Epic-K FAIL: revert all K-1..K-7 commits in reverse order; spec rev 4 stays as design artefact for future re-attempt.

## Risks & Mitigations (verbatim from spec rev 4 Sec 6 R1-R17; reproduced for plan completeness)

| # | Risk | Mitigation |
|---|---|---|
| R1 | Cell-compression bug — cell-mode bounds invariant violated | T3 cell ≡ obs to 1e-10 (K-4) |
| R2 | Algorithm 2 step-size adaptation unstable on ill-conditioned A | Fallback to Alg 1 if θ_k < 1e-15 (K-3); reset τ/σ to fixed values (not freeze) |
| R3 | u_max=Inf violates γ > 0 precondition | Auto-fallback to PDHG; verbose log; T7 verifies (K-3) |
| R4 | r_bridge dispatch arm wires wrong field set | Mirror newton_kl arm pattern at r_bridge.cpp:608-609 (K-1 step 7); spec reviewer audits |
| R5 | harvest.R match.arg whitelist drift | T1-T7 use method="cp"; whitelist drift → all tests fail |
| R6 | NEWS.md bullet under wrong section | Place under "## New features" (additive); reviewer audits |
| R7 | RK_ALG_CP enum collision | Verify 12 free at K-1 step 6 (currently RK_ALG_NEWTON_KL=11 last) |
| R8 | A_cell construction differs from obs-level A | T3 direct comparison + assertion sum(A_cell.x) == M_cell × K |
| R9 | Power-iter on cell A vs obs A different ‖A‖ | Self-consistent within each path; documented in result A_norm_estimate |
| R10 | accelerate=TRUE slower than Alg 1 (over-acceleration) | T6 verifies Alg 1 baseline; default TRUE because spike rate β=-1.05 R²=0.96 strong |
| R11 | T2 stepstone fixture not in CI | `skip_if_not_installed("arrow")` graceful skip pattern (K-5) |
| R12 | Algorithm bug not caught by spike sanity | Spec compliance reviewer audits cp_calib.cpp on stepstone trace |
| R13 | else-if dispatch order shifts existing fall-through | Place cp arm BEFORE catch-all else (K-1 step 7) |
| R14 | bounds_mode="unit" cell expansion incorrect | bounds_mode="unit" forces obs-level (K-4) |
| R15 | A_cell sparse vs dense storage ambiguous | Eigen::SparseMatrix<double> default; dense fallback Epic-K.2 future |
| R16 | OpenMP interaction with cell_table reuse | CP single-threaded; cell_table inherits ieppa parallelism; comment header K-1 step 4 |
| R17 | tools/check_research_isolation.R blocks pre-commit because cp_* now in src/leafblower.so | K-1 step 1 (FIRST step): UPDATE forbidden list to REMOVE cp_solve_R + cp_calibrate; KEEP ipm_solve_R + ipm_calibrate |

**Discontinuation triggers** (per spec rev 4 Sec 6):
- R5 fires (whitelist drift) → halt; verify match.arg before any further commit.
- R7 fires (enum collision) → halt; pick next free slot; audit downstream switch tables.
- R12 fires (algorithm bug) → halt; revert; re-spike before continuing.

## Success Criteria (Epic close)

- [ ] All 7 WU tickets closed via bd.
- [ ] `tests/testthat/test-cp.R` T1-T7 PASS via `devtools::test()`.
- [ ] `R CMD INSTALL --preclean .` clean.
- [ ] `R CMD check` no NOTES related to cp.
- [ ] AUTO routing regression tests PASS (`tests/testthat/test-algo-selection.R`).
- [ ] K-7 adversarial reviewers (3) APPROVED.
- [ ] `bd close leafblower-pcs9` with verdict comment.
- [ ] `bd update leafblower-ylsy --notes` records stepstone-CP productionization pointer (K-7 step 8).

## Post-merge action

Plan-review-gate iter 1 must PASS (3 adversarial reviewers: Feasibility, Completeness, Scope & Alignment per `metaswarm:plan-review-gate` convention) before any WU starts.
