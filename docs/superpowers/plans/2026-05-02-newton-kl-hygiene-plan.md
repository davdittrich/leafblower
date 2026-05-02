# Plan: Epic-H — Newton-KL post-Epic-Dβ hygiene + AUTO routing safety

- **Date:** 2026-05-02
- **Author:** Dennis Alexis Valin Dittrich
- **Parent epic:** `leafblower-8aex`
- **Master HEAD at planning time:** `3c8eabd` (Epic-Dβ verdict commit)

## Mechanism

| Path | Target Pattern / Lib                                                                                   |
|------|--------------------------------------------------------------------------------------------------------|
| c    | NEWS.md `## Newton-KL calibration` consolidation; spec erratum (markdown only)                          |
| d    | C++ struct field cleanup in `src/newton_calib.hpp` + producer-site removal in `src/newton_calib.cpp`    |
| e    | New `harvest()` parameter routed via `src/c_api.cpp` → `rk_options_t` → `src/newton_calib.cpp`; two `Rprintf`/`snprintf` diagnostic strings (I1/I2) |
| g    | AUTO routing: target-skew ratio gate added in `src/c_api.cpp:182` and `src/r_bridge.cpp` ~line 425       |

## Forbidden

- DO NOT modify Newton-KL inner solver semantics (LM scaling, TSVD ratio default, Steihaug-CG step). Algorithm is frozen post-Epic-Dβ.
- DO NOT loosen ANY existing test gate (T1, T4, T5, T7, stepstone basin floor T2 documented at 2.61e-4 — unchanged).
- DO NOT introduce a new convergence criterion or change the default `convergence` API.
- DO NOT delete `lm_mu_final` until/unless WH-d explicitly chooses the "delete" branch with user approval; default branch is "surface via r_bridge SEXP-pack" matching existing `n_projected_dims` precedent (see r_bridge.cpp:752).
- DO NOT add Python-side parameter plumbing for `newton_tsvd_ratio` in WH-e — Python harvest API mirrors R via the same C ABI; if the C `rk_options_t` field is added, Python parameter exposure is a separate ticket (out-of-scope here).
- DO NOT bundle multiple paths into one bead ticket (one ticket per task — per CLAUDE.md global rule).
- DO NOT skip pre-commit hooks via `--no-verify`.

## Audit (Spy/Mock Strategy)

- **WH-c:** No audit (markdown only). Verify by re-reading the committed NEWS.md and spec on a fresh checkout.
- **WH-d:** Compile gate (`R CMD INSTALL --preclean .`) + grep audit confirming zero remaining references to deleted field across `src/`, `R/`, `tests/`, `inst/`. Spy: `grep -rn "n_homotopy_levels_used"` post-edit must return zero hits.
- **WH-e:** New testthat unit test that calls `harvest(method = "newton_kl", newton_tsvd_ratio = 1e-6)` and asserts the diagnostic field round-trip differs from default `1e-8` run on a contrived rank-deficient fixture. I1/I2 audit: capture `Rprintf` via `capture.output(...)` and assert string contains the expected snippet ("dsy_info=" for I1; "λ_max" for I2).
- **WH-g:** Spy via `LEAFBLOWER_NEWTON_TRACE`/equivalent — write a unit test with K=6, max_T/min_T = 7 (severe-skew) confirming AUTO selects `RK_ALG_IEPPA` not `RK_ALG_NEWTON_KL`. Mirror test with max_T/min_T = 2 (moderate) confirming AUTO still selects `RK_ALG_NEWTON_KL`.

## Context

Five sequential Newton-KL epics ran from 2026-04-XX through 2026-05-01:

- **Epic-A (LM damping, 5k08)** — verdict NEEDS_HOMOTOPY. Shipped scale-invariant Levenberg-Marquardt + best-iterate fallback. T2 stepstone K=9 max_err improved 0.988 → 2.8e-4 but stalled at the 1e-4 gate. Field `lm_mu_final` added to `NewtonCalibResult`; R-side surfacing deferred.
- **Epic-B (target homotopy, 91u7)** — verdict BLOCKED. Implemented level-stepping over target margins; failed to break the basin floor. Field `n_homotopy_levels_used` shipped but the homotopy code path was reverted; field is now dead.
- **Epic-C (IEPPA warm-start, usg8)** — verdict BLOCKED. Spec subsection (`docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md` §IEPPA warm-start, line 162) was committed at 32fcee6 but the implementation never landed; field `warmstart_max_err_at_handoff` is absent from master (verified via grep).
- **Epic-Dα (cheap-first sweep, y6tw)** — verdict ESCALATE. Triggered Epic-Dβ.
- **Epic-Dβ (TSVD+Steihaug-CG, wkmq)** — verdict PARTIAL (HEAD `3c8eabd`). Stepstone K=9 max_err: 2.79e-4 → 2.61e-4 (6.5% improvement, gate 1e-4 unmet). kk1204 K=20 severe-skew: status=1 → status=0, gap=6.24e-2 (above 1e-3 PARTIAL threshold but eliminates pathological drift). Field `n_projected_dims` shipped and surfaced via `src/r_bridge.cpp:752`.

T2 stepstone basin floor at 2.61e-4 is now declared intrinsic to the dual landscape; further closure deferred to a future Epic-E (continuation methods or alternative algorithm). The kk1204 K=20 severe-skew K≥5 fixed point at 6.24e-2 is the user-visible failure mode that motivates **path g** (AUTO routing safety) — Epic-A spec routes K≥5 + cell-ratio≥0.9 unconditionally to Newton-KL, but Newton-KL is now known to converge to a high-error fixed point on severe-skew K≥5 problems. The AUTO router must guard against this regime.

Hygiene paths c/d/e are independent of the algorithmic verdict and represent the smallest-valuable closure before stopping the algorithmic chase.

## Plan Steps

### WH-c — NEWS.md consolidation + IEPPA warm-start spec erratum

- **Files:** `NEWS.md`, `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md`
- **Action:**
  1. Collapse the four stacked Newton-KL bullets in `NEWS.md` `## Newton-KL calibration` (lines 4–34) into a single coherent narrative bullet covering: LM damping (Epic-A), TSVD pseudoinverse (Epic-Dβ), Steihaug-CG trust-region (Epic-Dβ), best-iterate fallback, basin-floor disclosure. Drop redundant prose on the same diagnostic field appearing twice.
  2. Append a NEWS bullet under `## Breaking changes` (line 36) documenting that AUTO routing for K≥5 severe-skew (max_T/min_T > 5) will redirect to IEPPA — this NEWS bullet is added in WH-c so the `## Newton-KL calibration` consolidation lands first; **WH-g implements** the routing change.
  3. Append an erratum block to spec subsection `### IEPPA warm-start` (line 162) noting Epic-C was attempted at 32fcee6 but verdict BLOCKED; section retained for future revival, not current behavior.
- **Why first:** All other tickets touch NEWS.md; consolidation must land before per-ticket bullets are appended to avoid merge churn.

### WH-d — Stale field cleanup

- **Files:** `src/newton_calib.hpp`, `src/newton_calib.cpp`
- **Action:**
  1. Delete `n_homotopy_levels_used` from `NewtonCalibResult` (newton_calib.hpp:14). Verified zero consumers via grep across `src/`, `R/`, `tests/`, `inst/`.
  2. Delete its single producer assignment in `newton_calib.cpp` (if any — re-grep to confirm).
  3. Surface `lm_mu_final` via `r_bridge.cpp` SEXP-pack mirroring the existing `n_projected_dims` pattern at `src/r_bridge.cpp:752` (string slot + ScalarReal slot, increment `res_names` length by one).
  4. `warmstart_max_err_at_handoff` is absent from master (verified) — no action required; document this in commit message.
  5. Compile gate via `R CMD INSTALL --preclean .` after each sub-step.
- **Why second:** Deleting the dead field before adding new diagnostics in WH-e avoids unnecessary churn on the same struct.

### WH-e — `newton_tsvd_ratio` user parameter + I1/I2 diagnostics

- **Files:** `src/newton_calib.cpp`, `src/newton_calib.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`, `R/harvest.R`, `man/harvest.Rd`, `tests/testthat/test-newton-kl-tsvd-ratio.R` (new)
- **Action:**
  1. Add `double newton_tsvd_ratio` to `rk_options_t` (default `1e-8`, matching current hardcoded default in `newton_calib.cpp:454` or thereabouts — re-read to find exact line).
  2. Plumb through `c_api.cpp` AUTO dispatch and explicit `RK_ALG_NEWTON_KL` dispatch.
  3. Plumb through `r_bridge.cpp` argument decode.
  4. Add `newton_tsvd_ratio = 1e-8` argument to `R/harvest.R` and document in `man/harvest.Rd`.
  5. Add I1: `Rprintf("[newton_kl] dsy_info=%d (LAPACK dsyevd)\n", info)` on non-zero LAPACK return inside the eigendecomposition block.
  6. Add I2: `Rprintf("[newton_kl] degenerate λ_max=%.3e ≤ 0; skipping TSVD\n", lambda_max)` on the degenerate-eigenvalue path.
  7. New testthat: `tests/testthat/test-newton-kl-tsvd-ratio.R` round-trips two distinct ratios and asserts result diverges (or `n_projected_dims` differs) on a contrived rank-deficient fixture.
- **Why third:** Builds on WH-d struct cleanup. Backwards-compatible additive API change (default 1e-8 preserves existing behavior).

### WH-g — AUTO routing target-skew gate

- **Files:** `src/c_api.cpp` (around line 182), `src/r_bridge.cpp` (around line 425), `tests/testthat/test-auto-routing-severe-skew.R` (new), `NEWS.md` (append to the `## Breaking changes` bullet seeded in WH-c with full implementation date)
- **Action:**
  1. Compute `target_skew = max(targets) / max(min(targets), 1e-12)` at the AUTO branch.
  2. New routing rule:
     ```
     if      (K >= 5 && M_cell_ratio >= 0.9 && target_skew >  5.0) -> RK_ALG_IEPPA   (accelerate=TRUE)
     else if (K >= 5 && M_cell_ratio >= 0.9 && target_skew <= 5.0) -> RK_ALG_NEWTON_KL
     else                                                          -> RK_ALG_RAKING
     ```
  3. Mirror identical logic at `src/r_bridge.cpp` AUTO branch (~line 425).
  4. Two new unit tests: severe-skew K=6 → asserts IEPPA path; moderate-skew K=6 → asserts Newton-KL path.
  5. Update `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md` `## AUTO Routing` (line 240) to reflect the new threshold.
- **Why last:** Behavior change for users; depends on NEWS bullet seed (WH-c) for breaking-change disclosure and on diagnostics (WH-e) for post-merge debugging if a user complains about routing.

## Cost

- **WH-c:** ~50 LOC markdown edits, ~10 min, zero compile cost.
- **WH-d:** ~20 LOC C++ deletions + ~12 LOC r_bridge.cpp insertions, ~30 min including compile gate.
- **WH-e:** ~80 LOC across 6 files including testthat, ~90 min.
- **WH-g:** ~30 LOC C++ + ~40 LOC testthat + ~15 LOC docs, ~60 min including full regression FAIL=0 verification.
- **Total Epic-H:** ~3.5 h on a clean checkout, plus full regression run (~25 min).

## Risks

| Risk                                                                                    | Likelihood | Mitigation                                                                                                                |
|-----------------------------------------------------------------------------------------|------------|---------------------------------------------------------------------------------------------------------------------------|
| **WH-g** breaks user pipelines that currently rely on Newton-KL for severe-skew K≥5     | Medium     | NEWS.md `## Breaking changes` bullet (seeded in WH-c, finalized in WH-g); semver bump on next release.                     |
| **WH-d** `lm_mu_final` surface adds a SEXP slot that downstream code does not expect    | Low        | Mirrors existing `n_projected_dims` pattern; R consumers either read the named slot or ignore it.                          |
| **WH-e** I1/I2 diagnostics fire spuriously and pollute test output                       | Low        | Wrap in `if (verbose)` flag (existing pattern) or test capture via `capture.output()`.                                      |
| **WH-e** `newton_tsvd_ratio` default change breaks existing test fixtures                | Low        | Default `1e-8` preserves current hardcoded value; tests must hash-identical pre/post.                                      |
| Cumulative NEWS.md merge conflicts if tickets parallelize                                | High       | Sequential dep chain WH-c → WH-d → WH-e → WH-g (justified below).                                                          |
| **WH-g** `target_skew` calculation underflows on `min(targets) ≤ 0`                      | Medium     | Floor `min(targets)` at `1e-12`; test fixture covers zero-target edge.                                                     |

## Out of Scope

- Re-attempting Epic-B target homotopy or Epic-C IEPPA warm-start (both BLOCKED, deferred to Epic-E or beyond).
- Closing the T2 stepstone basin floor (2.61e-4 → <1e-4 gate); intrinsic to dual landscape per Epic-Dβ verdict.
- Python-side `newton_tsvd_ratio` parameter exposure (separate ticket if desired post-Epic-H).
- Refactoring `r_bridge.cpp` SEXP-pack to a generic helper (touched only minimally per Surgical Changes directive).
- Renaming `n_homotopy_levels_used` (deletion is simpler than renaming a dead field).

## Dep Chain — Sequential Justified

Chain: **WH-c → WH-d → WH-e → WH-g**

**Why sequential, not parallel:**
1. **NEWS.md merge conflict surface:** WH-c rewrites the entire `## Newton-KL calibration` block (lines 4–34) AND seeds a `## Breaking changes` bullet. WH-d, WH-e, WH-g each append additional NEWS bullets. Parallel execution causes guaranteed merge conflicts on overlapping line ranges.
2. **Struct ordering:** WH-d deletes a field from `NewtonCalibResult`; WH-e adds new diagnostic plumbing through the same struct + `rk_options_t`. Sequential ordering avoids two independent edits to the same hpp/cpp file pair.
3. **Behavioral safety on WH-g:** The breaking-change NEWS bullet must land before the routing change to ensure release notes ship coherent (avoid "behavior changed but NEWS missing" race).

Confidence: 90 (sequential is provably safe; parallel has identifiable conflict surfaces).

## Pre-empted Plan-Review-Gate Concerns

- **WH-g is BREAKING for K≥5 severe-skew users** — addressed via explicit `## Breaking changes` NEWS bullet seeded in WH-c, finalized in WH-g, and a semver-minor bump recommendation in the WH-g commit message. Confidence: 95.
- **WH-e is API addition only** — additive new parameter with backwards-compatible default (`1e-8` matching current hardcoded value); no existing-caller regression possible. Confidence: 95.
- **WH-d field deletion** — re-verified zero consumers across `src/` `R/` `tests/` `inst/` via grep at planning time (output: only producer site in `newton_calib.hpp` and zero `_used` references in `r_bridge.cpp`); WH-d implementation must re-grep before deletion as a final guard. Confidence: 90.
- **WH-c spec erratum** — markdown-only, additive; no code path affected. Confidence: 99.

## Verification Gates (per ticket)

- WH-c: `git diff` shows consolidated NEWS section + erratum block; manual re-read on fresh checkout.
- WH-d: `R CMD INSTALL --preclean .` succeeds; `grep -rn "n_homotopy_levels_used" .` returns zero hits; `lm_mu_final` accessible from R via `attr(result, "result")$lm_mu_final`.
- WH-e: New testthat passes; full testthat FAIL=0; manual test triggers I1 message via contrived rank-deficient fixture.
- WH-g: New AUTO-routing testthat passes; full regression FAIL=0; spec doc reflects new gate; NEWS breaking-change bullet finalized.

## Confidence Summary

- Plan correctness: 90 (paths c/d/e are mechanical; path g rests on the assumption that IEPPA outperforms Newton-KL on severe-skew K≥5, which Epic-Dβ data supports for kk1204 K=20 specifically — generalization to K=5,6,7 not yet measured).
- Sequential chain necessity: 90.
- Cost estimate: 80.
