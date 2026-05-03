# Epic A: r_bridge reliability + maintainability — Implementation Plan

**Date:** 2026-05-03
**Beads epic:** `leafblower-rnu2`
**Beads tasks:** `rnu2.3` (A1), `rnu2.1` (A2), `rnu2.2` (A3)
**Motivation:** Code review 2026-05-03 — REQUIRED severity issues in src/r_bridge.cpp.

## Mechanism

Three independent surgical edits to src/r_bridge.cpp + R/harvest.R. No algorithm changes. No ABI additions (one removal in A3).

**Mechanism per WU:**
- A1: `constexpr int N_RESULT_FIELDS = <count>;` inserted before result VECSXP/STRSXP allocations.
- A2: kAlgMap-based or static name table in C++ emits algorithm name string directly; R reads `calib_result$algorithm_used` as character — eliminates `alg_names` vector in harvest.R:518.
- A3: Remove `jacobi_sweep_sexp` from `.Call` signature (decreases arg count by 1); extract `CalibState init_calib_state(...)` helper replacing 47-line `st.X = p.X` copy block.

**Forbidden:** Algorithm changes; new .Call args (only removal for A3); new slots; any change to test-safety.R expected arg count without updating it simultaneously.

**Audit:**
- A1: `grep -n "VECSXP, [0-9]\+\|STRSXP, [0-9]\+" src/r_bridge.cpp` → empty after edit.
- A2: `devtools::test()` — all `attr(r,"algorithm")` assertions PASS for all methods.
- A3: `grep -rn "jacobi_sweep" src/ R/ tests/` → empty (or only removal-comment).

## Work Units

| WU | Bead | Title | Dep | Model | Wall |
|---|---|---|---|---|---|
| A1 | `leafblower-rnu2.3` | Replace magic `37` with N_RESULT_FIELDS | — | Haiku | ~20min |
| A2 | `leafblower-rnu2.1` | Eliminate alg_names two-language duplication | A1 | Gemini | ~1h |
| A3 | `leafblower-rnu2.2` | Remove dead jacobi_sweep param + extract CalibState helper | A2 | Gemini | ~1h |

**Dep chain:** A1 → A2 → A3 (linear).

## Decision Rule

Epic A closes PASS iff:
- All 3 WU tickets closed.
- `devtools::test()` FAIL=2 only (pre-existing T2 basin floor).
- `R CMD INSTALL --preclean .` clean.
- All 3 grep audits empty.

## Risks

| # | Risk | Mitigation |
|---|---|---|
| R1 | A2 breaks `attr(r,"algorithm")` for existing solvers | A2 DoD: run devtools::test() on all algorithm assertions; test-algo-selection.R must PASS |
| R2 | A3 .Call arg removal breaks test-safety.R (which counts args) | A3 step 4: update test-safety.R expected count simultaneously |
| R3 | A1 slot count miscounted | A1 step 2: count SET_VECTOR_ELT calls; confirm match with 37 |
