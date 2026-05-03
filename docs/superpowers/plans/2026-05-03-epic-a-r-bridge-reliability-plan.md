# Epic A: r_bridge reliability + maintainability — Implementation Plan (rev 2)

**Date:** 2026-05-03
**Beads epic:** `leafblower-rnu2`
**Beads tasks:** `rnu2.3` (A1), `rnu2.1` (A2), `rnu2.2` (A3)
**Motivation:** Code review 2026-05-03 — REQUIRED severity issues in src/r_bridge.cpp.
**Revision:** rev 2 closes 9 plan-review-gate iter 1 blockers (Feasibility: 5; Completeness: 4).

## Mechanism

Three surgical edits. No algorithm changes. A3 removes one `.Call` arg (not adds).

### A1 — constexpr N_RESULT_FIELDS

Current: `Rf_allocVector(VECSXP, 37)` at r_bridge.cpp:734 and `Rf_allocVector(STRSXP, 37)` at :735 — two bare literals.

Fix: insert `constexpr int N_RESULT_FIELDS = 37;` (or verified count from `grep -c "SET_VECTOR_ELT(res_list" src/r_bridge.cpp`) immediately before the `PROTECT(res_list = ...)` line. Replace both literals with `N_RESULT_FIELDS`. N_RESULT_FIELDS stays 37 — A2 will overwrite slot 3 in place (see A2 mechanism), NOT add a new slot.

**Audit (rev 2 — literal match):** `grep -n " 37)" src/r_bridge.cpp | grep -E "VECSXP|STRSXP"` → empty.

### A2 — eliminate alg_names two-language duplication

**Root state (verified):** r_bridge.cpp:749 currently emits `Rf_ScalarInteger(res_alg_used)` into slot 3 (`algorithm_used`). R/harvest.R:519 reads it as integer via `alg_names[calib_result$algorithm_used + 1L]`.

**Mechanism:**
1. In src/r_bridge.cpp: replace `Rf_ScalarInteger(res_alg_used)` at slot 3 with `Rf_mkString(alg_name_cstr)` — a character string. `alg_name_cstr` is derived from a `static const char* kAlgNames[]` table indexed by `rk_algorithm_t` int. This table is the **single source of truth**. It must cover ALL enum values including `RK_ALG_IEPPA_SOFT=8` (currently missing from the ternary chain). Include `""` for AUTO=0, GRAKE=7(deprecated), any undefined slots.
2. Pre-consumer scan: `grep -rn "algorithm_used" R/ tests/` before A2 edit — verify no other consumer reads slot 3 as integer.
3. In R/harvest.R:519: delete `alg_names <- c(...)` + `alg_used <- alg_names[calib_result$algorithm_used + 1L]`; replace with `alg_used <- calib_result$algorithm_used` (C++ now emits string).

N_RESULT_FIELDS stays 37 (overwrite slot 3 type, no new slot). Downstream: `attr(r, "algorithm")` still "ieppa", "newton_kl", etc.

**Audit:** `devtools::test()` — all algorithm-name assertions PASS; `grep -n "alg_names" R/harvest.R` → empty; ieppa_soft name correct.

### A3 — remove jacobi_sweep from public API + extract CalibState helper

**Root state (verified):** `jacobi_sweep` is a documented public `harvest()` parameter:
- `@param jacobi_sweep` at harvest.R:93
- `jacobi_sweep = FALSE` in signature at harvest.R:242
- Passed as `.Call` arg #34

It is LIVE at the R API, dead only at C++ (r_bridge.cpp:135 `(void)jacobi_sweep_sexp`).

**Mechanism:**
1. Remove `jacobi_sweep = FALSE` from harvest.R:242 signature.
2. Remove `@param jacobi_sweep` roxygen documentation.
3. Remove arg from `.Call(...)` invocation in harvest.R:~404 (arg #34).
4. Remove `SEXP jacobi_sweep_sexp` + `(void)jacobi_sweep_sexp` from r_bridge.cpp C++ function signature + body.
5. Update r_bridge.cpp:63 registration arity from 35 → 34.
6. Update tests/testthat/test-safety.R: delete line `as.integer(0L), # 34: jacobi_sweep`, renumber subsequent comment index `# 35: newton_tsvd_ratio` → `# 34: newton_tsvd_ratio`.
7. Extract `static void init_calib_state(lbw::CalibState& st, const rk_params_t& p, ...)` as a **static function within r_bridge.cpp** (NOT a new header — no new public surface per Surgical Changes rule). Replaces 47-line field-copy block with one call.

**Audit:** `grep -rn "jacobi_sweep" src/ R/ tests/` → empty except removal-comment; registration arity verified at 34; `devtools::test()` PASS.

## Work Units

| WU | Bead | Title | Dep | Model | Wall |
|---|---|---|---|---|---|
| A1 | `leafblower-rnu2.3` | Replace magic `37` with N_RESULT_FIELDS | — | Haiku | ~20min |
| A2 | `leafblower-rnu2.1` | Eliminate alg_names duplication — change slot 3 to STRSXP | A1 | Gemini | ~1.5h |
| A3 | `leafblower-rnu2.2` | Remove jacobi_sweep from public API + extract static CalibState helper | A2 | Gemini | ~1.5h |

**Dep chain:** A1 → A2 → A3 (linear).

## Decision Rule

Epic A closes PASS iff:
- All 3 WU tickets closed.
- `devtools::test()` FAIL=2 only (pre-existing T2 basin floor).
- `R CMD INSTALL --preclean .` clean.
- All 3 grep audits empty/correct.

## Risks

| # | Risk | Mitigation |
|---|---|---|
| R1 | A2 breaks attr(r,"algorithm") for existing solvers | Pre-A2: `grep -rn "algorithm_used" R/ tests/`; full devtools::test() post-A2 |
| R2 | A3 arg removal breaks test-safety.R | A3 DoD requires explicit renumbering of comment indices in test-safety.R |
| R3 | A1 count wrong (not 37) | A1 step: `grep -c "SET_VECTOR_ELT(res_list" src/r_bridge.cpp` and confirm = 37 |
| R4 | A2 kAlgNames table incomplete — ieppa_soft or future enum missing | A2 DoD: match table against all rk_algorithm_t values; add static_assert or R unit test asserting all 12 enum values map to non-empty string |
| R5 | A3 callers pass jacobi_sweep=TRUE; behavior change | Deprecation note in NEWS.md: "parameter removed; passing it now causes 'unused argument' error" |
