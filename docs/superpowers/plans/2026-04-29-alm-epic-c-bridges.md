# ALM ieppa_soft — Epic C: R-Bridge + C-API Wiring

> **For agentic workers:** Use `superpowers:subagent-driven-development` to implement this plan.

**Mechanism:** Insert `capacity_penalty` at `.Call` slot 9; resolve `st.capacity_mu` from cell_table auto or user value; pack 4 ALM diagnostic fields into R result; mirror resolution in c_api.cpp for direct-C callers.
**Forbidden:** Changing any solver logic in `src/ieppa.cpp`; touching `src/types.hpp`, `src/leafblower.h`, `src/cell_table.*`, or `R/harvest.R` — those belong to Epic B and Epic D.
**Audit:** Test outcome must be FAIL==3 after T5 (capacity_mu=0 means ALM block inactive — behavior identical to pre-Epic-C; the 3 failing tests are pre-existing from Epic B). No new failures permitted.

**Branch:** `fix/correctness-performance-2026-04-28`

**Prerequisite state (verify before starting):**
```bash
# cell_table.hpp must expose capacity_mu_auto (Epic B landed)
grep -n 'capacity_mu_auto' src/cell_table.hpp
# Expected: at least 1 hit (field declaration)

# leafblower.h must have alm_* fields in rk_result_t (Epic B landed)
grep -n 'alm_capacity_mu_final\|alm_n_growth_events\|alm_max_dual_norm\|alm_sum_drift' src/leafblower.h
# Expected: 4 hits

# leafblower.h must have capacity_penalty in rk_params_t (Epic B landed)
grep -n 'capacity_penalty' src/leafblower.h
# Expected: at least 1 hit

# Pre-condition: check current arity
grep 'C_rk_calibrate.*31' src/r_bridge.cpp
# Expected: 1 hit (the R_CallMethodDef registration)
```

If any prerequisite fails, **HALT** — Epic B has not landed. Do not proceed.

---

## Task 5 (T5): Wire `capacity_penalty` into `src/r_bridge.cpp`

**Ticket:** leafblower-T5
**File:** `src/r_bridge.cpp`

### Ground Truth (read before editing)

Current function signature (lines 89–107) has 31 named SEXP parameters. The forward declaration at lines 24–27 already has 32 SEXP placeholders. The result list (`res_list`) is a 30-element VECSXP (line 548), indices 0–29 occupied. The `R_CallMethodDef` arity at line 43 reads `31`.

Current slot 9 (0-indexed: slot index 8) = `tol_abs_sexp`. Inserting `capacity_penalty_sexp` at slot 9 (1-indexed; after `start_weights_sexp`) means every subsequent arg shifts by one position in the function signature. The body assignment (`p.tol_abs = REAL(tol_abs_sexp)[0]` etc.) uses named variables — no index arithmetic needed; only the signature line changes.

### Step 5.1 — Read the file

```bash
# Confirm current line count and key anchors
grep -n 'start_weights_sexp\|tol_abs_sexp\|accelerate_sexp\|capacity_penalty\|R_CallMethodDef\|res_list.*30' src/r_bridge.cpp
```

Expected lines:
- `92: SEXP inner_max_iter_sexp, SEXP start_weights_sexp,`
- `93: SEXP tol_abs_sexp, SEXP bounds_mode_sexp,`
- `43: {"C_rk_calibrate",       (DL_FUNC)&C_rk_calibrate,       31},`
- `107: SEXP accelerate_sexp) {`
- `548: SEXP res_list  = PROTECT(Rf_allocVector(VECSXP,  30));`
- `549: SEXP res_names = PROTECT(Rf_allocVector(STRSXP,  30));`

### Step 5.2 — Insert `capacity_penalty_sexp` at slot 9 in function signature

**Edit target** — lines 92–93 of `src/r_bridge.cpp`:

Old:
```cpp
                    SEXP inner_max_iter_sexp, SEXP start_weights_sexp,
                    SEXP tol_abs_sexp, SEXP bounds_mode_sexp,
```

New:
```cpp
                    SEXP inner_max_iter_sexp, SEXP start_weights_sexp,
                    SEXP capacity_penalty_sexp,
                    SEXP tol_abs_sexp, SEXP bounds_mode_sexp,
```

This inserts the new parameter between `start_weights_sexp` (slot 8) and `tol_abs_sexp` (slot 10 after insertion), exactly matching the spec's slot table.

### Step 5.3 — Resolve `st.capacity_mu` from cell_table auto or user value

**Location:** After the `CalibState st;` block setup (around line 314–318) and before the method dispatch at line 369. The cell_table must already be built by this point.

**Problem:** The current r_bridge.cpp bypass (`// WU-E: call C++ solvers directly`) does NOT call `build_cell_table` — it builds `st` directly. The `ct.capacity_mu_auto` is therefore not available without adding a `build_cell_table` call. Per the spec plumbing chain: `r_bridge.cpp calls build_cell_table → CellTable.capacity_mu_auto = M_cell/n`.

Add a `build_cell_table` call and capacity_mu resolution block. Insert immediately after the `CalibState st;` setup block (after line 318, `st.accelerate = ...`), before the try-block:

```cpp
    // Resolve capacity_mu for ieppa_soft: build cell table to obtain auto value.
    // Done here (not inside solver) so r_bridge.cpp controls the resolution contract.
    {
        lbw::CellTable ct_tmp;
        int ct_rc = lbw::build_cell_table(n, K,
            const_cast<const int32_t**>(group_ids.data()),
            cat_counts.data(),
            weights.data(),
            ct_tmp);
        if (ct_rc != 0) {
            // build_cell_table failure is non-fatal here — capacity_mu falls back to 1.0.
            // Infeasibility (empty cell with positive target) is caught by validate_inputs above.
            st.capacity_mu = 1.0;
        } else {
            const double cp_val = Rf_isNull(capacity_penalty_sexp)
                ? -1.0
                : (LENGTH(capacity_penalty_sexp) == 1 ? REAL(capacity_penalty_sexp)[0] : -1.0);
            st.capacity_mu = (cp_val <= 0.0) ? ct_tmp.capacity_mu_auto : cp_val;
        }
    }
```

**Why here:** The existing `ieppa_soft` dispatch at line 481–501 sets `st.use_admm_capacity = true` and then calls `lbw::ieppa_solve(st)`. `ieppa_solve` reads `st.capacity_mu` — it must be set before the dispatch. The resolution block above runs unconditionally (cheap: `build_cell_table` is O(n·K)); `st.capacity_mu` is zero-initialized and remains 0 for non-ieppa_soft callers (harmless since ALM block is gated by `st.use_admm_capacity`).

### Step 5.4 — Pack 4 ALM diagnostic fields into result list

**Location:** After line 620 (`SET_VECTOR_ELT(res_list, 29, Rf_ScalarInteger(res_conv_minimized_metric));`), before `Rf_setAttrib(res_list, R_NamesSymbol, res_names);` at line 621.

**Precondition:** The `ieppa_soft` dispatch block (lines 481–501) already calls `pack_solver_result(res)` and copies ieppa-specific fields. The `IeppaResult` struct must expose `alm_capacity_mu_final`, `alm_n_growth_events`, `alm_max_dual_norm`, `alm_sum_drift` (added in Epic B to `src/ieppa.hpp`). Add corresponding `res_alm_*` locals alongside the existing result locals:

Add 4 local variables after line 350 (`std::vector<double> res_best_weights;`):

```cpp
    double res_alm_capacity_mu_final = 0.0;
    int    res_alm_n_growth_events   = 0;
    double res_alm_max_dual_norm     = 0.0;
    double res_alm_sum_drift         = 0.0;
```

Extend `pack_solver_result` lambda — or add a separate capture block in the `ieppa_soft` dispatch only. The cleaner approach (no lambda change) is to extract them explicitly in the `ieppa_soft` branch. In the `ieppa_soft` dispatch (lines 481–501), add after `res_sor_n_damped = res.sor_n_damped;`:

```cpp
            res_alm_capacity_mu_final = res.alm_capacity_mu_final;
            res_alm_n_growth_events   = res.alm_n_growth_events;
            res_alm_max_dual_norm     = res.alm_max_dual_norm;
            res_alm_sum_drift         = res.alm_sum_drift;
```

Extend the result list: change `VECSXP, 30` → `VECSXP, 34` and `STRSXP, 30` → `STRSXP, 34` at lines 548–549. Then after index 29:

```cpp
    /* Elements 30-33: ALM diagnostics (non-zero only for ieppa_soft) */
    SET_STRING_ELT(res_names, 30, Rf_mkChar("alm_capacity_mu_final"));
    SET_STRING_ELT(res_names, 31, Rf_mkChar("alm_n_growth_events"));
    SET_STRING_ELT(res_names, 32, Rf_mkChar("alm_max_dual_norm"));
    SET_STRING_ELT(res_names, 33, Rf_mkChar("alm_sum_drift"));
    SET_VECTOR_ELT(res_list,  30, Rf_ScalarReal(res_alm_capacity_mu_final));
    SET_VECTOR_ELT(res_list,  31, Rf_ScalarInteger(res_alm_n_growth_events));
    SET_VECTOR_ELT(res_list,  32, Rf_ScalarReal(res_alm_max_dual_norm));
    SET_VECTOR_ELT(res_list,  33, Rf_ScalarReal(res_alm_sum_drift));
```

### Step 5.5 — Update `R_CallMethodDef` arity 31 → 32

**Edit target** — line 43 of `src/r_bridge.cpp`:

Old:
```cpp
        {"C_rk_calibrate",       (DL_FUNC)&C_rk_calibrate,       31},
```

New:
```cpp
        {"C_rk_calibrate",       (DL_FUNC)&C_rk_calibrate,       32},
```

### Step 5.6 — Fix test-safety.R B7 .Call arity

> **Note:** This overlaps with Epic E T12 — Epic E T12 should be a no-op (or removed) if Epic C already applied this fix. Add a note to Epic E T12: "Epic C applies this; T12 verifies only."

In `tests/testthat/test-safety.R`, find the B7 `.Call("C_rk_calibrate", ...)` block with 31 args. Insert `NULL` at slot 9 (after `bad_sw=` at slot 8):

```r
  NULL,               # 9: capacity_penalty (NULL=auto; not under test here)
```

Update the comment for what was slot 9 (`tol_abs`) to read slot 10.

Verify:
```bash
grep -c 'C_rk_calibrate' tests/testthat/test-safety.R
# Expected: 1
```

### Step 5.7 — Compile gate

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected last line: `* DONE (leafblower)`. If compile fails with "no matching function for call to `build_cell_table`", verify `cell_table.hpp` is included (already at line 11) and `lbw::build_cell_table` signature matches. If it fails with "has no member named `alm_capacity_mu_final`", Epic B has not landed — **HALT**.

### Step 5.8 — Test gate

```bash
Rscript -e "devtools::test()" 2>&1 | tail -10
```

Expected: exactly **3 FAILs** (pre-existing from Epic B; capacity_mu=0 means ALM block inactive). Zero new failures. If there are more than 3 failures, **HALT** and diagnose.

### Step 5.9 — grep audit (spec acceptance criterion 10)

```bash
grep -rn 'C_rk_calibrate' R/ tests/ benchmarks/ python/ src/r_bridge.cpp
```

For each call site found: verify that `capacity_penalty_sexp` (or `NULL`) is now at slot 9. The primary call site is `R/harvest.R`. Any site with a positional `.Call(...)` must have `NULL` inserted at position 9 if `capacity_penalty` is not yet wired in that caller (Epic D). Do NOT modify `R/harvest.R` here — that is Epic D scope. Simply confirm the count of call sites and record it.

### Step 5.10 — Commit

```bash
git add src/r_bridge.cpp tests/testthat/test-safety.R
git commit -m "$(cat <<'EOF'
feat(r_bridge): wire capacity_penalty slot 9; resolve auto/manual capacity_mu; pack 4 alm diagnostics; arity 31→32

- Insert capacity_penalty_sexp at C_rk_calibrate slot 9 (after start_weights_sexp,
  before tol_abs_sexp); forward decl already had 32 placeholders.
- Call build_cell_table before dispatch to obtain ct.capacity_mu_auto (M_cell/n);
  resolve st.capacity_mu: NULL or ≤0 → auto, positive → user value.
- Extend result list 30→34 elements; pack alm_capacity_mu_final, alm_n_growth_events,
  alm_max_dual_norm, alm_sum_drift (non-zero only for ieppa_soft dispatch).
- Update R_CallMethodDef arity 31→32; static guard catches any missed .Call site at
  runtime (R errors "Incorrect number of arguments to .Call" on arity mismatch).
- Fix test-safety.R B7 .Call arity: insert NULL at slot 9 (capacity_penalty=auto)
  to match new 32-arg signature; Epic E T12 may verify only (no-op if already fixed).

Epic C bridge, §Plumbing chain. Spec: docs/superpowers/specs/2026-04-29-ieppa-alm-soft-capacity.md
EOF
)"
```

---

## Task 6 (T6): Wire `capacity_penalty` into `src/c_api.cpp`

**Ticket:** leafblower-T6
**File:** `src/c_api.cpp`

### Ground Truth (read before editing)

The `RK_ALG_IEPPA_SOFT` dispatch block at lines 288–320 of `src/c_api.cpp` sets `st.use_admm_capacity = true` and calls `lbw::ieppa_solve(st)`. It does NOT set `st.capacity_mu`. The CalibState `st` is zero-initialized so `st.capacity_mu = 0.0`, meaning the ALM Newton step denominator `(1 + rho)` = `(1 + 0 * X_tilde_c)` = `1`, which degrades the ALM update to `X_alm = X_tilde_c * (1 - lambda_cell[c])` — pure dual-driven, no penalty signal. This is a functional bug for direct-C callers.

The `build_cell_table` call in `c_api.cpp`: check whether it is called and where. It must be called before the dispatch block (needed for `ct.capacity_mu_auto`). If it is already present, simply add the resolution block. If absent, add it.

### Step 6.1 — Read and orient

```bash
grep -n 'build_cell_table\|capacity_mu\|capacity_penalty\|RK_ALG_IEPPA_SOFT\|st\.capacity_mu' src/c_api.cpp
```

Use output to confirm:
1. Where `build_cell_table` is called (obtain `ct` variable name and line).
2. Where `RK_ALG_IEPPA_SOFT` dispatch begins (confirmed: line 288).
3. Whether `st.capacity_mu` is already assigned (expected: NOT present — confirming the bug).

### Step 6.2 — Add capacity_mu resolution after `build_cell_table`

Find the `build_cell_table` call in `c_api.cpp` (produces a `CellTable` named `ct` or similar). Insert the resolution block immediately after:

```cpp
    /* capacity_penalty: direct C API callers bypass R-layer validation.
       Contract: p->capacity_penalty <= 0.0 selects auto (M_cell/n from cell_table);
       positive value is used directly without range-checking — callers must validate
       externally (out-of-range values produce degraded convergence, not UB). */
    st.capacity_mu = (p->capacity_penalty <= 0.0)
                     ? ct.capacity_mu_auto
                     : p->capacity_penalty;
```

**Where to place:** Immediately after the `lbw::build_cell_table(...)` call returns, before the `alg == RK_ALG_LBFGSB` dispatch block at line 241. This ensures `st.capacity_mu` is set for all algorithms; non-ieppa_soft algorithms ignore it (harmless, since `st.use_admm_capacity` remains false).

**If `build_cell_table` is not called in `c_api.cpp`:** The function uses `cat_counts`, `group_ids`, etc. directly. In that case, add a local cell_table build:

```cpp
    // Build cell table for capacity_mu auto-resolution (also needed by ALM block).
    lbw::CellTable ct_for_capacity;
    {
        int ct_rc = lbw::build_cell_table(n, K, group_ids, cat_counts, weights, ct_for_capacity);
        // Non-zero rc means degenerate input; validate_inputs above should have caught it.
        // Fall back to safe default rather than error (direct callers own their validation).
        st.capacity_mu = (ct_rc != 0)
                         ? 1.0
                         : (p->capacity_penalty <= 0.0)
                             ? ct_for_capacity.capacity_mu_auto
                             : p->capacity_penalty;
    }
```

**Contract:** `p->capacity_penalty <= 0.0` → auto (`capacity_mu_auto`). Positive → used directly, no clamping. Do NOT use `std::max(p->capacity_penalty, 1.0)` — that silently corrupts user values in (0, 1). The build failure fallback is `1.0` (safe sentinel), not the user value (which may be ≤0 meaning auto).

Adapt the variable name (`ct_for_capacity` vs. `ct`) to whatever `build_cell_table` produces in the existing code.

### Step 6.3 — Compile gate

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected: `* DONE (leafblower)`.

### Step 6.4 — Test gate

```bash
Rscript -e "devtools::test()" 2>&1 | tail -10
```

Expected: same **3 FAILs** as after T5 (count must not increase). Zero new failures.

### Step 6.5 — Commit

```bash
git add src/c_api.cpp
git commit -m "$(cat <<'EOF'
feat(c_api): route capacity_penalty to st.capacity_mu; document direct-caller contract

After build_cell_table, resolve st.capacity_mu from p->capacity_penalty:
  ≤ 0.0 → use ct.capacity_mu_auto (M_cell/n, same default as R-layer)
  > 0.0 → use directly (direct-C callers own range validation; no UB).

The RK_ALG_IEPPA_SOFT dispatch already sets st.use_admm_capacity=true; with
st.capacity_mu=0 the ALM Newton step was rho=0 (pure dual-driven, no penalty
enforcement). This fix restores the intended ALM behavior for C-API callers.

Epic C bridge, §Plumbing chain. Spec: docs/superpowers/specs/2026-04-29-ieppa-alm-soft-capacity.md
EOF
)"
```

---

## Acceptance Criteria

| # | Criterion | Verification |
|---|-----------|-------------|
| 1 | `C_rk_calibrate` arity = 32 | `grep 'C_rk_calibrate.*32' src/r_bridge.cpp` → 1 hit |
| 2 | `capacity_penalty_sexp` is slot 9 in function body | `grep -n 'capacity_penalty_sexp' src/r_bridge.cpp` shows it between `start_weights_sexp` and `tol_abs_sexp` |
| 3 | `st.capacity_mu` resolved before dispatch in r_bridge | `grep -n 'st.capacity_mu' src/r_bridge.cpp` → ≥1 hit |
| 4 | 4 ALM fields in result list at indices 30–33 | `grep -n 'alm_capacity_mu_final\|alm_n_growth_events\|alm_max_dual_norm\|alm_sum_drift' src/r_bridge.cpp` → 8 hits (4 name + 4 value) |
| 5 | result list size = 34 | `grep 'VECSXP.*34\|STRSXP.*34' src/r_bridge.cpp` → 2 hits |
| 6 | `st.capacity_mu` resolved in c_api before dispatch | `grep -n 'st.capacity_mu' src/c_api.cpp` → ≥1 hit before line 241 |
| 7 | Compile clean | `R CMD INSTALL --preclean . 2>&1 | tail -1` = `* DONE (leafblower)` |
| 8 | Test count ≤ pre-Epic-C failures | FAIL count = 3 (not 4+) |
| 9 | No silent .Call arity mismatch | `grep -rn 'C_rk_calibrate' R/ tests/ benchmarks/ python/` → all sites account for 32 args or await Epic D |

---

## Notes for Epic D (harvest.R wiring — not in scope here)

Epic D must insert `capacity_penalty_sexp` at slot 9 of the `.Call("C_rk_calibrate", ...)` call in `R/harvest.R`. Until Epic D lands, the `.Call` in `harvest.R` will pass 31 args and R will error at runtime for any harvest() call. This is expected — Epic C wires the C layer; Epic D wires the R layer. Do not attempt to fix harvest.R here.

The test suite's `Rscript -e "devtools::test()"` gate (Step 5.7 and 6.4) only fails with 3 pre-existing failures because those tests do NOT call `harvest()` (they call `.Call` directly or use fixtures from Epic A).
