# Validation Target-Sum Dedup Implementation Plan

**Goal:** Eliminate the duplicate target-sum-to-1 check that runs once in `validate_calibrate_inputs` and again in `calib_validate_preentry`.
**Architecture:** `src/validation.hpp:74-86` already validates `|sum(targets[k]) - 1| <= 1e-6` for every margin during entry validation. `src/calib_validate.cpp:78-88` repeats the same check inside `calib_validate_preentry` (called by cell-table solvers). Both run on every R-bridge call.
**Tech Stack:** C++17, R CMD INSTALL --preclean ., testthat.

**Mechanism:** delete the duplicate target-sum block in `calib_validate.cpp`; rely on the upstream check in `validation.hpp`.
**Forbidden:** changing the entry-validation error message or error code (would break tests pinned to exact strings); altering the cell-bound checks at `calib_validate.cpp:32-75`; introducing a new shared helper file (over-abstraction for a 10-line block).
**Audit:** verify `validate_calibrate_inputs` is called *before* `calib_validate_preentry` on every code path that reaches the cell-table solvers.

---

## Task T1: Audit call ordering

Steps:

1. Read `src/r_bridge.cpp` to map: `harvest_cpp` → `validate_calibrate_inputs` → solver dispatch → `calib_validate_preentry`. Confirm validation runs before pre-entry on every R bridge entry point.
2. Read `src/c_api.cpp` for the C ABI path — same ordering must hold for `rk_calibrate` callers.
3. Grep all callers of `calib_validate_preentry`:
   ```bash
   grep -n "calib_validate_preentry" src/*.cpp
   ```
4. For each caller, trace upward to confirm `validate_calibrate_inputs` ran first. If any caller reaches `calib_validate_preentry` *without* prior validation, this plan is **blocked** — log finding to ticket and halt.

Confidence: 70 — call ordering is conventional but must be verified, not assumed.

---

## Task T2: Delete the duplicate block

Edit `src/calib_validate.cpp` lines 77-88 (the `// 5. Target sum validation` block):

```cpp
// REMOVE these lines:
    // 5. Target sum validation
    for (int k = 0; k < st.K; k++) {
        double s = 0.0;
        for (int j = 0; j < st.cat_counts[k]; j++) s += st.targets[k][j];
        if (std::fabs(s - 1.0) > 1e-6) {
            char msg[256];
            std::snprintf(msg, sizeof(msg),
                "margin %d targets sum to %.8f (expected 1.0±1e-6); "
                "normalize targets before calling", k, s);
            return fail(RK_ERR_BADARG, msg);
        }
    }
```

Replace with a one-line comment:

```cpp
// Target-sum validation lives in validate_calibrate_inputs (validation.hpp).
// It always runs before this function on both R-bridge and C-ABI paths.
```

Update the contract block in `src/calib_validate.hpp:12-23` — drop list item 5 and renumber. Final list:

```
* 1. n_cats_total <= kNCatsTotalMax        → RK_ERR_BADARG
* 2. L_c <= U_c for all cells              → RK_ERR_BADARG with cell index
* 3. X_init[c]==0 && L_c>0 for any cell    → RK_ERR_INFEAS
* 4. sum(L_c) <= n <= sum(U_c)             → RK_ERR_INFEAS
```

Confidence: 90 — exact source quoted from `calib_validate.cpp:77-88` and `validation.hpp:74-86`.

---

## Task T3: Rebuild + test

```bash
R CMD INSTALL --preclean .
Rscript -e 'devtools::test(filter = "calib")'
Rscript -e 'devtools::test(filter = "validation")'
```

Pass criteria:
- All `test-calib-*.R` and existing validation tests pass.
- The error message from `validate_calibrate_inputs` (line 85: `"targets[k] does not sum to 1 (within 1e-6)"`) — not the deleted message — is what fires when a test feeds bad targets. Update any test pinned to the deleted message string.

---

## Task T4: Verify removed-message has no callers in tests

```bash
grep -rn "margin.*targets sum to" tests/ 2>/dev/null
```

If hits exist: update those tests to expect the upstream `"targets[k] does not sum to 1"` message.

---

## Task T5: Commit

`refactor(validation): drop duplicate target-sum check in calib_validate_preentry`

Body cites the duplication: validation.hpp:74-86 + calib_validate.cpp:77-88, single source of truth at entry validation (leafblower-c8dg).
