# Greenkhorn + Logit — Epic D: Bridge Wiring

**Mechanism:** Add dispatch blocks to `r_bridge.cpp` and `c_api.cpp` — mirror the raking pattern for Greenkhorn, the greg/sinkhorn pattern for Logit  
**Forbidden:** Modifying solver logic inside `greenkhorn.cpp` / `logit_calib.cpp` (those are Epics B/C); adding new fields to `rk_result_t` or `rk_params_t`; using `Rf_error` inside the `try {}` block for non-argument errors (use `throw std::runtime_error` instead)  
**Audit:** After each file edit: `R CMD INSTALL --preclean . 2>&1 | tail -5`; after both edits: `devtools::test()` FAIL count must remain 3; T1 greens if Epics B+C landed first

**Branch:** `fix/correctness-performance-2026-04-28`  
**Depends on:** Epic A (enum values 9, 10 in leafblower.h), Epics B+C (greenkhorn.hpp, greenkhorn.cpp, logit_calib.hpp, logit_calib.cpp)

---

## D1: `src/r_bridge.cpp` — add greenkhorn and logit dispatch blocks

**File:** `src/r_bridge.cpp`

**CRITICAL — TWO separate chains must both be extended (plan-review-gate finding):**

**Chain 1** (~lines 236-245): `p.algorithm = ...` if/else-if cascade — sets the algorithm for `validate_calibrate_inputs`. Currently falls through to `RK_ALG_IEPPA` for unknown methods — this silently routes greenkhorn/logit to IEPPA if missed.

**Chain 2** (~lines 252-260): `alg_for_validation` ternary — second dispatch for validation. Falls through to `RK_ALG_RAKING`.

Both chains must include GREENKHORN and LOGIT entries. See Step 1b below.

### Step 1 — Read and locate structure

Read lines 1–30 to see the `#include` block. Current includes:
```cpp
#include "logit.hpp"
#include "cell_table.hpp"
#include "types.hpp"
#include "ieppa.hpp"
#include "raking.hpp"
#include "sinkhorn.hpp"
#include "greg.hpp"
#include "chebyshev.hpp"
#include "lbfgsb_solver.hpp"
```

Read lines 465–510 to find the `greg` dispatch block (template for logit) and the `sinkhorn` block that immediately precedes it — this region is inside the `try {}` and after `pack_solver_result` is defined.

Confirmed pattern for `sinkhorn` and `greg` (lines ~471–492):
```cpp
    } else if (strcmp(method_str, "sinkhorn") == 0) {
        auto res = lbw::sinkhorn_solve(st);
        pack_solver_result(res);
        res_status     = res.status;
        res_iterations = res.iterations;
        res_max_error  = res.max_error;
        res_alg_used   = static_cast<int>(RK_ALG_SINKHORN);
        if (!res.best_weights.empty())
            res_best_weights = std::move(res.best_weights);
        else
            res_best_weights.assign(st.n, 0.0);
    } else if (strcmp(method_str, "greg") == 0) {
        auto res = lbw::greg_solve(st);
        pack_solver_result(res);
        res_status     = res.status;
        res_iterations = res.iterations;
        res_max_error  = res.max_error;
        res_alg_used   = static_cast<int>(RK_ALG_GREG);
        if (!res.best_weights.empty())
            res_best_weights = std::move(res.best_weights);
        else
            res_best_weights.assign(st.n, 0.0);
    } else {
```

The new dispatch blocks go before the `} else {` that opens the chebyshev/grake/ieppa_soft/ieppa catch-all.

Read lines 245–265 to locate the `alg_for_validation` chain and confirm it does NOT reference `greenkhorn` or `logit` yet (they must be added here).

### Step 2 — Add `#include` directives

The two new headers (`greenkhorn.hpp`, `logit_calib.hpp`) do not exist yet (Epic B/C create them). Add the includes now so the file is wired; the compiler will error until Epic B/C land — that is expected.

Replace the include block ending at `"lbfgsb_solver.hpp"`:
```cpp
#include "lbfgsb_solver.hpp"
```
with:
```cpp
#include "lbfgsb_solver.hpp"
#include "greenkhorn.hpp"
#include "logit_calib.hpp"
```

### Step 1b — Extend `p.algorithm` assignment chain (Chain 1, lines ~236-245)

**This step is mandatory and was added by plan-review-gate.** Without it, `method="greenkhorn"` silently routes to IEPPA.

Find the if/else-if cascade ending in `else p.algorithm = RK_ALG_IEPPA;`. Add two new entries before the final `else`:

```cpp
    else if (strcmp(method_str, "greenkhorn") == 0) p.algorithm = RK_ALG_GREENKHORN;
    else if (strcmp(method_str, "logit")      == 0) p.algorithm = RK_ALG_LOGIT;
    else                                             p.algorithm = RK_ALG_IEPPA;
```

Exact location (lines 236-245 in current file):
```
if      (strcmp(method_str, "ieppa")      == 0) p.algorithm = RK_ALG_IEPPA;
else if (strcmp(method_str, "ieppa_soft") == 0) p.algorithm = RK_ALG_IEPPA_SOFT;
...
else if (strcmp(method_str, "auto")       == 0) p.algorithm = RK_ALG_AUTO;
+else if (strcmp(method_str, "greenkhorn")== 0) p.algorithm = RK_ALG_GREENKHORN;   // ADD
+else if (strcmp(method_str, "logit")     == 0) p.algorithm = RK_ALG_LOGIT;        // ADD
else                                             p.algorithm = RK_ALG_IEPPA;
```

### Step 3 — Extend `alg_for_validation` chain (Chain 2)

Read lines 251–261 to find the `alg_for_validation` ternary chain:
```cpp
        rk_algorithm_t alg_for_validation =
            (strcmp(method_str, "ieppa")      == 0) ? RK_ALG_IEPPA :
            (strcmp(method_str, "ieppa_soft") == 0) ? RK_ALG_IEPPA_SOFT :
            (strcmp(method_str, "lbfgsb")     == 0) ? RK_ALG_LBFGSB :
            (strcmp(method_str, "auto")       == 0) ? RK_ALG_AUTO :
            (strcmp(method_str, "sinkhorn")   == 0) ? RK_ALG_SINKHORN :
            (strcmp(method_str, "greg")       == 0) ? RK_ALG_GREG :
            (strcmp(method_str, "chebyshev")  == 0) ? RK_ALG_CHEBYSHEV :
            (strcmp(method_str, "grake")      == 0) ? RK_ALG_GRAKE :
                                                       RK_ALG_RAKING;
```

Replace the last two lines of the chain:
```cpp
            (strcmp(method_str, "grake")      == 0) ? RK_ALG_GRAKE :
                                                       RK_ALG_RAKING;
```
with:
```cpp
            (strcmp(method_str, "grake")      == 0) ? RK_ALG_GRAKE :
            (strcmp(method_str, "greenkhorn") == 0) ? RK_ALG_GREENKHORN :
            (strcmp(method_str, "logit")      == 0) ? RK_ALG_LOGIT :
                                                       RK_ALG_RAKING;
```

### Step 4 — Add input validation for greenkhorn

The spec requires an `Rf_error` guard before `greenkhorn_solve` when `min_weight >= max_weight && max_weight > 0.0`. This prevents UB in `std::clamp(x, L, U)` when `L > U`. Insert it inside the `greenkhorn` dispatch block (see Step 5 below).

### Step 5 — Add greenkhorn dispatch block

Locate the exact text to find the insertion point:
```cpp
    } else if (strcmp(method_str, "greg") == 0) {
        auto res = lbw::greg_solve(st);
```

Insert the greenkhorn block BEFORE the greg block. The full replacement:

Replace:
```cpp
    } else if (strcmp(method_str, "greg") == 0) {
        auto res = lbw::greg_solve(st);
        pack_solver_result(res);
        res_status     = res.status;
        res_iterations = res.iterations;
        res_max_error  = res.max_error;
        res_alg_used   = static_cast<int>(RK_ALG_GREG);
        if (!res.best_weights.empty())
            res_best_weights = std::move(res.best_weights);
        else
            res_best_weights.assign(st.n, 0.0);
```
with:
```cpp
    } else if (strcmp(method_str, "greenkhorn") == 0) {
        // Input guard: clamp(x, L, U) is UB in C++17 when L > U ([alg.clamp]).
        if (p.min_weight >= p.max_weight && p.max_weight > 0.0)
            Rf_error("leafblower: min_weight (%.4g) >= max_weight (%.4g)",
                     p.min_weight, p.max_weight);
        auto res = lbw::greenkhorn_solve(st);
        pack_solver_result(res);
        res_status     = res.status;
        res_iterations = res.iterations;
        res_max_error  = res.max_error;
        res_alg_used   = (int)RK_ALG_GREENKHORN;
        if (!res.best_weights.empty())
            res_best_weights = std::move(res.best_weights);
        else
            res_best_weights.assign(st.n, 0.0);
    } else if (strcmp(method_str, "greg") == 0) {
        auto res = lbw::greg_solve(st);
        pack_solver_result(res);
        res_status     = res.status;
        res_iterations = res.iterations;
        res_max_error  = res.max_error;
        res_alg_used   = static_cast<int>(RK_ALG_GREG);
        if (!res.best_weights.empty())
            res_best_weights = std::move(res.best_weights);
        else
            res_best_weights.assign(st.n, 0.0);
```

### Step 6 — Add logit dispatch block

Locate the greg block closing (the `else` before `res_best_weights.assign(st.n, 0.0);` for greg). The logit block goes AFTER greg and BEFORE the `} else {` catch-all. 

Locate exact insertion boundary — the closing of the greg block:
```cpp
        else
            res_best_weights.assign(st.n, 0.0);
    } else {
        // Dispatch for chebyshev, grake (shared solver), ieppa_soft, and default ieppa.
```

Replace the greg block's closing and the `} else {` line:
```cpp
        else
            res_best_weights.assign(st.n, 0.0);
    } else {
        // Dispatch for chebyshev, grake (shared solver), ieppa_soft, and default ieppa.
```
with:
```cpp
        else
            res_best_weights.assign(st.n, 0.0);
    } else if (strcmp(method_str, "logit") == 0) {
        auto res = lbw::logit_calibrate(st);
        pack_solver_result(res);
        res_status     = res.status;
        res_iterations = res.iterations;
        res_max_error  = res.max_error;
        res_alg_used   = (int)RK_ALG_LOGIT;
        if (!res.best_weights.empty())
            res_best_weights = std::move(res.best_weights);
        else
            res_best_weights.assign(st.n, 0.0);
    } else {
        // Dispatch for chebyshev, grake (shared solver), ieppa_soft, and default ieppa.
```

Note: the logit dispatch uses the same `pack_solver_result` lambda used by all other solvers. `LogitCalibResult` has identical fields to `GregResult` (per spec §LogitCalibResult struct), so no custom packing is needed.

### Step 7 — Compile gate

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected states:
- **If Epics B+C are NOT landed yet:** compile fails with `greenkhorn.hpp: No such file or directory` — expected, commit anyway.  
- **If Epics B+C ARE landed:** must show `* DONE (leafblower)`.

If the compile fails for any other reason (syntax error, missing symbol), halt and diagnose before committing.

### Step 8 — Test gate (if compile succeeded)

```bash
Rscript -e "devtools::test()" 2>&1 | tail -3
```
Expected: FAIL count remains 3. T1 and T5 should GREEN if Epics B+C landed.

### Step 9 — Smoke test (if compile succeeded)

```bash
Rscript -e "
  library(leafblower)
  df  <- data.frame(sex=factor(c('M','F','M','F','M')))
  tgt <- list(sex=c(M=0.5, F=0.5))
  r   <- harvest(df, tgt, method='greenkhorn', max_iterations=50L)
  cat('status:', attr(r,'result')\$status, '\n')
  cat('alg_used:', attr(r,'result')\$algorithm_used, '\n')
"
```
Expected: `status: 0`, `alg_used: 9`.

### Step 10 — Commit

```bash
git add src/r_bridge.cpp
git commit -m "feat(r_bridge): add greenkhorn and logit dispatch blocks"
```

---

## D2: `src/c_api.cpp` — add GREENKHORN and LOGIT cases

**File:** `src/c_api.cpp`

### Step 1 — Read and locate structure

Read lines 1–17 to see the current `#include` block:
```cpp
#include "leafblower.h"
#include "validation.hpp"
#include "types.hpp"
#include "lbfgsb_solver.hpp"
#include "ieppa.hpp"
#include "raking.hpp"
#include "sinkhorn.hpp"
#include "greg.hpp"
#include "chebyshev.hpp"
#include "cell_table.hpp"
```

Read lines 270–290 to locate the dispatch chain — specifically the `RK_ALG_GREG` and `RK_ALG_CHEBYSHEV` cases:
```cpp
    } else if (alg == RK_ALG_SINKHORN) {
        auto sres = lbw::sinkhorn_solve(st);
        pack_solver_result(result, sres, alg);
        return sres.status;
    } else if (alg == RK_ALG_GREG) {
        auto gres = lbw::greg_solve(st);
        pack_solver_result(result, gres, alg);
        return gres.status;
    } else {
        if (alg == RK_ALG_CHEBYSHEV) { ...
```

### Step 2 — Add `#include` directives

Replace:
```cpp
#include "chebyshev.hpp"
#include "cell_table.hpp"
```
with:
```cpp
#include "chebyshev.hpp"
#include "cell_table.hpp"
#include "greenkhorn.hpp"
#include "logit_calib.hpp"
```

### Step 3 — Add RK_ALG_GREENKHORN case

The Greenkhorn case mirrors `RK_ALG_GREG` (both use `pack_solver_result` with three args and return immediately). Insert AFTER the GREG case and BEFORE the `} else {` that opens the chebyshev/grake/ieppa_soft/ieppa catch-all.

Replace:
```cpp
    } else if (alg == RK_ALG_GREG) {
        auto gres = lbw::greg_solve(st);
        pack_solver_result(result, gres, alg);
        return gres.status;
    } else {
```
with:
```cpp
    } else if (alg == RK_ALG_GREG) {
        auto gres = lbw::greg_solve(st);
        pack_solver_result(result, gres, alg);
        return gres.status;
    } else if (alg == RK_ALG_GREENKHORN) {
        auto res = lbw::greenkhorn_solve(st);
        pack_solver_result(result, res, RK_ALG_GREENKHORN);
        return res.status;
    } else if (alg == RK_ALG_LOGIT) {
        /* Direct C API callers bypass R-layer validation.
           Caller must validate p.capacity_penalty range externally. */
        st.capacity_mu = (p->capacity_penalty <= 0.0) ? ct_for_logit_mu : p->capacity_penalty;
        auto res = lbw::logit_calibrate(st);
        pack_solver_result(result, res, RK_ALG_LOGIT);
        return res.status;
    } else {
```

**Note on `capacity_mu` for logit:** The logit solver does not use `capacity_mu` — the logit link enforces bounds analytically with no penalty parameter. The `ct_for_logit_mu` fallback above is only needed if `CalibState.capacity_mu` is read by shared infrastructure (like `calib_dispatch.hpp`). If `logit_calibrate` does not read `capacity_mu`, omit the `st.capacity_mu` line entirely. Read `src/logit_calib.cpp` (once it exists in Epic C) to verify before finalizing this step.

**Simplified version** (if `logit_calibrate` does not read `capacity_mu`):
```cpp
    } else if (alg == RK_ALG_LOGIT) {
        /* Direct C API callers bypass R-layer validation.
           Caller must validate p.capacity_penalty range externally. */
        auto res = lbw::logit_calibrate(st);
        pack_solver_result(result, res, RK_ALG_LOGIT);
        return res.status;
    } else {
```

Prefer the simplified version unless Epic C's logit_calib.cpp reads `st.capacity_mu`.

### Step 4 — Compile gate

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected: `* DONE (leafblower)` (all Epics A, B, C, D1 must be landed for a clean compile).

If compile fails:
- Missing `greenkhorn.hpp` → Epic B not landed; expected until B lands.
- Missing `logit_calib.hpp` → Epic C not landed; expected until C lands.
- Syntax error → halt, read the error, fix before committing.

### Step 5 — Test gate

```bash
Rscript -e "devtools::test()" 2>&1 | tail -3
```
Expected: FAIL count remains 3 (T1–T8 are GREEN if all epics landed).

### Step 6 — Commit

```bash
git add src/c_api.cpp
git commit -m "feat(c_api): add RK_ALG_GREENKHORN and RK_ALG_LOGIT dispatch cases"
```

---

## Acceptance Criteria — Epic D

| Check | Expected |
|-------|----------|
| `grep greenkhorn_solve src/r_bridge.cpp` | 1 match in dispatch block |
| `grep logit_calibrate src/r_bridge.cpp` | 1 match in dispatch block |
| `grep greenkhorn_solve src/c_api.cpp` | 1 match in dispatch block |
| `grep logit_calibrate src/c_api.cpp` | 1 match in dispatch block |
| `grep greenkhorn.hpp src/r_bridge.cpp` | 1 match |
| `grep logit_calib.hpp src/r_bridge.cpp` | 1 match |
| `grep greenkhorn.hpp src/c_api.cpp` | 1 match |
| `grep logit_calib.hpp src/c_api.cpp` | 1 match |
| `R CMD INSTALL --preclean .` (all epics) | `* DONE (leafblower)` |
| `devtools::test()` FAIL count | remains 3 |
| `harvest(df, tgt, method="greenkhorn")$algorithm_used` | `"greenkhorn"` |
| `harvest(df, tgt, method="logit")$algorithm_used` | `"logit"` |

---

## Dependency Order

Epic D can be committed in two sub-steps:
1. `src/r_bridge.cpp` (D1) — can commit with compile failure if B/C not landed.
2. `src/c_api.cpp` (D2) — can commit with compile failure if B/C not landed.

The full clean compile (`* DONE`) requires A + B + C + D1 + D2 to all be present.
