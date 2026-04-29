# Greenkhorn + Logit — Epic A: Infrastructure

**Mechanism:** Enum extension + build wiring + R-layer routing stubs  
**Forbidden:** Touching solver .cpp/.hpp files (those are Epic B/C); any field additions to `rk_params_t` or `rk_result_t` (ABI-stable — spec confirms zero new fields)  
**Audit:** `R CMD INSTALL --preclean .` compile gate after every step; `devtools::test()` FAIL count must remain at 3 throughout

**Branch:** `fix/correctness-performance-2026-04-28`  
**Spec:** `docs/superpowers/specs/2026-04-29-greenkhorn-solver.md`

---

## A1: `src/leafblower.h` — extend `rk_algorithm_t` enum

**File:** `src/leafblower.h`  
**Verified baseline:** enum ends at `RK_ALG_IEPPA_SOFT = 8` (line 49). No static_assert on enum size. `EXPECTED_RK_PARAMS_BYTES` / `EXPECTED_RK_RESULT_BYTES` are unaffected (no new struct fields).

### Steps

**1. Read lines 39–51** to confirm `RK_ALG_IEPPA_SOFT = 8` is the last value and the closing `} rk_algorithm_t;` follows immediately.

Expected text:
```c
typedef enum {
    RK_ALG_AUTO   = 0,
    RK_ALG_IEPPA  = 1,
    RK_ALG_LBFGSB = 2,
    RK_ALG_RAKING    = 3,
    RK_ALG_SINKHORN  = 4,
    RK_ALG_CHEBYSHEV = 5,
    RK_ALG_GREG      = 6,
    RK_ALG_GRAKE      = 7,
    RK_ALG_IEPPA_SOFT = 8    /* ieppa + ADMM soft capacity enforcement */
} rk_algorithm_t;
```

**2. Edit** — replace the closing brace line only:

```c
    RK_ALG_IEPPA_SOFT = 8,   /* ieppa + ADMM soft capacity enforcement */
    RK_ALG_GREENKHORN = 9,   /* greedy coordinate-descent IPF (autumn::harvest style) */
    RK_ALG_LOGIT      = 10   /* Deville-Sarndal 1992 logit Newton calibration */
} rk_algorithm_t;
```

Exact edit — replace:
```c
    RK_ALG_IEPPA_SOFT = 8    /* ieppa + ADMM soft capacity enforcement */
} rk_algorithm_t;
```
with:
```c
    RK_ALG_IEPPA_SOFT = 8,   /* ieppa + ADMM soft capacity enforcement */
    RK_ALG_GREENKHORN = 9,   /* greedy coordinate-descent IPF (autumn::harvest style) */
    RK_ALG_LOGIT      = 10   /* Deville-Sarndal 1992 logit Newton calibration */
} rk_algorithm_t;
```

**3. Compile gate:**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```
Expected: `* DONE (leafblower)`. Failure here means a syntax error — halt and diagnose.

**4. Test gate:**
```bash
Rscript -e "devtools::test()" 2>&1 | tail -3
```
Expected: FAIL count unchanged at 3. A regression here means an existing test checks the enum size — halt and investigate.

**5. Commit:**
```bash
git add src/leafblower.h
git commit -m "feat(api): add RK_ALG_GREENKHORN=9, RK_ALG_LOGIT=10 to rk_algorithm_t enum"
```

---

## A2: `src/Makevars.in` — register new source files

**File:** `src/Makevars.in`  
**Verified baseline** (line 3):
```
PKG_SOURCES = c_api.cpp logit.cpp lbfgsb_solver.cpp ieppa.cpp cell_table.cpp r_bridge.cpp calib_validate.cpp sinkhorn.cpp calib_linalg.cpp greg.cpp chebyshev.cpp grake.cpp
```
Note: `raking.cpp` is NOT in this list — it is compiled via another mechanism. Do not add it.

### Steps

**1. Read** `src/Makevars.in` to confirm the exact `PKG_SOURCES` line shown above.

**2. Edit** — append the two new files to `PKG_SOURCES`:

Replace:
```
PKG_SOURCES = c_api.cpp logit.cpp lbfgsb_solver.cpp ieppa.cpp cell_table.cpp r_bridge.cpp calib_validate.cpp sinkhorn.cpp calib_linalg.cpp greg.cpp chebyshev.cpp grake.cpp
```
with:
```
PKG_SOURCES = c_api.cpp logit.cpp lbfgsb_solver.cpp ieppa.cpp cell_table.cpp r_bridge.cpp calib_validate.cpp sinkhorn.cpp calib_linalg.cpp greg.cpp chebyshev.cpp grake.cpp greenkhorn.cpp logit_calib.cpp
```

**3. Verify:**
```bash
grep PKG_SOURCES src/Makevars.in
```
Must show both `greenkhorn.cpp` and `logit_calib.cpp` on the same line.

**4. Compile gate:**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```
Expected outcome: compile will FAIL with "greenkhorn.cpp: No such file or directory" — this is correct and expected. Epic A only registers the files; Epics B and C create them. The install failure is not a regression; it confirms the build system now knows about the files.

If the failure is anything other than a missing-file error, halt and diagnose.

**5. Commit** (even with known compile failure — the plan explicitly permits this):
```bash
git add src/Makevars.in
git commit -m "build: add greenkhorn.cpp and logit_calib.cpp to PKG_SOURCES"
```

---

## A3: `R/harvest.R` — R-layer routing stubs for both methods

**File:** `R/harvest.R`  
**Three surgical edits. Read the file first; verify line numbers before editing.**

### Step 1 — Read target regions

Read lines 375–385 (status==2 stop), 432–440 (alg_names), and 488–496 (map_method) to confirm exact current text before editing.

Verified current text at each location:

**Location 1 — status==2 stop (line ~379):**
```r
  if (calib_result$status == 2L)
    stop("leafblower: infeasible problem — persistent empty cell with positive target (detected after 5 consecutive outer iterations).")
```

**Location 2 — alg_names (line ~437):**
```r
  # Enum: RK_ALG_AUTO=0, IEPPA=1, LBFGSB=2, RAKING=3, SINKHORN=4, CHEBYSHEV=5, GREG=6, GRAKE=7, IEPPA_SOFT=8
  alg_names <- c("", "ieppa", "lbfgsb", "raking", "sinkhorn", "chebyshev", "greg", "grake", "ieppa_soft")
```

**Location 3 — map_method() (line ~494):**
```r
  match.arg(method, c("auto", "ieppa", "ieppa_soft", "lbfgsb", "raking", "sinkhorn", "chebyshev", "greg", "grake"))
```

### Step 2 — Fix 1: `stop()` for `status == 2L`

The hardcoded message discards the C++ solver's actual diagnostic. Replace with `res$message` passthrough:

Replace:
```r
  if (calib_result$status == 2L)
    stop("leafblower: infeasible problem — persistent empty cell with positive target (detected after 5 consecutive outer iterations).")
```
with:
```r
  if (calib_result$status == 2L)
    stop("leafblower: ", if (nchar(calib_result$message) > 0) calib_result$message else "infeasible problem")
```

Rationale: the logit solver returns `RK_ERR_INFEAS=2` with `res$message = "logit: singular normal equations"`. A hardcoded override discards that diagnostic and confuses T5/T8 failures.

### Step 3 — Fix 2: extend `alg_names` to 11 elements

Replace:
```r
  # Enum: RK_ALG_AUTO=0, IEPPA=1, LBFGSB=2, RAKING=3, SINKHORN=4, CHEBYSHEV=5, GREG=6, GRAKE=7, IEPPA_SOFT=8
  alg_names <- c("", "ieppa", "lbfgsb", "raking", "sinkhorn", "chebyshev", "greg", "grake", "ieppa_soft")
```
with:
```r
  # Enum: RK_ALG_AUTO=0, IEPPA=1, LBFGSB=2, RAKING=3, SINKHORN=4, CHEBYSHEV=5, GREG=6, GRAKE=7, IEPPA_SOFT=8, GREENKHORN=9, LOGIT=10
  alg_names <- c("", "ieppa", "lbfgsb", "raking", "sinkhorn", "chebyshev", "greg", "grake", "ieppa_soft",
                  "greenkhorn", "logit")
```

Without this fix, `algorithm_used` returns `""` for any result with `algorithm_used == 9` or `== 10`, causing T1/T5 `algorithm_used` assertions to fail.

### Step 4 — Fix 3: `map_method()` — add both methods to `match.arg()`

Replace:
```r
  match.arg(method, c("auto", "ieppa", "ieppa_soft", "lbfgsb", "raking", "sinkhorn", "chebyshev", "greg", "grake"))
```
with:
```r
  match.arg(method, c("auto", "ieppa", "ieppa_soft", "lbfgsb", "raking", "sinkhorn", "chebyshev", "greg", "grake", "greenkhorn", "logit"))
```

Without this fix, `harvest()` rejects `method="greenkhorn"` and `method="logit"` at the R layer with a wrong-method error, before dispatch ever reaches C++.

### Step 5 — Test gate

```bash
Rscript -e "devtools::test()" 2>&1 | tail -3
```
Expected: FAIL count unchanged at 3. A new failure here means one of the three edits has a typo — halt and read the error.

### Step 6 — Verify routing reaches C++ (not rejected by R layer)

```bash
Rscript -e "
  library(leafblower)
  tryCatch(
    harvest(data.frame(x=factor(c('a','b'))), list(x=c(a=0.5, b=0.5)), method='greenkhorn'),
    error = function(e) cat('error:', conditionMessage(e), '\n')
  )
"
```
Expected: error message contains "internal solver" or "greenkhorn" (a C-level error from dispatch), NOT "should be one of" (which would indicate the R-layer match.arg is still rejecting it).

### Step 7 — Commit

```bash
git add R/harvest.R
git commit -m "feat(harvest): add greenkhorn+logit to map_method; extend alg_names to 11; fix stop() for status==2"
```

---

## Acceptance Criteria — Epic A

| Check | Expected |
|-------|----------|
| `grep RK_ALG_GREENKHORN src/leafblower.h` | `= 9` |
| `grep RK_ALG_LOGIT src/leafblower.h` | `= 10` |
| `grep PKG_SOURCES src/Makevars.in` | contains `greenkhorn.cpp logit_calib.cpp` |
| `grep '"greenkhorn"' R/harvest.R` | appears in both `match.arg()` and `alg_names` |
| `grep '"logit"' R/harvest.R` | appears in both `match.arg()` and `alg_names` |
| `devtools::test()` FAIL count | remains 3 |
| `harvest(..., method="greenkhorn")` error | C-level error, NOT `match.arg` error |
| `harvest(..., method="logit")` error | C-level error, NOT `match.arg` error |

Epic A does not require a passing compile (`R CMD INSTALL`) after A2 — the build fails with "file not found" until Epics B and C land. This is expected and documented.
