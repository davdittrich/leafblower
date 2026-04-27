# iEPPA Overflow Fix (T1.B) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Prevent linear-space overflow in iEPPA K=20 benchmark by replacing the broken T1.A geometric-mean renorm with T1.B (per-cell product tracking via `cell_lf`), bringing kk1204 from 31s → ~8.6s.

**Architecture:** Reuse existing `lf[]` vector (line 135) as a shadow log of `f_lin`. Add `cell_lf[c] = Σ_k lf[k][g_k(c)]` (per-cell log product). Track overflow via a high-water mark; correct before X_tilde exceeds `kLinearOverflowThreshold`. T2.A (log-path acceleration) falls out for free since `cell_lf` is already maintained.

**Tech Stack:** C++17, `src/ieppa.cpp` only (plus test update). Build: `R CMD INSTALL --preclean . 2>&1 | tail -3`. Test: `Rscript -e 'devtools::test()' 2>&1 | tail -4`.

---

## File Map

| File | Role |
|------|------|
| `src/ieppa.cpp` | All C++ changes |
| `tests/testthat/test-calibration-solvers.R` | Replace n=1000 T1.A test with n=100000 T1.B test |

---

## Task 0: Remove T1.A block (leafblower-mmre, step 1)

**Files:** Modify `src/ieppa.cpp:599-651`

- [ ] **Step 0.1: Find and delete T1.A block**

Search for the comment at line 599:
```
// T1.A: Renormalize f_lin per margin when the worst-case X_tilde product
```
Delete from this line through the closing `}` of the outer `if (!overflow_trip)` block at line 651. After deletion, line 597 (`}` closing the BCD else-branch) is immediately followed by line 653 (`if (overflow_trip && !linear_fallback_used)`).

- [ ] **Step 0.2: Compile gate**

```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -3
```
Expected: `* DONE (leafblower)`

- [ ] **Step 0.3: Commit**

```bash
git add src/ieppa.cpp
git commit -m "fix(ieppa): remove broken T1.A block (product trigger did not fire)

Geometric mean trigger (avg log f_lin per margin) stays near 3 while
threshold is 17.7 when one category grows and others shrink. X_cur /= c_k
direction was also wrong. T1.B replaces this with exact per-cell product
tracking via cell_lf."
```

---

## Task 1: TDD RED — Write failing T1.B test (new ticket)

**Files:** Modify `tests/testthat/test-calibration-solvers.R` (replace last test, ~lines 297-340)

- [ ] **Step 1.1: Replace the n=1000 T1.A test with the n=100000 T1.B test**

Find and replace the entire last `test_that(...)` block (starting at the `# ───` comment before `test_that("ieppa: no linear overflow trip on K=20 skewed targets (T1.A)"`) with:

```r
# ──────────────────────────────────────────────────────────────────────────────
# T1.B regression: no linear overflow on skewed multi-margin problems.
# Uses K=20, n=100000 (large enough to avoid S_lin collapse), max_weight=3.
# Math: kLinearOverflowTrip = (DBL_MAX/2)^(1/20) ≈ 2.1e15.
# Each BCD sweep: f_lin for 0.3-target category grows ~1.5x/sweep.
# Product across K=20 margins after N sweeps: (1.5)^(20*N).
# At N=5:  (1.5)^100 ≈ 6e17 >> trip ≈ 2.1e15 → overflow in ~4 sweeps.
# T1.B renorm fires at sqrt(trip) ≈ 4.6e7, preventing overflow entirely.
# Before T1.B: overflow trip fires, status != 0.
# After  T1.B: renorm fires (verbose=2 shows "T1.B renorm"), status == 0.
# ──────────────────────────────────────────────────────────────────────────────
test_that("ieppa: no linear overflow trip on K=20 skewed targets (T1.B)", {
  set.seed(42); K <- 20L; n <- 100000L
  cols <- paste0("v", seq_len(K))
  data <- as.data.frame(lapply(seq_len(K),
    function(k) factor(sample(5L, n, TRUE), levels = 1:5)))
  names(data) <- cols
  skewed <- c("1" = 0.3, "2" = 0.175, "3" = 0.175, "4" = 0.175, "5" = 0.175)
  target <- setNames(lapply(seq_len(K), function(.) skewed), cols)

  r <- leafblower::harvest(data, target, method = "ieppa",
                           min_weight = 0.2, max_weight = 3,
                           max_iterations = 50,
                           attach_weights = FALSE, verbose = 0)
  res <- attr(r, "result")

  # Before T1.B: overflow fires, solver falls back to log-space, fails to
  # converge in 50 iters (status != 0).
  # After T1.B: overflow prevented, linear path maintained, status == 0.
  expect_equal(res$status, 0L,
               label = "ieppa converges without linear overflow with T1.B")
})
```

- [ ] **Step 1.2: Verify RED — test must FAIL before T1.B**

```bash
cd /home/dd/Gemini/leafblower
Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | grep -E "FAIL|PASS|overflow|T1" | head -10
```
Expected: `FAIL` — overflow fires, status ≠ 0.

If test PASSES (status already 0), the problem is too easy. This should not happen since T1.A was removed in Task 0.

- [ ] **Step 1.3: Commit RED test**

```bash
git add tests/testthat/test-calibration-solvers.R
git commit -m "test(ieppa): RED — T1.B overflow prevention (K=20, n=100000)"
```

---

## Task 2: T1.B — cell_lf declaration + linear path maintenance (leafblower-mmre)

**Files:** Modify `src/ieppa.cpp`

### Step 2.1: Declare cell_lf and cell_lf_hwm

After line 139 (`std::vector<double> X_tilde(ct.M_cell);`), insert:

```cpp
    // T1.B: per-cell log-product shadow. cell_lf[c] = Σ_k lf[k][g_k(c)].
    // Reuses lf[] (already at line 135) as log(f_lin) in the linear path.
    // cell_lf is also used by T2.A in the log path (free side-effect).
    std::vector<double> cell_lf(ct.M_cell, 0.0);
    // High-water mark: max_c(log_X_init[c] + cell_lf[c]) ≈ max_c log(X_tilde[c]).
    // Monotone-nondecreasing between corrections (stale-high on negative deltas
    // is intentional: triggers early-but-never-missed correction, O(1) cost).
    double cell_lf_hwm = std::numeric_limits<double>::lowest();
```

- [ ] **Compile gate:**
```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -3
```

### Step 2.2: Maintain lf + cell_lf inside apply_single_margin_linear

The per-j loop in `apply_single_margin_linear` runs from line 413. Inside the loop, after line 457 (`f_lin[off + j] = new_f;`) and before line 458 (`}`), insert:

```cpp
                // T1.B: update lf shadow and propagate delta to cell_lf.
                // new_f is guaranteed finite and <= kLinearOverflowTrip by line 422.
                // Guard for subnormal f_lin (< 1e-300): skip to avoid log(0).
                if (new_f > 1e-300) {
                    double lf_new = std::log(new_f);
                    double delta = lf_new - lf[off + j];
                    lf[off + j] = lf_new;
                    if (std::fabs(delta) > 1e-12) {
                        for (int c : cells_by_margin_cat[off + j]) {
                            cell_lf[c] += delta;
                            if (delta > 0.0) {
                                // Only growing updates can raise the HWM.
                                double val = cell_lf[c] + log_X_init[c];
                                if (std::isfinite(val) && val > cell_lf_hwm)
                                    cell_lf_hwm = val;
                            }
                        }
                    }
                }
```

Full context of the modified per-j block (lines 453-465 after edit):
```cpp
                if (!std::isfinite(new_f) || new_f > kLinearOverflowTrip) {
                    return true;
                }
                rescale_lin[j] = new_f * inv_f_old_lin[j];
                f_lin[off + j] = new_f;
                // T1.B: update lf shadow and propagate delta to cell_lf.
                if (new_f > 1e-300) {
                    double lf_new = std::log(new_f);
                    double delta = lf_new - lf[off + j];
                    lf[off + j] = lf_new;
                    if (std::fabs(delta) > 1e-12) {
                        for (int c : cells_by_margin_cat[off + j]) {
                            cell_lf[c] += delta;
                            if (delta > 0.0) {
                                double val = cell_lf[c] + log_X_init[c];
                                if (std::isfinite(val) && val > cell_lf_hwm)
                                    cell_lf_hwm = val;
                            }
                        }
                    }
                }
            }
            for (int c = 0; c < ct.M_cell; c++) {
                int j = gk[c];
                if (j < 0 || j >= nj) continue;
                if (X_init[c] <= 0.0) continue;
                X_cur[c] *= rescale_lin[j];
            }
            return false;
        };
```

- [ ] **Compile gate:**
```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -3
```

### Step 2.3: Add post-sweep correction block

After the BCD loop closes at line 597 (`}` ending the `else { for (int k...) }` block), and before the `if (overflow_trip && !linear_fallback_used)` at line 653, insert:

```cpp
            // T1.B: correct before X_tilde product overflows.
            // Fires when max_c log(X_init[c] * Π_k f_lin[k][g_k(c)]) reaches log(threshold).
            // shift > 0 guaranteed by >=; x_scale in (0,1]; distributes correction
            // evenly across K margins to preserve relative calibration state.
            if (!overflow_trip && cell_lf_hwm >= std::log(kLinearOverflowThreshold)) {
                double shift = cell_lf_hwm - std::log(kLinearOverflowThreshold);
                double lf_correction = -shift / static_cast<double>(st.K);
                double x_scale = std::exp(-shift);
                for (int k2 = 0; k2 < st.K; k2++) {
                    for (int j = 0; j < st.cat_counts[k2]; j++) {
                        lf[cat_offset[k2] + j] += lf_correction;
                        f_lin[cat_offset[k2] + j] = std::exp(lf[cat_offset[k2] + j]);
                    }
                }
                // X_cur /= exp(shift): maintains X_cur = X_init × W × Π f_lin.
                for (int c = 0; c < ct.M_cell; c++) {
                    cell_lf[c] -= shift;
                    X_cur[c] *= x_scale;
                }
                cell_lf_hwm = std::log(kLinearOverflowThreshold);
                if (st.verbose >= 2) {
                    char msg[128];
                    std::snprintf(msg, sizeof(msg), "iEPPA T1.B renorm shift=%.2e", shift);
                    st.log(msg);
                }
            }
```

### Step 2.4: Add cell_lf + hwm reset to BOTH fallback blocks

**Fallback block 1** (line ~653, `if (overflow_trip && !linear_fallback_used)`):
After `std::fill(lf.begin(), lf.end(), 0.0);` (line 658), add:
```cpp
                std::fill(cell_lf.begin(), cell_lf.end(), 0.0);
                cell_lf_hwm = std::numeric_limits<double>::lowest();
```

**Fallback block 2** (line ~766, `if (overflow_detected)`):
After `std::fill(lf.begin(), lf.end(), 0.0);` (line 772), add:
```cpp
                std::fill(cell_lf.begin(), cell_lf.end(), 0.0);
                cell_lf_hwm = std::numeric_limits<double>::lowest();
```

- [ ] **Compile gate:**
```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -3
```

### Step 2.5: Verify GREEN — T1.B test must PASS

```bash
cd /home/dd/Gemini/leafblower
Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | grep -E "FAIL|PASS|T1|overflow" | head -10
```
Expected: PASS. Also run with verbose=2 to confirm renorm fires:
```bash
Rscript -e '
suppressPackageStartupMessages(library(leafblower))
set.seed(42); K <- 20L; n <- 100000L
cols <- paste0("v", seq_len(K))
data <- as.data.frame(lapply(seq_len(K), function(k) factor(sample(5L, n, TRUE), levels=1:5)))
names(data) <- cols
skewed <- c("1"=0.3,"2"=0.175,"3"=0.175,"4"=0.175,"5"=0.175)
target <- setNames(lapply(seq_len(K), function(.) skewed), cols)
r <- leafblower::harvest(data, target, method="ieppa",
  min_weight=0.2, max_weight=3, max_iterations=50,
  attach_weights=FALSE, verbose=2)
cat("status:", attr(r,"result")$status, "\n")
' 2>&1 | grep -E "T1.B|overflow|status"
```
Expected: "iEPPA T1.B renorm" messages visible, NO "overflow trip", status=0.

### Step 2.6: Run full test suite

```bash
cd /home/dd/Gemini/leafblower && Rscript -e 'devtools::test()' 2>&1 | tail -4
```
Expected: FAIL 0, PASS ≥ 330.

### Step 2.7: Commit T1.B

```bash
git add src/ieppa.cpp
git commit -m "fix(ieppa): T1.B — cell_lf linear-path overflow prevention

Maintain lf[k][j] = log(f_lin[k][j]) and cell_lf[c] = Σ_k lf[k][g_k(c)]
incrementally in apply_single_margin_linear. High-water mark detects
when max X_tilde product approaches kLinearOverflowThreshold. Post-sweep
correction shifts all lf by -shift/K and scales X_cur by exp(-shift),
maintaining invariant X_cur = X_init × W × Π f_lin exactly.

~1.5ms/iter overhead at n=1M (O(M_cell) additions, bucket-sequential).

Closes: leafblower-mmre"
```

---

## Task 3: T4.B — Deferred X_tilde allocation (leafblower-5hru)

**Files:** Modify `src/ieppa.cpp:139, ~662, ~685, ~787`

- [ ] **Step 3.1: Change X_tilde to deferred at line 139**

Find:
```cpp
    std::vector<double> X_tilde(ct.M_cell);
```
Replace with:
```cpp
    std::vector<double> X_tilde;  // deferred: allocated at first log-path/fallback use
```

- [ ] **Step 3.2: Add guard at fallback block 1 (site 1, ~line 662)**

Find the start of the `for` loop that writes `X_tilde[c] = X_init[c]` inside `if (overflow_trip && !linear_fallback_used)`:
```cpp
                for (int c = 0; c < ct.M_cell; c++) {
                    X_tilde[c] = X_init[c];
```
Insert immediately before this `for` loop:
```cpp
                if (X_tilde.empty()) X_tilde.assign(ct.M_cell, 0.0);
```

- [ ] **Step 3.3: Add guard at greedy X_tilde rebuild (site 2, ~line 685)**

Find the greedy X_tilde rebuild loop inside `if (use_greedy)`:
```cpp
                for (int c = 0; c < ct.M_cell; c++) {
                    if (X_init[c] <= 0.0) { X_tilde[c] = 0.0; continue; }
                    double s = log_X_init[c];
                    for (int m = 0; m < st.K; m++) {
```
Insert immediately before this `for` loop:
```cpp
                if (X_tilde.empty()) X_tilde.assign(ct.M_cell, 0.0);
```

- [ ] **Step 3.4: Add guard at main log X_tilde rebuild (site 3, ~line 787)**

Find the main log-path X_tilde rebuild loop inside the `else` of `if (use_greedy)` (after line ~784 comment `// Log-path: X_tilde + capacity`):
```cpp
            for (int c = 0; c < ct.M_cell; c++) {
                if (X_init[c] <= 0.0) { X_tilde[c] = 0.0; continue; }
                double s = log_X_init[c];
                for (int m = 0; m < st.K; m++) {
```
Insert immediately before this `for` loop:
```cpp
            if (X_tilde.empty()) X_tilde.assign(ct.M_cell, 0.0);
```

- [ ] **Step 3.5: Compile gate**
```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -3
```

- [ ] **Step 3.6: Run full test suite**
```bash
cd /home/dd/Gemini/leafblower && Rscript -e 'devtools::test()' 2>&1 | tail -4
```
Expected: same pass count as after Task 2.

- [ ] **Step 3.7: Commit T4.B**
```bash
git add src/ieppa.cpp
git commit -m "fix(ieppa): T4.B — defer X_tilde allocation to log-path entry

Saves 8MB allocation when T1.B prevents overflow (linear path succeeds).
Guards at 3 sites: overflow_trip fallback block, greedy X_tilde rebuild,
main log-path X_tilde rebuild.

Closes: leafblower-5hru"
```

---

## Task 4: T2.A — Log-path cell_lf maintenance (leafblower-9cee)

**Files:** Modify `src/ieppa.cpp:512-519, ~685-693, ~787-800`

### Step 4.1: Update apply_single_margin_log to maintain cell_lf

Find the two-branch lf assignment in `apply_single_margin_log` (lines 512-519):
```cpp
                if (net_log == 1.0) {
                    lf[cat_offset[k] + j] = log_target - log_S_kj;
                } else {
                    double lf_old = lf[cat_offset[k] + j];
                    lf[cat_offset[k] + j] =
                        (1.0 - net_log) * lf_old
                        + net_log * (log_target - log_S_kj);
                }
```
Replace with:
```cpp
                {
                    double lf_old = lf[cat_offset[k] + j];
                    double lf_new = (net_log == 1.0)
                        ? (log_target - log_S_kj)
                        : ((1.0 - net_log) * lf_old + net_log * (log_target - log_S_kj));
                    lf[cat_offset[k] + j] = lf_new;
                    // T2.A: maintain cell_lf incrementally (eliminates K=20 DRAM streams).
                    double delta = lf_new - lf_old;
                    if (std::fabs(delta) > 1e-12) {
                        for (int c : cells_by_margin_cat[cat_offset[k] + j])
                            cell_lf[c] += delta;
                    }
                }
```

### Step 4.2: Replace greedy X_tilde rebuild with cell_lf single-pass

Find the greedy K-loop X_tilde rebuild (~lines 686-693, inside `if (use_greedy)`):
```cpp
                    if (X_init[c] <= 0.0) { X_tilde[c] = 0.0; continue; }
                    double s = log_X_init[c];
                    for (int m = 0; m < st.K; m++) {
                        int gm = ct.g_per_cell[m][c];
                        s += lf[cat_offset[m] + gm];
                    }
                    double s_clip = (s > kLogClip) ? kLogClip : s;
                    X_tilde[c] = std::exp(s_clip);
```
Replace with:
```cpp
                    if (X_init[c] <= 0.0) { X_tilde[c] = 0.0; continue; }
                    // T2.A: single-stream exp via cell_lf (was K=20 DRAM streams)
                    double s = log_X_init[c] + cell_lf[c];
                    double s_clip = (s > kLogClip) ? kLogClip : s;
                    X_tilde[c] = std::exp(s_clip);
```

### Step 4.3: Replace main log-path X_tilde rebuild with cell_lf single-pass

Find the main log-path K-loop rebuild (~lines 787-799):
```cpp
                if (X_init[c] <= 0.0) { X_tilde[c] = 0.0; continue; }
                double s = log_X_init[c];
                for (int m = 0; m < st.K; m++) {
                    int gm = ct.g_per_cell[m][c];
                    s += lf[cat_offset[m] + gm];
                }
                if (s > max_log_X_tilde) max_log_X_tilde = s;
                double s_clip = (s > kLogClip) ? kLogClip : s;
                if (s > kLogClip && U_cell[c] >= 1e299) {
                    overflow_detected = true;
                }
                X_tilde[c] = std::exp(s_clip);
```
Replace with:
```cpp
                if (X_init[c] <= 0.0) { X_tilde[c] = 0.0; continue; }
                // T2.A: single-stream exp via cell_lf (was K=20 DRAM streams)
                double s = log_X_init[c] + cell_lf[c];
                if (s > max_log_X_tilde) max_log_X_tilde = s;
                double s_clip = (s > kLogClip) ? kLogClip : s;
                if (s > kLogClip && U_cell[c] >= 1e299) {
                    overflow_detected = true;
                }
                X_tilde[c] = std::exp(s_clip);
```

- [ ] **Step 4.4: Compile gate**
```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -3
```

- [ ] **Step 4.5: Run full test suite**
```bash
cd /home/dd/Gemini/leafblower && Rscript -e 'devtools::test()' 2>&1 | tail -4
```
Expected: FAIL 0, PASS ≥ 330.

- [ ] **Step 4.6: Commit T2.A**
```bash
git add src/ieppa.cpp
git commit -m "fix(ieppa): T2.A — eliminate K=20 DRAM streams in log-path X_tilde rebuild

Replace K-inner-loop X_tilde rebuild with single sequential exp pass using
cell_lf[c] (maintained incrementally by apply_single_margin_log). Eliminates
K=20 simultaneous g_per_cell[] DRAM streams that caused hardware prefetcher
saturation (~400ms/iter → ~15ms). Applied to both greedy and main log paths.

Closes: leafblower-9cee"
```

---

## Task 5: T3.A — Verify SIMD exp vectorization (leafblower-4nwa)

- [ ] **Step 5.1: Check for libmvec vectorized exp**
```bash
objdump -d /home/dd/R/x86_64-pc-linux-gnu-library/4.6/leafblower/libs/leafblower.so \
  2>/dev/null | grep -c "_ZGVdN4v_exp"
```
If output > 0: already vectorized. Proceed to Step 5.3.
If output == 0: proceed to Step 5.2.

- [ ] **Step 5.2: (only if not vectorized) Add ivdep pragma**

Find the T2.A single-pass exp loop in the main log-path rebuild:
```cpp
            for (int c = 0; c < ct.M_cell; c++) {
                if (X_init[c] <= 0.0) { X_tilde[c] = 0.0; continue; }
                // T2.A: single-stream exp via cell_lf
```
Add before the `for` loop:
```cpp
            #pragma GCC ivdep
```

Rebuild and re-check:
```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -3
objdump -d /home/dd/R/x86_64-pc-linux-gnu-library/4.6/leafblower/libs/leafblower.so \
  2>/dev/null | grep -c "_ZGVdN4v_exp"
```

- [ ] **Step 5.3: Commit**
```bash
git add src/ieppa.cpp
git commit -m "fix(ieppa): T3.A — verify/enable SIMD exp in log-path rebuild

objdump confirms _ZGVdN4v_exp (libmvec 4-double AVX2) in X_tilde exp loop.
GCC -O3 -mavx2 -lmvec auto-vectorizes with or without ivdep pragma.

Closes: leafblower-4nwa"
```
If no code change was needed: `git commit --allow-empty -m "fix(ieppa): T3.A ..."` or just close the ticket with a note.

---

## Task 6: Benchmark verification (leafblower-d9zs)

- [ ] **Step 6.1: Run kk1204 benchmark**

```bash
cd /home/dd/Gemini/leafblower
timeout 120 Rscript -e '
suppressPackageStartupMessages(library(leafblower))
set.seed(42); K <- 20L; n <- 1000000L
cols <- paste0("v", seq_len(K))
data <- as.data.frame(lapply(seq_len(K), function(k) factor(sample(5L, n, TRUE), levels=1:5)))
names(data) <- cols
skewed <- c("1"=0.3,"2"=0.175,"3"=0.175,"4"=0.175,"5"=0.175)
target <- setNames(lapply(seq_len(K), function(.) skewed), cols)
t0 <- proc.time()["elapsed"]
r <- leafblower::harvest(data, target, method="ieppa",
  min_weight=0.2, max_weight=3,
  max_iterations=500, attach_weights=FALSE, verbose=2)
wall <- proc.time()["elapsed"] - t0
res <- attr(r, "result")
cat(sprintf("wall=%.1fs iters=%d max_err=%.4e status=%d\n",
  wall, res$iterations, res$max_error, res$status))
' 2>&1 | grep -E "^wall|T1.B|overflow"
```

Expected:
- "iEPPA T1.B renorm" messages visible (confirms T1.B firing)
- NO "overflow trip" message
- `wall < 30s` for 80 iters (target is P1 gate)

- [ ] **Step 6.2: Close benchmark ticket**

```bash
bd close leafblower-d9zs --reason="T1.B prevents overflow; linear path maintained throughout; wall < 30s confirmed"
```

---

## Task 7: Final verification + ticket closure

- [ ] **Step 7.1: Full test suite**
```bash
cd /home/dd/Gemini/leafblower && Rscript -e 'devtools::test()' 2>&1 | tail -4
```
Expected: FAIL 0, PASS ≥ 330.

- [ ] **Step 7.2: Verify Makevars has no debug flags**
```bash
cat /home/dd/Gemini/leafblower/src/Makevars
```
Must show: `PKG_CXXFLAGS = -std=c++17 -O3 -fopenmp-simd -mavx2 -I. -DSTRICT_R_HEADERS`
No `-pg`, no `-g`.

- [ ] **Step 7.3: Close all tickets**
```bash
bd close leafblower-mmre leafblower-5hru leafblower-9cee leafblower-4nwa 2>&1 | tail -5
bd update leafblower-svbx --notes="Resolved by T1.B+T4.B+T2.A+T3.A in plan 2026-04-27-ieppa-overflow-fix (rev5 spec)"
bd close leafblower-svbx
```

---

## Self-Review

**Spec coverage:**

| Spec item | Task |
|-----------|------|
| Step 0: Remove T1.A (lines 599-651) | Task 0 |
| T1.B: cell_lf + cell_lf_hwm declaration | Task 2.1 |
| T1.B: lf+cell_lf maintenance in apply_single_margin_linear | Task 2.2 |
| T1.B: post-sweep correction block | Task 2.3 |
| T1.B: cell_lf resets in BOTH fallback blocks (lines 658+772) | Task 2.4 |
| T4.B: deferred X_tilde at line 139 | Task 3.1 |
| T4.B: guard site 1 (overflow_trip fallback, ~line 662) | Task 3.2 |
| T4.B: guard site 2 (greedy rebuild, ~line 685) | Task 3.3 |
| T4.B: guard site 3 (main log rebuild, ~line 787) | Task 3.4 |
| T2.A: apply_single_margin_log lf delta + cell_lf update | Task 4.1 |
| T2.A: greedy X_tilde rebuild → cell_lf single-pass | Task 4.2 |
| T2.A: main X_tilde rebuild → cell_lf single-pass | Task 4.3 |
| T3.A: SIMD exp verification | Task 5 |
| Benchmark < 30s | Task 6 |
| Regression test K=20 n=100000 | Task 1 |

**Placeholder scan:** None found.

**Type consistency:** `cell_lf` as `std::vector<double>` used consistently. `cell_lf_hwm` as `double`. `lf[]` reused (no new declaration). `X_tilde` changed to default-constructed (empty).
