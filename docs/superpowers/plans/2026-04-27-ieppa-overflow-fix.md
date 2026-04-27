# iEPPA Linear Overflow Fix Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix iEPPA K=20 benchmark (n=1M) from 118s to <30s by preventing linear-space overflow that permanently falls back to slow log-space path.

**Architecture:** T1.A (f_lin renormalization) keeps the fast linear path active indefinitely by detecting when multiplicative factors approach overflow and rescaling them. T4.B defers the X_tilde allocation to save 8MB when the linear path succeeds. T2.A+T2.B (log-path acceleration) make the rare fallback 10-20× faster via incremental cell_lf maintenance.

**Tech Stack:** C++17, R (devtools::test), ieppa.cpp only. Build: `R CMD INSTALL --preclean . 2>&1 | tail -3`.

---

## File Map

| File | Role |
|------|------|
| `src/ieppa.cpp` | All changes — single file |
| `tests/testthat/test-calibration-solvers.R` | New regression test (Task 1) |

---

## Task 1: TDD RED — Write failing overflow test (leafblower-61l4)

**Files:**
- Modify: `tests/testthat/test-calibration-solvers.R` (append at end)

- [ ] **Step 1: Append the RED test to test-calibration-solvers.R**

```r
# ──────────────────────────────────────────────────────────────────────────────
# T1.A regression: no linear overflow on skewed multi-margin problems
# ──────────────────────────────────────────────────────────────────────────────
test_that("ieppa: no linear overflow trip on K=3 highly skewed (T1.A)", {
  set.seed(42); n <- 2000L
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, TRUE)),
    b = factor(sample(c("1","2","3"), n, TRUE)),
    c = factor(sample(c("1","2","3"), n, TRUE))
  )
  target <- list(
    a = c("1"=0.7, "2"=0.2, "3"=0.1),
    b = c("1"=0.7, "2"=0.2, "3"=0.1),
    c = c("1"=0.7, "2"=0.2, "3"=0.1)
  )
  # Capture verbose log — overflow message goes to stderr via REprintf
  msg_file <- tempfile(fileext = ".txt")
  sink(file(msg_file, "w"), type = "message")
  r <- tryCatch(
    leafblower::harvest(data, target, method = "ieppa",
                        min_weight = 0.1, max_weight = 10,
                        max_iterations = 200,
                        attach_weights = FALSE, verbose = 1),
    finally = sink(type = "message")
  )
  log_lines <- readLines(msg_file)
  overflow_msgs <- grep("overflow trip", log_lines, value = TRUE)

  # Before T1.A: this test FAILS (overflow fires on iter ~30-80 with max_weight=10)
  # After T1.A:  this test PASSES (renorm prevents overflow)
  expect_length(overflow_msgs, 0L,
                label = "no linear overflow trip with T1.A renormalization")
  expect_equal(attr(r, "result")$status, 0L,
               label = "ieppa converges without overflow")
})
```

- [ ] **Step 2: Verify RED**

```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | grep -E "overflow|FAIL|PASS" | head -5
```
Expected: FAIL with "overflow trip" message appearing, `length(overflow_msgs) > 0`.

- [ ] **Step 3: Commit RED test**

```bash
git add tests/testthat/test-calibration-solvers.R
git commit -m "test(ieppa): RED — no overflow trip on K=3 skewed (T1.A)"
```

---

## Task 2: T1.A — f_lin renormalization (leafblower-mmre)

**Files:**
- Modify: `src/ieppa.cpp:212-215` (add threshold), `src/ieppa.cpp:592-593` (insert renorm block)

### Step 2.1: Add kLinearOverflowThreshold after kLinearOverflowTrip

- [ ] **Read** `src/ieppa.cpp` lines 212-215 to confirm:
```cpp
    const double kLinearOverflowTrip = std::pow(
        std::numeric_limits<double>::max() / (2.0 * max_X_init_val),
        1.0 / static_cast<double>(st.K));
```

- [ ] **Add** immediately after line 214 (the closing `;`):
```cpp
    // Renorm threshold: halfway in log-space before the trip.
    // At this point, f_lin geometric mean = sqrt(trip); product f_lin^K = trip^(K/2),
    // leaving sqrt(trip) multiplicative headroom before the trip fires.
    const double kLinearOverflowThreshold = std::sqrt(kLinearOverflowTrip);
```

- [ ] **Compile gate:**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
# Expected: * DONE (leafblower)
```

### Step 2.2: Insert renorm block between K-loop and fallback check

- [ ] **Read** `src/ieppa.cpp` lines 589-594 to confirm this exact text exists:
```cpp
            } else {
                for (int k = 0; k < st.K && !overflow_trip; k++) {
                    if (apply_single_margin_linear(k)) overflow_trip = true;
                }
            }
            if (overflow_trip && !linear_fallback_used) {
```

- [ ] **Insert** between line 593 (`}`) and line 594 (`if (overflow_trip`) — add the renorm block. The full context after the edit:
```cpp
            } else {
                for (int k = 0; k < st.K && !overflow_trip; k++) {
                    if (apply_single_margin_linear(k)) overflow_trip = true;
                }
            }

            // T1.A: Renormalize f_lin per margin when factors approach overflow.
            // Identity: dividing f_lin[k][j] by c_k and multiplying X_cur by c_k
            // cancels per-cell (each cell maps to exactly one j per margin k).
            // Inserted BEFORE the P1.1 capacity block — X_cur still equals
            // X_init * W_prev * Π_k f_lin[k][g_k(c)] at this point.
            if (!overflow_trip) {
                double log_cumul = 0.0;
                for (int k = 0; k < st.K && !overflow_trip; k++) {
                    double log_sum = 0.0; int cnt = 0;
                    for (int j = 0; j < st.cat_counts[k]; j++) {
                        double fkj = f_lin[cat_offset[k] + j];
                        if (std::isfinite(fkj) && fkj > 1e-300) {
                            log_sum += std::log(fkj);
                            cnt++;
                        }
                    }
                    if (cnt == 0) continue;
                    double log_c = log_sum / cnt;
                    if (std::fabs(log_c) < std::log(kLinearOverflowThreshold)) continue;
                    log_cumul += log_c;
                    if (std::fabs(log_cumul) > 700.0) {
                        // Cumulative renorm would overflow X_cur — trip to log-space.
                        overflow_trip = true; break;
                    }
                    double c_k = std::exp(log_c);
                    for (int j = 0; j < st.cat_counts[k]; j++)
                        f_lin[cat_offset[k] + j] /= c_k;
                    for (int c = 0; c < ct.M_cell; c++) X_cur[c] *= c_k;
                    if (st.verbose >= 2) {
                        char msg[128];
                        std::snprintf(msg, sizeof(msg),
                            "iEPPA linear renorm k=%d c_k=%.2e", k, c_k);
                        st.log(msg);
                    }
                }
            }

            if (overflow_trip && !linear_fallback_used) {
```

- [ ] **Compile gate:**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
# Expected: * DONE (leafblower)
```

### Step 2.3: Run RED test — verify it goes GREEN

```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | grep -E "overflow|FAIL|PASS" | head -5
```
Expected: PASS — no overflow trip message, status==0.

- [ ] **Run full test suite:**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -4
# Expected: FAIL 0 (or ≤ existing), PASS ≥ 330, E1+E2 GREEN
```

- [ ] **Commit T1.A:**
```bash
git add src/ieppa.cpp
git commit -m "fix(ieppa): T1.A — renormalize f_lin to prevent linear overflow

Inserts renorm block between K-margin BCD sweep (line 593) and P1.1
capacity block. When f_lin geometric mean per margin exceeds sqrt(trip),
divides f_lin by c_k and compensates X_cur. Identity holds exactly.
Cumulative guard (log_cumul > 700) prevents K simultaneous renorms from
overflowing X_cur. Verbose>=2 logs renorm events for debugging.

Closes: leafblower-mmre"
```

---

## Task 3: T4.B — Defer X_tilde allocation (leafblower-5hru)

**Files:**
- Modify: `src/ieppa.cpp:139` (change allocation), `:594`, `:711`, `:728` (add guards)

- [ ] **Step 3.1: Change line 139 to deferred allocation**

Find:
```cpp
    std::vector<double> X_tilde(ct.M_cell);
```
Replace with:
```cpp
    std::vector<double> X_tilde;  // deferred: allocated at first fallback/log-path use
```

- [ ] **Step 3.2: Add guard at line ~594 (overflow_trip block, site 1)**

Find the start of:
```cpp
            if (overflow_trip && !linear_fallback_used) {
                // One-shot fallback: reset all solver state, switch to log-space,
```
Insert immediately after the `{`:
```cpp
                if (X_tilde.empty()) X_tilde.assign(ct.M_cell, 0.0);
```

- [ ] **Step 3.3: Add guard at line ~711 (overflow_detected block, site 2)**

Find:
```cpp
            if (overflow_detected) {
                // Full state reset on mid-loop break; partial writes to W/X/X_cur undone.
                std::fill(X_cur.begin(),   X_cur.end(),   0.0);
                std::fill(W.begin(),       W.end(),       1.0);
                std::fill(X.begin(),       X.end(),       0.0);
                std::fill(X_tilde.begin(), X_tilde.end(), 0.0);
```
Insert before `std::fill(X_tilde.begin()...`:
```cpp
                if (X_tilde.empty()) X_tilde.assign(ct.M_cell, 0.0);
```

- [ ] **Step 3.4: Add guard at log-path entry (site 3)**

Find the log-path rebuild loop (around line 728):
```cpp
        } else {
            // Log-path: X_tilde + capacity + X_cur unchanged from current implementation.
            bool overflow_detected = false;
            double max_log_X_tilde = -std::numeric_limits<double>::infinity();
            for (int c = 0; c < ct.M_cell; c++) {
                if (X_init[c] <= 0.0) { X_tilde[c] = 0.0; continue; }
```
Insert before the for loop:
```cpp
            if (X_tilde.empty()) X_tilde.assign(ct.M_cell, 0.0);
```

- [ ] **Compile gate:**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
# Expected: * DONE (leafblower)
```

- [ ] **Run full test suite:**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -4
# Expected: same pass count, all GREEN
```

- [ ] **Commit T4.B:**
```bash
git add src/ieppa.cpp
git commit -m "fix(ieppa): T4.B — defer X_tilde allocation to log-path entry

Saves 8MB allocation on every ieppa_solve call when linear path succeeds.
X_tilde lazily allocated at 3 sites: overflow_trip fallback block (~line
594), overflow_detected fallback block (~line 711), and log-path rebuild
entry (~line 728). With T1.A preventing overflow, X_tilde is never
allocated on the kk1204 K=20 benchmark.

Closes: leafblower-5hru"
```

---

## Task 4: Benchmark verification (leafblower-d9zs)

- [ ] **Step 4.1: Run kk1204 benchmark with verbose=1**

```bash
Rscript -e '
suppressPackageStartupMessages(library(leafblower))
set.seed(42); K <- 20L; n <- 1000000L
cols <- paste0("v", seq_len(K))
data <- as.data.frame(lapply(seq_len(K), function(k) factor(sample(5L, n, TRUE), levels=1:5)))
names(data) <- cols
skewed <- c("1"=0.3,"2"=0.175,"3"=0.175,"4"=0.175,"5"=0.175)
target <- setNames(lapply(seq_len(K), function(.) skewed), cols)
t0 <- proc.time()["elapsed"]
r <- harvest(data, target, method="ieppa", min_weight=0.2, max_weight=3,
             max_iterations=500, attach_weights=FALSE, verbose=1)
wall <- proc.time()["elapsed"] - t0
res <- attr(r, "result")
cat(sprintf("wall=%.1fs iters=%d max_err=%.4e status=%d\n",
  wall, res$iterations, res$max_error, res$status))
' 2>&1 | grep -E "^wall|overflow|renorm"
```
Expected:
- NO "overflow trip" line
- wall < 30s
- per_iter stays ~0.105s throughout (no step-change after iter 45)

- [ ] **Step 4.2: Verify per-iteration cost is stable**

```bash
Rscript -e '
suppressPackageStartupMessages(library(leafblower))
set.seed(42); K <- 20L; n <- 1000000L
cols <- paste0("v", seq_len(K))
data <- as.data.frame(lapply(seq_len(K), function(k) factor(sample(5L, n, TRUE), levels=1:5)))
names(data) <- cols
skewed <- c("1"=0.3,"2"=0.175,"3"=0.175,"4"=0.175,"5"=0.175)
target <- setNames(lapply(seq_len(K), function(.) skewed), cols)
# Compare 50-iter vs 80-iter: if ratio stays ~1.6x (not 7x), linear path maintained
for (max_it in c(50L, 80L)) {
  t0 <- proc.time()["elapsed"]
  r <- harvest(data, target, method="ieppa", min_weight=0.2, max_weight=3,
               max_iterations=max_it, attach_weights=FALSE, verbose=0)
  wall <- proc.time()["elapsed"] - t0
  cat(sprintf("max_it=%d wall=%.2fs per_iter=%.3fs\n",
    max_it, wall, wall/attr(r,"result")$iterations))
}
' 2>&1 | grep "max_it"
```
Expected: `per_iter` should be similar for both (both ~0.105-0.120s/iter), NOT 0.105s vs 0.75s.

- [ ] **Step 4.3: Close benchmark ticket**

```bash
bd close leafblower-d9zs 2>&1 | tail -1
```

---

## Task 5: T2.A + T2.B — Incremental cell_lf for log fallback (leafblower-9cee)

**Context**: Only needed when log fallback triggers (rare with T1.A). Implements faster X_tilde rebuild.

**Files:**
- Modify: `src/ieppa.cpp` — add cell_lf vector + incremental maintenance in log path

- [ ] **Step 5.1: Add cell_lf declaration near X_tilde declaration**

After the `std::vector<double> X_tilde;` line (modified in Task 3), add:
```cpp
    std::vector<double> cell_lf;  // deferred: Σ_k lf[k][g_k(c)] per cell; allocated at log-path entry
```

- [ ] **Step 5.2: Initialize cell_lf at log-path entry (row-major population)**

Find the log-path entry guard added in Task 3 (site 3):
```cpp
            if (X_tilde.empty()) X_tilde.assign(ct.M_cell, 0.0);
```
After it, add cell_lf initialization using row-major loop (T2.B):
```cpp
            if (cell_lf.empty()) {
                cell_lf.assign(ct.M_cell, 0.0);
                // T2.B: row-major population — one sequential pass per margin
                // avoids K=20 simultaneous DRAM streams of current column-major rebuild
                for (int m = 0; m < st.K; m++) {
                    for (int c = 0; c < ct.M_cell; c++)
                        cell_lf[c] += lf[cat_offset[m] + ct.g_per_cell[m][c]];
                }
            }
```

- [ ] **Step 5.3: Replace X_tilde rebuild loop with cell_lf-based rebuild**

Find the current X_tilde rebuild loop (the K-inner-loop version, lines ~728-739):
```cpp
            for (int c = 0; c < ct.M_cell; c++) {
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
            // T2.A: X_tilde rebuild via maintained cell_lf — ONE sequential exp pass
            // instead of K=20 simultaneous DRAM streams (eliminates K-inner-loop DRAM thrash)
            for (int c = 0; c < ct.M_cell; c++) {
                if (X_init[c] <= 0.0) { X_tilde[c] = 0.0; continue; }
                double s = log_X_init[c] + cell_lf[c];
                if (s > max_log_X_tilde) max_log_X_tilde = s;
                double s_clip = (s > kLogClip) ? kLogClip : s;
                if (s > kLogClip && U_cell[c] >= 1e299) {
                    overflow_detected = true;
                }
                X_tilde[c] = std::exp(s_clip);
```

- [ ] **Step 5.4: Update cell_lf in apply_single_margin_log**

Find `apply_single_margin_log` function (search for its definition). Inside it, after `lf[cat_offset[k] + j] += delta_lf;` (the lf update), add:
```cpp
                // T2.A: maintain cell_lf incrementally
                if (!cell_lf.empty()) {
                    for (int c : cells_by_margin_cat[cat_offset[k] + j])
                        cell_lf[c] += delta_lf;
                }
```
Note: `delta_lf` is the change applied to `lf[cat_offset[k]+j]`. Find the exact variable name in `apply_single_margin_log` before implementing.

- [ ] **Compile gate:**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
# Expected: * DONE (leafblower)
```

- [ ] **Run tests:**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -4
# Expected: FAIL 0, PASS ≥ 330
```

- [ ] **Commit T2.A+T2.B:**
```bash
git add src/ieppa.cpp
git commit -m "fix(ieppa): T2.A+T2.B — incremental cell_lf eliminates K-stream DRAM thrash

Maintain cell_lf[c] = Σ_k lf[k][g_k(c)] incrementally. X_tilde rebuild
reduces from O(K*M_cell) with K=20 simultaneous streams to O(M_cell)
single sequential exp pass. Initial population uses row-major loop (one
4MB pass per margin). Only allocated at log-fallback entry.

Closes: leafblower-9cee"
```

---

## Task 6: T3.A — Verify SIMD exp vectorization (leafblower-4nwa)

- [ ] **Step 6.1: Check if exp loop is vectorized**

```bash
# Look for libmvec vectorized exp calls in the .so
objdump -d /home/dd/R/x86_64-pc-linux-gnu-library/4.6/leafblower/libs/leafblower.so \
  2>/dev/null | grep -c "_ZGVdN4v_exp"
```
If output > 0: already vectorized. Done — skip Step 6.2.

- [ ] **Step 6.2: (Only if not vectorized) Add ivdep pragma**

Find the X_tilde rebuild loop in the log path (the one now using cell_lf):
```cpp
            for (int c = 0; c < ct.M_cell; c++) {
                if (X_init[c] <= 0.0) { X_tilde[c] = 0.0; continue; }
                double s = log_X_init[c] + cell_lf[c];
```
Add before the loop:
```cpp
            #pragma GCC ivdep
```

- [ ] **Step 6.3: If changed, compile and check again:**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
objdump -d /home/dd/R/x86_64-pc-linux-gnu-library/4.6/leafblower/libs/leafblower.so \
  2>/dev/null | grep -c "_ZGVdN4v_exp"
```

- [ ] **Step 6.4: Commit (even if no change needed):**
```bash
git add src/ieppa.cpp
git commit -m "fix(ieppa): T3.A — verify SIMD exp vectorization in log-path rebuild

objdump confirms _ZGVdN4v_exp (libmvec) called in X_tilde exp loop.
GCC -O3 -mavx2 -lmvec auto-vectorizes 4 doubles/cycle.

Closes: leafblower-4nwa"
```
If no change was needed, `git commit --allow-empty -m "..."` or simply close the ticket with a note.

---

## Task 7: Final verification + close all tickets

- [ ] **Step 7.1: Full test suite**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -4
# Expected: FAIL 0 (or ≤ existing pre-session baseline), PASS ≥ 330
```

- [ ] **Step 7.2: Verify no -pg or -g debug flags remain in Makevars**
```bash
cat src/Makevars
# Must show: PKG_CXXFLAGS = -std=c++17 -O3 -fopenmp-simd -mavx2 -I. -DSTRICT_R_HEADERS
# No -pg, no -g
```

- [ ] **Step 7.3: Close all tickets**
```bash
bd close leafblower-61l4 leafblower-mmre leafblower-5hru leafblower-d9zs leafblower-9cee leafblower-4nwa 2>&1 | tail -6
```

- [ ] **Step 7.4: Delete cheb-renorm worktree (already done earlier) and update svbx ticket**
```bash
bd update leafblower-svbx --notes="Resolved by T1.A+T4.B+T2.A+T2.B+T3.A in plan 2026-04-27-ieppa-overflow-fix" 2>&1 | tail -1
bd close leafblower-svbx 2>&1 | tail -1
```

---

## Self-Review

**Spec coverage:**
| Spec item | Task |
|-----------|------|
| T1.A kLinearOverflowThreshold declaration | Task 2 Step 2.1 |
| T1.A renorm block insertion after line 593 | Task 2 Step 2.2 |
| T1.A isfinite + 1e-300 guard | Task 2 Step 2.2 |
| T1.A cumulative log_cumul ≤ 700 guard | Task 2 Step 2.2 |
| T1.A verbose≥2 renorm log | Task 2 Step 2.2 |
| T4.B deferred X_tilde at line 139 | Task 3 Step 3.1 |
| T4.B guard at line ~594 | Task 3 Step 3.2 |
| T4.B guard at line ~711 | Task 3 Step 3.3 |
| T4.B guard at line ~728 | Task 3 Step 3.4 |
| kk1204 <30s benchmark | Task 4 |
| K=3 skewed regression test | Task 1 |
| T2.A cell_lf incremental | Task 5 |
| T2.B row-major population | Task 5 |
| T3.A SIMD exp verification | Task 6 |

**Placeholder scan:** None found.

**Type consistency:** `cell_lf` used consistently as `std::vector<double>`. `kLinearOverflowThreshold` used in the comparison. `kLinearOverflowTrip` existing name preserved.
