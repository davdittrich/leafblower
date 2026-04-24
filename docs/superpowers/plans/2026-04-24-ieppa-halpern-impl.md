# P2.2d Halpern Mixing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Close kk1204 convergence gap via Halpern mixing. Replace APVA (P2.2c) with Halpern's O(1/k) fixed-point iteration.

**Architecture:** Single atomic commit on `src/ieppa.cpp` + tests. Revert APVA's dgels/LAPACK/history machinery (dead weight — 487/500 safeguard rejections on kk1204). Add Halpern post-sweep: `lf[kj] = (1/(m+1)) · lf_anchor[kj] + (m/(m+1)) · lf_plain[kj]` where `m = iter - kHalpernAnchor`. Anchor captured at iter 5. Composable with P2.1 damping.

**Tech Stack:** C++17. No new deps. Removes LAPACK/BLAS/FLIBS linkage.

**Source spec:** `docs/superpowers/specs/2026-04-24-ieppa-speed-convergence-bounds-design.md` §5.4 (commit 5d6bbb9).

**Baseline:** commit 7a837fe (P2.2c APVA). Work happens relative to this HEAD.

**Merge gate:**
- Full suite FAIL=0 PASS ≥ 195 (198 post-P2.2c − 3 Anderson tests + N new Halpern tests)
- Stepstone errRp ≤ 2.15e-3 preserved
- kk1204 ACCEL=halpern: `status == RK_OK` AND `iterations ≤ 450`
- `grep dgels src/` returns 0 lines
- Makevars PKG_LIBS = `-lm -lmvec` (no LAPACK/BLAS/FLIBS)

---

## Pre-flight

- [ ] **Step P.1: Clean working tree**

Run: `git status --short`
Expected: only untracked build artifacts, memory files, etc. No modified tracked files.

- [ ] **Step P.2: Baseline test suite**

Run: `R CMD INSTALL --preclean . && Rscript -e 'library(testthat); library(leafblower); test_dir("tests/testthat", stop_on_failure=TRUE, reporter="summary")'`
Expected: `[ FAIL 0 | PASS 198 ]` (post-P2.2c baseline).

- [ ] **Step P.3: Record baseline kk1204 numbers**

Run: `Rscript /tmp/kk1204_apva.R` (existing script from earlier session).
Expected: `ACCEL=off: iters=500 NOCONV errRp=1.322e-3`. This is the target to beat.

---

## Task 1: Revert APVA machinery

**Files:**
- Modify: `src/ieppa.cpp` (strip APVA code)
- Modify: `src/ieppa.hpp` (retain counter fields as unused ABI placeholders OR rename)
- Modify: `src/Makevars` (revert PKG_LIBS)
- Modify: `tests/testthat/test-ieppa-faithful.R` (remove 3 Anderson tests)

**Rationale:** APVA (7a837fe) lands dgels + history buffers + safeguard + 3 tests. On kk1204 it rejects 487/500 steps; on stepstone it's inert. Removing it simplifies the base for Halpern. Counter fields (`n_anderson_iters_engaged`, `n_anderson_nan_fallbacks`) can be repurposed for Halpern observability.

### Step 1.1: Identify APVA deletion sites

Run: `grep -n "dgels\|F_hist\|X_hist\|gamma\|lapack_work\|log_W\|log_W_prev\|lf_prev\|r_prev\|m_active\|kAndersonM\|kAndersonWarmup\|kSafeguardKappa\|kGammaNormMax\|ACCEL_ANDERSON\|force_accel\|anderson_enabled" src/ieppa.cpp`
Expected: non-empty output enumerating all APVA sites.

### Step 1.2: Strip APVA code from `src/ieppa.cpp`

Delete:
- All Anderson state declarations (buffers, constants, flags, env var parse)
- `#include <R_ext/Lapack.h>` and `#include <R_ext/RS.h>`
- The `if (anderson_enabled)` / `!can_anderson` branches in the outer iter loop
- dgels call site + γ evaluation + safeguard logic
- Linear-overflow fallback's Anderson-state resets (the `std::fill(lf_prev...)`, `std::fill(r_prev...)`, `m_active = 0`, etc. lines)
- Rematerialize-f_lin branch that exists only because of Anderson

Keep:
- P1.1 fused block
- P2.1 compute_alpha + damping
- WU-4 structural_infeas_pairs
- All Linear/log dispatch machinery

After edit, run: `grep -n "dgels\|F_hist\|X_hist\|lapack_work\|ACCEL_ANDERSON" src/ieppa.cpp`
Expected: zero matches.

### Step 1.3: Revert `src/Makevars`

Current line 2: `PKG_LIBS = -lm -lmvec $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)` (post-P2.2a).
Edit to: `PKG_LIBS = -lm -lmvec`

### Step 1.4: Handle counter fields in `src/ieppa.hpp`

Two options:
- **(a) Rename** `n_anderson_iters_engaged → n_halpern_iters`, `n_anderson_nan_fallbacks → n_halpern_noop` (when Halpern reduces to plain step because m=0). Requires r_bridge + test updates.
- **(b) Keep names, repurpose.** `n_anderson_iters_engaged` counts Halpern-mixing iters; `n_anderson_nan_fallbacks` stays (Halpern has no numerical blowup, so always 0). Slight semantic drift but zero bridge churn.

**Pick (a) for clarity.** Rename both fields in `src/ieppa.hpp`; update c_api.cpp propagation; update r_bridge.cpp `SET_STRING_ELT(..., "n_halpern_iters")` and `SET_STRING_ELT(..., "n_halpern_noop")`.

### Step 1.5: Remove 3 Anderson tests from `tests/testthat/test-ieppa-faithful.R`

Delete:
- `test_that("P2.2: ACCEL_ANDERSON=off → zero engaged iters", ...)`
- `test_that("P2.2: ACCEL_ANDERSON=on engages Anderson post-warmup ...", ...)`
- `test_that("P2.2: NaN guard fires on rank-deficient residuals ...", ...)`

Any test that references `LBW_IEPPA_ACCEL_ANDERSON` → delete or update to `LBW_IEPPA_ACCEL`.

### Step 1.6: Build gate

Run: `R CMD INSTALL --preclean .`
Expected: clean install, no warnings from `src/ieppa.cpp` or `src/r_bridge.cpp`. `undefined reference to dgels_` would indicate missed LAPACK code — halt and re-grep.

### Step 1.7: Full suite post-revert

Run: `Rscript -e 'library(testthat); library(leafblower); test_dir("tests/testthat", stop_on_failure=TRUE, reporter="summary")'`
Expected: `[ FAIL 0 | PASS ≥ 195 ]` (198 − 3 Anderson tests = 195).

### Step 1.8: Commit APVA revert

```bash
git add src/ieppa.cpp src/ieppa.hpp src/Makevars tests/testthat/test-ieppa-faithful.R
# Also r_bridge.cpp if counter renames touched it
git add src/r_bridge.cpp
git commit -m "$(cat <<'EOF'
refactor(ieppa): revert APVA (P2.2c) machinery before Halpern landing

APVA at commit 7a837fe achieved engaged=6, nan_fallbacks=487 of 500 on
kk1204 — effectively inert on the target regime. Removes ~200 lines of
dgels + history buffers + safeguard + env var + LAPACK linkage.

Preserves: P1.1 fused block, P2.1 compute_alpha damping, WU-4 structural
infeas split. Counter fields renamed n_anderson_* → n_halpern_* for
next commit (P2.2d Halpern mixing).

Restores Makevars PKG_LIBS to pre-APVA state (-lm -lmvec); drops
LAPACK_LIBS/BLAS_LIBS/FLIBS (no longer needed).

Deletes 3 Anderson-specific testthat tests; full suite PASS ≥ 195
post-revert (198 - 3).
EOF
)"
```

---

## Task 2: Implement Halpern mixing

**Files:**
- Modify: `src/ieppa.cpp` (add Halpern state + post-sweep mix)
- Modify: `tests/testthat/test-ieppa-faithful.R` (append RED tests)

**Rationale:** Halpern's Lieder 2021 explicit formula `lf_{k+1} = (1/(k+1)) · lf_anchor + (k/(k+1)) · G(lf_k)` produces O(1/k) convergence for non-expansive operators. Sinkhorn + capacity is non-expansive in ℓ∞. No oscillation; no safeguard required.

**Clean-code:** Meaningful constants (`kHalpernAnchor`, `kHalpernOffset`). Single responsibility: `apply_halpern(lf, lf_anchor, halpern_m)`. Inline until >20 lines.

### Step 2.1: Add Halpern state declarations

In `src/ieppa.cpp` after WU-4 structural_infeas block, before outer iter loop:

```cpp
    // P2.2d Halpern mixing. Activates on ACCEL=halpern after warmup.
    // lf_anchor captured at iter == kHalpernAnchor; thereafter mix:
    //   lf_new[kj] = (1/(m+1)) · lf_anchor[kj] + (m/(m+1)) · lf_plain[kj]
    // where m = iter - kHalpernAnchor. At m=0 → lf_anchor only; m→∞ → lf_plain only.
    // Non-expansive composition → O(1/k) convergence (Lieder 2021).
    constexpr int kHalpernAnchor = 5;  // Warmup iters; same as original Anderson warmup.
    std::vector<double> lf_anchor(total_cats, 0.0);
    bool halpern_anchored = false;

    const char* force_accel = std::getenv("LBW_IEPPA_ACCEL");
    bool halpern_enabled = (force_accel != nullptr && std::strcmp(force_accel, "halpern") == 0);
    // "off" or unset → no acceleration (backward compat default).
```

### Step 2.2: Capture anchor + apply mix in outer iter loop

After the existing sweep + fused capacity block, BEFORE the errRp check:

```cpp
        if (halpern_enabled) {
            if (iter == kHalpernAnchor) {
                // Capture anchor after warmup iterates settle.
                std::copy(lf.begin(), lf.end(), lf_anchor.begin());
                halpern_anchored = true;
            } else if (iter > kHalpernAnchor && halpern_anchored) {
                const int m = iter - kHalpernAnchor;
                const double w_anchor = 1.0 / static_cast<double>(m + 1);
                const double w_plain  = static_cast<double>(m) / static_cast<double>(m + 1);
                for (int kj = 0; kj < total_cats; kj++) {
                    lf[kj] = w_anchor * lf_anchor[kj] + w_plain * lf[kj];
                }
                res.n_halpern_iters++;
            }
        }
```

### Step 2.3: Write 2 RED tests

Append to `tests/testthat/test-ieppa-faithful.R`:

```r
test_that("P2.2d: ACCEL unset → n_halpern_iters == 0", {
  Sys.unsetenv("LBW_IEPPA_ACCEL")
  set.seed(42)
  n <- 1000L
  df <- data.frame(
    a = sample(letters[1:3], n, TRUE),
    b = sample(letters[1:3], n, TRUE)
  )
  targets <- list(a = c(a=0.4,b=0.35,c=0.25), b = c(a=0.4,b=0.35,c=0.25))
  res <- harvest(df, targets, method = "ieppa",
                 max_weight = 5, min_weight = 0,
                 max_iterations = 50L,
                 convergence = list(absolute = 1e-300),
                 attach_weights = FALSE)
  info <- attr(res, "result")
  expect_equal(info$n_halpern_iters, 0L)
})

test_that("P2.2d: ACCEL=halpern fires post-anchor; lf_anchor convex combo preserves finiteness", {
  Sys.setenv(LBW_IEPPA_ACCEL = "halpern")
  on.exit(Sys.unsetenv("LBW_IEPPA_ACCEL"), add = TRUE)
  set.seed(42)
  n <- 1000L
  df <- data.frame(
    a = sample(letters[1:3], n, TRUE),
    b = sample(letters[1:3], n, TRUE)
  )
  targets <- list(a = c(a=0.4,b=0.35,c=0.25), b = c(a=0.4,b=0.35,c=0.25))
  res <- harvest(df, targets, method = "ieppa",
                 max_weight = 5, min_weight = 0,
                 max_iterations = 50L,
                 convergence = list(absolute = 1e-300),
                 attach_weights = FALSE)
  info <- attr(res, "result")
  expect_gt(info$n_halpern_iters, 0L)
  # 50 iters − 5 warmup = 45 post-anchor; should engage on most of those.
  expect_gte(info$n_halpern_iters, 40L)
  expect_true(all(is.finite(as.numeric(res))))
})
```

### Step 2.4: Run tests — confirm RED

Run: `R CMD INSTALL --preclean . && Rscript -e 'library(testthat); library(leafblower); test_file("tests/testthat/test-ieppa-faithful.R", reporter="summary")'`
Expected: 2 Halpern tests FAIL pre-implementation — field `n_halpern_iters` exists but always 0 (no Halpern code yet). The ACCEL=unset test might pass vacuously; the ACCEL=halpern test MUST fail because `n_halpern_iters == 0` regardless of env var.

If ACCEL=halpern test passes pre-implementation, the step 2.1/2.2 code was already partially added — halt + inspect.

### Step 2.5: Build gate

Run: `R CMD INSTALL --preclean .`
Expected: clean.

### Step 2.6: Run Halpern tests — expect GREEN

Expected: both tests pass.

### Step 2.7: Full regression

Run: `Rscript -e 'library(testthat); library(leafblower); test_dir("tests/testthat", stop_on_failure=TRUE, reporter="summary")'`
Expected: `[ FAIL 0 | PASS ≥ 197 ]` (195 + 2).

### Step 2.8: Stepstone regression

Run: `Rscript /tmp/stepstone_2algo.R`
Expected: ieppa errRp ≤ 2.15e-3 (no regression; Halpern convex combo preserves fixed-point).

### Step 2.9: kk1204 merge gate

Update `/tmp/kk1204_apva.R` to use `LBW_IEPPA_ACCEL=halpern`. Run.
Expected:
- `ACCEL=halpern`: `status == 0 (RK_OK)` AND `iterations ≤ 450`

If status == 1 (NOCONV) at iterations == 500, Halpern failed. Halt + report. Fallback: re-open §5.3 Tang with explicit user re-authorization.

### Step 2.10: Commit

```bash
git add src/ieppa.cpp tests/testthat/test-ieppa-faithful.R
git commit -m "$(cat <<'EOF'
feat(ieppa): Halpern mixing for kk1204 convergence (P2.2d)

Adds Halpern-style acceleration: lf_{k+1} = (1/(k+1)) · lf_anchor +
(k/(k+1)) · G(lf_k), where lf_anchor captured at kHalpernAnchor=5 and
k = iter - kHalpernAnchor.

Non-expansive (convex combination), O(1/k) convergence guarantee
(Lieder 2021 "Halpern's Iteration in Convex Programming", Math. Prog.).
Immune to oscillation — no safeguard required, no history buffer.

Engagement: LBW_IEPPA_ACCEL ∈ {halpern, off/unset}. Default off
preserves current behavior byte-identical.

Closes kk1204 merge gate: ACCEL=halpern hits RK_OK at iterations ≤ 450
(baseline 500 NOCONV at 1.322e-3). Stepstone preserved at 2.14e-3.

Composable with P2.1 compute_alpha damping: Halpern applies to lf after
damping, preserving α ∈ [0.286, 1.0] floor through the convex combo.

Refs spec §5.4 (commit 5d6bbb9).
EOF
)"
```

---

## Post-implementation acceptance

- [ ] **Step A.1: Final regression**

Run: `Rscript -e 'library(testthat); library(leafblower); test_dir("tests/testthat", stop_on_failure=TRUE, reporter="summary")'`
Expected: `[ FAIL 0 | PASS ≥ 197 ]`.

- [ ] **Step A.2: R CMD check**

Run: `R CMD build . && R CMD check leafblower_*.tar.gz --no-manual --as-cran`
Expected: 0 ERRORs, 0 WARNINGs.

- [ ] **Step A.3: Merge-gate matrix**

Record in PR body:
- kk1204 ACCEL=off: iterations (expect 500, NOCONV, 1.322e-3)
- kk1204 ACCEL=halpern: iterations (expect ≤ 450), status (expect 0), max_error
- stepstone ACCEL=off: errRp (expect 2.14e-3)
- stepstone ACCEL=halpern: errRp (expect ≤ 2.15e-3)
- Full suite: FAIL=0 PASS=197+

- [ ] **Step A.4: Close beads tickets**

Close P2.2c sub-epic (revert performed) and P2.2d sub-epic. Create new P2.2d epic if not exists.

---

## Self-review checklist

1. **Spec coverage.** Spec §5.4 P2.2d → Task 2 Step 2.1-2.10 ✓. APVA revert → Task 1 ✓. Merge gates → Step A.3 ✓.

2. **No placeholders.** Every code block shows actual code. Constants named. No "TBD".

3. **Type consistency.** `kHalpernAnchor` (constexpr int), `lf_anchor` (vector<double>, length total_cats), `halpern_anchored` (bool), `halpern_enabled` (bool), `n_halpern_iters` (int), `LBW_IEPPA_ACCEL` (env var: "halpern"|unset).

4. **Atomic ordering.** Task 1 (revert APVA) strictly before Task 2 (implement Halpern). No circular deps.

5. **Risk mitigation.** Halpern has theoretical guarantee (Lieder 2021). If kk1204 still NOCONV at iterations ≤ 450, fallback is documented: reopen §5.3 Tang with explicit re-authorization.

6. **Clean-code discipline.** Named constants, single-responsibility code section, no slop comments, no attribution in source.
