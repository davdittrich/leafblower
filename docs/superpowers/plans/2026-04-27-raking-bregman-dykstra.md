# Raking: Bregman Dykstra + SOR + Greedy Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace Euclidean Dykstra corrections in raking with Bregman/KL Dykstra to guarantee convergence; wire SOR and Greedy scheduler; A/B test on stepstone-fulldata.

**Architecture:** Three independent changes to `src/raking.cpp`: (1) Bregman Dykstra replaces additive `p[c]`/`q_hyp` corrections with multiplicative ratio corrections — unifies all projections under KL geometry. (2) SOR wires `st.sor_cfg` into the IPF scaling step. (3) Greedy sorts the K-margin loop by per-margin residual. A/B test: Approach A (Bregman+SOR) vs C (Bregman+SOR+Greedy).

**Tech Stack:** C++17 (`src/raking.cpp`), R (`R/harvest.R`, `data-raw/gen_raking_obs_ref.R`). Build: `R CMD INSTALL --preclean . 2>&1 | tail -3`. Test: `Rscript -e 'devtools::test()' 2>&1 | tail -3`.

**Ticket**: leafblower-9tua
**Spec**: `docs/superpowers/specs/2026-04-27-raking-bregman-dykstra-design.md`

---

## File Map

| File | Change |
|------|--------|
| `tests/testthat/test-calibration-solvers.R` | RED test (Task 1) |
| `src/raking.cpp:104-105,141-147,177-183,334-346,156-175` | Bregman Dykstra + SOR + Greedy (Tasks 2-4) |
| `R/harvest.R:48` | sor docstring (Task 5) |
| `data-raw/gen_raking_obs_ref.R` | Fixture regen (Task 6) |
| `tests/testthat/fixtures/raking_obs_reference_stepstone.rds` | Regenerated (Task 6) |

---

## Task 1: RED Test (leafblower-9tua, step 1)

**Files:** `tests/testthat/test-calibration-solvers.R` (append)

- [ ] **Step 1.1: Append the RED test**

```r
# ── Raking Bregman Dykstra RED test ─────────────────────────────────────────
# Before fix: Euclidean hyperplane correction changes fixed point vs pure IPF.
# After fix:  multiplicative hyperplane = KL projection → same fixed point as ieppa.
# RED: expect_equal(wkl_raking, wkl_ieppa, tol=1e-4) FAILS before Bregman Dykstra.
# solver_objective field confirmed: raking.cpp:349, harvest.R:280.
# ────────────────────────────────────────────────────────────────────────────
test_that("raking-bregman: unconstrained raking matches ieppa weight_kl (unified KL fixed point)", {
  set.seed(1L); n <- 2000L
  df <- data.frame(
    v1 = factor(sample(3L, n, TRUE)),
    v2 = factor(sample(2L, n, TRUE))
  )
  tgt <- list(v1 = c("1"=0.5, "2"=0.3, "3"=0.2),
              v2 = c("1"=0.6, "2"=0.4))

  r_raking <- leafblower::harvest(df, tgt, method="raking",
    min_weight=0, max_weight=Inf, max_iterations=500L,
    attach_weights=FALSE)
  r_ieppa  <- leafblower::harvest(df, tgt, method="ieppa",
    min_weight=0, max_weight=Inf, max_iterations=500L,
    attach_weights=FALSE)

  wkl_raking <- attr(r_raking, "result")$convergence_used$solver_objective
  wkl_ieppa  <- attr(r_ieppa,  "result")$convergence_used$solver_objective

  # Both are unconstrained IPF → same KL minimum → same weight_kl.
  # Before Bregman: Euclidean hyperplane shifts raking fixed point → FAIL.
  # After  Bregman: unified KL geometry → PASS.
  expect_equal(wkl_raking, wkl_ieppa, tolerance=1e-4,
               label="unconstrained raking must reach same weight_kl as ieppa")
})
```

- [ ] **Step 1.2: Verify RED (must FAIL, not error)**

```bash
cd /home/dd/Gemini/leafblower
Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | grep -E "raking-bregman|FAIL|PASS" | head -5
```
Expected: test FAILS (wkl values differ, e.g., 0.XX vs 0.YY).

- [ ] **Step 1.3: Commit RED test**

```bash
git add tests/testthat/test-calibration-solvers.R
git commit -m "test(raking): RED — Bregman Dykstra KL fixed point

Unconstrained raking must match ieppa weight_kl when Bregman Dykstra
unifies KL geometry. Before fix: Euclidean hyperplane shifts fixed point.

Closes: leafblower-9tua (step 1)"
```

---

## Task 2: Bregman Dykstra — p[], q_hyp, hyperplane_step, finalizer (leafblower-9tua, step 2)

**Files:** `src/raking.cpp:104-105,141-147,177-183,334-346`

### Step 2.1: Change p[] and q_hyp initialization (lines 104-105)

Find:
```cpp
    std::vector<double> p(ct.M_cell, 0.0);
    double q_hyp = 0.0;
```
Replace with:
```cpp
    std::vector<double> p(ct.M_cell, 1.0);   // Bregman: multiplicative identity
    double q_hyp = 1.0;                        // Bregman: multiplicative identity
```

### Step 2.2: Replace hyperplane_step lambda (lines 141-147)

Find the entire lambda:
```cpp
    auto hyperplane_step = [&]() {
        double s = 0.0;
        for (int c = 0; c < ct.M_cell; c++) { X[c] += q_hyp; s += X[c]; }
        double shift = (static_cast<double>(st.n) - s) / static_cast<double>(ct.M_cell);
        for (int c = 0; c < ct.M_cell; c++) X[c] += shift;
        q_hyp = -shift;
    };
```
Replace with:
```cpp
    // Bregman/KL Dykstra hyperplane: multiplicative scale instead of additive shift.
    // KL projection onto {sum(X)=n}: scale all cells uniformly.
    auto hyperplane_step = [&]() {
        double s = 0.0;
        for (int c = 0; c < ct.M_cell; c++) { X[c] *= q_hyp; s += X[c]; }
        const double scale_hp = static_cast<double>(st.n) / s;
        for (int c = 0; c < ct.M_cell; c++) X[c] *= scale_hp;
        q_hyp = (scale_hp > 0.0) ? 1.0 / scale_hp : 1.0;
    };
```

### Step 2.3: Replace box correction (lines 177-183)

Find:
```cpp
        // Dykstra box: X[c] = clamp(X[c] + p[c], L_cell[c], U_cell[c])
        for (int c = 0; c < ct.M_cell; c++) {
            double yc = X[c] + p[c];
            double Xc = std::clamp(yc, L_cell[c], U_cell[c]);
            p[c] = yc - Xc;
            X[c] = Xc;
        }
```
Replace with:
```cpp
        // Bregman/KL Dykstra box: multiplicative correction p[c] (init=1.0).
        // Guard: if Xc=0 (min_weight=0 and structural zero cell), p[c]=1.0 (no correction).
        for (int c = 0; c < ct.M_cell; c++) {
            double yc = X[c] * p[c];
            double Xc = std::clamp(yc, L_cell[c], U_cell[c]);
            p[c] = (Xc > 0.0) ? yc / Xc : 1.0;
            X[c] = Xc;
        }
```

### Step 2.4: Replace post-loop finalizer (lines 334-346)

Find:
```cpp
    // Post-loop Dykstra finalizer at cell level.
    for (int fixup = 0; fixup < 20; fixup++) {
        bool box_ok = true;
        for (int c = 0; c < ct.M_cell; c++) {
            double yc = X[c] + p[c];
            double Xc = std::clamp(yc, L_cell[c], U_cell[c]);
            p[c] = yc - Xc;
            if (yc != Xc) box_ok = false;
            X[c] = Xc;
        }
        hyperplane_step();
        if (box_ok) break;
    }
```
Replace with:
```cpp
    // Post-loop Bregman/KL Dykstra finalizer (multiplicative pattern).
    for (int fixup = 0; fixup < 20; fixup++) {
        bool box_ok = true;
        for (int c = 0; c < ct.M_cell; c++) {
            double yc = X[c] * p[c];
            double Xc = std::clamp(yc, L_cell[c], U_cell[c]);
            p[c] = (Xc > 0.0) ? yc / Xc : 1.0;
            if (yc != Xc) box_ok = false;
            X[c] = Xc;
        }
        hyperplane_step();
        if (box_ok) break;
    }
```

### Step 2.5: Compile gate

```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -3
```
Must show `* DONE (leafblower)`.

### Step 2.6: Verify RED test goes GREEN

```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | grep -E "raking-bregman|PASS|FAIL" | head -3
```
Expected: `raking-bregman` test PASSES.

### Step 2.7: Full test suite

```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL ≤ 2 (pre-existing: test-sor.R/lhs, test-ieppa-nonuniform-d.R).

### Step 2.8: Commit Bregman Dykstra

```bash
git add src/raking.cpp
git commit -m "fix(raking): Bregman/KL Dykstra replaces Euclidean Dykstra

Unifies all projections under KL geometry (Bauschke & Lewis 2000).
- p[c]: init 1.0, ratio correction yc/Xc (was 0.0, difference yc-Xc)
- q_hyp: init 1.0, scale factor 1/scale (was 0.0, additive shift)
- hyperplane_step: multiplicative X*=scale (was additive X+=shift)
- Post-loop finalizer: same multiplicative pattern
Guards: Xc=0 → p[c]=1.0; scale=0 → q_hyp=1.0

Closes: leafblower-9tua (step 2)"
```

---

## Task 3: SOR wiring (leafblower-9tua, step 3)

**Files:** `src/raking.cpp:149-175` (inner loop), `R/harvest.R:48`

SOR is already in `CalibState::sor_cfg` (accessible as `st.sor_cfg`). Need to:
1. Read SOR state from CalibState
2. Apply per-margin omega to IPF scaling step
3. Update harvest.R docstring

### Step 3.1: Add SOR state variables after the convergence tracking declarations (~line 125)

After `std::vector<double> W_best(ct.M_cell, 0.0);`, add:
```cpp
    // SOR: wire st.sor_cfg into raking's IPF step (same API as ieppa).
    // Apply only when bounds are active (oscillation risk).
    const bool sor_active  = st.sor_cfg.enabled &&
                              (st.min_weight > 0.0 || hi < 1e300);
    const bool sor_auto    = st.sor_cfg.auto_adapt;
    const double omega_min = st.sor_cfg.omega_min;    // default 0.3
    const double omega_init = st.sor_cfg.omega_init;  // default 1.0
    std::vector<double> sor_omega(st.K, omega_init);
    std::vector<double> sor_prev_errRp(st.K, std::numeric_limits<double>::infinity());
```

### Step 3.2: Replace IPF scaling application (lines 171-174) to apply SOR

Find:
```cpp
            for (int c = 0; c < ct.M_cell; c++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k]) X[c] *= scale[g];
            }
```
Replace with:
```cpp
            // Apply IPF scaling with optional SOR damping.
            const double eff_omega = sor_active ? sor_omega[k] : 1.0;
            for (int c = 0; c < ct.M_cell; c++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k]) {
                    const double s_g = scale[g];
                    if (s_g > 0.0)  // guard: S=0 bucket → skip
                        X[c] *= (eff_omega == 1.0) ? s_g : std::pow(s_g, eff_omega);
                }
            }

            // SOR auto-adaptation: sign flip in per-margin errRp → damp omega.
            if (sor_active && sor_auto) {
                double ek = 0.0;
                for (int j = 0; j < st.cat_counts[k]; j++) {
                    double e = std::fabs(scale[j] - 1.0);
                    if (e > ek) ek = e;
                }
                if (sor_prev_errRp[k] < ek) {
                    sor_omega[k] = std::max(omega_min, sor_omega[k] * 0.7);
                } else {
                    sor_omega[k] = std::min(1.0, sor_omega[k] * 1.05);
                }
                sor_prev_errRp[k] = ek;
            }
```

### Step 3.3: Update harvest.R sor docstring (line 48)

Find:
```r
#' @param sor Named list for SOR adaptive under-relaxation (iEPPA only).
```
Replace with:
```r
#' @param sor Named list for SOR adaptive under-relaxation (iEPPA and raking).
```

### Step 3.4: Compile + full test suite

```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -3
Rscript -e 'devtools::test()' 2>&1 | tail -3
```

### Step 3.5: Commit SOR

```bash
git add src/raking.cpp R/harvest.R
git commit -m "feat(raking): wire SOR adaptive under-relaxation

Uses CalibState::sor_cfg (same struct as ieppa). Applied only when
bounds are active (st.min_weight>0 || hi<1e300). S=0 bucket guard
prevents div-by-zero in pow(). Auto-adaptation mirrors ieppa logic.

Closes: leafblower-9tua (step 3)"
```

---

## Task 4: Greedy margin ordering (leafblower-9tua, step 4)

**Files:** `src/raking.cpp:149-175`

### Step 4.1: Add per-margin errRp tracking and order vector

After the SOR state variables (end of Task 3's additions), add:
```cpp
    // Greedy scheduler: per-margin residuals + sort order.
    // errRp_k[k] updated each iter during sweep using bucket[] sums.
    std::vector<double> errRp_k(st.K, 1.0 / st.K);  // init uniform
    std::vector<int> margin_order(st.K);
    std::iota(margin_order.begin(), margin_order.end(), 0);
    const bool use_greedy = (st.scheduler.mode == SchedulerMode::GREEDY);
```

### Step 4.2: Add Greedy sort before K-margin sweep and update errRp_k during sweep

Find the K-margin sweep header (line 156):
```cpp
        for (int k = 0; k < st.K; k++) {
```

Replace with:
```cpp
        // Greedy: sort margins descending by per-margin errRp from previous iter.
        if (use_greedy)
            std::sort(margin_order.begin(), margin_order.end(),
                      [&](int a, int b){ return errRp_k[a] > errRp_k[b]; });

        for (int ki = 0; ki < st.K; ki++) {
            const int k = use_greedy ? margin_order[ki] : ki;
```
(update closing `}` at the end of the k-loop to close `ki` instead — this is purely a rename from `k` to `ki`/`k`)

Also update the close of the per-margin block to compute `errRp_k[k]` using `bucket[]` (already filled during the sweep):
```cpp
            // Track per-margin residual for next iter's Greedy sort.
            // bucket[j] filled above; W_total computed at top of iter loop.
            {
                double ek = 0.0;
                for (int j = 0; j < st.cat_counts[k]; j++) {
                    double e = std::fabs(bucket[j] / W_total - st.targets[k][j]);
                    if (e > ek) ek = e;
                }
                errRp_k[k] = ek;
            }
```
Insert this block after the SOR adaptation block (end of per-margin body), before the closing `}` of the `ki`-loop.

### Step 4.3: Add `#include <numeric>` if needed

```bash
grep -n "iota\|#include.*numeric" src/raking.cpp | head -3
```
If `<numeric>` not included, add to the include block at the top of raking.cpp.

### Step 4.4: Compile + full test suite

```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -3
Rscript -e 'devtools::test()' 2>&1 | tail -3
```

### Step 4.5: Commit Greedy

```bash
git add src/raking.cpp
git commit -m "feat(raking): Greedy margin scheduler

Per-margin errRp_k computed during sweep using existing bucket[] sums.
When scheduler='greedy', margins sorted descending by errRp_k each iter.
Same SchedulerMode::GREEDY API as ieppa.

Closes: leafblower-9tua (step 4)"
```

---

## Task 5: A/B Benchmark (leafblower-9tua, step 5)

Run stepstone-fulldata with Approach A and C. Compare max_err, marg_kl, weight_kl, wall.

### Step 5.1: Run A/B benchmark

```bash
cd /home/dd/Gemini/leafblower
OMP_NUM_THREADS=1 Rscript -e '
suppressPackageStartupMessages({library(arrow); library(leafblower); library(jsonlite)})
df <- read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet"); df$uuid <- NULL
tgt <- lapply(fromJSON("benchmarks/stepstone_fulldata_bench_targets.json"),function(t){t<-unlist(t);t/sum(t)})
for(nm in names(tgt)) df[[nm]] <- factor(df[[nm]])

run <- function(label, ...) {
  t0 <- proc.time()["elapsed"]
  r  <- suppressWarnings(leafblower::harvest(df, tgt, method="raking",
    max_weight=5, min_weight=0, max_iterations=500L, attach_weights=FALSE, ...))
  wall <- proc.time()["elapsed"] - t0
  res  <- attr(r,"result")
  W    <- sum(as.numeric(r)); w <- as.numeric(r)
  max_err <- 0; marg_kl <- 0
  for(nm in names(tgt)){
    S <- tapply(w, df[[nm]], sum, default=0)[names(tgt[[nm]])]; S[is.na(S)] <- 0
    tk <- unname(tgt[[nm]]); Sr <- S/W
    max_err <- max(max_err, max(abs(Sr-tk)))
    safe <- tk>0 & Sr>0
    marg_kl <- marg_kl + sum(ifelse(safe, tk*log(tk/pmax(Sr,1e-300)), 0))
  }
  cat(sprintf("%-30s wall=%5.1fs iters=%3d max_err=%.4e marg_kl=%.3e weight_kl=%.3e\n",
    label, wall, res$iterations, max_err, marg_kl,
    res$convergence_used$solver_objective))
}

run("baseline (no SOR/Greedy)")
run("Approach A (Bregman+SOR)",   sor=list(auto=TRUE))
run("Approach C (Bregman+SOR+G)", sor=list(auto=TRUE), scheduler="greedy")
' 2>&1 | grep -v "^$"
```

### Step 5.2: Verify AC4 — winner must beat baseline

Expected: at least one of A or C achieves max_err < 5.44e-3 (current raking baseline).

If both fail AC4, **HALT** and do NOT proceed to fixture regeneration (Task 6). File a blocker ticket.

### Step 5.3: Document winner in commit

```bash
git commit --allow-empty -m "bench(raking): A/B test results — [winner] chosen

Approach A (Bregman+SOR): max_err=X marg_kl=Y weight_kl=Z
Approach C (Bregman+SOR+G): max_err=X marg_kl=Y weight_kl=Z
Winner: [A or C] — lower max_err and marg_kl

Closes: leafblower-9tua (step 5)"
```

---

## Task 6: Fixture Regeneration (leafblower-9tua, step 6)

**PREREQUISITE**: Task 5 (AC4) must pass before this task. If AC4 failed, HALT.

### Step 6.1: Regenerate obs-level reference

```bash
cd /home/dd/Gemini/leafblower
OMP_NUM_THREADS=1 Rscript data-raw/gen_raking_obs_ref.R 2>&1 | tail -5
```

### Step 6.2: Verify fixture quality

```bash
Rscript -e '
ref <- readRDS("tests/testthat/fixtures/raking_obs_reference_stepstone.rds")
cat("max_error:", ref$max_error, "\n")
cat("A5 gate (< 5.44e-3):", ref$max_error < 5.44e-3, "\n")
'
```
Expected: max_error < 5.44e-3.

### Step 6.3: Run A5 test

```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | grep "A5\|PASS\|FAIL" | head -5
```

### Step 6.4: Full suite

```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```

### Step 6.5: Commit fixture

```bash
git add data-raw/gen_raking_obs_ref.R tests/testthat/fixtures/raking_obs_reference_stepstone.rds
git commit -m "fix(raking): regenerate A5 fixture with Bregman Dykstra baseline

Old fixture used Euclidean Dykstra fixed point. New fixture matches
the Bregman KL-minimum. AC4 confirmed (max_err < 5.44e-3) before regen.

Closes: leafblower-9tua (step 6)"
```

---

## Self-Review

**Spec coverage:**

| Spec item | Task |
|-----------|------|
| Bregman p[c] box correction (ratio, Xc>0 guard) | Task 2.3 |
| Bregman hyperplane_step (multiplicative, scale>0 guard) | Task 2.2 |
| Bregman finalizer (multiplicative) | Task 2.4 |
| RED test provably RED before fix | Task 1 |
| SOR wired via st.sor_cfg, bounds-active condition | Task 3.1-3.2 |
| S=0 guard in IPF scaling | Task 3.2 |
| Greedy: per-margin errRp_k in sweep, sort descending | Task 4.1-4.2 |
| harvest.R @param sor docstring updated | Task 3.3 |
| A/B benchmark: A vs C on stepstone-fulldata | Task 5 |
| AC4 gate before fixture regen | Task 5.2 |
| Fixture regeneration command | Task 6.1 |
| Full test suite ≤ 2 pre-existing fails | Tasks 2.7, 3.4, 4.4, 6.4 |

**Placeholder scan:** None.

**Type consistency:** `sor_omega[k]` is `double`, `margin_order` is `vector<int>`, `errRp_k[k]` is `double`. Consistent throughout.
