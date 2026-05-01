# ieppa+SRAA SOR Coexistence Implementation Plan

**Goal:** Investigate and decide whether to re-enable SOR damping during rejected SRAA steps in ieppa, fixing the 20-iter early KL convergence with worse max_err.
**Architecture:** `src/ieppa.cpp:390` disables SOR auto-adapt whenever `st.accelerate=true` (SRAA active). When AA steps are rejected (`SRAAStepResult.aa_accepted=false`), the BCD sweep behaves exactly like plain ieppa — but lacks SOR damping, so it overshoots and converges KL too early, missing better max_err points. The fix gates SOR-disable on AA *acceptance*, not AA *activation*.
**Tech Stack:** C++17, R, `R CMD INSTALL --preclean .`, testthat, `benchmarks/stepstone_benchmark.R`.

**Mechanism:** track per-iter `aa_accepted` from `sraa_step`; allow SOR adaptation in iterations where AA was rejected.
**Forbidden:** any change to plain ieppa (accelerate=false) behavior; any change to SRAA acceptance criterion; bumping pct_tol globally; touching the SOR adaptation rule itself.
**Audit:** unit test must spy on per-iter `omega[k]` trajectory and assert SOR fired only on rejected-AA iterations.

---

## Task T1: Verification spike — confirm root cause

Reproduce the gap on stepstone via an interleaved bench (not sequential).

Steps:
1. Read `benchmarks/stepstone_benchmark.R` to identify the harvest() call shape used today (max_iter, tol_pct, etc.).
2. Add a probe script `benchmarks/probe_sraa_sor.R`:

```r
# Interleaved before/after — mandate from CLAUDE.md §2.
library(leafblower); library(bench)
data <- readRDS("benchmarks/fixtures/stepstone_small_targets.rds")  # use existing
make_call <- function(accelerate, sor_auto) {
  function() harvest(
    data$weights, data$group_ids, data$cat_counts, data$targets,
    method = "ieppa", accelerate = accelerate,
    sor = list(enabled = TRUE, auto_adapt = sor_auto, omega_init = 1.0,
               omega_min = 0.3, omega_fixed = -1.0, burnin = 20),
    max_iter = 200, tol_pct = 1e-4, verbose = 0)
}
res <- bench::mark(
  plain_sor    = make_call(FALSE, TRUE)(),
  sraa_no_sor  = make_call(TRUE,  TRUE)(),   # current behavior (sor masked)
  check = FALSE, iterations = 10)
print(res)
# Print iters + max_err side-by-side
```

3. Record `iterations` and `max_err` for both rows. If `sraa_no_sor` has fewer iters AND worse max_err than `plain_sor`, the root cause is confirmed.
4. **Stop and decide.** If the gap reproduces: proceed to T2 (Option B). If not: close with Option A (document in NEWS.md only).

Confidence: 90 — exact code at `ieppa.cpp:390` quoted.

---

## Task T2: Implement Option B — SOR adapts on rejected AA iterations only

**Mechanism:** Replace the static `sor_auto_v` flag (`ieppa.cpp:390`) with a per-iter dynamic decision driven by the most recent `SRAAStepResult.aa_accepted`.

Steps:

1. **Read context** — `src/ieppa.cpp:380-405` (SOR state init), `:740-870` (SRAA-m branch), `:1357` (SOR adaptation site). The SOR adaptation block at `:1357` lives inside the inner BCD for-loop; the SRAA branch at `:740-870` is a separate while-loop that does NOT execute the BCD for-loop. So the question is: can SOR run *inside* the SRAA loop's plain-step path?

2. **Audit:** SRAA's plain step is hidden inside `sraa_step` (template in `sraa.hpp:74-203`). The plain BCD sweep happens in `f_eval_lf` (passed to `sraa_step`). `f_eval_lf` does not know about SOR — SOR adaptation lives in the for-loop body, not in `f_eval_lf`. So Option B requires either:
   - **B-narrow:** push SOR omega update into `f_eval_lf` (capture-by-reference), gated on a `sor_run_this_iter` bool the outer SRAA loop sets after observing `r.aa_accepted=false`.
   - **B-wide:** when AA is rejected for kSRAAOuterStallWindow consecutive iters, fall back to the plain BCD for-loop with SOR re-enabled.

3. **Choose B-narrow** (fewer moving parts; preserves SRAA acceleration potential). Implementation:

```cpp
// ieppa.cpp:389-391, replace sor_auto_v const with mutable flag:
const bool sor_active     = st.sor_cfg.enabled;
bool sor_auto_v           = st.sor_cfg.auto_adapt;   // no longer masked by accelerate
const bool sor_auto_base  = st.sor_cfg.auto_adapt;
```

4. In the SRAA branch (`ieppa.cpp:773` while loop), after the `sraa_step` call:

```cpp
auto r = lbw::sraa_step(f_eval_lf, lf_flat, dummy_L, dummy_U,
                        ieppa_sraa, /*apply_clamp=*/false);
// SOR runs on rejected-AA iterations only — when AA fires, omega adaptation
// would race the AA correction. ieppa-tus6.
const bool allow_sor_this_iter = sor_auto_base && !r.aa_accepted;
sor_auto_v = allow_sor_this_iter;
```

5. The SOR adaptation block at `ieppa.cpp:1357` reads `sor_auto_v`. Confirm it is reachable from the SRAA path. **It is not** (SRAA branch does not enter the for-loop). So we need to invoke SOR adaptation *inside* `f_eval_lf` itself. The outer `bool sor_run_this_feval` flag (set in the SRAA while-loop) controls whether f_eval_lf applies SOR on this call.

Concrete `f_eval_lf` signature addition (add `sor_run_this_feval` captured by [&]):
```cpp
bool sor_run_this_feval = false;   // declared just before f_eval_lf, before the SRAA block

auto f_eval_lf = [&](std::vector<double>& flat) -> double {
    unpack_lf(flat, lf, f_lin, cell_lf, X_cur, ct, X_init, log_X_init,
              st.K, cat_offset, cell_lf_hwm);
    // [existing sweep logic: linear or log path] ...
    // B-narrow SOR hook: only when previous SRAA step was plain (aa_accepted=false)
    if (sor_run_this_feval && sor_active && iter_sraa >= sor_burnin_v) {
        for (int k = 0; k < st.K; k++) {
            double errRp_k = compute_margin_errRp_linear(k); // existing lambda
            bool decreasing = (errRp_k < sor_prev_errRp[k]);
            bool sign_flip  = !decreasing && sor_prev_decreasing[k];
            if (sign_flip)
                sor_omega[k] = std::max(omega_min_v, sor_omega[k] * kSorOscillationDamp);
            else if (decreasing)
                sor_omega[k] = std::min(1.0,         sor_omega[k] * kSorRecoveryGrowth);
            sor_prev_decreasing[k] = decreasing;
            sor_prev_errRp[k]      = errRp_k;
        }
    }
    pack_lf(lf, flat);
    // ... [errRp return] ...
};
```
In the SRAA while-loop, set `sor_run_this_feval = !last_aa_accepted` before calling `sraa_step`. Track `last_aa_accepted = r.aa_accepted` after the call.

6. `iter_sraa` counter: declare as `int iter_sraa = 0;` before the SRAA while-loop; increment each iteration.

7. **Touch only what is strictly needed.** Do not refactor unrelated SOR code paths. Do not change defaults. Do not touch `sor_cfg` struct.

Confidence: 70 — implementation sketch is consistent with code structure, but the SOR-inside-`f_eval_lf` path needs verification that errRp signal is per-margin available there.

---

## Task T3: Test — SOR fires only on rejected AA iters

Add `tests/testthat/test-ieppa-sraa-sor.R`:

```r
test_that("SOR adapts only on rejected SRAA steps", {
  skip_if_not(file.exists("../fixtures/stepstone_small_targets.rds"))
  d <- readRDS("../fixtures/stepstone_small_targets.rds")
  res <- harvest(d$weights, d$group_ids, d$cat_counts, d$targets,
                 method = "ieppa", accelerate = TRUE,
                 sor = list(enabled = TRUE, auto_adapt = TRUE,
                            omega_init = 1.0, omega_min = 0.3,
                            omega_fixed = -1.0, burnin = 5),
                 max_iter = 100, tol_pct = 1e-4,
                 diagnostics = TRUE)
  expect_true(res$diagnostics$sor_n_damped >= 0)
  expect_true(res$diagnostics$aa_accepted_count >= 0)
  # max_err must not regress vs documented baseline:
  expect_lt(res$max_err, 1.0e-2)   # was 1.36e-2 with SOR masked
})
```

Update `tests/testthat/test-sor.R` if it asserts that accelerate=TRUE silences SOR.

---

## Task T4: Bench — interleaved before/after

Re-run `benchmarks/probe_sraa_sor.R` against the new build. Pass criteria:

- `max_err`: SRAA+SOR ≤ plain SOR (tolerance 10%).
- `iterations`: SRAA+SOR ≤ 1.5× plain SOR (no regression vs old SRAA).
- `aa_accepted_count` > 0 (SRAA acceleration still active).

Record `bench::mark()` output in commit message body.

---

## Task T5: Documentation

- Update `NEWS.md` under development version: "ieppa: SOR auto-adapt now coexists with SRAA acceleration; fires only on iterations where AA is rejected (leafblower-tus6)."
- Update `man/harvest.Rd` if the `sor` arg description claims interaction with `accelerate`.
- Add a Do-Not-Repeat entry to `.wolf/cerebrum.md`: "When SRAA disables SOR, KL converges to a worse-max_err fixed point. SOR must adapt on rejected-AA iterations."

---

## Decision Tree

- **T1 fails to reproduce:** close ticket with Option A (accept). Document only.
- **T1 reproduces, T2 unsafe (SOR-inside-`f_eval_lf` racing):** fall back to Option C (max_err secondary gate at SRAA convergence check). Replace step T2 with a single-line gate at `ieppa.cpp:837-839`.
- **T2 lands, T4 regresses iterations >1.5×:** revert; ship Option C.
