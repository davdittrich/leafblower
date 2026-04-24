# Move Weight Normalization from Wrappers to Solvers (rev 3)

**Goal:** Enforce strict bounds on the final user-visible weight vector by moving sum-to-n normalization inside each solver, so `bounds_mode="unit"` water-filling sees post-normalize weights and can clamp them correctly.

**User directive (2026-04-24):** "normalization should never be in the wrapper."

**Current bug:** `bounds_mode="unit"` doesn't enforce per-obs bounds. `R/harvest.R:118` runs `weights <- weights / mean(weights)` AFTER the C solver returns. ieppa's solver-level water-fill clamps to [min, max] BEFORE the wrapper normalize — normalize can push weights back over the cap. On stepstone: 36,745 weights > 5.0 in both cell- AND unit-mode (identical output).

## Rev history
- **rev 2 → rev 3** (plan-review-gate iter 2):
  - Feasibility: expanded Step 3.2 to explicitly state max_err scale-invariance rationale when updating lbfgsb comment
  - Completeness: Task 6 now says "add contract docstring" (c_api.cpp has no existing docstring at rk_calibrate)
  - Completeness: Step P.3 tightened — r_bridge.cpp:150 comment documents INPUT-side (start_weights pre-normalized via R/harvest.R:218); rev 2 changes OUTPUT-side only, so that comment stays. Grep verified no output-side normalize in r_bridge.cpp or c_api.cpp.
  - Scope: Step 7.3 demoted to optional test-hygiene; skip if numerically bitwise identical to existing reference

- **rev 1 → rev 2** (plan-review-gate iter 1):
  - Feasibility: removed `X[c] *= norm` in ieppa (dead code, false justification)
  - Feasibility: corrected raking rationale (hyperplane finalizer already enforces sum=n unconditionally — normalize is defensive, not NOCONV fix)
  - Feasibility: made lbfgsb normalize placement explicit (wrap around `compute_final_weights_and_error`; max_err = S[j]/Wn is scale-invariant so either placement correct, but outer placement avoids ambiguity)
  - Completeness: explicit total_w == 0 contract (no-op, leaves weights unchanged; documented)
  - Completeness: widened grep pattern to `weight.*\b(mean|sum)\(|\b(mean|sum)\([^)]*weight`
  - Completeness: added explicit r_bridge.cpp audit step (Scope warn)
  - Completeness: stepstone compare script will be committed at `benchmarks/stepstone_compare_current.R` (benchmarks/ is already .Rbuildignored; tests/manual/ is NOT, so R CMD check would pick it up)
  - Scope: added C-API contract shift documentation step
  - Line-ref fix: Python block actually at lines 192-194 (plan cited 189-193)

## Scope

One atomic commit. Files:
- `R/harvest.R` — remove `weights / mean(weights)` line
- `python/leafblower/_harvest.py` — remove `weights_out / w_mean` equivalent
- `src/ieppa.cpp` — normalize AFTER expansion, BEFORE bounds scan/water-fill; operates on st.weights only
- `src/raking.cpp` — append defensive normalize at solver exit (current finalizer already enforces sum=n; normalize is a no-op at convergence)
- `src/lbfgsb_solver.cpp` — wrap `compute_final_weights_and_error` call with post-normalize (max_err is scale-invariant)
- `src/r_bridge.cpp` — audit for any post-solve normalization, remove if present
- `tests/testthat/` — update assertions asserting `mean(weights) == 1` → `sum(weights) ≈ n`
- `tests/testthat/fixtures/stepstone_reference.rds` — regenerate (weights mathematically identical but algebraically via different path)

## Pre-flight

- [ ] **Step P.1:** `git status --short` shows only `.wolf/anatomy.md` dirty (tool artifact).
- [ ] **Step P.2:** Baseline test count: 205 pass, 0 fail.
- [ ] **Step P.3:** Audit `src/r_bridge.cpp` and `src/c_api.cpp` for any OUTPUT-side post-solve normalization. Verified 2026-04-24: only hit is the INPUT-side comment at `r_bridge.cpp:150` describing `normalize_start_weights` (R layer pre-normalizes start_weights via `R/harvest.R:218` — this invariant is UNCHANGED by rev 3 and the comment stays). No output-side normalize exists in either file. If a future audit surfaces one, remove in the same commit.
- [ ] **Step P.4:** Copy `/tmp/stepstone_compare_current.R` → `benchmarks/stepstone_compare_current.R` (benchmarks/ is already `.Rbuildignore`d per current manifest).

## Task 1: ieppa.cpp — normalize between expansion and bounds handling

- [ ] **Step 1.1:** After the expansion loop at `src/ieppa.cpp:535-537` (writes `st.weights[i]`), BEFORE the `if (st.bounds_mode == RK_BOUNDS_CELL)` branch at line 539, insert:

```cpp
// Normalize to sum = n (moved in-solver per user directive 2026-04-24).
// Applied AFTER expansion and BEFORE bounds_mode post-processing so unit-mode
// water-fill sees final-scale weights and can guarantee [min_weight, max_weight] strictly.
// Contract: if total_w == 0 (degenerate: all X_init[c] == 0 or all mult[c] == 0),
// weights are left as-is (all zero). Solver status remains whatever upstream set it to.
// X[c] is NOT rescaled: it is dead after the expansion loop (water-fill uses
// (void)target_sum at line 561 and redistributes within each cell, so cell-aggregate
// scale is irrelevant to the bounds post-processing).
double total_w = 0.0;
for (int i = 0; i < st.n; i++) total_w += st.weights[i];
if (total_w > 0.0) {
    const double norm = static_cast<double>(st.n) / total_w;
    for (int i = 0; i < st.n; i++) st.weights[i] *= norm;
}
```

- [ ] **Step 1.2:** No change to cell-mode scan or unit-mode water-fill — they operate on `st.weights[i]` directly with absolute `min_weight`/`max_weight` bounds. Post-normalize scale is now the correct frame for those bounds.
- [ ] **Step 1.3:** Build gate: `R CMD INSTALL --preclean .` clean.

## Task 2: raking.cpp — defensive normalize at exit

- [ ] **Step 2.1:** The Dykstra hyperplane finalizer at `src/raking.cpp:261-277` already enforces `sum(w) = n` on every iteration unconditionally (NOT just at convergence). Adding a normalize here is **defensive** — strictly a no-op at both CONV and NOCONV. Rationale: solver contract moves from wrapper to solver, so every solver self-enforces even if the invariant is currently guaranteed by prior logic. After the Dykstra finalizer, before final diagnostic logging:

```cpp
// Solver-owned normalization (moved from wrapper per user directive 2026-04-24).
// Defensive: the preceding hyperplane finalizer already enforces sum(w) = n
// unconditionally; this is a belt-and-braces guard against future refactors that
// might weaken that invariant. total_w == 0 is pathological (caller passed all-zero
// design weights) — leave unchanged and let downstream surface the anomaly.
double total_w = 0.0;
for (int i = 0; i < st.n; i++) total_w += st.weights[i];
if (total_w > 0.0) {
    const double norm = static_cast<double>(st.n) / total_w;
    for (int i = 0; i < st.n; i++) st.weights[i] *= norm;
}
```

- [ ] **Step 2.2:** Build gate.

## Task 3: lbfgsb_solver.cpp — wrap the compute step

- [ ] **Step 3.1:** `lbfgsb_solve` at `src/lbfgsb_solver.cpp:564-582` returns directly from `lbfgsb_solve_inner` (which tail-calls `compute_final_weights_and_error`). No BADARG/early-exit branches exist at this layer (argument validation is upstream). Refactor to capture result and normalize before return:

```cpp
LBFGSResult lbfgsb_solve(CalibState& st) {
    auto off = build_offsets(st);
    std::vector<double> T(off[st.K]);
    double W_sum = compute_targets_abs(st, T);
    std::vector<double> d(st.n);
    for (int i = 0; i < st.n; i++) d[i] = st.weights[i];
    st.alm_lambda = 0.0;
    st.alm_mu     = 0.0;
    LBFGSResult res = lbfgsb_solve_inner(st, off, T, d, W_sum);

    // Solver-owned normalization (moved from wrapper per user directive 2026-04-24).
    // Placement after inner solve is safe because max_err (res.max_error) is computed
    // inside compute_final_weights_and_error as S[j]/Wn — scale-invariant under uniform
    // weight rescaling. Therefore normalizing st.weights here does NOT invalidate
    // res.max_error. Contract: total_w == 0 leaves weights untouched.
    double total_w = 0.0;
    for (int i = 0; i < st.n; i++) total_w += st.weights[i];
    if (total_w > 0.0) {
        const double norm = static_cast<double>(st.n) / total_w;
        for (int i = 0; i < st.n; i++) st.weights[i] *= norm;
    }
    return res;
}
```

- [ ] **Step 3.2:** Update the internal comments at `lbfgsb_solver.cpp:174-178` and `576-578` to reflect the contract shift. Replace the "Caller (harvest.R) normalises to mean=1 post-call" language with:

```cpp
// Solver self-normalizes at exit: after lbfgsb_solve_inner returns, we scale
// st.weights so Σ w_i = n. max_error is computed inside
// compute_final_weights_and_error as max_{k,j} |S_kj/Wn - target_kj|; both S_kj
// and Wn scale identically under uniform weight normalization, so
// res.max_error is scale-invariant and remains valid after the outer normalize.
// Caller contract (harvest.R, _harvest.py): solver returns sum(w) = n.
```
- [ ] **Step 3.3:** Build gate.

## Task 4: R wrapper

- [ ] **Step 4.1:** `R/harvest.R:117-118` — delete:

```r
# Normalize to mean=1 (preserves calibration constraints which are proportional).
weights <- weights / mean(weights)
```

Replace with comment:
```r
# Solver returns sum(weights) = n (enforced in src/ieppa.cpp, src/raking.cpp,
# src/lbfgsb_solver.cpp per user directive 2026-04-24). No wrapper-level
# normalization — removing it preserves bounds_mode="unit" strict-bounds guarantee.
```

- [ ] **Step 4.2:** Leave `enforce_mean` param entry as-is (TODO comment at line 28 already says "Ignored (retained for compatibility)"). Out-of-scope cleanup.

## Task 5: Python wrapper

- [ ] **Step 5.1:** `python/leafblower/_harvest.py:192-194` — delete the `w_mean = weights_out.mean(); if w_mean > 0: weights_out = weights_out / w_mean` block. Replace with comment matching R wrapper.

## Task 6: C-API contract documentation

- [ ] **Step 6.1:** `src/c_api.cpp` has no existing docstring at `rk_calibrate` (verified: line 138 declaration has no preceding comment). Add a contract docstring immediately above the function declaration:

```cpp
// Return contract: on RK_OK or RK_ERR_NOCONV, output weight vector
// satisfies Σ weights[i] = n (where n = number of observations). Solvers
// enforce this invariant internally (see ieppa.cpp, raking.cpp,
// lbfgsb_solver.cpp post-exit normalize blocks). Third-party callers
// should NOT apply their own sum/mean normalization — doing so will
// silently invalidate bounds_mode="unit" strict-bounds guarantees.
```

## Task 7: Test fixture regeneration

- [ ] **Step 7.1:** Grep for mean/sum assertions on weight vectors with a broad pattern:
  ```bash
  grep -rnE 'weight.*\b(mean|sum)\(|\b(mean|sum)\([^)]*weight' tests/
  ```
  For each hit: if the assertion is `mean(weights) ≈ 1`, tighten to `sum(weights) ≈ n` (mathematically equivalent under the new contract; no looser tolerance needed since solvers compute the normalization exactly rather than averaging floating-point-noisy outputs).
- [ ] **Step 7.2:** Known hits from prior audit: `tests/testthat/test-lbfgsb.R:60,73,85` assert `mean(res$weights) == 1.0`. Convert to `expect_equal(sum(res$weights), st$n, tolerance = 1e-10)`.
- [ ] **Step 7.3 (optional):** `tests/testthat/fixtures/stepstone_reference.rds` — compute cell-mode weights with new code and compare bitwise to existing fixture. If `max(abs(w_new - ref$lb_weights)) == 0`, SKIP regeneration (pure test hygiene, adds binary diff noise). If nonzero delta exists (should not — normalize sequence is identical), regenerate and document the delta source. Autumn weights unchanged (autumn has its own normalize path, untouched).
- [ ] **Step 7.4:** Full suite: `[FAIL 0 | PASS ≥ 205]`. Any regression: halt + diagnose.

## Task 8: Stepstone re-verification

- [ ] **Step 8.1:** Run `benchmarks/stepstone_compare_current.R` — assertions to verify:
  - cell-mode: `max(w) ≤ 5.0 + 1e-9` — still may be violated without water-fill on cell-mode (violation counter should now surface the 36,745)
  - unit-mode: `max(w) ≤ 5.0 + 1e-9` — MUST hold; `n_bounds_clamped > 0` expected (~36,745 obs over 5.0 pre-normalize → water-fill fires)
  - Stepstone errRp preserved ≤ 2.22e-3 (cell-mode byte-identical to reference; unit-mode ≈ reference within clamp redistribution)

## Task 9: Commit

- [ ] **Step 9.1:**

```bash
git add R/harvest.R python/leafblower/_harvest.py \
        src/ieppa.cpp src/raking.cpp src/lbfgsb_solver.cpp \
        src/c_api.cpp src/r_bridge.cpp \
        tests/testthat/test-lbfgsb.R \
        benchmarks/stepstone_compare_current.R
# Only include tests/testthat/fixtures/stepstone_reference.rds if Step 7.3 found a nonzero delta.
git commit -m "$(cat <<'EOF'
refactor(normalize): move weight normalization from wrappers into solvers

User directive (2026-04-24): "normalization should never be in the wrapper."

Each solver (ieppa, raking, lbfgsb) now owns its sum-to-n normalization at
exit. In ieppa, normalize runs AFTER expansion and BEFORE bounds_mode
post-processing, so unit-mode water-fill now sees final-scale weights and
can strictly enforce [min_weight, max_weight]. Previously, R/harvest.R
post-solver `weights / mean(weights)` re-pushed weights over max_weight
AFTER water-fill had clamped them — silently invalidating the strict-bounds
guarantee.

raking finalizer already enforced sum=n unconditionally via the hyperplane
step; the added normalize is defensive (no-op at both CONV and NOCONV),
kept so each solver is self-contained against future refactors.

lbfgsb: normalize is placed after compute_final_weights_and_error because
max_error = S[j]/Wn is scale-invariant, so moving the normalization does
not invalidate the reported error.

C-API contract shift: output weights now satisfy sum(w) = n from the
solver (formerly enforced in wrappers). Documented in src/c_api.cpp.

Stepstone verification:
- cell-mode: 36,745 weights > 5.0 surfaces as n_bounds_violated + warning
- unit-mode: 0 weights > 5.0 after water-fill (previously also 36,745)
EOF
)"
```

## Self-review

1. Net complexity: 3 solver-local normalize blocks added, 2 wrapper normalize lines removed. Small net positive, but each solver now self-contains its contract — no cross-file invariant.
2. User contract `sum(weights) = n` preserved end-to-end; `bounds_mode="unit"` gains teeth.
3. `bounds_mode="cell"` numerically byte-identical (normalize under cell-mode ran in wrapper before; now runs in solver — same sequence of operations, same result).
4. Backward-compat: any caller relying on `mean(weights) == 1` still sees exactly that (algebraically equivalent to `sum == n`).
5. Reviewer-addressed fixes folded in — rev 2 addresses all Feasibility blockers, all Completeness critical/moderate gaps, both Scope warns.
