# Move Weight Normalization from Wrappers to Solvers

**Goal:** Enforce strict bounds on the final user-visible weight vector by moving sum-to-n normalization inside each solver, so `bounds_mode="unit"` water-filling sees post-normalize weights and can clamp them correctly.

**User directive (2026-04-24):** "normalization should never be in the wrapper."

**Current bug:** `bounds_mode="unit"` doesn't enforce per-obs bounds. `R/harvest.R:118` runs `weights <- weights / mean(weights)` AFTER the C solver returns. ieppa's solver-level water-fill (per plan §6.2) clamps to [min, max] BEFORE the wrapper normalize — normalize can push weights back over the cap. On stepstone: 36,745 weights > 5.0 in both cell- AND unit-mode (identical output).

## Scope

One atomic commit. Files:
- `R/harvest.R` — remove `weights / mean(weights)` line
- `python/leafblower/_harvest.py` — remove `weights_out / w_mean` equivalent
- `src/ieppa.cpp` — normalize + water-fill ordering change: (a) compute mult, (b) expand, (c) normalize to sum=n, (d) if unit-mode, water-fill with sum-preserving redistribution
- `src/raking.cpp` — verify sum=n post-finalizer (comment claims no-op at convergence); add explicit normalize for NOCONV case so caller contract holds
- `src/lbfgsb_solver.cpp` — same verification; add explicit normalize for NOCONV
- `tests/testthat/fixtures/stepstone_reference.rds` — note: REFERENCE WILL CHANGE if current behavior was wrapper-normalized. Must recapture.
- Tests that assert `mean(weights) == 1` — update to use solver-internal invariants

## Pre-flight

- [ ] **Step P.1:** `git status --short` shows only `.wolf/anatomy.md` dirty (tool artifact).
- [ ] **Step P.2:** Baseline test count: 205 pass, 0 fail.
- [ ] **Step P.3:** Capture CURRENT solver output (pre-normalize-by-wrapper) on a tiny test input via verbose, so we know what each solver returns natively.

## Task 1: ieppa.cpp — normalize + water-fill

- [ ] **Step 1.1:** Current end of `ieppa_solve` computes `mult[c] = X[c]/X_init[c]`, then `st.weights[i] = st.weights[i] * mult[ct.cell_of[i]]`. After this, add:

```cpp
// Normalize to sum = n (match wrapper contract, moved in-solver per user directive 2026-04-24).
// Applied BEFORE bounds_mode-specific post-processing so unit-mode water-fill sees
// final weights and can guarantee [min_weight, max_weight] strictly.
double total_w = 0.0;
for (int i = 0; i < st.n; i++) total_w += st.weights[i];
if (total_w > 0.0) {
    const double norm = static_cast<double>(st.n) / total_w;
    for (int i = 0; i < st.n; i++) st.weights[i] *= norm;
    // X[c] scales too, so n_bounds_clamped / water-fill use consistent cell totals:
    for (int c = 0; c < ct.M_cell; c++) X[c] *= norm;
}
```

- [ ] **Step 1.2:** Move the existing cell-mode violation scan AND unit-mode water-fill block to AFTER the normalize. Per-cell bounds for water-fill are now `[min_weight, max_weight]` directly on `st.weights[i]` — the cell-aggregate cap at `X[c]` already reflects the post-normalize scale.
- [ ] **Step 1.3:** Build gate: `R CMD INSTALL --preclean .` clean.

## Task 2: raking.cpp

- [ ] **Step 2.1:** Read current `raking_solve` exit. Comment at line 250 claims sum(w)=n at convergence. For NOCONV exits (max_iter hit), sum may drift. Append explicit normalize at solver exit:

```cpp
// Solver-owned normalization (moved from wrapper 2026-04-24).
double total_w = 0.0;
for (int i = 0; i < st.n; i++) total_w += st.weights[i];
if (total_w > 0.0) {
    const double norm = static_cast<double>(st.n) / total_w;
    for (int i = 0; i < st.n; i++) st.weights[i] *= norm;
}
```

Place after the Dykstra finalizer, before any final diagnostic logging.

- [ ] **Step 2.2:** Build gate.

## Task 3: lbfgsb_solver.cpp

- [ ] **Step 3.1:** Same pattern as raking — append explicit normalize at solver exit (handles NOCONV cases where the invariant drifts).
- [ ] **Step 3.2:** Build gate.

## Task 4: R wrapper

- [ ] **Step 4.1:** `R/harvest.R:117-118` — delete:

```r
# Normalize to mean=1 (preserves calibration constraints which are proportional).
weights <- weights / mean(weights)
```

Replace with comment:
```r
# Solver returns sum(weights) = n. Normalization is solver-internal
# (see src/ieppa.cpp, src/raking.cpp, src/lbfgsb_solver.cpp).
# No-op here preserves bounds_mode="unit" strict-bounds guarantee.
```

- [ ] **Step 4.2:** `R/harvest.R` — remove `enforce_mean` parameter entry from signature (or keep as deprecated/ignored; existing TODO comment at line 28 says "Ignored (retained for compatibility)"). Leave as-is.

## Task 5: Python wrapper

- [ ] **Step 5.1:** `python/leafblower/_harvest.py:189-193` — delete the `w_mean = weights_out.mean(); if w_mean > 0: weights_out = weights_out / w_mean` block. Replace with comment matching R wrapper.

## Task 6: Test fixture regeneration

- [ ] **Step 6.1:** Existing tests may assert `mean(weights) ≈ 1` post-harvest. Grep for patterns: `grep -rn "mean(.*weight" tests/`. For each hit: tighten assertion to `sum(weights) ≈ n` (mathematically equivalent under the new contract) OR loosen bound to allow solver-internal normalization precision.
- [ ] **Step 6.2:** `tests/testthat/fixtures/stepstone_reference.rds` — stored leafblower weights WILL CHANGE (they were post-wrapper-normalize; new weights are post-solver-normalize which is the same VALUE but via a different path). Worth regenerating for cleanliness. Autumn weights unchanged (autumn package has its own normalize, untouched).
- [ ] **Step 6.3:** Full suite: `[FAIL 0 | PASS ≥ 205]`. Any regression: halt + diagnose. Expected new-behavior delta: weights where water-fill now fires may differ slightly from pre-change.

## Task 7: Stepstone re-verification

- [ ] **Step 7.1:** Run `/tmp/stepstone_compare_current.R` — assertions to verify:
  - cell-mode: `max(w) ≤ 5.0 + 1e-9` — still may be violated without water-fill on cell-mode (the violation count SHOULD now be non-zero if violations exist)
  - unit-mode: `max(w) ≤ 5.0 + 1e-9` — MUST hold; `n_bounds_clamped > 0` expected on stepstone (36,745 obs currently over 5.0 → water-fill fires)
  - Stepstone errRp preserved ≤ 2.22e-3

## Task 8: Commit

- [ ] **Step 8.1:**

```bash
git add R/harvest.R python/leafblower/_harvest.py \
        src/ieppa.cpp src/raking.cpp src/lbfgsb_solver.cpp
git commit -m "$(cat <<'EOF'
refactor(normalize): move weight normalization from wrappers into solvers

User directive (2026-04-24): "normalization should never be in the wrapper."

Each solver (ieppa, raking, lbfgsb) now owns its sum-to-n normalization,
applied at exit before any bounds_mode post-processing. This fixes a
silent bounds_mode='unit' bug where R/harvest.R's post-solver
`weights <- weights / mean(weights)` re-pushed weights over max_weight
AFTER the solver's water-fill had clamped them.

Post-change: bounds_mode='unit' delivers strict per-observation bounds
as documented. cell-mode unchanged (existing violations still surface via
result$n_bounds_violated + warning).

Stepstone cell-mode: 36,745 weights > 5.0 → surfaces as warning + counter
(previously silent); unit-mode: 0 weights > 5.0 after water-fill.
EOF
)"
```

## Self-review

1. Moves normalization from 2 wrappers + 0 solvers → 0 wrappers + 3 solvers. Net complexity up by 1 fn (3 additions - 2 deletions).
2. User contract `sum(weights) = n` preserved.
3. `bounds_mode="unit"` gains teeth (unit-mode water-fill now effective).
4. `bounds_mode="cell"` behavior: weights identical BEFORE wrapper normalize; after normalize, identical. No user-visible change.
5. Backward-compat: any caller relying on `mean(weights) == 1` still sees exactly that (algebraically equivalent to `sum == n`).
