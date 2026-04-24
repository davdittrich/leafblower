# Fix max_err_of() name-alignment bug + re-run WU-6 benchmark

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans or superpowers:subagent-driven-development.

**Goal:** Fix the positional-subtraction bug in `max_err_of()` and re-run the stepstone-fulldata benchmark to get correct external errRp values for all configs, enabling accurate assessment of the WU-6 merge gate.

**Architecture:** Single-file bug fix in `benchmarks/stepstone_fulldata_homotopy.R`. No source code changes. Re-run the existing benchmark script. Re-evaluate merge gates.

**Root cause (confirmed, commit d9dc5f9 diagnostic):**
`tapply(w, col, sum)` returns alphabetically-sorted names. `target[[k]]` preserves JSON order. R subtraction `tab - target[[k]]` is POSITIONAL, not name-aligned. For `rk_gender_time` (German compound category names), positional mismatch inflates error from true ~0.001 to phantom 0.218 for ALL weight vectors including autumn's converged solution.

**Fix:** `tab[names(target[[k]])] - target[[k]]` — index `tab` by target's names before subtracting.

**Reference data (commit d9dc5f9, internal solver errRp):**
- Baseline (overlays off): errRp = 2.223e-3
- Autumn reference (memory, commit 8146894): max_err ≤ 1.60e-3 (using autumn's own error=mean metric)
- AB (homotopy + greedy): errRp = 6.567e-3 (solver-reported, expected to be validated by corrected external metric)

**Tech Stack:** R, `benchmarks/stepstone_fulldata_homotopy.R`, existing parquet fixtures.

---

## Task 1: Fix max_err_of() + verify on known-good case

**Files:**
- Modify: `benchmarks/stepstone_fulldata_homotopy.R`
- No other files.

- [ ] **Step 1: Read current `max_err_of()` in the benchmark script.**

Current (wrong):
```r
max_err_of <- function(w, data, target) {
  errs <- numeric(length(target))
  for (k in seq_along(target)) {
    tab <- tapply(w, data[[names(target)[k]]], sum) / sum(w)
    errs[k] <- max(abs(tab - target[[k]]))
  }
  max(errs)
}
```

- [ ] **Step 2: Apply the name-alignment fix.**

```r
max_err_of <- function(w, data, target) {
  errs <- numeric(length(target))
  for (k in seq_along(target)) {
    tab <- tapply(w, data[[names(target)[k]]], sum) / sum(w)
    # Index tab by target's names to force name-aligned subtraction.
    # R's positional subtraction of differently-ordered named vectors
    # produces phantom errors; German compound category names in
    # rk_gender_time are ordered differently by tapply vs JSON.
    tab_aligned <- tab[names(target[[k]])]
    errs[k] <- max(abs(tab_aligned - target[[k]]), na.rm = TRUE)
  }
  max(errs)
}
```

Note `na.rm = TRUE` guards against categories present in target but absent from data column (would produce NA from `tab[name]` lookup).

- [ ] **Step 3: Verify fix on uniform weights.**

Add a diagnostic block ABOVE the benchmark runs (remove after verification):

```r
# Sanity check: uniform weights should show max_err = max departure from uniform target
w_uniform <- rep(1, nrow(data))
err_uniform <- max_err_of(w_uniform, data, target)
cat(sprintf("Uniform weight max_err: %.4f (should reflect true marginal deviation)\n",
            err_uniform))
```

Run: `Rscript -e 'source("benchmarks/stepstone_fulldata_homotopy.R")' 2>&1 | head -5`
Verify `err_uniform` is a plausible departure from uniform distribution (e.g., 0.1–0.5), NOT the phantom 0.218 seen before.

Remove the diagnostic block before committing.

---

## Task 2: Re-run benchmark + assess merge gates

**Files:**
- Read: `benchmarks/stepstone_fulldata_homotopy_report.rds` (re-saved after re-run)
- Create: `benchmarks/stepstone_fulldata_homotopy_report_v2.rds` if prior RDS needs preservation

**Note:** The benchmark takes ~3–5 min per config (5 configs + autumn = ~25 min total). Use `Rscript benchmarks/stepstone_fulldata_homotopy.R` and let it run to completion.

- [ ] **Step 1: Run the corrected benchmark.**

```bash
Rscript benchmarks/stepstone_fulldata_homotopy.R 2>&1 | tee benchmarks/stepstone_fulldata_homotopy_run_v2.log
```

Expected with corrected max_err_of():
- Autumn `max_err` should now be ≈ 1.60e-3 or lower (confirming fix works)
- Leafblower baseline (overlays off) `max_err` should be ≈ 2.22e-3 (matching internal errRp)
- AB config `max_err` tells us whether the gate is actually met

- [ ] **Step 2: Assess merge gate.**

Read the run output and check:
- Autumn `max_err_corrected` — confirms diagnostic is correct
- AB `max_err_corrected` vs threshold 1.60e-3
- If AB passes: update `test-bench-gate.R` to use `max_err` column (now correct) and commit with GATE PASS
- If AB fails: update `leafblower-b25a` HUMAN ticket with correct external errRp values and proceed to user decision

- [ ] **Step 3: Save updated report + commit.**

```bash
git add benchmarks/stepstone_fulldata_homotopy.R \
        benchmarks/stepstone_fulldata_homotopy_report.rds \
        benchmarks/stepstone_fulldata_homotopy_run_v2.log
git commit -m "$(cat <<'EOF'
fix(bench): correct max_err_of() positional-subtraction bug + re-run WU-6 gate

tapply() alphabetizes category names; R subtraction aligns by position not
name, inflating error for rk_gender_time (German compound names reordered).
Fix: tab[names(target[[k]])] forces name-aligned lookup. All prior max_err
values (0.163–0.168) were phantom errors; corrected values pending re-run.
EOF
)"
```

---

## Self-review checklist

- `grep -n "tab - target" benchmarks/stepstone_fulldata_homotopy.R` → zero hits (replaced by `tab_aligned - target[[k]]`)
- Autumn `max_err` after fix ≈ 1.60e-3 (confirms real convergence)
- Merge gate re-assessed with correct numbers
- No changes to `src/`, `R/`, or `tests/` — benchmark fix only

---

## Gate dependency

After this plan, EITHER:

**A. AB max_err ≤ 1.60e-3:** WU-6 passes. Proceed to WU-7 (Python parity + docs). Close `leafblower-b25a` HUMAN ticket.

**B. AB max_err > 1.60e-3:** Hypothesis partially falsified. Report correct errRp to user. `leafblower-b25a` HUMAN ticket stays open with updated numbers. User decides next step.
