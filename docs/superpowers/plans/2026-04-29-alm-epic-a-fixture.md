# ALM ieppa_soft — Epic A: Fixture Capture

> **For agentic workers:** Use `superpowers:subagent-driven-development` to implement this plan.

**Goal:** Capture pre-ALM ieppa reference fixture before any code change.
**Mechanism:** `Rscript data-raw/gen_ieppa_pre_alm_ref.R`
**Forbidden:** Any modification to `src/` before this task completes
**Audit:** `grep -c 'lambda_cell\|capacity_mu' src/ieppa.cpp` must return 0

---

### Task 0: Capture pre-ALM ieppa reference fixture

**Branch:** `fix/correctness-performance-2026-04-28`
**Files:** `data-raw/gen_ieppa_pre_alm_ref.R` (create), `tests/testthat/fixtures/ieppa_pre_alm_ref.rds` (generate)

#### Step 0 — Prereq guard (STOP if non-zero)

```bash
grep -c 'lambda_cell\|capacity_mu' src/ieppa.cpp
# Expected: 0
# If non-zero: ALM already landed — fixture cannot be captured now. HALT.
```

#### Step 1 — Create generator script

Create `data-raw/gen_ieppa_pre_alm_ref.R`:

```r
library(leafblower)
set.seed(3); n <- 5000L
df  <- data.frame(v1 = factor(sample(5, n, TRUE)))
tgt <- list(v1 = setNames(c(0.4, 0.3, 0.15, 0.1, 0.05), as.character(1:5)))
r   <- harvest(df, tgt, method = "ieppa",
               max_weight = 1.8, min_weight = 0, max_iterations = 500,
               convergence = list(improvement = 1e-4),
               attach_weights = FALSE)
saveRDS(list(df           = df,
             tgt          = tgt,
             max_weight   = 1.8,
             min_weight   = 0,
             max_iterations = 500L,
             convergence  = list(improvement = 1e-4),
             weights      = as.numeric(r),
             result       = attr(r, "result")),
        "tests/testthat/fixtures/ieppa_pre_alm_ref.rds")
cat("Fixture written. n_weights:", length(as.numeric(r)),
    "max_error:", attr(r, "result")$max_error, "\n")
```

#### Step 2 — Run generator

```bash
Rscript data-raw/gen_ieppa_pre_alm_ref.R
# Expected output: "Fixture written. n_weights: 5000 max_error: <value>"
```

#### Step 3 — Verify fixture integrity

```bash
Rscript -e "f <- readRDS('tests/testthat/fixtures/ieppa_pre_alm_ref.rds'); stopifnot(length(f\$weights) == 5000); cat('OK\n')"
# Expected: OK
```

#### Step 4 — Confirm test suite baseline unchanged

```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
# Expected: same FAIL count as before this task (should be 3)
```

#### Step 5 — Commit both files

```bash
git add data-raw/gen_ieppa_pre_alm_ref.R tests/testthat/fixtures/ieppa_pre_alm_ref.rds
git commit -m "$(cat <<'EOF'
test(ieppa): capture pre-ALM reference fixture for backward compat regression

Canonical ieppa weights on set.seed(3) 5-category tight-bounds problem
(max_weight=1.8) captured before any ALM code lands. Used by T9 backward
compat test to assert method='ieppa' produces bit-identical weights pre/post
ALM merge. Step 0 prerequisite per spec §Implementation Order.
EOF
)"
```

#### Success criteria

| Check | Expected |
|---|---|
| `grep -c 'lambda_cell\|capacity_mu' src/ieppa.cpp` | 0 |
| `length(f$weights)` | 5000 |
| `Rscript -e "…stopifnot…"` | `OK` |
| `devtools::test()` FAIL count | unchanged (3) |
| `git log --oneline -1` | commit message starts with `test(ieppa):` |
