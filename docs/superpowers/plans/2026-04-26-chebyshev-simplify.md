# Chebyshev/Grake Simplification Fixes

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix three copy-paste and duplication issues found by /critical-code-reviewer in the Plan E chebyshev+grake implementation.

**Architecture:** Two C++ fixes (c_api.cpp, r_bridge.cpp), one in-file constant dedup (chebyshev.cpp). No new files. No behavior change — purely structural.

**Tickets:** leafblower-t3oj, leafblower-bryp, leafblower-ri13

**Baseline:** FAIL 0 | PASS 374 | SKIP 5

---

## File Structure

| File | Task | Ticket |
|---|---|---|
| `src/c_api.cpp` | 1 | leafblower-t3oj |
| `src/r_bridge.cpp` | 2 | leafblower-bryp |
| `src/chebyshev.cpp` | 3 | leafblower-ri13 |

---

## Task 1 — c_api.cpp: extract CHEBYSHEV/GRAKE 17-field pack (leafblower-t3oj)

**Problem:** `RK_ALG_CHEBYSHEV` and `RK_ALG_GRAKE` dispatch blocks are identical 17-field copy-pastes from `ChebyshevResult` → `rk_result_t`. Any new field added to one silently omits from the other.

### Step 1.1: Read the two blocks

```bash
grep -n "alg == RK_ALG_CHEBYSHEV\|alg == RK_ALG_GRAKE\|gres2\|cres\b" src/c_api.cpp | head -10
```

### Step 1.2: Replace both blocks with a shared lambda

Find the first block (GRAKE, using `gres`) and the second (CHEBYSHEV, using `cres`) and third (GRAKE again, using `gres2`). Read exact lines:

```bash
grep -n "alg == RK_ALG_CHEBYSHEV\|alg == RK_ALG_GRAKE\|return gres\.\|return cres\.\|return gres2\." src/c_api.cpp | head -10
```

Replace ALL three with a single lambda + two call sites:

```cpp
    // Shared packer for ChebyshevResult (used by CHEBYSHEV and GRAKE).
    auto pack_cheb = [&](const lbw::ChebyshevResult& r, rk_algorithm_t a) -> int {
        if (result) {
            result->status                       = r.status;
            result->iterations                   = r.iterations;
            result->max_error                    = r.max_error;
            result->convergence_metric           = r.convergence_metric;
            result->convergence_rule             = r.convergence_rule;
            result->convergence_tol              = r.convergence_tol;
            result->convergence_iter             = r.convergence_iter;
            result->convergence_objective        = r.convergence_objective;
            result->convergence_minimized_metric = r.convergence_minimized_metric;
            result->best_error                   = r.best_error;
            result->best_iter                    = r.best_iter;
            result->mean_error                   = r.mean_error;
            result->kl                           = r.kl;
            result->chi2                         = r.chi2;
            result->grake_norm                   = r.grake_norm;
            result->l1_weight_change             = r.l1_weight_change;
            result->algorithm_used               = static_cast<int>(a);
            std::strncpy(result->message, r.message, sizeof(result->message) - 1);
        }
        return r.status;
    };
    if (alg == RK_ALG_CHEBYSHEV) {
        return pack_cheb(lbw::chebyshev_ipm(st, lbw::LpVariant::CHEBYSHEV), alg);
    } else if (alg == RK_ALG_GRAKE) {
        return pack_cheb(lbw::chebyshev_ipm(st, lbw::LpVariant::GRAKE), alg);
    } else {
```

**IMPORTANT:** Place the lambda BEFORE the `if (alg == RK_ALG_CHEBYSHEV)` block. The lambda must be in scope for both calls. Read the exact surrounding structure before editing.

### Step 1.3: Build gate
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```

### Step 1.4: Regression
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 374.

### Step 1.5: Commit
```bash
git add src/c_api.cpp
git commit -m "refactor(c_api): extract CHEBYSHEV/GRAKE 17-field pack into shared lambda

Both dispatch branches were identical copy-pastes. Single lambda
pack_cheb() eliminates sync hazard on future ChebyshevResult field additions."
bd close leafblower-t3oj 2>/dev/null || true
```

---

## Task 2 — r_bridge.cpp: deduplicate chebyshev+grake dispatch (leafblower-bryp)

**Problem:** The chebyshev and grake dispatch branches in `r_bridge.cpp` are identical except for `res_alg_used`. Extract to a shared lambda.

### Step 2.1: Read the two blocks

Lines ~458-490 in r_bridge.cpp. Both have:
```
} else if (strcmp(method_str, "chebyshev/grake") == 0) {
    auto res = lbw::chebyshev_ipm(st, LpVariant::CHEBYSHEV/GRAKE);
    pack_solver_result(res);
    res_status = res.status; res_iterations = ...; res_max_error = ...;
    res_alg_used = RK_ALG_CHEBYSHEV/GRAKE;
    res_mean_error = ...; res_kl = ...; res_chi2 = ...;
    res_best_error = ...; res_best_iter = ...;
    if (!res.best_weights.empty()) res_best_weights = ...; else assign(0);
```

### Step 2.2: Extract shared lambda before the chebyshev branch

**Place the lambda BEFORE the `} else if (strcmp(method_str, "chebyshev") == 0)` line:**

```cpp
    // Shared dispatch for both chebyshev and grake (same solver, different variant).
    auto dispatch_cheb = [&](lbw::LpVariant variant, int alg_code) {
        auto res = lbw::chebyshev_ipm(st, variant);
        pack_solver_result(res);
        res_status     = res.status;
        res_iterations = res.iterations;
        res_max_error  = res.max_error;
        res_alg_used   = alg_code;
        res_mean_error = res.mean_error;
        res_kl         = res.kl;
        res_chi2       = res.chi2;
        res_best_error = res.best_error;
        res_best_iter  = res.best_iter;
        if (!res.best_weights.empty())
            res_best_weights = std::move(res.best_weights);
        else
            res_best_weights.assign(st.n, 0.0);
    };
    } else if (strcmp(method_str, "chebyshev") == 0) {
        dispatch_cheb(lbw::LpVariant::CHEBYSHEV, static_cast<int>(RK_ALG_CHEBYSHEV));
    } else if (strcmp(method_str, "grake") == 0) {
        dispatch_cheb(lbw::LpVariant::GRAKE, static_cast<int>(RK_ALG_GRAKE));
    } else {
```

**Note:** The lambda captures `[&]` — it reads `st`, `pack_solver_result`, `res_status`, `res_iterations`, etc. which are all local variables in scope. The `std::move(res.best_weights)` is correct since `res` is local to the lambda.

### Step 2.3: Build + regression
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
Rscript -e 'devtools::test()' 2>&1 | tail -3
```

### Step 2.4: Commit
```bash
git add src/r_bridge.cpp
git commit -m "refactor(r_bridge): deduplicate chebyshev+grake dispatch branches

Lambda dispatch_cheb() replaces two identical copy-paste blocks.
Only alg_code differs between the two callers."
bd close leafblower-bryp 2>/dev/null || true
```

---

## Task 3 — chebyshev.cpp: hoist kMetricEps/kChi2Floor (leafblower-ri13)

**Problem:** Two pairs of identical constants in the same function:
- Line 167 (in-loop): `constexpr double kMetricEps2 = 1e-10, kChi2Floor2 = 1.0;`
- Line 407 (post-loop): `constexpr double kMetricEps = 1e-10, kChi2Floor = 1.0;`

Same values, different names. Hoist to function-level.

### Step 3.1: Add static constexpr at function entry

Find the constants at the top of `sinkhorn_solve`... wait, these are in `chebyshev_ipm`. Find the top of `chebyshev_ipm`:

```bash
grep -n "ChebyshevResult chebyshev_ipm\|kMaxIpm\|static constexpr" src/chebyshev.cpp | head -10
```

Add after existing `static constexpr` block (after `kStepScale`):
```cpp
static constexpr double kMetricEps = 1e-10;
static constexpr double kChi2Floor = 1.0;
```

### Step 3.2: Replace in-loop usage

Find `kMetricEps2` and `kChi2Floor2` and replace with `kMetricEps` and `kChi2Floor`:
```bash
grep -n "kMetricEps2\|kChi2Floor2" src/chebyshev.cpp
```

Delete the in-loop `constexpr double kMetricEps2 = 1e-10, kChi2Floor2 = 1.0;` line and replace `kMetricEps2` → `kMetricEps`, `kChi2Floor2` → `kChi2Floor` in the in-loop block.

Also delete the post-loop `constexpr double kMetricEps = 1e-10, kChi2Floor = 1.0;` line (now redundant since defined at function level).

### Step 3.3: Build + regression
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
Rscript -e 'devtools::test()' 2>&1 | tail -3
```

### Step 3.4: Commit
```bash
git add src/chebyshev.cpp
git commit -m "refactor(chebyshev): hoist kMetricEps/kChi2Floor to function-level

Was defined twice with different names (kMetricEps2/kChi2Floor2 in loop,
kMetricEps/kChi2Floor post-loop). Single function-level definition."
bd close leafblower-ri13 2>/dev/null || true
```

---

## Final Verification

- [ ] `Rscript -e 'devtools::test()' 2>&1 | tail -3` → FAIL 0, PASS ≥ 374
- [ ] `grep "gres2\|cres\b\|kMetricEps2\|kChi2Floor2" src/c_api.cpp src/chebyshev.cpp` → 0 output
- [ ] `grep "pack_cheb\|dispatch_cheb" src/c_api.cpp src/r_bridge.cpp` → 2 definitions present
