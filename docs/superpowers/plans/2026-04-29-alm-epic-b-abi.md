# ALM ieppa_soft — Epic B: ABI & Data Structures

> **For agentic workers:** Use `superpowers:subagent-driven-development`

**Branch:** `fix/correctness-performance-2026-04-28`
**Spec:** `docs/superpowers/specs/2026-04-29-ieppa-alm-soft-capacity.md`
**Goal:** Add `capacity_mu` / `capacity_mu_auto` / `rk_params_t` / `rk_result_t` fields + ABI bumps
**Mechanism:** Struct field additions, `n > 0` guard, `static_assert` tripwires
**Forbidden:** Changing algorithm logic; modifying fields used by lbfgsb (`alm_lambda`, `alm_mu` in `CalibState`); reordering existing fields
**Audit:** `EXPECTED_RK_RESULT_BYTES=480` and `EXPECTED_RK_PARAMS_BYTES=232` must `static_assert` clean after Task 3; compile gate (`R CMD INSTALL --preclean .`) after every task

---

## Baseline (read before touching anything)

Current struct sizes (verified from source 2026-04-29):

| Symbol | Current value | After Epic B |
|--------|--------------|-------------|
| `EXPECTED_RK_PARAMS_BYTES` | 224 | **232** |
| `EXPECTED_RK_RESULT_BYTES` | 448 | **480** |

`rk_params_t` last field before Epic B: `sor_burnin` (int, 4B) + 4B pad → total 224B.
Adding `capacity_penalty` (double, 8B) → 224 + 8 = **232B** (no new pad needed; double is 8-aligned, ends block on 8-byte boundary).

`rk_result_t` last field before Epic B: `convergence_minimized_metric` (int, 4B) + 4B trailing pad → 448B.
Adding 4 new fields: `alm_capacity_mu_final` (double, 8B) + `alm_n_growth_events` (int, 4B) + 4B pad + `alm_max_dual_norm` (double, 8B) + `alm_sum_drift` (double, 8B) = 32B → 448 + 32 = **480B**.

---

## Task 1: `src/types.hpp` — add `capacity_mu` to `CalibState`

**File:** `src/types.hpp`
**Ticket:** one independent bead per task (do not bundle)

### Step 1 — Read the file and locate the insertion point

```bash
grep -n "alm_mu" /home/dd/Gemini/leafblower/src/types.hpp
```

Expected hit (verified):
```
88:    double alm_mu     = 0.0;  // penalty coefficient; 0.0 = ALM inactive
```

### Step 2 — Add `capacity_mu` immediately after `alm_mu`

Edit `src/types.hpp`. After line 88 (`double alm_mu = 0.0;`), insert:

```cpp
    double capacity_mu = 0.0;  // ieppa_soft ALM penalty (capacity box constraint); 0.0 = inactive
```

Full context for the edit — `old_string`:
```cpp
    double alm_lambda = 0.0;  // dual variable for sum(w)=n; only read when alm_mu > 0
    double alm_mu     = 0.0;  // penalty coefficient; 0.0 = ALM inactive
    rk_bounds_mode_t bounds_mode = RK_BOUNDS_CELL;  /* P3.1: per-obs vs cell-aggregate bounds */
```

`new_string`:
```cpp
    double alm_lambda = 0.0;  // dual variable for sum(w)=n; only read when alm_mu > 0
    double alm_mu     = 0.0;  // penalty coefficient; 0.0 = ALM inactive
    double capacity_mu = 0.0;  // ieppa_soft ALM penalty (capacity box constraint); 0.0 = inactive
    rk_bounds_mode_t bounds_mode = RK_BOUNDS_CELL;  /* P3.1: per-obs vs cell-aggregate bounds */
```

### Step 3 — Compile gate

```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```

Expected: last line is `* DONE (leafblower)`.
If it fails: read the compiler error; do NOT guess a fix — stop and diagnose.

### Step 4 — Run tests

```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
```

Expected: `FAIL == 3` (pre-existing failures from Epic A / unimplemented ALM logic), all other tests pass. Zero new failures introduced by this task.

### Step 5 — Commit

```bash
git add src/types.hpp
git commit -m "feat(types): add capacity_mu field to CalibState for ieppa_soft ALM"
```

---

## Task 2: `src/cell_table.hpp` + `src/cell_table.cpp` — add `capacity_mu_auto`

**Files:** `src/cell_table.hpp`, `src/cell_table.cpp`
**Ticket:** one independent bead per task

### Step 1 — Read both files

Both files are already read in this session. Key facts:

- `CellTable` struct (cell_table.hpp, lines 7–13) currently ends with `double W_input;`
- `build_cell_table` in cell_table.cpp sets `out.M_cell` at line 90 (packed path) and line 134 (general path), then sets `out.W_input` at line 137–138 and returns 0 at line 139.

### Step 2A — Add `capacity_mu_auto` to `CellTable` struct in `src/cell_table.hpp`

Edit `src/cell_table.hpp`. `old_string`:
```cpp
    double W_input;                              // sum of input weights
};
```

`new_string`:
```cpp
    double W_input;                              // sum of input weights
    double capacity_mu_auto = 0.0;               // auto-computed ALM default: M_cell/n
};
```

### Step 2B — Compute `capacity_mu_auto` in `build_cell_table` in `src/cell_table.cpp`

Insert after `out.W_input` accumulation, before the `return 0`. `old_string`:
```cpp
    out.W_input = 0.0;
    for (int i = 0; i < n; i++) out.W_input += weights[i];
    return 0;
```

`new_string`:
```cpp
    out.W_input = 0.0;
    for (int i = 0; i < n; i++) out.W_input += weights[i];
    out.capacity_mu_auto = (n > 0 && out.M_cell > 0)
        ? static_cast<double>(out.M_cell) / static_cast<double>(n)
        : 1.0;
    return 0;
```

**Rationale for guard:** `n > 0 && out.M_cell > 0` — both are checked at the top of `build_cell_table` (`if (K <= 0 || n <= 0) return -1;`) so guard is defensive but harmless. `M_cell` is always `>= 1` after the sort loop if `n > 0`, but the guard makes the division unconditionally safe for future callers.

### Step 3 — Compile gate

```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```

Expected: `* DONE (leafblower)`.

### Step 4 — Run tests

```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
```

Expected: same FAIL count as after Task 1 (no new failures).

### Step 5 — Commit both files in one commit

```bash
git add src/cell_table.hpp src/cell_table.cpp
git commit -m "feat(cell_table): add capacity_mu_auto = M_cell/n auto-computed ALM default"
```

---

## Task 3: `src/leafblower.h` — ABI bumps

**File:** `src/leafblower.h`
**Ticket:** one independent bead per task

This task has the highest failure risk: the `static_assert` tripwires will reject a wrong byte count at compile time with an exact error message showing the computed size. Read that message before guessing.

### Step 1 — Read the file

Already read in this session. Key structural facts:

`rk_params_t` (lines 53–92):
- Last field before Epic B: `sor_burnin` (int) at line 90, then closing `}` at line 92.
- Current `EXPECTED_RK_PARAMS_BYTES = 224` (line 189).

`rk_result_t` (lines 95–129):
- Last field before Epic B: `convergence_minimized_metric` (int) at line 128, then closing `}` at line 129.
- Current `EXPECTED_RK_RESULT_BYTES = 448` (line 171).

### Step 2A — Add `capacity_penalty` to `rk_params_t`

Insert one field after `sor_burnin` and before the closing brace. `old_string`:
```c
    int    sor_burnin;
    /* ── End convergence/SOR config ── */
} rk_params_t;
```

`new_string`:
```c
    int    sor_burnin;
    /* ── End convergence/SOR config ── */
    double capacity_penalty;  /* ieppa_soft ALM penalty; <=0.0 = use auto (M_cell/n) */
} rk_params_t;
```

**ABI arithmetic:** `sor_burnin` (int, 4B) currently has 4B trailing pad to reach 224B. Adding `capacity_penalty` (double, 8B) consumes that 4B pad + adds 4B = 224 - 4 + 8 = **228B**... wait.

Re-derive: current 224B layout ends with `sor_burnin` int (4B) + **4B implicit pad** = 224B total. Inserting `double capacity_penalty` (8B) after `sor_burnin` before the closing brace: the int `sor_burnin` is followed by 4B pad (to align the double to 8B), then the double (8B). So the new total = 224 (without old trailing pad) - 0 + 4B pad + 8B = 224 + 4 + 8 - 4 = **232B**.

More precisely: old struct = 220B of fields + 4B pad = 224B. New struct = 220B + 4B align-pad + 8B `capacity_penalty` = 232B. No trailing pad needed (struct ends on 8-byte boundary). ✓

### Step 2B — Add 4 diagnostic fields to `rk_result_t`

Insert before the closing brace of `rk_result_t`. `old_string`:
```c
    double convergence_solver_objective;  /* solver's mathematical objective at best_iter */
    int    convergence_minimized_metric; /* CalibMetric: which metric was minimized */
    /* ── End extended quality metrics ── */
} rk_result_t;
```

`new_string`:
```c
    double convergence_solver_objective;  /* solver's mathematical objective at best_iter */
    int    convergence_minimized_metric; /* CalibMetric: which metric was minimized */
    /* ── End extended quality metrics ── */
    /* ── ieppa_soft ALM diagnostics ── */
    double alm_capacity_mu_final;   /* final capacity_mu after adaptive scaling; 0 if not ieppa_soft */
    int    alm_n_growth_events;     /* adaptive growth fire count; 0 if not ieppa_soft */
    double alm_max_dual_norm;       /* max |lambda_cell[c]| at exit */
    double alm_sum_drift;           /* |sum(X) - n| after final projection */
    /* ── End ieppa_soft ALM diagnostics ── */
} rk_result_t;
```

**ABI arithmetic:** Current `rk_result_t` = 448B. Last field before insert: `convergence_minimized_metric` (int, 4B) + 4B trailing pad = 448B total. Adding:
- `alm_capacity_mu_final` (double, 8B): consumes existing 4B pad + new 4B → 448 + 4 + 4 = 456B? No.

Re-derive: old struct = 444B of declared fields + 4B trailing pad = 448B. Inserting after `convergence_minimized_metric` (int):
- 4B pad to align `alm_capacity_mu_final` (double) → 4B
- `alm_capacity_mu_final` (double) → 8B
- `alm_n_growth_events` (int) → 4B
- 4B pad to align `alm_max_dual_norm` (double)
- `alm_max_dual_norm` (double) → 8B
- `alm_sum_drift` (double) → 8B

New total = 444 (old fields) + 4 (align pad) + 8 + 4 + 4 + 8 + 8 = 444 + 36 = 480B. No trailing pad (struct ends on 8-byte boundary). ✓

### Step 2C — Update `EXPECTED_RK_PARAMS_BYTES`

`old_string`:
```c
/* ABI layout (2026-04-24): added overlay fields after bounds_mode.
 *   rk_homotopy_cfg_t (n_levels int + 4B pad + 3 doubles + enabled int + 4B pad = 40B)
 *   rk_scheduler_t (int, 4B) + rk_eta_mode_t (int, 4B)
 *   eta_start (double, 8B) + eta_end (double, 8B) + eta_schedule_power (double, 8B)
 * WU-A (2026-04-25): convergence redesign — criterion→metric+rule in rk_params_t:
 *   pct_tol (double, 8B) + absolute_tol (double, 8B)
 *   metric (int, 4B) + rule (int, 4B) + stop_when (int, 4B)
 *   sor_enabled (int, 4B) + sor_auto (int, 4B)
 *   sor_omega_init (double, 8B) + sor_omega_min (double, 8B) + sor_omega_fixed (double, 8B)
 *   sor_burnin (int, 4B) + 4B pad
 * Total: 224B. Verified 2026-04-25 Linux x86_64. */
#define EXPECTED_RK_PARAMS_BYTES 224
```

`new_string`:
```c
/* ABI layout (2026-04-24): added overlay fields after bounds_mode.
 *   rk_homotopy_cfg_t (n_levels int + 4B pad + 3 doubles + enabled int + 4B pad = 40B)
 *   rk_scheduler_t (int, 4B) + rk_eta_mode_t (int, 4B)
 *   eta_start (double, 8B) + eta_end (double, 8B) + eta_schedule_power (double, 8B)
 * WU-A (2026-04-25): convergence redesign — criterion→metric+rule in rk_params_t:
 *   pct_tol (double, 8B) + absolute_tol (double, 8B)
 *   metric (int, 4B) + rule (int, 4B) + stop_when (int, 4B)
 *   sor_enabled (int, 4B) + sor_auto (int, 4B)
 *   sor_omega_init (double, 8B) + sor_omega_min (double, 8B) + sor_omega_fixed (double, 8B)
 *   sor_burnin (int, 4B) + 4B align-pad
 * Epic B (2026-04-29): ieppa_soft ABI:
 *   capacity_penalty (double, 8B) — consumes prior 4B pad + adds 4B = +8B
 * Total: 232B. Verified 2026-04-29 Linux x86_64. */
#define EXPECTED_RK_PARAMS_BYTES 232
```

### Step 2D — Update `EXPECTED_RK_RESULT_BYTES`

`old_string`:
```c
/* rk_result_t tripwire. Linux x86_64, verified 2026-04-24: 448 bytes. */
#define EXPECTED_RK_RESULT_BYTES 448
```

`new_string`:
```c
/* rk_result_t tripwire. Linux x86_64, verified 2026-04-29: 480 bytes.
 * Epic B additions: alm_capacity_mu_final (double,8B) + 4B align-pad after
 * convergence_minimized_metric (int) + alm_n_growth_events (int,4B) +
 * 4B align-pad + alm_max_dual_norm (double,8B) + alm_sum_drift (double,8B) = +32B. */
#define EXPECTED_RK_RESULT_BYTES 480
```

### Step 3 — Compile gate (tripwire check)

```bash
R CMD INSTALL --preclean . 2>&1 | tail -10
```

**If `static_assert` fires:** the error message contains the compiler-computed sizeof. Read that value and update the macro to match it — do NOT guess. Example error:
```
src/leafblower.h:172:1: error: static assertion failed: rk_result_t size changed; update EXPECTED_RK_RESULT_BYTES
note: sizeof(rk_result_t) evaluates to 480
```
If the message shows a value other than 480, update `EXPECTED_RK_RESULT_BYTES` to that value and recompile. Do not silently accept a wrong value.

Expected success: `* DONE (leafblower)`.

### Step 4 — Run tests

```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
```

Expected: same FAIL count as before (no new failures). The new fields are zero-initialized by `rk_result_init` (memset), so no behavioral change.

### Step 5 — Commit

```bash
git add src/leafblower.h
git commit -m "feat(api): add capacity_penalty/alm_* fields; bump ABI tripwires to 232/480"
```

---

## Execution order

All 3 tasks are independent. Suggested serial order: Task 1 → Task 2 → Task 3.
Parallel execution is safe (each touches a distinct file set) if the worker supports it.

## Definition of Done

- [ ] `src/types.hpp` has `capacity_mu` after `alm_mu`
- [ ] `src/cell_table.hpp` has `capacity_mu_auto` in `CellTable`
- [ ] `src/cell_table.cpp` computes `capacity_mu_auto = M_cell/n` (guarded) after `W_input`
- [ ] `src/leafblower.h` has `capacity_penalty` in `rk_params_t`
- [ ] `src/leafblower.h` has 4 `alm_*` fields at end of `rk_result_t`
- [ ] `EXPECTED_RK_PARAMS_BYTES == 232` and `static_assert` passes
- [ ] `EXPECTED_RK_RESULT_BYTES == 480` and `static_assert` passes
- [ ] `R CMD INSTALL --preclean .` → `* DONE` after each task
- [ ] No new test failures beyond pre-existing FAIL==3
- [ ] 3 commits pushed to `fix/correctness-performance-2026-04-28`
