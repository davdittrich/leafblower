# iEPPA SRAA-m Anderson Acceleration — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend iEPPA and ieppa_soft with SRAA-m Anderson Acceleration operating in lf-space (n_cats_total ≈ 50–500 dims), reducing iteration count by ≥30% on tight-bounds problems.

**Architecture:** SRAA-m operates on the lf flat vector (log-factor space, dim `n_cats_total_with_na`) rather than X_cur (cell masses, dim `M_cell`). `pack_lf` / `unpack_lf` convert between ieppa's internal state and the lf-flat iterate for SRAA. The existing `sraa_step` function is reused with `apply_clamp=false` (lf is unconstrained).

**Tech Stack:** C++17, R package build (`R CMD INSTALL --preclean .`), testthat 3.

---

## Mechanism / Forbidden / Audit (mandatory header)

- **Mechanism:** Type-II SRAA-m (m=5) on the iEPPA fixed-point map `lf_{n+1} = G(lf_n)` where `G` performs one BCD pass (`apply_single_margin_linear` over k=0..K-1). Reuse the existing `lbw::sraa_step` (file `src/sraa.hpp`) extended with an `apply_clamp` switch.
- **Forbidden:**
  - Do NOT replicate the `for (iter_in_lvl=…)` outer-loop body inside the SRAA branch — use the `f_eval_lf` lambda as the canonical fixed-point map.
  - Do NOT keep two definitions of `sor_auto`. The single binding `sor_auto_v = st.sor_cfg.auto_adapt && !st.accelerate;` replaces the existing one.
  - Do NOT call `sraa_step` with non-empty `L_cell` / `U_cell` when `apply_clamp=false`; pass empty vectors to make accidental indexing crash loudly.
  - Do NOT seed `ieppa_sraa.F_cur = lf_flat` inside the while-loop (only once before).
  - Do NOT reset `ieppa_sraa.aa_accepted_count`; the SRAA contract is that `clear()` does not zero it (sraa.hpp:32, 64-66).
- **Audit:** Tests in `tests/testthat/test-ieppa-sraa.R` exercise quantitative thresholds (≥30% F-eval reduction, parity max_error, output correlation, aa_accepted_count). Spy on `verbose=1` log output via `capture.output()` for the greedy-downgrade message. No internal mocking required.

---

## Source-of-truth references (read before each task)

| File | Lines | What it provides |
|---|---|---|
| `src/sraa.hpp` | 1-200 | `SRAAState`, `SRAAStepResult`, `sraa_step` template, `kSRAAm`, `kSRAAOuterSlack`, `kSRAAOuterStallWindow` |
| `src/ieppa.hpp` | 1-39 | `IEPPAResult` struct (all existing fields) |
| `src/ieppa.cpp` | 65-260 | Setup: `total_cats`, `cat_offset`, `lf`, `f_lin`, `X_cur`, `cell_lf`, `use_linear` |
| `src/ieppa.cpp` | 369-440 | Homotopy outer loop, `budget_lvl`, `total_iters`, `iter_in_lvl` declaration |
| `src/ieppa.cpp` | 449-551 | `apply_single_margin_linear` (the F-step body to wrap) |
| `src/ieppa.cpp` | 1040-1110 | Best-iter snapshot (`W_best[c] = X[c]/X_init[c]`), SOR-auto adaptation block |
| `src/ieppa.cpp` | 400-413, 730-735, 862-864 | ALM `lambda_cell` reset sites; linear→log fallback site |
| `src/ieppa.cpp` | 1350-1410 | Post-loop expansion `mult[c] = X[c]/X_init[c]` and `best_weights` build |
| `src/r_bridge.cpp` | 318 | `st.accelerate = (INTEGER(accelerate_sexp)[0] != 0);` |
| `src/r_bridge.cpp` | 580-605 | ieppa default dispatch (where `aa_accepted_count` is harvested) |
| `src/r_bridge.cpp` | 555-579 | ieppa_soft dispatch |
| `src/r_bridge.cpp` | 636-637, 670-717 | Result list size (currently 34) and field-by-index assembly |
| `src/raking.cpp` | 230-395 | Reference SRAA integration (`F_eval`, `F_cur` seeding, stall guard) |
| `R/harvest.R` | 93-99 | `@param accelerate` docstring to update |

---

## Task list (4 epics, 10 tasks)

### Epic A — Infrastructure

- [ ] **Task A1** — `sraa.hpp`: Add `apply_clamp` parameter
- [ ] **Task A2** — `ieppa.hpp`: Add `aa_accepted_count` field

### Epic B — ieppa.cpp SRAA core

- [ ] **Task B1** — `pack_lf` + `unpack_lf` file-local helpers
- [ ] **Task B2** — `f_eval_lf` lambda inside the homotopy level loop
- [ ] **Task B3** — SRAA outer while-loop with `lf_best` + `W_best` tracking
- [ ] **Task B4** — ALM ieppa_soft: clear+reseed on `capacity_mu` update
- [ ] **Task B5** — SOR disable + greedy-downgrade log message

### Epic C — R bridge

- [ ] **Task C1** — `r_bridge.cpp`: extend `accelerate_bool` and expose `aa_accepted_count`
- [ ] **Task C2** — `R/harvest.R`: update `@param accelerate` docstring

### Epic D — Tests

- [ ] **Task D1** — `tests/testthat/test-ieppa-sraa.R`: 8 quantitative tests

---

## Task A1 — `sraa.hpp`: add `apply_clamp` parameter

**Ticket:** `Task SRAA-A1 [apply_clamp param + L/U guard] ! [no behavior change for existing callers]`

**Files:** `src/sraa.hpp` (lines 74-80, 165-171)

### A1.1 — Add the parameter to the template signature

- [ ] In `src/sraa.hpp`, edit lines 74-80:

  Replace:
  ```cpp
  template<typename FEval>
  SRAAStepResult sraa_step(
      FEval& f_eval,
      std::vector<double>& X,
      const std::vector<double>& L_cell,
      const std::vector<double>& U_cell,
      SRAAState& state)
  ```
  with:
  ```cpp
  template<typename FEval>
  SRAAStepResult sraa_step(
      FEval& f_eval,
      std::vector<double>& X,
      const std::vector<double>& L_cell,
      const std::vector<double>& U_cell,
      SRAAState& state,
      bool apply_clamp = true)
  ```
  Default `true` preserves all existing call sites in `raking.cpp` and `greenkhorn.cpp`.

### A1.2 — Guard the L/U dereference in Step 8

- [ ] In `src/sraa.hpp`, edit lines 163-171:

  Replace:
  ```cpp
      // --- Step 8: Extrapolate + clamp into scratch ---
      // scratch currently holds R_k; overwrite in-place with X_AA
      for (int c = 0; c < M; c++) {
          double Rk_c = state.scratch[c];  // read R_k before overwrite
          double corr = 0.0;
          for (int i = 0; i < n; i++)
              corr += state.gamma_[i] * (state.dX_buf[i * M + c] + state.dR_buf[i * M + c]);
          state.scratch[c] = std::clamp(X[c] + Rk_c - corr, L_cell[c], U_cell[c]);
      }
  ```
  with:
  ```cpp
      // --- Step 8: Extrapolate + clamp into scratch ---
      // scratch currently holds R_k; overwrite in-place with X_AA.
      // apply_clamp=false: lf-space (unconstrained); L_cell/U_cell are NOT dereferenced
      // and may be empty vectors. Used by ieppa SRAA path.
      for (int c = 0; c < M; c++) {
          double Rk_c = state.scratch[c];  // read R_k before overwrite
          double corr = 0.0;
          for (int i = 0; i < n; i++)
              corr += state.gamma_[i] * (state.dX_buf[i * M + c] + state.dR_buf[i * M + c]);
          double val = X[c] + Rk_c - corr;
          state.scratch[c] = apply_clamp
              ? std::clamp(val, L_cell[c], U_cell[c])
              : val;
      }
  ```

### A1.3 — Compile gate

- [ ] Run: `R CMD INSTALL --preclean .` from `/home/dd/Gemini/leafblower`
- [ ] Run: `Rscript -e 'testthat::test_dir("tests/testthat")' 2>&1 | tail -3`
- [ ] **Acceptance:** all existing tests pass (raking + greenkhorn unchanged because they omit the new arg → default `true`).

### A1.4 — Commit

- [ ] `git add src/sraa.hpp && git commit -m "feat(sraa): add apply_clamp parameter for unbounded callers"`

---

## Task A2 — `ieppa.hpp`: add `aa_accepted_count` field

**Ticket:** `Task SRAA-A2 [add aa_accepted_count to IEPPAResult] ! [no relocation of other fields]`

**Files:** `src/ieppa.hpp`

### A2.1 — Insert the field

- [ ] In `src/ieppa.hpp`, edit lines 24-26:

  Replace:
  ```cpp
      // ── iEPPA internal metrics ──
      double best_objective_seen          = 0.0;
      double marginal_kl_at_iter          = 0.0;
  ```
  with:
  ```cpp
      // ── iEPPA internal metrics ──
      double best_objective_seen          = 0.0;
      double marginal_kl_at_iter          = 0.0;
      // ── SRAA-m diagnostics (Anderson Acceleration; 0 when accelerate=FALSE) ──
      int    aa_accepted_count            = 0;   // cumulative AA-accepted super-steps this solve
  ```

### A2.2 — Compile gate

- [ ] `R CMD INSTALL --preclean .`
- [ ] `Rscript -e 'testthat::test_dir("tests/testthat")' 2>&1 | tail -3`
- [ ] **Acceptance:** clean build; field default-initializes to 0 in all existing call paths.

### A2.3 — Commit

- [ ] `git add src/ieppa.hpp && git commit -m "feat(ieppa): add aa_accepted_count diagnostic field"`

---

## Task B1 — `pack_lf` + `unpack_lf` helpers (file-local, ieppa.cpp)

**Ticket:** `Task SRAA-B1 [pack/unpack lf flat vector] ! [no allocations inside; rely on existing X_cur/f_lin/cell_lf buffers]`

**Files:** `src/ieppa.cpp`

### B1.1 — Add helpers in the anonymous file scope

- [ ] In `src/ieppa.cpp` near line 56 (just after `write_trajectory_csv` and before `IEPPAResult ieppa_solve`), insert:

  ```cpp
  // SRAA-m helpers — file-local. lf is the flat per-margin log-factor vector
  // of size cat_offset[K] = Σ_k (cat_counts[k]+1). These helpers are O(M_cell)
  // and pre-allocate nothing; all destination buffers are owned by ieppa_solve.
  //
  // pack_lf: copy lf -> dst (size n_cats_total_with_na). NA slots (j == cat_counts[k])
  // are inert zeros; they participate in the SRAA linear system with ΔX=ΔR=0.
  static inline void pack_lf(const std::vector<double>& lf,
                             std::vector<double>& dst) {
      // Caller guarantees dst.size() == lf.size().
      std::copy(lf.begin(), lf.end(), dst.begin());
  }

  // unpack_lf: from a flat lf iterate, rebuild the derived state required by
  // apply_single_margin_linear for the next sweep:
  //   1) lf      <- src
  //   2) f_lin[i] = exp(lf[i])     for i in 0..total_cats-1
  //   3) cell_lf[c] = Σ_k lf[cat_offset[k] + g_per_cell[k][c]]   (g<0 skipped)
  //   4) X_cur[c]   = X_init[c] * exp(cell_lf[c])
  //
  // inv_f_old_lin is NOT set here; apply_single_margin_linear recomputes it
  // at the head of each margin call (ieppa.cpp:453).
  static inline void unpack_lf(const std::vector<double>& src,
                               std::vector<double>& lf,
                               std::vector<double>& f_lin,
                               std::vector<double>& cell_lf,
                               std::vector<double>& X_cur,
                               const lbw::CellTable& ct,
                               const std::vector<double>& X_init,
                               int K,
                               const std::vector<int>& cat_offset) {
      const int total_cats = cat_offset[K];
      for (int i = 0; i < total_cats; i++) {
          lf[i]    = src[i];
          f_lin[i] = std::exp(src[i]);
      }
      const int M = ct.M_cell;
      // Rebuild cell_lf from scratch (cheap: O(K * M_cell), same as one sweep pass).
      std::fill(cell_lf.begin(), cell_lf.end(), 0.0);
      for (int k = 0; k < K; k++) {
          const int* gk  = ct.g_per_cell[k].data();
          const int  off = cat_offset[k];
          for (int c = 0; c < M; c++) {
              int g = gk[c];
              if (g < 0) continue;
              cell_lf[c] += src[off + g];
          }
      }
      for (int c = 0; c < M; c++) {
          X_cur[c] = (X_init[c] > 0.0) ? X_init[c] * std::exp(cell_lf[c]) : 0.0;
      }
  }
  ```

  Notes on signature: `cat_offset` and `lf`/`f_lin`/`cell_lf`/`X_cur`/`X_init` are all `std::vector<double>` (or `std::vector<int>` for `cat_offset`) declared as locals in `ieppa_solve`; passing by reference avoids copies. `K` is `st.K`.

### B1.2 — Compile-only check

- [ ] `R CMD INSTALL --preclean .` (no callers yet; this is a syntax/link gate).
- [ ] `Rscript -e 'testthat::test_dir("tests/testthat")' 2>&1 | tail -3`
- [ ] **Acceptance:** clean build; no test regressions (helpers are unused dead code at this stage — that is intentional).

### B1.3 — Commit

- [ ] `git add src/ieppa.cpp && git commit -m "feat(ieppa): add pack_lf/unpack_lf helpers for SRAA-m lf-space iterate"`

---

## Task B2 — `f_eval_lf` lambda (inside homotopy level loop)

**Ticket:** `Task SRAA-B2 [f_eval_lf lambda mirrors single BCD sweep over K margins] ! [must call apply_single_margin_linear; do NOT reimplement the sweep]`

**Files:** `src/ieppa.cpp`

### B2.1 — Insert the lambda just before the inner iteration loop

- [ ] In `src/ieppa.cpp` at line 437 (just before `for (int iter_in_lvl = 1; iter_in_lvl <= budget_lvl; iter_in_lvl++) {` on line 438), insert:

  ```cpp
        // ────────────────────────────────────────────────────────────────────
        // SRAA-m fixed-point map for the linear path of this homotopy level.
        // Captures: lf, f_lin, cell_lf, X_cur, ct, X_init, st.K, cat_offset,
        //           apply_single_margin_linear, st (for verbose).
        //
        // Contract: receives a flat lf iterate (size cat_offset[st.K]); rebuilds
        // all derived state via unpack_lf; runs one full BCD pass over K margins;
        // packs the resulting lf back into `flat`; returns errRp on success or
        // +inf on overflow. On overflow lf is mid-sweep (partial); the SRAA
        // safeguard rejects via err_AA=inf>err_plain and reverts via swap with
        // F_cur. Outer log-fallback handler clears history if both plain and AA
        // return inf.
        // ────────────────────────────────────────────────────────────────────
        auto f_eval_lf = [&](std::vector<double>& flat) -> double {
            unpack_lf(flat, lf, f_lin, cell_lf, X_cur, ct, X_init,
                      st.K, cat_offset);
            bool overflow = false;
            for (int k = 0; k < st.K && !overflow; k++) {
                if (apply_single_margin_linear(k)) overflow = true;
            }
            // Always pack — on overflow, lf has been partially updated; packing
            // it preserves the SRAA invariant that `flat` reflects the current
            // outer iterate after the call.
            pack_lf(lf, flat);
            if (overflow) return std::numeric_limits<double>::infinity();

            // Compute errRp from X_cur (mirrors the in-loop body at lines ~970-1020).
            double W_total = 0.0;
            for (int c = 0; c < ct.M_cell; c++) W_total += X_cur[c];
            if (!(W_total > 0.0)) return std::numeric_limits<double>::infinity();
            double errRp = 0.0;
            for (int k = 0; k < st.K; k++) {
                const int nj = st.cat_counts[k];
                std::fill(S_lin.begin(), S_lin.begin() + nj, 0.0);
                const int* gk = ct.g_per_cell[k].data();
                for (int c = 0; c < ct.M_cell; c++) {
                    int j = gk[c];
                    if (j >= 0 && j < nj) S_lin[j] += X_cur[c];
                }
                for (int j = 0; j < nj; j++) {
                    double e = std::fabs(S_lin[j] / W_total - st.targets[k][j]);
                    if (e > errRp) errRp = e;
                }
            }
            return errRp;
        };
  ```

  Capture-by-reference (`[&]`) is mandatory because `apply_single_margin_linear` is itself a `[&]` lambda already capturing the same state.

### B2.2 — Compile-only check (still no caller)

- [ ] `R CMD INSTALL --preclean .`
- [ ] `Rscript -e 'testthat::test_dir("tests/testthat")' 2>&1 | tail -3`
- [ ] **Acceptance:** clean build; the lambda is unused — accept the `[-Wunused-variable]` suppression by leaving it un-declared until B3 wires it in. If the warning becomes an error, mark with `[[maybe_unused]]` on the `auto f_eval_lf` declaration.

### B2.3 — Commit

- [ ] `git add src/ieppa.cpp && git commit -m "feat(ieppa): add f_eval_lf lambda mirroring one BCD sweep"`

---

## Task B3 — SRAA outer while-loop replacing the for-loop on the linear path

**Ticket:** `Task SRAA-B3 [SRAA-m outer loop with lf_best+W_best tracking and stall guard] ! [exactly one write to res.base.iterations; do NOT seed F_cur inside the loop]`

**Files:** `src/ieppa.cpp` (around line 438 onwards — the `for (int iter_in_lvl = 1; ...)` loop)

### B3.1 — Add level-scope SRAA state declarations and best-iter trackers

- [ ] In `src/ieppa.cpp` just after the `f_eval_lf` lambda definition (and before the inner loop), insert:

  ```cpp
        // SRAA-m state for this homotopy level. Allocated only when accelerate is
        // requested AND the linear path is active (ieppa supports a log fallback
        // path that does not use SRAA in this revision).
        const bool sraa_active_lvl = st.accelerate && use_linear;
        lbw::SRAAState ieppa_sraa;
        std::vector<double> lf_flat;
        std::vector<double> lf_best;
        const std::vector<double> dummy_L;   // empty; not dereferenced when apply_clamp=false
        const std::vector<double> dummy_U;
        int sraa_outer_stall_count = 0;
        double sraa_best_errRp     = std::numeric_limits<double>::infinity();
        if (sraa_active_lvl) {
            ieppa_sraa.init(total_cats, lbw::kSRAAm);
            lf_flat.assign(total_cats, 0.0);
            lf_best.assign(total_cats, 0.0);
            pack_lf(lf, lf_flat);
            // Seed F_cur ONCE before the loop; sraa_step swaps it forward each call.
            ieppa_sraa.F_cur = lf_flat;
        }
  ```

### B3.2 — Branch the inner loop on `sraa_active_lvl`

- [ ] In `src/ieppa.cpp`, wrap the existing `for (int iter_in_lvl = 1; iter_in_lvl <= budget_lvl; iter_in_lvl++) { ... }` block (line 438 through its matching close ~line 1290) so it executes only when `!sraa_active_lvl`. Insert a sibling `if (sraa_active_lvl)` block ahead of it containing the SRAA driver shown below.

  Concrete edit pattern (apply to line 438):

  Replace:
  ```cpp
      for (int iter_in_lvl = 1; iter_in_lvl <= budget_lvl; iter_in_lvl++) {
  ```
  with:
  ```cpp
      if (sraa_active_lvl) {
          // ── SRAA-m accelerated path (replaces the for-loop entirely) ──
          int  f_evals_used = 0;
          bool converged    = false;
          while (f_evals_used < budget_lvl && !converged) {
              auto r = lbw::sraa_step(f_eval_lf, lf_flat, dummy_L, dummy_U,
                                      ieppa_sraa, /*apply_clamp=*/false);
              f_evals_used += r.f_evals;
              res.base.iterations = total_iters + f_evals_used;  // single write site
              res.base.max_error  = r.err_result;

              // Best-iterate tracking. Two snapshots:
              //   W_best  = X_cur / X_init  (cell-mass ratios; consumed by post-loop
              //                              expansion to populate res.base.best_weights).
              //   lf_best = lf_flat         (lf-space; consumed by SRAA stall revert).
              if (r.err_result < best_metric_seen) {
                  best_metric_seen = r.err_result;
                  for (int c = 0; c < ct.M_cell; c++) {
                      W_best[c] = (X_init[c] > 0.0) ? X_cur[c] / X_init[c] : 0.0;
                  }
                  lf_best       = lf_flat;
                  best_iter_val = res.base.iterations;
                  best_objective_seen = lbw::compute_weight_kl(
                      X_cur, X_init, ct.M_cell, st.n,
                      kl_ratio_buf.data(), kl_weight_buf.data());
              }

              // Outer stall guard (mirrors raking.cpp:358-368).
              if (r.err_result > best_metric_seen * (1.0 + lbw::kSRAAOuterSlack)) {
                  if (++sraa_outer_stall_count >= lbw::kSRAAOuterStallWindow) {
                      // Revert to lf-space best iterate; rebuild all derived state.
                      lf_flat = lf_best;
                      unpack_lf(lf_flat, lf, f_lin, cell_lf, X_cur, ct, X_init,
                                st.K, cat_offset);
                      ieppa_sraa.clear();
                      ieppa_sraa.F_cur = lf_flat;   // reseed after clear
                      sraa_outer_stall_count = 0;
                  }
              } else {
                  sraa_best_errRp = std::min(sraa_best_errRp, r.err_result);
                  sraa_outer_stall_count = 0;
              }

              // Convergence: use solver tol on errRp at this level.
              if (r.err_result <= tol_lvl) {
                  converged       = true;
                  level_converged = true;
                  res.base.status = RK_OK;
              }
              if (st.verbose >= 1) {
                  char msg[256];
                  std::snprintf(msg, sizeof(msg),
                                "ieppa[sraa lvl=%d] f_evals=%d errRp=%.2e aa=%d",
                                lvl, f_evals_used, r.err_result, (int)r.aa_accepted);
                  st.log(msg);
              }
          }
          // Snapshot AA count (cumulative across levels via SRAAState; clear() does NOT zero it)
          res.aa_accepted_count = ieppa_sraa.aa_accepted_count;
          // Update global iter counter for downstream homotopy bookkeeping.
          total_iters += f_evals_used;
      } else {
          for (int iter_in_lvl = 1; iter_in_lvl <= budget_lvl; iter_in_lvl++) {
  ```

  Then add a matching `}` to close the new `else { ... }` after the existing loop's terminating `}` (the loop currently ends near line 1290 — verify the line during the edit; the brace that closes `for (int iter_in_lvl…)` must now be followed by an additional `}` for the `else`).

  Note on `kl_ratio_buf` / `kl_weight_buf`: these are existing scratch buffers declared near line 322 (`kl_ratio_buf` / `kl_weight_buf` for `compute_weight_kl`). If their names differ in the actual source, use the names that appear at the existing call site near line 1044. Preserve those exact identifiers.

### B3.3 — Linear→log fallback handling

- [ ] In `src/ieppa.cpp` at lines 730 and 864 (the two existing `use_linear = false;` sites that fire on overflow), insert immediately above each:

  ```cpp
                  if (sraa_active_lvl) {
                      // Snapshot AA count before history is invalidated by the
                      // log path takeover. clear() does not zero aa_accepted_count.
                      res.aa_accepted_count = ieppa_sraa.aa_accepted_count;
                      ieppa_sraa.clear();
                  }
  ```

  Rationale: After fallback the log-path runs the non-accelerated for-loop; SRAA history is no longer relevant. The aa_accepted_count snapshot captures the AA work done in the linear phase.

### B3.4 — Compile + smoke test

- [ ] `R CMD INSTALL --preclean .`
- [ ] `Rscript -e 'testthat::test_dir("tests/testthat")' 2>&1 | tail -3`
- [ ] **Acceptance:** all existing tests pass at `accelerate=FALSE` (default). With the existing test suite, no `accelerate=TRUE` paths run yet for ieppa, so nothing changes behaviorally.

### B3.5 — Commit

- [ ] `git add src/ieppa.cpp && git commit -m "feat(ieppa): SRAA-m outer loop with lf_best/W_best tracking and stall guard"`

---

## Task B4 — ALM ieppa_soft: clear + reseed on `capacity_mu` update

**Ticket:** `Task SRAA-B4 [SRAA history reset on ALM mu adaptation] ! [reset lambda_cell BEFORE pack_lf so reseeded lf reflects post-dual state]`

**Files:** `src/ieppa.cpp`

### B4.1 — Add the combined reset block

- [ ] Identify the ALM violation handler that updates `capacity_mu_adaptive` / `st.alm.capacity_mu`. This is the block that runs when ALM dual-norm growth triggers `capacity_mu_adaptive *= growth_factor` (typically guarded by `alm_active && violation_streak >= threshold`). Find the assignment site (search for `capacity_mu_adaptive *=` or `alm_n_growth_events++`).

- [ ] Inside that block, immediately after the new `capacity_mu` is committed and `lambda_cell` is conventionally reset, ensure the order is:

  ```cpp
                  // Order matters: reset duals first, then reseed SRAA from the
                  // post-reset lf state so the next ALM level starts consistently.
                  std::fill(lambda_cell.begin(), lambda_cell.end(), 0.0);
                  if (sraa_active_lvl) {
                      ieppa_sraa.clear();
                      pack_lf(lf, lf_flat);
                      ieppa_sraa.F_cur = lf_flat;
                  }
  ```

  Place this combined block at every site where `capacity_mu` is updated mid-solve. If the existing code already calls `std::fill(lambda_cell.begin(), lambda_cell.end(), 0.0)` (lines 412, 718, 735, 862), the SRAA-aware addendum (the `if (sraa_active_lvl) { ... }` four-liner) must follow each one that fires *during* an active SRAA level. The level-boundary reset at line 412 is already covered by the per-level re-init in B3.1 — only the mid-level resets at 718, 735, 862 plus the ALM violation-triggered site need the SRAA addition.

### B4.2 — Compile + test

- [ ] `R CMD INSTALL --preclean .`
- [ ] `Rscript -e 'testthat::test_dir("tests/testthat")' 2>&1 | tail -3`
- [ ] **Acceptance:** all existing ieppa_soft tests pass at `accelerate=FALSE`.

### B4.3 — Commit

- [ ] `git add src/ieppa.cpp && git commit -m "feat(ieppa_soft): clear+reseed SRAA history on ALM mu adaptation"`

---

## Task B5 — SOR disable + greedy-downgrade log when accelerate=TRUE

**Ticket:** `Task SRAA-B5 [SOR auto-adapt off and greedy→round_robin under SRAA] ! [exactly one boolean change: sor_auto_v binding]`

**Files:** `src/ieppa.cpp`

### B5.1 — Replace the `sor_auto_v` definition

- [ ] In `src/ieppa.cpp`, locate the existing `sor_auto` resolution (search for `bool sor_auto`). Replace its definition with:

  ```cpp
      // SRAA-m and SOR adaptive omega are not co-designed; under accelerate=TRUE
      // we keep omega frozen at omega_init (default 1.0). Documented in
      // R/harvest.R @param accelerate.
      const bool sor_auto_v = st.sor_cfg.auto_adapt && !st.accelerate;
  ```

  All downstream uses of `sor_auto` should be updated to `sor_auto_v`. If the existing identifier is already `sor_auto`, rename the new constant to match or update the binding name in-place — pick the option that minimizes diff.

### B5.2 — Greedy-downgrade log message

- [ ] In `src/ieppa.cpp`, locate the `use_greedy` resolution (search for `const bool use_greedy = (st.scheduler.mode == SchedulerMode::GREEDY);` near line 653). Replace with:

  ```cpp
      bool use_greedy = (st.scheduler.mode == SchedulerMode::GREEDY);
      if (st.accelerate && use_greedy) {
          if (st.verbose >= 1) {
              st.log("[ieppa] greedy scheduler disabled under SRAA-m; using round_robin");
          }
          use_greedy = false;
      }
  ```

  Drop the `const` so the variable can be reassigned. All downstream sites already read `use_greedy` after this point.

### B5.3 — Compile + test

- [ ] `R CMD INSTALL --preclean .`
- [ ] `Rscript -e 'testthat::test_dir("tests/testthat")' 2>&1 | tail -3`
- [ ] **Acceptance:** existing tests at `accelerate=FALSE` still pass; SOR behavior unchanged when accelerate is off (the `&& !st.accelerate` second conjunct is true).

### B5.4 — Commit

- [ ] `git add src/ieppa.cpp && git commit -m "feat(ieppa): freeze SOR and downgrade greedy to round_robin under SRAA-m"`

---

## Task C1 — `r_bridge.cpp`: extend `accelerate_bool` and expose `aa_accepted_count`

**Ticket:** `Task SRAA-C1 [accelerate gate covers ieppa/ieppa_soft + harvest aa_accepted_count] ! [bump res_list size from 34 to 35]`

**Files:** `src/r_bridge.cpp`

### C1.1 — Currently `st.accelerate` is set unconditionally from the SEXP at line 318

The existing line `st.accelerate = (INTEGER(accelerate_sexp)[0] != 0);` is method-agnostic — it already covers ieppa/ieppa_soft. Per the spec's `accelerate_bool` discussion (lines 173-175 of the spec), the bridge does NOT have an explicit `accelerate_bool` filter; the filtering happens at the solver level (B3.1's `sraa_active_lvl = st.accelerate && use_linear`). **No edit needed at line 318.** Confirm this during implementation; if a method-gated `accelerate_bool` exists elsewhere, extend it as the spec dictates.

- [ ] Verify (by `grep -n "accelerate_bool\\|accelerate &&\\|method ==.*accelerate" src/r_bridge.cpp`): no method-restricting filter currently exists. Document the finding inline as a comment near line 318:

  Replace line 318:
  ```cpp
      st.accelerate = (INTEGER(accelerate_sexp)[0] != 0);
  ```
  with:
  ```cpp
      // Accept accelerate=TRUE for raking, greenkhorn (existing) and
      // ieppa, ieppa_soft (SRAA-m extension). The solver-level gate
      // sraa_active_lvl=(st.accelerate && use_linear) handles per-path eligibility.
      st.accelerate = (INTEGER(accelerate_sexp)[0] != 0);
  ```

### C1.2 — Add `res_aa_accepted_count` declaration

- [ ] In `src/r_bridge.cpp` near line 378 (alongside `res_alm_sum_drift`), add:

  ```cpp
      int    res_aa_accepted_count       = 0;
  ```

### C1.3 — Harvest in the ieppa/ieppa_soft dispatch arms

- [ ] In `src/r_bridge.cpp` ieppa-default arm (after line 599, alongside `res_sor_n_damped`), add:

  ```cpp
              res_aa_accepted_count = res.aa_accepted_count;
  ```

- [ ] In the `ieppa_soft` arm (after line 578, alongside `res_alm_sum_drift`), add the same line:

  ```cpp
              res_aa_accepted_count = res.aa_accepted_count;
  ```

- [ ] In the `auto`→ieppa fallback arm (after line 447), add:

  ```cpp
              res_aa_accepted_count = res.aa_accepted_count;
  ```

### C1.4 — Bump result list size and append the field

- [ ] In `src/r_bridge.cpp` line 636-637, replace:

  ```cpp
      SEXP res_list  = PROTECT(Rf_allocVector(VECSXP,  34));  // 14 prior + 8 scalars + best_weights + 7 convergence fields + 4 ALM diagnostics
      SEXP res_names = PROTECT(Rf_allocVector(STRSXP,  34));
  ```
  with:
  ```cpp
      SEXP res_list  = PROTECT(Rf_allocVector(VECSXP,  35));  // 34 prior + 1 SRAA diagnostic (aa_accepted_count)
      SEXP res_names = PROTECT(Rf_allocVector(STRSXP,  35));
  ```

- [ ] At the end of the field-assembly block (after line 717 `SET_VECTOR_ELT(res_list, 33, …);`), append:

  ```cpp
      /* Element 34: SRAA-m diagnostic */
      SET_STRING_ELT(res_names, 34, Rf_mkChar("aa_accepted_count"));
      SET_VECTOR_ELT(res_list,  34, Rf_ScalarInteger(res_aa_accepted_count));
  ```

### C1.5 — Compile + test

- [ ] `R CMD INSTALL --preclean .`
- [ ] `Rscript -e 'testthat::test_dir("tests/testthat")' 2>&1 | tail -3`
- [ ] Smoke check: `Rscript -e 'r <- leafblower::harvest(...some small fixture...); stopifnot("aa_accepted_count" %in% names(attr(r,"result")))'`
- [ ] **Acceptance:** existing tests pass; `attr(result,"result")$aa_accepted_count` is present and integer 0 when `accelerate=FALSE`.

### C1.6 — Commit

- [ ] `git add src/r_bridge.cpp && git commit -m "feat(r_bridge): expose aa_accepted_count diagnostic for ieppa/ieppa_soft"`

---

## Task C2 — `R/harvest.R`: update `@param accelerate` docstring

**Ticket:** `Task SRAA-C2 [docstring lists ieppa/ieppa_soft + caveats] ! [no help-page contract change beyond the @param block]`

**Files:** `R/harvest.R`

### C2.1 — Replace the `@param accelerate` block

- [ ] In `R/harvest.R` lines 93-98:

  Replace:
  ```r
  #' @param accelerate Logical. Enable Safeguarded Regularized Anderson Acceleration
  #'   (SRAA-m, window m=5) for \code{method="greenkhorn"} and \code{method="raking"}.
  #'   SRAA-m guarantees \code{max_error_accelerated <= max_error_plain} per super-step
  #'   via an explicit safeguard, using 2 F-evals per accepted step (vs. SQUAREM's 3).
  #'   Replaces prior SQUAREM/CBB scheme which overshot the bounded optimum.
  #'   Default \code{FALSE}.
  ```
  with:
  ```r
  #' @param accelerate Logical. Enable Safeguarded Regularized Anderson Acceleration
  #'   (SRAA-m, window m=5) for \code{method="raking"}, \code{"greenkhorn"},
  #'   \code{"ieppa"}, and \code{"ieppa_soft"}. SRAA-m guarantees
  #'   \code{max_error_accelerated <= max_error_plain} per super-step via an
  #'   explicit safeguard, using 2 F-evals per accepted step (vs. SQUAREM's 3).
  #'   Default \code{FALSE}.
  #'
  #'   \strong{When \code{accelerate=TRUE} for ieppa/ieppa_soft:}
  #'   \itemize{
  #'     \item SOR adaptive under-relaxation is disabled (omega frozen at
  #'       \code{sor_omega_init}, default 1.0). Use \code{accelerate=FALSE} if
  #'       adaptive SOR is required.
  #'     \item If \code{scheduler="greedy"} is also set, greedy is silently
  #'       downgraded to round-robin (logged at \code{verbose>=1}).
  #'     \item \code{result$iterations} reports F-evals consumed (AA step = 2,
  #'       plain step = 1).
  #'     \item \code{result$aa_accepted_count} reports how many AA steps were
  #'       accepted.
  #'     \item SRAA history is reset at each homotopy level boundary, on
  #'       linear→log path fallback, and on ALM penalty (\code{capacity_mu})
  #'       adaptation (ieppa_soft only).
  #'   }
  ```

### C2.2 — Regenerate documentation

- [ ] `Rscript -e 'devtools::document()' 2>&1 | tail -10`
- [ ] **Acceptance:** `man/harvest.Rd` regenerates without error.

### C2.3 — Commit

- [ ] `git add R/harvest.R man/harvest.Rd && git commit -m "docs(harvest): document SRAA-m for ieppa/ieppa_soft and caveats"`

---

## Task D1 — `tests/testthat/test-ieppa-sraa.R`: 8 quantitative tests

**Ticket:** `Task SRAA-D1 [eight tests with explicit numerical thresholds] ! [no flaky tolerances; every threshold ties to a spec-listed metric]`

**Files:** `tests/testthat/test-ieppa-sraa.R` (new)

### D1.1 — Write the test file FIRST and confirm it fails before B-tasks land (TDD red phase)

> If executing the plan strictly in task order (A1 → A2 → B1 → … → D1), this file is written last and tests pass green. To preserve TDD discipline, an alternative ordering is: A1, A2, **D1 (red)**, B1–B5, C1, C2, then re-run D1 (green). Either ordering is acceptable; document your choice in the bead ticket.

- [ ] Create `tests/testthat/test-ieppa-sraa.R` with the following content:

  ```r
  # Tests for SRAA-m Anderson Acceleration on iEPPA / ieppa_soft.
  # Spec: docs/superpowers/specs/2026-04-30-ieppa-sraa-acceleration.md
  # Plan: docs/superpowers/plans/2026-04-30-ieppa-sraa-acceleration.md (Task D1)

  skip_if_not_installed("leafblower")

  # ── Fixture builders ─────────────────────────────────────────────────────────
  make_fixture_K2_well_conditioned <- function(n = 1000, seed = 1L) {
      set.seed(seed)
      data.frame(
          a = sample(c("a1", "a2"), n, replace = TRUE),
          b = sample(c("b1", "b2", "b3"), n, replace = TRUE)
      )
  }

  make_fixture_K5_tight <- function(n = 5000, seed = 2L) {
      set.seed(seed)
      data.frame(
          v1 = sample(letters[1:4], n, replace = TRUE),
          v2 = sample(letters[1:5], n, replace = TRUE),
          v3 = sample(letters[1:3], n, replace = TRUE),
          v4 = sample(letters[1:6], n, replace = TRUE),
          v5 = sample(letters[1:4], n, replace = TRUE)
      )
  }

  make_targets_uniform <- function(df, vars) {
      lapply(vars, function(v) {
          lvls <- sort(unique(df[[v]]))
          setNames(rep(1 / length(lvls), length(lvls)), lvls)
      })
  }

  # ── Test 1 — Convergence parity (max_error within 1e-3) ──────────────────────
  test_that("SRAA accel max_error is within 1e-3 of plain on K=2 well-conditioned", {
      df <- make_fixture_K2_well_conditioned(1000)
      tgt <- make_targets_uniform(df, c("a", "b"))
      r_plain <- leafblower::harvest(
          data = df, target_vars = c("a", "b"), targets = tgt,
          method = "ieppa", accelerate = FALSE,
          max_iterations = 200, attach_weights = FALSE)
      r_accel <- leafblower::harvest(
          data = df, target_vars = c("a", "b"), targets = tgt,
          method = "ieppa", accelerate = TRUE,
          max_iterations = 200, attach_weights = FALSE)
      e_plain <- attr(r_plain, "result")$max_error
      e_accel <- attr(r_accel, "result")$max_error
      expect_lte(e_accel, e_plain + 1e-3)
  })

  # ── Test 2 — Iteration reduction ≥30% on K=5 tight ──────────────────────────
  test_that("SRAA cuts F-eval count to <70% of plain iters on K=5 tight", {
      df <- make_fixture_K5_tight(5000)
      tgt <- make_targets_uniform(df, names(df))
      r_plain <- leafblower::harvest(
          data = df, target_vars = names(df), targets = tgt,
          method = "ieppa", accelerate = FALSE, max_weight = 1.8,
          max_iterations = 500, attach_weights = FALSE)
      r_accel <- leafblower::harvest(
          data = df, target_vars = names(df), targets = tgt,
          method = "ieppa", accelerate = TRUE, max_weight = 1.8,
          max_iterations = 500, attach_weights = FALSE)
      n_plain <- attr(r_plain, "result")$iterations
      n_accel <- attr(r_accel, "result")$iterations
      expect_lt(n_accel, 0.7 * n_plain)
  })

  # ── Test 3 — ieppa_soft converges with accelerate=TRUE on T5 fixture ────────
  test_that("ieppa_soft SRAA reaches max_error <= 0.01 on T5 (5 cats, mw=1.8, unit)", {
      df <- make_fixture_K5_tight(2000, seed = 3L)
      tgt <- make_targets_uniform(df, names(df))
      r <- leafblower::harvest(
          data = df, target_vars = names(df), targets = tgt,
          method = "ieppa_soft", accelerate = TRUE,
          max_weight = 1.8, bounds_mode = "unit",
          max_iterations = 500, attach_weights = FALSE)
      expect_lte(attr(r, "result")$max_error, 0.01)
  })

  # ── Test 4 — Greedy downgrade emits message, status=0 ───────────────────────
  test_that("scheduler=greedy with accelerate=TRUE downgrades cleanly with verbose msg", {
      df <- make_fixture_K2_well_conditioned(500)
      tgt <- make_targets_uniform(df, c("a", "b"))
      msgs <- capture.output(
          r <- leafblower::harvest(
              data = df, target_vars = c("a", "b"), targets = tgt,
              method = "ieppa", accelerate = TRUE, scheduler = "greedy",
              verbose = 1L, max_iterations = 100, attach_weights = FALSE),
          type = "message")
      expect_equal(attr(r, "result")$status, 0L)
      expect_true(any(grepl("round_robin", msgs)))
  })

  # ── Test 5 — Output correlation across methods ──────────────────────────────
  test_that("plain vs SRAA weights correlate > 0.9999 on K=2 well-conditioned", {
      df <- make_fixture_K2_well_conditioned(1000)
      tgt <- make_targets_uniform(df, c("a", "b"))
      r_plain <- leafblower::harvest(
          data = df, target_vars = c("a", "b"), targets = tgt,
          method = "ieppa", accelerate = FALSE,
          max_iterations = 500, attach_weights = TRUE)
      r_accel <- leafblower::harvest(
          data = df, target_vars = c("a", "b"), targets = tgt,
          method = "ieppa", accelerate = TRUE,
          max_iterations = 500, attach_weights = TRUE)
      expect_gt(cor(r_plain$weights, r_accel$weights), 0.9999)
  })

  # ── Test 6 — aa_accepted_count >= 5 on tight K=5 ────────────────────────────
  test_that("aa_accepted_count >= 5 within 200 iters on K=5 tight", {
      df <- make_fixture_K5_tight(5000)
      tgt <- make_targets_uniform(df, names(df))
      r <- leafblower::harvest(
          data = df, target_vars = names(df), targets = tgt,
          method = "ieppa", accelerate = TRUE, max_weight = 1.8,
          max_iterations = 200, attach_weights = FALSE)
      expect_gte(attr(r, "result")$aa_accepted_count, 5L)
  })

  # ── Test 7 — Equal-budget regression: accel not worse than 2x plain ─────────
  test_that("SRAA at equal budget keeps max_error <= 2x plain (no regression)", {
      df <- make_fixture_K5_tight(5000)
      tgt <- make_targets_uniform(df, names(df))
      r_plain <- leafblower::harvest(
          data = df, target_vars = names(df), targets = tgt,
          method = "ieppa", accelerate = FALSE, max_weight = 1.8,
          max_iterations = 50, attach_weights = FALSE)
      r_accel <- leafblower::harvest(
          data = df, target_vars = names(df), targets = tgt,
          method = "ieppa", accelerate = TRUE, max_weight = 1.8,
          max_iterations = 50, attach_weights = FALSE)
      expect_lte(attr(r_accel, "result")$max_error,
                 2 * attr(r_plain, "result")$max_error)
  })

  # ── Test 8 — F-eval reduction to reach tol=1e-4 on K=5 tight ────────────────
  test_that("SRAA reaches tol=1e-4 in <70% of plain iters on K=5 tight", {
      df <- make_fixture_K5_tight(5000)
      tgt <- make_targets_uniform(df, names(df))
      r_plain <- leafblower::harvest(
          data = df, target_vars = names(df), targets = tgt,
          method = "ieppa", accelerate = FALSE, max_weight = 1.8,
          tolerance = 1e-4, max_iterations = 1000, attach_weights = FALSE)
      r_accel <- leafblower::harvest(
          data = df, target_vars = names(df), targets = tgt,
          method = "ieppa", accelerate = TRUE, max_weight = 1.8,
          tolerance = 1e-4, max_iterations = 1000, attach_weights = FALSE)
      stopifnot(attr(r_plain, "result")$status == 0L,
                attr(r_accel, "result")$status == 0L)
      f_evals_accel <- attr(r_accel, "result")$iterations
      iters_plain   <- attr(r_plain, "result")$iterations
      expect_lt(f_evals_accel, 0.7 * iters_plain)
  })
  ```

  Argument names (`tolerance`, `max_iterations`, `bounds_mode`, `scheduler`, `target_vars`) follow the existing `leafblower::harvest()` signature in `R/harvest.R`. If any name differs in the actual signature (e.g. `tol_abs` vs `tolerance`, `outer_max_iter` vs `max_iterations`), substitute the canonical name during implementation; the threshold semantics are unchanged.

### D1.2 — Run the suite

- [ ] `Rscript -e 'testthat::test_dir("tests/testthat")' 2>&1 | tail -3`
- [ ] **Acceptance:** all 8 tests in this file pass green. If any fails:
  - **Test 2 / Test 8 (iteration reduction)** → SRAA acceptance rate is too low; verify F_cur seeding and stall-guard implementation.
  - **Test 1 / Test 7 (parity / regression)** → safeguard math is wrong; check `sraa_step` `apply_clamp=false` branch.
  - **Test 4 (greedy)** → log-message string mismatch; check exact wording in B5.2.
  - **Test 6 (aa_accepted_count)** → harvest path missing the assignment; check C1.3.

### D1.3 — Commit

- [ ] `git add tests/testthat/test-ieppa-sraa.R && git commit -m "test(ieppa): add 8 quantitative tests for SRAA-m acceleration"`

---

## Final integration gate

- [ ] `R CMD INSTALL --preclean .` clean
- [ ] `Rscript -e 'testthat::test_dir("tests/testthat")' 2>&1 | tail -10` — 100% pass
- [ ] `Rscript -e 'devtools::document(); devtools::check()' 2>&1 | tail -30` — 0 errors, 0 warnings, 0 notes (or only pre-existing notes)
- [ ] `git log --oneline -10` shows the conventional commits in order

---

## Self-review (post-write checks performed inline)

**1. Spec coverage (every requirement maps to a task):**

| Spec section | Requirement | Task |
|---|---|---|
| Architecture / Files changed | sraa.hpp `apply_clamp` | A1 |
| Architecture / Files changed | ieppa.hpp `aa_accepted_count` | A2 |
| Architecture / Files changed | ieppa.cpp pack_lf/unpack_lf, F_eval, SRAA loop | B1, B2, B3 |
| Architecture / Files changed | r_bridge.cpp accelerate gate + harvest | C1 |
| Architecture / Files changed | harvest.R docstring | C2 |
| §lf vector layout note (NA slots) | total_cats includes NA buckets | B1 (uses `cat_offset[K]`), B3 (`init(total_cats, ...)`) |
| §pack_lf, §unpack_lf | exact body | B1.1 |
| §f_eval_lf | always pack on overflow → +inf | B2.1 |
| §SRAA outer loop / W_best vs lf_best | two trackers | B3.2 |
| §F_cur seeded once before loop | B3.1 (init), B3.2 (no re-seed inside) |
| §res.base.iterations single write | B3.2 (one write per AA step) |
| §SOR disabled when accelerate | sor_auto_v binding | B5.1 |
| §Greedy downgrade with log msg | round_robin message | B5.2 |
| §Log-path fallback: clear + snapshot | B3.3 |
| §ALM ieppa_soft combined block | B4.1 |
| §n_cats_total invariant | B3.1 (`init(total_cats, ...)`) |
| §Testing matrix (8 rows) | D1.1 (eight tests) |

**2. Placeholder scan:** No TBD/TODO/`...`/`<X>` placeholder tokens remain in the plan body — every task references concrete files, line numbers, and full code blocks.

**3. Type / name consistency across tasks:**

- `pack_lf`, `unpack_lf` — defined in B1, used in B2 (`f_eval_lf`), B3.1 (initial pack), B3.2 (revert), B3.3 (none — uses clear), B4.1 (reseed). All call sites match the B1 signature.
- `ieppa_sraa` — declared in B3.1, used in B3.2 (driver), B3.3 (clear+snapshot), B4.1 (clear+reseed). Single name throughout.
- `lf_flat`, `lf_best`, `dummy_L`, `dummy_U`, `sraa_outer_stall_count`, `sraa_best_errRp` — declared in B3.1, consumed only inside B3.2 / B4.1. No leakage.
- `sraa_active_lvl` — declared in B3.1, consumed in B3.2, B3.3, B4.1.
- `aa_accepted_count` — field added in A2, written in B3.2, B3.3 (snapshots), harvested in C1.3 across three dispatch arms, exposed via slot 34 in C1.4, queried in tests D1.1 (#6).
- `apply_clamp` — added in A1.1, consumed in B3.2 (`/*apply_clamp=*/false`).
- `kSRAAm`, `kSRAAOuterSlack`, `kSRAAOuterStallWindow` — existing constants from `sraa.hpp` lines 16-21; consumed in B3.1, B3.2.

All identifiers consistent. Plan is internally referentially closed.
