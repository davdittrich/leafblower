# Implementation Plan: Greenkhorn Logit — Epic B (Greenkhorn Solver)

- **Date:** 2026-04-29
- **Branch:** `fix/correctness-performance-2026-04-28`
- **Spec:** `docs/superpowers/specs/2026-04-29-greenkhorn-solver.md` (Part 1 algorithm section)
- **Scope:** Epic B — implement the Greenkhorn solver as a parallel artifact to the existing Raking solver, mirroring its result struct and lifecycle.

---

## Plan Header

- **Mechanism:** Greedy coordinate-descent IPF (Greenkhorn) over the cell table `CellTable`, with per-margin residual tracking, incremental total mass `W`, bucket-sum table `S_flat`, capacity clamping against `[L_cell, U_cell]`, best-iterate tracking, and convergence dispatch via `lbw::compute_cell_metrics` + `lbw::check_convergence`. Result struct `GreenkornResult` mirrors `lbw::RakingResult` field-for-field.
- **Forbidden:**
  - Forbidden: rebuilding `S_flat` from scratch each step (must be incremental, except per-margin re-sum on the chosen margin to absorb clamp drift).
  - Forbidden: full O(n) sweep over observations inside the Greenkhorn loop (work is over cells, not obs).
  - Forbidden: introducing a new `RakingResult`-like definition that drifts from `src/raking.hpp` — fields must be copied verbatim.
  - Forbidden: skipping the `R CMD INSTALL --preclean .` compile gate after each task.
  - Forbidden: `--no-verify` on commits, `--amend`, or bundling B1 + B2 in one commit.
  - Forbidden: silently picking a single interpretation if the spec is ambiguous; halt and ask.
- **Audit:** No tests are written in this epic. Compile gates only (`R CMD INSTALL --preclean .` → `* DONE`). Behavioral validation is owned by later epics (E: tests, G: benchmarks). The compile gate is necessary but not sufficient — completion of B2 requires the install to succeed, but does NOT certify correctness.

---

## Pre-flight

Before starting, confirm:

1. Working tree is on `fix/correctness-performance-2026-04-28`:
   ```bash
   git rev-parse --abbrev-ref HEAD
   ```
   Must print `fix/correctness-performance-2026-04-28`.
2. Reference files exist and are unchanged in this session:
   - `src/raking.hpp` (RakingResult definition)
   - `src/raking.cpp` lines 75-120 (cells_per_cat construction pattern)
   - `src/calib_dispatch.hpp` lines 140-180 (`compute_cell_metrics`, `check_convergence` signatures)
   - `src/cell_table.hpp` (`build_cell_table`, `CellTable`)
   - `src/types.hpp` (`CalibState`, `CalibConvergenceCfg`)
   - `src/leafblower.h` (status codes `RK_OK`, `RK_ERR_NOCONV`, `RK_ERR_INFEAS`, `RK_ERR_BUDGET`, `RK_ERR_STALL`)
3. No uncommitted edits in `src/greenkhorn.{hpp,cpp}` (these files do not exist yet).

If any of the above fails, **halt** and report.

---

## Task B1 — `src/greenkhorn.hpp`

**Goal:** Declare `GreenkornResult` (mirroring `RakingResult`) and the `greenkhorn_solve` entry point.

### Step B1.1 — Re-read `src/raking.hpp`

Read `src/raking.hpp` and copy `RakingResult` fields verbatim into a new struct `GreenkornResult`. The two structs must have identical fields, identical default initializers, identical declaration order, and identical comments. Only the type name differs.

### Step B1.2 — Write `src/greenkhorn.hpp`

Create the file with exactly this content:

```cpp
#pragma once
#include "lbw_config.h"
#include "types.hpp"
#include "leafblower.h"
#include <limits>
#include <vector>

namespace lbw {

struct GreenkornResult {
    int status;
    int iterations;
    double max_error;
    // ── Extended quality metrics (WU-A scaffold; populated in WU-B+) ──
    double mean_error          = 0.0;
    double kl                  = 0.0;
    double chi2                = 0.0;
    double l1_weight_change    = 0.0;  // WU-A: renamed from pct_change; computation in WU-B
    double grake_norm          = 0.0;  // WU-A stub; computation in WU-D
    int    convergence_metric  = 0;    // WU-A stub; CalibMetric at exit
    int    convergence_rule    = 1;    // WU-A stub; CalibRule at exit (IMPROVEMENT)
    double convergence_tol     = 0.001; // WU-A stub; threshold that fired
    int    convergence_iter    = -1;   // WU-A stub; iteration at convergence (-1=max_iter)
    double best_error          = std::numeric_limits<double>::infinity();
    int    best_iter           = 0;
    std::vector<double> best_weights;  // obs-level; length n; sum-normalized to n; empty if never checked
    double best_objective_seen          = 0.0;   // internal: weight KL at best_iter
    double convergence_solver_objective = 0.0;   // exposed: solver's mathematical objective
    int    convergence_minimized_metric = 0;     // CalibMetric: which metric was minimized
    // ── End extended quality metrics ──

    char message[256] = {0};  // status string (see greenkhorn.cpp)
};

GreenkornResult greenkhorn_solve(CalibState& st);

}  // namespace lbw
```

> **Note on `message[256]`:** the spec block in this plan's directives uses `res.message` via `std::snprintf(res.message, sizeof(res.message), ...)`. Verify that `RakingResult` does or does not carry an analogous `message` field. If `RakingResult` does NOT have it, `GreenkornResult` adds a `message[256]` buffer (above) — this is the only authorized deviation from "verbatim copy of fields", and is required by the B2 implementation's `snprintf` calls. If `RakingResult` does have a `message` field, copy it verbatim and drop the addendum above.

### Step B1.3 — Compile gate (B1)

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Required tail: `* DONE (leafblower)`. If anything else, **halt** and diagnose; do not proceed to B2.

### Step B1.4 — Commit (B1)

```bash
git add src/greenkhorn.hpp
git commit -m "feat(greenkhorn): add greenkhorn.hpp struct definition mirroring RakingResult"
```

No `--amend`, no `--no-verify`. If pre-commit hooks fail, fix the underlying issue and create a NEW commit.

---

## Task B2 — `src/greenkhorn.cpp`

**Goal:** Implement `greenkhorn_solve` per the algorithm in the spec — greedy coordinate-descent IPF with incremental `W`, bucket-sum table `S_flat`, capacity clamping, periodic convergence dispatch, and best-iterate tracking.

### Step B2.1 — Pre-read

Re-read (only if not yet read this session, per OpenWolf token discipline):
- `src/raking.cpp` lines 75-120 — `cells_per_cat` construction pattern. The B2 code below mirrors it exactly.
- `src/calib_dispatch.hpp` lines 140-180 — exact signatures of `compute_cell_metrics(CalibState&, CellTable&, std::vector<double>& X, double W, std::vector<double>& bucket_scratch)` and `check_convergence(CalibConvergenceCfg&, CellMetrics&, double& prev_metric, double tol_abs)`. Adjust the call sites in B2.4 step (5) if the actual signatures differ from what the spec block assumes.
- `src/cell_table.hpp` — confirm `CellTable` exposes `cell_of`, `n_per_cell`, `g_per_cell` (vector-of-vectors indexed `[k][c]`), and `M_cell`. Confirm `build_cell_table` signature matches the call site in B2.3.

If any signature differs, **halt** and surface the divergence before writing code.

### Step B2.2 — Write `src/greenkhorn.cpp` — file scaffold and includes

Create the file beginning with:

```cpp
#include "lbw_config.h"
#include "greenkhorn.hpp"
#include "calib_dispatch.hpp"
#include "cell_table.hpp"
#include "calib_validate.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <numeric>

namespace lbw {

GreenkornResult greenkhorn_solve(CalibState& st) {
    GreenkornResult res;
    res.status = RK_ERR_NOCONV;
```

> Note: `<cstdio>` is added beyond the spec block to make `std::snprintf` available portably.

### Step B2.3 — Cell-table build + `cells_per_cat`

Append to `greenkhorn_solve`:

```cpp
    // Build cell table
    CellTable ct;
    {
        std::vector<const int32_t*> gids(st.K);
        for (int k = 0; k < st.K; k++) gids[k] = st.group_ids[k];
        if (build_cell_table(st.n, st.K, gids.data(), st.cat_counts, st.weights, ct) != RK_OK) {
            res.status = RK_ERR_INFEAS;
            std::snprintf(res.message, sizeof(res.message),
                "greenkhorn: cell table build failed (structural infeasibility)");
            return res;
        }
    }

    const int M = ct.M_cell;
    const int K = st.K;

    // Build local bucket-to-cell lists from g_per_cell
    int max_cats = 0;
    for (int k = 0; k < K; k++) max_cats = std::max(max_cats, st.cat_counts[k]);

    std::vector<std::vector<std::vector<int>>> cells_per_cat(K);
    for (int k = 0; k < K; k++) {
        cells_per_cat[k].assign(st.cat_counts[k], {});
        for (int c = 0; c < M; c++) {
            int g = ct.g_per_cell[k][c];
            if (g >= 0 && g < st.cat_counts[k])
                cells_per_cat[k][g].push_back(c);
        }
    }
```

This mirrors `src/raking.cpp` lines 75-120 exactly.

### Step B2.4 — State init: `X`, bounds, `W`, `S_flat`, `errRp`

Append:

```cpp
    // Cell masses from design weights
    std::vector<double> X(M, 0.0);
    for (int i = 0; i < st.n; i++) X[ct.cell_of[i]] += st.weights[i];
    const std::vector<double> X_init = X;

    // Capacity bounds
    std::vector<double> L_cell(M), U_cell(M);
    for (int c = 0; c < M; c++) {
        L_cell[c] = st.min_weight * ct.n_per_cell[c];
        U_cell[c] = st.max_weight * ct.n_per_cell[c];
        X[c] = std::clamp(X[c], L_cell[c], U_cell[c]);
    }

    // Total mass W (maintained incrementally)
    double W = 0.0;
    for (int c = 0; c < M; c++) W += X[c];

    // Bucket sums S_flat[k * max_cats + j]
    const int S_stride = max_cats;
    std::vector<double> S_flat(K * S_stride, 0.0);
    for (int k = 0; k < K; k++)
        for (int j = 0; j < st.cat_counts[k]; j++)
            for (int c : cells_per_cat[k][j]) S_flat[k*S_stride+j] += X[c];

    // Per-margin residuals
    std::vector<double> errRp(K, 0.0);
    auto compute_errRp_k = [&](int k) -> double {
        double e = 0.0;
        for (int j = 0; j < st.cat_counts[k]; j++) {
            double achieved = (W > 0.0) ? S_flat[k*S_stride+j] / W : 0.0;
            e = std::max(e, std::abs(achieved - st.targets[k][j]));
        }
        return e;
    };
    for (int k = 0; k < K; k++) errRp[k] = compute_errRp_k(k);

    // Best-iterate tracking
    double best_errRp = *std::max_element(errRp.begin(), errRp.end());
    res.best_error = best_errRp;
    res.best_iter  = 0;
    std::vector<double> X_best = X;

    // Scratch for compute_cell_metrics
    std::vector<double> bucket_scratch(max_cats, 0.0);
    double prev_metric = std::numeric_limits<double>::infinity();
    const CalibConvergenceCfg& cfg = st.convergence_cfg;
    constexpr int kErrCheckInterval = 10;
    constexpr double kEmptyBucketThreshold = 1e-15;
```

### Step B2.5 — Greenkhorn main loop

Append:

```cpp
    for (int iter = 0; iter < st.inner_max_iter; iter++) {
        // W <= 0 guard
        if (W <= 0.0) {
            res.status = RK_ERR_INFEAS;
            std::snprintf(res.message, sizeof(res.message),
                "greenkhorn: total mass W<=0 (all cells at zero bound)");
            break;
        }

        // (1) Pick k* = argmax errRp
        int k_star = (int)(std::max_element(errRp.begin(), errRp.end()) - errRp.begin());

        // (2) Single-margin IPF update for margin k*
        for (int j = 0; j < st.cat_counts[k_star]; j++) {
            double S_kj = S_flat[k_star * S_stride + j];
            if (S_kj < kEmptyBucketThreshold * W) continue;
            double target_mass = st.targets[k_star][j] * W;
            double f = target_mass / S_kj;

            for (int c : cells_per_cat[k_star][j]) {
                double X_old = X[c];
                double X_new = std::clamp(X[c] * f, L_cell[c], U_cell[c]);
                double delta = X_new - X_old;
                X[c] = X_new;
                W += delta;

                // Incremental bucket-sum updates for other margins
                if (std::abs(delta) < 1e-300) continue;
                for (int k2 = 0; k2 < K; k2++) {
                    if (k2 == k_star) continue;
                    int g2 = ct.g_per_cell[k2][c];
                    if (g2 >= 0 && g2 < st.cat_counts[k2])
                        S_flat[k2 * S_stride + g2] += delta;
                }
            }
            // Recompute S[k*][j] exactly (clamp may cause drift from f-scaling)
            S_flat[k_star * S_stride + j] = 0.0;
            for (int c : cells_per_cat[k_star][j]) S_flat[k_star * S_stride + j] += X[c];
        }

        // (3) Recompute errRp for all margins
        for (int k = 0; k < K; k++) errRp[k] = compute_errRp_k(k);

        // (4) Best-iterate tracking
        double curr_max = *std::max_element(errRp.begin(), errRp.end());
        if (curr_max < best_errRp) {
            best_errRp = curr_max;
            res.best_error = best_errRp;
            res.best_iter  = iter + 1;
            X_best = X;
        }

        // (5) Convergence check every kErrCheckInterval steps
        res.iterations = iter + 1;
        if ((iter + 1) % kErrCheckInterval == 0 || iter == st.inner_max_iter - 1) {
            lbw::CellMetrics m = lbw::compute_cell_metrics(st, ct, X, W, bucket_scratch);
            bool converged = lbw::check_convergence(cfg, m, prev_metric, st.tol_abs);
            if (converged) {
                res.status = RK_OK;
                res.convergence_iter = iter + 1;
                X_best = X;  // final is best
                break;
            }
        }
    }
```

### Step B2.6 — Post-loop classification + weight reconstruction

Append:

```cpp
    // Post-loop status classification
    if (res.status == RK_ERR_NOCONV) {
        double final_errRp = *std::max_element(errRp.begin(), errRp.end());
        res.status = (final_errRp < prev_metric * 0.999) ? RK_ERR_BUDGET : RK_ERR_STALL;
        std::snprintf(res.message, sizeof(res.message),
            "greenkhorn: %s after %d steps; best max_err=%.4e",
            res.status == RK_ERR_BUDGET ? "budget exhausted" : "stall",
            res.iterations, res.best_error);
    }

    // Weight reconstruction from X_best
    res.best_weights.resize(st.n);
    for (int i = 0; i < st.n; i++) {
        int c = ct.cell_of[i];
        res.best_weights[i] = (X_init[c] > 0.0)
            ? st.weights[i] * X_best[c] / X_init[c]
            : st.weights[i];
    }
    res.max_error = best_errRp;

    return res;
}

}  // namespace lbw
```

### Step B2.7 — Compile gate (B2)

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Required tail: `* DONE (leafblower)`.

If the install fails:
1. **Halt** — do not commit.
2. Read the actual error.
3. Most likely failure modes (in order):
   - Signature mismatch on `compute_cell_metrics` / `check_convergence` (read `src/calib_dispatch.hpp` and adjust the call site).
   - `CalibState` field-name drift (`group_ids`, `cat_counts`, `targets`, `weights`, `min_weight`, `max_weight`, `inner_max_iter`, `convergence_cfg`, `tol_abs`, `n`, `K`) — re-read `src/types.hpp`.
   - `CellTable` field-name drift (`cell_of`, `n_per_cell`, `g_per_cell`, `M_cell`) — re-read `src/cell_table.hpp`.
   - `build_cell_table` parameter list — re-read `src/cell_table.hpp`.
   - `RK_OK` / `RK_ERR_*` codes — re-read `src/leafblower.h`.
   - Missing `message[256]` field on `GreenkornResult` — see B1.2 note; add if needed and re-run B1 compile gate.
4. Fix the root cause (per CLAUDE.md §3 "diagnose root cause before retrying"), do **not** patch over the symptom.
5. Re-run the compile gate.

### Step B2.8 — Commit (B2)

```bash
git add src/greenkhorn.cpp
git commit -m "feat(greenkhorn): implement greenkhorn_solve — greedy coordinate-descent IPF with W-incremental and best-iterate tracking"
```

If the B2 step had to amend `greenkhorn.hpp` (e.g., the `message` field decision flipped), `git add src/greenkhorn.hpp` as well in this commit and adjust the message to `"feat(greenkhorn): implement greenkhorn_solve and finalize header"`.

---

## Done Criteria

Epic B is complete when ALL of the following hold:

1. `src/greenkhorn.hpp` exists, declares `lbw::GreenkornResult` (field-equivalent to `lbw::RakingResult`) and `lbw::greenkhorn_solve(CalibState&)`.
2. `src/greenkhorn.cpp` exists, defines `greenkhorn_solve` per B2.3-B2.6.
3. `R CMD INSTALL --preclean .` ends with `* DONE (leafblower)` after B1 and again after B2.
4. Two commits exist on `fix/correctness-performance-2026-04-28`:
   - `feat(greenkhorn): add greenkhorn.hpp struct definition mirroring RakingResult`
   - `feat(greenkhorn): implement greenkhorn_solve — greedy coordinate-descent IPF with W-incremental and best-iterate tracking`
5. No tests are added in this epic (deferred to Epic E).
6. No call site is wired up to `greenkhorn_solve` in this epic (deferred to the dispatch epic). The function is unreferenced; this is expected and not a defect.
7. OpenWolf bookkeeping:
   - `.wolf/anatomy.md` updated with entries for `src/greenkhorn.hpp` and `src/greenkhorn.cpp`.
   - `.wolf/memory.md` appended with B1 and B2 outcome lines.

## Out of Scope (explicit)

- Tests for `greenkhorn_solve` (Epic E).
- R-level dispatch / user-facing API (later epic).
- Benchmarks vs. raking (Epic G).
- Logit-domain variant (separate epic).
- SIMD or threading inside the inner loop.
- Reformatting or refactoring `raking.cpp` / `raking.hpp` to share code with `greenkhorn.*`. Per CLAUDE.md §2 "Surgical Changes", duplication is acceptable here; deduplication is a separate ticket if/when justified.

## Halt Codes

If any of the following occur, output `SPEC_FAILURE` and stop:
- `R CMD INSTALL --preclean .` cannot be made to print `* DONE (leafblower)` after a single root-cause diagnosis pass.
- `compute_cell_metrics` or `check_convergence` does not exist in `src/calib_dispatch.hpp` (the spec assumed they do).
- `RakingResult` field set materially diverges from what is copied into `GreenkornResult` (e.g., new fields added between spec authoring and execution).
- `build_cell_table` returns a status code other than `RK_OK` on the project's existing test inputs (would indicate environmental breakage unrelated to this epic).

---

## Amendment: SQUAREM Acceleration (in scope per user confirmation)

**Background**: autumn's greedy rake uses SQUAREM with greedy scheduling. This is valid because autumn sorts margins **once at round entry** — F_eval runs K steps in that fixed order, making F_eval stationary. Contrast with leafblower's R8 fix (greedy disabled under raking SQUAREM) where greedy re-sorted inside F_eval per call.

**When `st.accelerate == true`**: implement SQUAREM where **F_eval = K greedy Greenkhorn steps, margins sorted once at F_eval entry**.

### B2 Amendment — SQUAREM F_eval

Define F_eval as a lambda inside greenkhorn_solve:

```cpp
auto F_eval = [&](const std::vector<double>& X_in,
                  std::vector<double>& S_in,
                  double& W_in) {
    // Sort margins by errRp ONCE at round entry (stationary F_eval guarantee)
    std::vector<int> order(K);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
        [&](int a, int b) { return errRp[a] > errRp[b]; });

    // K greedy steps in this fixed order
    for (int ki = 0; ki < K; ki++) {
        int k_step = order[ki];
        if (W_in <= 0.0) break;  // W<=0 guard
        for (int j = 0; j < st.cat_counts[k_step]; j++) {
            double S_kj = S_in[k_step * S_stride + j];
            if (S_kj < kEmptyBucketThreshold * W_in) continue;
            double f = st.targets[k_step][j] * W_in / S_kj;
            for (int c : cells_per_cat[k_step][j]) {
                double X_old = X_in[c];
                double X_new = std::clamp(X_in[c] * f, L_cell[c], U_cell[c]);
                double delta = X_new - X_old;
                const_cast<std::vector<double>&>(X_in)[c] = X_new;
                W_in += delta;
                if (std::abs(delta) < 1e-300) continue;
                for (int k2 = 0; k2 < K; k2++) {
                    if (k2 == k_step) continue;
                    int g2 = ct.g_per_cell[k2][c];
                    if (g2 >= 0 && g2 < st.cat_counts[k2]) S_in[k2*S_stride+g2] += delta;
                }
            }
            S_in[k_step*S_stride+j] = 0.0;
            for (int c : cells_per_cat[k_step][j]) S_in[k_step*S_stride+j] += X_in[c];
        }
        for (int k = 0; k < K; k++) errRp[k] = compute_errRp_k_from(X_in, S_in, W_in, k);
    }
};
```

**Note**: F_eval modifies X, S, W in-place via copies in the SQUAREM super-step. Pre-allocate scratch buffers before the main loop:

```cpp
std::vector<double> sq_w1, sq_w2, sq_X_snap;
std::vector<double> sq_S1(K*S_stride), sq_S2(K*S_stride), sq_Ssnap(K*S_stride);
double sq_W1, sq_W2, sq_Wsnap;
if (st.accelerate) {
    sq_w1.assign(M, 0.0); sq_w2.assign(M, 0.0); sq_X_snap.assign(M, 0.0);
}
```

### B2 Amendment — SQUAREM super-step (mirrors raking.cpp lines 351-480)

Inside the main loop, after every K basic Greenkhorn steps (one "round"), when `st.accelerate == true`:

```cpp
if (st.accelerate && (iter % K == K-1)) {  // End of each K-step round
    // w1 = F_eval(X)
    sq_X_snap = X; sq_Ssnap = S_flat; sq_Wsnap = W;
    sq_w1 = X; sq_S1 = S_flat; sq_W1 = W;
    F_eval(sq_w1, sq_S1, sq_W1);

    // w2 = F_eval(w1)
    sq_w2 = sq_w1; sq_S2 = sq_S1; sq_W2 = sq_W1;
    F_eval(sq_w2, sq_S2, sq_W2);

    // CBB alpha = -||r||_obs / ||v||_obs (obs-level norms)
    double r_sq = 0.0, v_sq = 0.0;
    for (int c = 0; c < M; c++) {
        double inv_nc = 1.0 / ct.n_per_cell[c];
        double r_c = sq_w1[c] - sq_X_snap[c];
        double v_c = sq_w2[c] - 2.0*sq_w1[c] + sq_X_snap[c];
        r_sq += r_c * r_c * inv_nc;
        v_sq += v_c * v_c * inv_nc;
    }
    double alpha = (v_sq > 0.0) ? -std::sqrt(r_sq / v_sq) : -1.0;
    alpha = std::max(alpha, -4.0);  // cap

    // Extrapolate: X_star = X_snap - 2*alpha*r + alpha^2*v
    std::vector<double> sq_X_star = sq_X_snap;
    for (int c = 0; c < M; c++) {
        double r_c = sq_w1[c] - sq_X_snap[c];
        double v_c = sq_w2[c] - 2.0*sq_w1[c] + sq_X_snap[c];
        sq_X_star[c] = std::clamp(sq_X_snap[c] - 2.0*alpha*r_c + alpha*alpha*v_c,
                                   L_cell[c], U_cell[c]);
    }
    // Recompute W and S_flat for X_star
    std::vector<double> sq_Sstar(K*S_stride, 0.0);
    double sq_Wstar = 0.0;
    for (int c = 0; c < M; c++) sq_Wstar += sq_X_star[c];
    // rebuild S_star from scratch (accept the O(M*K) cost at super-step)
    for (int k = 0; k < K; k++)
        for (int j = 0; j < st.cat_counts[k]; j++)
            for (int c : cells_per_cat[k][j]) sq_Sstar[k*S_stride+j] += sq_X_star[c];

    // Stabilization: X_star = F_eval(X_star)
    F_eval(sq_X_star, sq_Sstar, sq_Wstar);

    // Halving: if X_star is worse than w2, fall back
    double err_star = 0.0;
    for (int k = 0; k < K; k++) {
        for (int j = 0; j < st.cat_counts[k]; j++) {
            double ach = (sq_Wstar > 0.0) ? sq_Sstar[k*S_stride+j] / sq_Wstar : 0.0;
            err_star = std::max(err_star, std::abs(ach - st.targets[k][j]));
        }
    }
    double err_w2 = *std::max_element(errRp.begin(), errRp.end());
    bool accept_star = (err_star <= err_w2 * 1.01);  // 1% slack for numerical noise

    if (accept_star) {
        X = sq_X_star; S_flat = sq_Sstar; W = sq_Wstar;
    } else {
        X = sq_w2; S_flat = sq_S2; W = sq_W2;  // fall back to w2
    }
    // Update errRp from accepted state
    for (int k = 0; k < K; k++) errRp[k] = compute_errRp_k(k);
}
```

### Implementation order

B2 base (pure Greenkhorn, accelerate=false) first → compile gate → then add SQUAREM block → compile gate → commit single file. Do NOT split into two commits for this amendment.

### Test addition (amend E2 ticket)

Add T_acc to E2: verify `accelerate=TRUE` with greenkhorn works without error and produces reasonable output:

```r
test_that("Tacc: greenkhorn with accelerate=TRUE runs without error", {
  set.seed(99); n <- 2000L
  df  <- data.frame(x=factor(sample(letters[1:3],n,TRUE)),
                    y=factor(sample(c("M","F"),n,TRUE)))
  tgt <- list(x=c(a=0.3,b=0.4,c=0.3), y=c(M=0.5,F=0.5))
  r_acc <- harvest(df, tgt, method="greenkhorn", accelerate=TRUE, max_iterations=500L)
  r_plain <- harvest(df, tgt, method="greenkhorn", accelerate=FALSE, max_iterations=500L)
  # Both must converge to similar quality
  me_acc   <- attr(r_acc,   "result")$max_error
  me_plain <- attr(r_plain, "result")$max_error
  expect_lt(me_acc, 1e-3)
  # Accelerated should use fewer rounds (or same)
  iters_acc   <- attr(r_acc,   "result")$iterations
  iters_plain <- attr(r_plain, "result")$iterations
  # Not strict (problem-dependent) but shouldn't massively regress
  expect_lt(iters_acc, iters_plain * 2L)
})
```
