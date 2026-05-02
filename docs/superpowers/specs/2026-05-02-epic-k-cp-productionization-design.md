# Epic-K: stepstone-CP Productionization — Design (rev 2)

**Date:** 2026-05-02
**Status:** Design — pre-plan, gate iteration 2
**Predecessor:** Epic-J (`leafblower-y2ls`) FAIL verdict on kk1204 + side-finding CP wins stepstone (parity 0.45 vs `ieppa+sraa` baseline 1.13e-4)
**Investigation report:** `docs/investigations/2026-05-02-ylsy-cp-ipm-spike-result.md` Sec 2 + Sec 6 (Side-finding)
**Revision history:**
- rev 1 → rev 2: closes 8 blockers from design-review-gate iter 1 (Architect: dispatch line refs, harvest.R line refs, accelerate default mechanism. Designer: same accelerate mechanism, alg_name ternary chain, kAlgMap entry, result field SEXP-pack indices, status_code R-side semantics).
- rev 2 → rev 3: closes 4 blockers from design-review-gate iter 2 (Designer: warning EMIT SITE in harvest.R not C++; status_code=3 returns init-clamped weights not all-1.0. CTO: rate-fit exponents are spike-only diagnostic, not production-exposed; production CP stores NO per-iter trace — diagnostics are final-iteration scalars only). Folds Architect/Designer/Security/CTO/PM suggestions where high-value (T6 threshold tighten, T7 fallback test, NEWS.md bullet draft, K-7 reviewer set, R16 comment header step, include-guard rename).

## 1. Problem & Goals

### Problem

Epic-J spike confirmed Chambolle-Pock primal-dual (CP) outperforms `ieppa+sraa` on stepstone_K9 by 2.2× (5.083e-05 vs 4.388e-04) with last-iterate rate exponent $\beta = -1.05$, $R^2 = 0.96$. CP currently lives in `research/cp_calib.{hpp,cpp}` under `.Rbuildignore` — analyst-inaccessible.

### Goal

Ship `harvest(method="cp", ...)` as production solver:
- Opt-in (no AUTO routing change).
- Algorithm 1 (vanilla PDHG, O(1/k)) and Algorithm 2 (accelerated PDHG, O(1/k²)) per Chambolle-Pock 2011.
- `accelerate=TRUE` default; auto-fallback to Alg 1 when `u_max = Inf` (γ-strong-convexity precondition violated).
- Cell-compressed by default; obs-level fallback when `M_cell/n > 0.9` or `bounds_mode = "unit"`.
- Three user-facing knobs: `accelerate`, `max_iterations`, `convergence` (standard `harvest()` API).
- `cp_safety_factor = 1.05` hardcoded (no user override; deferred to Epic-K.2).

### Out of scope (deferred)

- AUTO routing carve-out (Epic-K.2, requires fixture-class sweep).
- Python bindings, vignette, SIMD optimization, warm-start, predictor-corrector PDHG.

### Decision criteria for Epic-K close

- All 6 tests in `tests/testthat/test-cp.R` PASS.
- `R CMD INSTALL --preclean .` clean.
- `harvest.Rd` + `NEWS.md` updated.
- AUTO routing untouched (existing `test-algo-selection.R` regression PASS).

## 2. Architecture & ABI

### File layout

```
src/
├── cp_calib.hpp                      ← move from research/, refactor to newton_calib.hpp pattern
├── cp_calib.cpp                      ← move from research/, extend with Alg 2 + cell compression
├── leafblower.h                      ← add RK_ALG_CP = 12 to rk_algorithm_t enum
└── r_bridge.cpp                      ← see "r_bridge.cpp edit points" below

R/
└── harvest.R                         ← see "harvest.R edit points" below

tools/
└── check_research_isolation.R        ← UPDATE forbidden-symbol list: REMOVE cp_solve_R, cp_calibrate
                                        (now legitimately in src/leafblower.so post-K-1).
                                        Keep ipm_solve_R, ipm_calibrate (still research-only — Epic-J FAIL).

tests/testthat/
└── test-cp.R                         ← NEW (6 tests: T1-T6)

man/harvest.Rd                        ← regenerate via devtools::document()
NEWS.md                               ← additive entry under "## New features"

research/                             ← stays under .Rbuildignore (Epic-J FAIL artefact)
├── cp_calib.{hpp,cpp}                ← STAYS as fossil; ADD header comment "MOVED TO src/cp_calib.{hpp,cpp};
                                        this copy retained for Epic-J spike traceability only — do not edit"
└── ipm_calib.{hpp,cpp}               ← STAYS, ABANDONED (Epic-J FAIL)
```

### r_bridge.cpp edit points (verified line numbers, master @ 7ad92f7)

| Line | Current state | Edit |
|---|---|---|
| 25–37 | `kAlgMap` table missing `cp` | Add `{"cp", RK_ALG_CP}` after `{"newton_kl", RK_ALG_NEWTON_KL}` (line 36) |
| 599–617 | `newton_kl` dispatch arm (existing) | NO CHANGE |
| 618 | catch-all `else { if (chebyshev) ... }` block start | Insert NEW `else if (strcmp(method_str, "cp") == 0) { ... }` arm immediately BEFORE this line |
| 709–718 | `alg_name` ternary chain | Insert `(res_alg_used == static_cast<int>(RK_ALG_CP)) ? "cp"` before the fallback `: "iEPPA"` |
| 819–822 | NewtonCalibResult SEXP-pack precedent (`lm_mu_final`, `n_projected_dims` at fixed slots) | Add CP-specific slots — see "Result SEXP-pack" below |

### harvest.R edit points (verified line numbers, master @ 7ad92f7)

| Line | Current state | Edit |
|---|---|---|
| ~16–30 | `@param method` docstring (lists 11 methods) | Add `\code{"cp"} (Chambolle-Pock primal-dual; for moderate-skew K≥5 problems — outperforms ieppa+sraa on stepstone-class fixtures by ~2x. Cell-compressed by default; obs-level fallback when M_cell/n > 0.9 or bounds_mode='unit'. Supports \code{accelerate=TRUE} for accelerated PDHG)`. Insertion position: adjacent to `\code{"chebyshev"}` (primal-dual splitting paradigm). |
| 236 | `accelerate = FALSE` default in signature | NO CHANGE to signature default |
| 319–322 | `if (isTRUE(accelerate) && !method %in% c("raking", "greenkhorn", "ieppa", "ieppa_soft"))` warning + `accelerate_bool <- isTRUE(accelerate) && method %in% c(...)` | Add `"cp"` to BOTH method whitelists (warning predicate + accelerate_bool predicate). |
| ~318 (NEW LINE) | — | INSERT BEFORE line 319: `accelerate_explicit <- !missing(accelerate)`. THEN INSERT BEFORE line 322: `if (method == "cp" && !accelerate_explicit) accelerate <- TRUE`. This makes CP default to `accelerate=TRUE` while preserving the global `accelerate=FALSE` signature default for all OTHER methods (no breaking change). |
| 590 | `match.arg(method, c("auto", "ieppa", ..., "newton_kl"))` (11 entries) | Add `"cp"` after `"newton_kl"` → 12-entry whitelist. |

### ABI changes (single-line)

`src/leafblower.h` `rk_algorithm_t`:
```cpp
enum rk_algorithm_t {
  RK_ALG_AUTO = 0, RK_ALG_IEPPA = 1, RK_ALG_LBFGSB = 2,
  RK_ALG_RAKING = 3, RK_ALG_SINKHORN = 4, RK_ALG_CHEBYSHEV = 5,
  RK_ALG_GREG = 6, /* 7 reserved */ RK_ALG_IEPPA_SOFT = 8,
  RK_ALG_GREENKHORN = 9, RK_ALG_LOGIT = 10, RK_ALG_NEWTON_KL = 11,
  RK_ALG_CP = 12        // NEW (Epic-K)
};
```

**No `rk_params_t` field changes.** `cp_safety_factor` hardcoded; `accelerate` reuses existing `st.accelerate` field. `EXPECTED_RK_PARAMS_BYTES` unchanged.

### Result struct

`CpCalibResult` in `cp_calib.hpp` (POD-compatible — no virtual, no inheritance other than `CalibResult` aggregate; matches `NewtonCalibResult` precedent):
- `CalibResult base` (status, iterations, max_error, weights, best_weights)
- `int n_cells` — `M_cell` if cell-compressed, `n` if obs-level
- `std::string algorithm_requested` — `"pdhg"` or `"accelerated_pdhg"` (user-asked, before fallback)
- `std::string algorithm_used` — `"pdhg"` or `"accelerated_pdhg"` (actual, after possible fallback)
- `double A_norm_estimate` — power-iter `‖A‖`
- `int n_power_iter` — power-iter convergence count
- `double final_theta` — accelerated_pdhg only, last θ_k (NaN if pdhg)
- `double final_tau`, `final_sigma` — accelerated_pdhg only (NaN if pdhg)
- `bool fell_back_to_pdhg` — true iff `algorithm_requested != algorithm_used`

### Result SEXP-pack (R-side surfacing)

Mirror NewtonCalibResult precedent (Epic-Dβ + Epic-H WH-d): r_bridge.cpp surfaces CP-specific fields by adding new slots to the result list. Slot indices append after existing newton_kl slots (which currently end at slot 36 = `lm_mu_final` per Epic-H WH-e). New CP slots:

| Slot | Field | Type |
|---|---|---|
| 37 | `n_cells` | INTSXP scalar |
| 38 | `algorithm_requested` | STRSXP scalar |
| 39 | `algorithm_used` | STRSXP scalar |
| 40 | `A_norm_estimate` | REALSXP scalar |
| 41 | `n_power_iter` | INTSXP scalar |
| 42 | `final_theta` | REALSXP scalar (NA_real_ if pdhg) |
| 43 | `final_tau` | REALSXP scalar (NA_real_ if pdhg) |
| 44 | `final_sigma` | REALSXP scalar (NA_real_ if pdhg) |
| 45 | `fell_back_to_pdhg` | LGLSXP scalar |

Result list size grows from 37 elements (slots 0-36) to 46 elements (slots 0-45). Update both the `Rf_allocVector(VECSXP, ...)` call site and the names STRSXP allocation accordingly.

**Population mechanism**: `pack_solver_result` is a generic lambda in r_bridge.cpp (~line 407) that touches only `res.base.*` shared fields — it does NOT dispatch on result type. CP-specific scalars (`n_cells`, `algorithm_used`, etc.) are captured in the CP dispatch arm (mirroring how the newton_kl arm captures `res_n_projected_dims` and `res_lm_mu_final` at r_bridge.cpp:608-609), then SEXP-pack writes slots 37-45 unconditionally. Non-CP solvers leave the C++ scalars at default-init values: `int = 0`, `double = NA_REAL`, `bool = FALSE`, `std::string = ""`. R-side these surface as `attr(result, "result")$n_cells`, etc. Empty/NA values for non-CP solvers are documented in `harvest.Rd` `@return` section.

### R-side warning/error semantics for CP status_codes

Match newton_kl precedent: warnings emit from R/harvest.R post-call (lines 479-488 show status-code → R `warning()` dispatch). NOT C++ `Rf_warning()`. Add CP-specific dispatch arm in harvest.R after the existing newton_kl arm:

| status_code | R behavior (emitted from harvest.R post-call) |
|---|---|
| 0 (converged) | Silent — return result |
| 1 (max_iterations) | `warning(sprintf("leafblower cp: max_iterations reached at iter=%d, max_error=%.2e; consider increasing max_iterations or trying method='ieppa'", iter, max_err))` |
| 2 (NaN/Inf detected) | `warning(sprintf("leafblower cp: NaN/Inf detected at iter=%d; weights truncated; consider tightening max_weight or trying method='ieppa'", iter))` |
| 3 (infeasible bounds) | `warning(sprintf("leafblower cp: infeasible bounds (lo[i] >= hi[i] - 2e-8 for some i); returning init-clamp weights"))` |
| 4 (power-iter divergence) | `warning("leafblower cp: power-iteration on ||A|| did not converge; aborted before solver started; report a bug")` |

**status_code=3 weight contract**: returned weights are the init-clamp values $w_i = \mathrm{clamp}(d_i, \ell_i + 10^{-8}, u_i - 10^{-8})$ (per spec Sec 3 common form init). NOT all-1.0. This matches the spec's init invariant and lets the caller see the design-weighted starting point even on infeasible bounds.

NEVER `stop()` on solver status — match existing solver patterns. Caller can inspect `attr(result, "result")$status` for programmatic dispatch.

### r_bridge.cpp dispatch arm

Place BEFORE the catch-all `else` (chebyshev/ieppa_soft/ieppa default), AFTER the `newton_kl` arm:

```cpp
} else if (strcmp(method_str, "cp") == 0) {
    auto res = lbw::cp_calibrate(st);
    res_status     = res.base.status;
    res_iterations = res.base.iterations;
    res_max_error  = res.base.max_error;
    res_alg_used   = (int)RK_ALG_CP;
    res_n_cells    = res.n_cells;
    pack_solver_result(res);
    if (!res.base.best_weights.empty())
        res_best_weights = std::move(res.base.best_weights);
    else
        res_best_weights.assign(st.n, 0.0);
}
```

## 3. Algorithm Specifications

### 3.1 Common framework

**Problem:**
$$\min_w \; \sum_i d_i \, \mathrm{KL}(w_i / d_i \| 1) \quad \text{s.t.} \quad A^\top w = b, \; \ell \le w \le u$$

**Init:** $w^0_i = \mathrm{clamp}(d_i, \ell_i + 10^{-8}, u_i - 10^{-8})$.

**Power iteration:** 50-iter on $A^\top A$ (or cell-compressed $A_c^\top A_c$); stop on $|\lambda_{k+1} - \lambda_k|/\lambda_k < 10^{-6}$; status_code=4 if relative delta > $10^{-3}$ at k=50.

**Stop criteria** (both algorithms):
- $\|A^\top w - b\|_\infty / \max(Z, 1) < 10^{-7}$ → status_code=0
- max_iter → status_code=1
- NaN/Inf → status_code=2
- Power-iter divergence → status_code=4

### 3.2 Cell compression (default; reverts to obs-level when M_cell/n > 0.9)

**Pre-flight**: build cell table via `lbw::cell_table` (existing pattern in `src/ieppa.cpp` + `src/cell_table.hpp`).

| State | Source |
|---|---|
| `M_cell` | unique `(g_1, ..., g_K)` tuples |
| `cell_d_sum` | $\sum_{i \in \mathrm{cell}_c} d_i$ (length M_cell) |
| `cell_to_obs` | obs index list per cell |
| `A_cell` | block-incidence on cells (M_cell × ΣJ_k, sparse, K nonzeros/row) |

**Decision**: if `static_cast<double>(M_cell) * 10.0 > static_cast<double>(n) * 9.0` (M_cell/n > 0.9), skip cell-compression and run obs-level. Else run cell-compressed. Use double cast to avoid 32-bit integer overflow at n > 2e8 (matches `src/ieppa.cpp` line ~168 `estimate_M_cell` precedent).

**Bounds_mode constraint**: cell-compressed CP requires `bounds_mode = "cell"` (per-obs bounds aren't cell-aggregable). For `bounds_mode = "unit"`, force obs-level path regardless of M_cell/n.

**Weight expansion**: at termination, $w_i = w_c[\text{cell}(i)]$ — cell-mode bounds invariant.

### 3.3 Algorithm 1 (vanilla PDHG, `accelerate = FALSE`)

Verbatim spike implementation:

```
σ = τ = 1 / (||A|| · 1.05)
Init: w⁰, w̄⁰ = w⁰, y⁰ = 0
For k = 0 .. max_iter:
  y^{k+1} = y^k + σ (A^T w̄^k − b)
  w^{k+1} = prox_{τ f}(w^k − τ A y^{k+1})
  w̄^{k+1} = 2 w^{k+1} − w^k
```

Convergence: O(1/k) ergodic on primal-dual gap.

### 3.4 Algorithm 2 (accelerated PDHG, `accelerate = TRUE` default)

Requires γ-strong convexity in f. KL on box $[\ell, u]$ with $u_{\max} < \infty$ has $\gamma = 1/u_{\max}$.

**Fallback (γ = 0 case)**: γ is computed from the SCALAR upper bound `st.max_weight` (single value broadcast to all i; matches existing newton_kl/ieppa convention — `harvest()` exposes only scalar `max_weight`/`min_weight`, never per-obs vectors). If `st.max_weight = Inf`, γ = 0 → set `accelerate = FALSE` automatically + log to `verbose ≥ 1`. `algorithm_used = "pdhg"` in result struct (and `algorithm_requested = "accelerated_pdhg"`, `fell_back_to_pdhg = true`).

```
γ = 1 / max(u_i)
σ_0 = τ_0 = 1 / (||A|| · 1.05)
Init: w⁰, w̄⁰ = w⁰, y⁰ = 0, θ_0 = 1
For k = 0 .. max_iter:
  y^{k+1} = y^k + σ_k (A^T w̄^k − b)
  w^{k+1} = prox_{τ_k f}(w^k − τ_k A y^{k+1})
  θ_{k+1} = 1 / sqrt(1 + 2 γ τ_k)
  τ_{k+1} = θ_{k+1} τ_k
  σ_{k+1} = σ_k / θ_{k+1}
  w̄^{k+1} = w^{k+1} + θ_{k+1} (w^{k+1} − w^k)
```

Step-size invariant maintained: $\sigma_k \tau_k \|A\|^2 \le 1$ for all k. Convergence: O(1/k²) on primal-dual gap.

**Stability fallback (θ_k underflow)**: if $\theta_k < 10^{-15}$ ($\tau_k$ approaches denormal) at any iter k_fallback, RESET adaptive step sizes to fixed $\tau = \sigma = 1/(\|A\| \cdot 1.05)$ AND continue under Algorithm 1 update from iter k_fallback onward. Record `final_theta = θ_{k_fallback}` (last adaptive value), set `algorithm_used = "pdhg"`, `fell_back_to_pdhg = true`. Do NOT freeze τ at the underflowed value (would stall progress).

### 3.5 prox_{τ f} per component

Direct Newton (NO Lambert-W); per spike:

```
FOC: τ d_i log(w/d_i) + (w - z_i) = 0
Newton: w_{n+1} = w_n - (τ d_i log(w_n/d_i) + (w_n - z_i)) / (τ d_i / w_n + 1)
Init: w_0 = clamp(d_i, ℓ + δ, u - δ)
Stop: |f'(w_n)| < 1e-12 OR n=5
Overflow: |z/(τd) - 1| > 700 → w ≈ z (asymptote) then clamp
Final: clamp to [ℓ + δ, u - δ]
```

### 3.6 Diagnostic outputs

`CpCalibResult` populates (matches Sec 2 SEXP-pack slots 37-45):
- `n_cells` — int (M_cell or n)
- `algorithm_requested` — string (`"pdhg"` or `"accelerated_pdhg"`, user-asked)
- `algorithm_used` — string (after possible fallback)
- `A_norm_estimate` — double (power-iter ‖A‖)
- `n_power_iter` — int
- `final_theta`, `final_tau`, `final_sigma` — accelerated_pdhg only (NA_real_ if pdhg)
- `fell_back_to_pdhg` — bool

**No per-iter trace stored in production CpCalibResult.** Memory bounded by O(n + ΣJ_k) — final-iteration scalars only. Spike-style trace (1000-row CSV with per-iter `max_err_last`, `max_err_ergodic`, `primal_resid`, `primal_stationarity_proxy`) is RESEARCH-ONLY (`benchmarks/research/cp_*_trace_rep*.csv` from Epic-J). Production analysts who want trace must use `verbose=2` for stderr logging, OR re-run via the spike harness.

**Rate exponents (β_last, β_ergodic, R²) are spike-only.** Sec 1 cites β=-1.05, R²=0.96 as the empirical motivation/headline of why we're productionizing CP — they are NOT runtime-exposed in CpCalibResult or surface via `attr(r, "result")`. Production users get convergence quality via `max_error` and final residuals; rate-fit is a research-tool concept (`lm(log(max_err) ~ log(iter))` on captured trace).

## 4. Test Suite

`tests/testthat/test-cp.R` — 6 tests:

| Test | Purpose |
|---|---|
| **T1** | K=3 small (n=1000): `status=0`, `max_err < 1e-6` (machine-precision sanity). |
| **T2** | stepstone K=9: `cp` `max_err ≤ 0.7 × ieppa+sraa max_err` (LOCKS IN headline 2.2× win with safety margin; spike showed parity ratio 0.45). |
| **T3** | Cell-compressed CP (`bounds_mode="cell"`) ≡ obs-level CP (`bounds_mode="unit"`) within 1e-10 weight diff on K=2 fixture. |
| **T4** | Bounds-active fallback: tight `max_weight=1.3` on 95/5 target → finite `max_error` (no NaN propagation). |
| **T5** | KL-form vs chi2 (`greg`): cp weights distinct from greg by >1% rel diff; all cp weights `> 0`. |
| **T6** | `accelerate=FALSE` (Algorithm 1) reaches stepstone parity (`max_err < 1e-4`) — covers Alg 1 code path. Tightened from rev 1's 1e-3 to actually exercise Alg 1 quality (spike showed Alg 1 reaches 5.08e-5; 1e-4 leaves 50% headroom). |
| **T7** | `accelerate=TRUE` with `max_weight=Inf` triggers γ=0 fallback. Asserts `attr(r, "result")$fell_back_to_pdhg == TRUE` AND `algorithm_used == "pdhg"` AND `algorithm_requested == "accelerated_pdhg"`. Verifies R3 mitigation (Sec 3.4 fallback path). |

kk1204 NOT a CP test fixture (spike showed CP fails kk1204; out of scope).

## 5. Out of Scope

(See Sec 5 of brainstorm transcript or Sec 7 of spec for full list.)

| Item | Reason |
|---|---|
| AUTO routing change (stepstone-class default) | Epic-K.2 deferred — fixture-class sweep needed |
| Python bindings | Standard sequencing |
| SIMD/OpenMP, warm-start, predictor-corrector | Tuning phase |
| `cp_safety_factor` user knob | Hardcoded 1.05; Epic-K.2 |
| Vignette docs | Separate doc epic |
| Cell compression for `bounds_mode="unit"` | Force obs-level for unit mode |
| ABI breaking changes | Hardcoded knobs; no rk_params_t fields |

## 6. Risks & Mitigations

| # | Risk | L | I | Mitigation |
|---|---|---|---|---|
| R1 | Cell compression bug — cell-mode bounds invariant violated | M | H | T3 test: cell vs obs path produce same weights to 1e-10 |
| R2 | Algorithm 2 step-size adaptation unstable on ill-conditioned A | M | M | Fallback to Alg 1 if θ_k underflows; diagnostic `final_theta` |
| R3 | `u_max = Inf` violates γ > 0 precondition | L | M | Auto-fallback to PDHG; verbose log; `algorithm_used="pdhg"`, `fell_back_to_pdhg=true` |
| R4 | r_bridge dispatch arm wires wrong field set | L | M | Mirror `newton_kl` dispatch arm; spec reviewer audits |
| R5 | harvest.R `match.arg` whitelist drift | L | H | T1-T6 use `method="cp"`; whitelist drift → all tests fail |
| R6 | NEWS.md bullet under wrong section | L | L | Place under `## New features`; reviewer audits |
| R7 | rk_algorithm_t enum collision | L | H | Verify before adding; 12 currently free |
| R8 | A_cell construction differs from obs-level A | M | H | T3 direct comparison + assertion `sum(A_cell.x) == M_cell * K` |
| R9 | Power-iter on cell A vs obs A different `‖A‖` estimates | L | M | Self-consistent within each path; documented in result |
| R10 | accelerate=TRUE slower than Alg 1 (over-acceleration) | M | M | T6 verifies Alg 1 stepstone parity; default TRUE because spike rate already strong |
| R11 | T2 stepstone fixture not in CI | M | M | `skip_if_not_installed("arrow")` graceful skip |
| R12 | Algorithm bug not caught by spike sanity (max_err=0 only proves trivial recovery) | M | M | Spec compliance reviewer audits cp_calib.cpp on stepstone trace inspection |
| R13 | else-if dispatch order shifts existing fall-through | L | H | Place `cp` arm BEFORE catch-all `else` |
| R14 | bounds_mode="unit" cell expansion incorrect | L | H | bounds_mode="unit" forces obs-level; comment in cp_calib.cpp |
| R15 | A_cell sparse vs dense storage choice ambiguous | L | M | Sparse `Eigen::SparseMatrix<double>` for A_cell (matches obs-level A storage); for power-iter use sparse matvec on $A^\top A v$. If profiling shows dense beats sparse on M_cell·ΣJ_k < 1e6, switch in Epic-K.2 (out of scope here). |
| R16 | OpenMP interaction with cell_table reuse | L | L | CP solver runs single-threaded. cell_table.cpp build inherits ieppa OpenMP behavior unchanged. Document in cp_calib.cpp comment header: "single-threaded; cell-table build inherits cell_table.cpp parallelism". |
| R17 | tools/check_research_isolation.R blocks pre-commit because cp_solve_R + cp_calibrate now legitimately in src/leafblower.so | H (will fire on every K-1 commit attempt) | H | WU-K-1 step 1: edit `tools/check_research_isolation.R` forbidden-symbol list to REMOVE `cp_solve_R` and `cp_calibrate`. Keep `ipm_solve_R` and `ipm_calibrate` (still research-only per Epic-J FAIL). |

**Discontinuation triggers:**
- R5 (whitelist drift) → halt; verify match.arg before commit.
- R7 (enum collision) → halt; pick next free slot; audit downstream switch tables.
- R12 (algorithm bug) → halt; revert; re-spike before continuing.

## 7. Implementation Phases (Atomic Work Units)

One bd ticket per WU. Sequential per `superpowers:subagent-driven-development`.

| WU | Title | Hard deps | Model | Wall | Decision Gate |
|---|---|---|---|---|---|
| **K-1** | Move + adapt cp_calib to src/ + isolation gate update + R wiring | — | Gemini | ~2.5h | All sub-steps below complete in single atomic commit. |

**K-1 sub-step list** (single ticket, single commit per global rule 9.1 atomicity; revertible as one SHA):
1. Edit `tools/check_research_isolation.R`: REMOVE `cp_solve_R` and `cp_calibrate` from forbidden-symbol list. KEEP `ipm_solve_R` and `ipm_calibrate` (still research-only per Epic-J FAIL).
2. `git mv research/cp_calib.hpp src/cp_calib.hpp` and `git mv research/cp_calib.cpp src/cp_calib.cpp` (preserves history). Then add fossil-pointer header comment to `research/cp_calib.hpp`: `// Epic-K: MOVED TO src/cp_calib.hpp; this copy retained for Epic-J spike traceability only — DO NOT EDIT`.
3. Refactor `src/cp_calib.hpp`: rename include guard `LEAFBLOWER_RESEARCH_CP_CALIB_HPP_` → `LEAFBLOWER_CP_CALIB_HPP_`. Replace standalone signature `cp_calibrate(int n_row, int n_col, ...)` with `lbw::CpCalibResult lbw::cp_calibrate(lbw::CalibState& st)` matching newton_calib.hpp pattern.
4. Add OpenMP/thread-safety comment header to top of `src/cp_calib.cpp`: `// Single-threaded solver; cell-table build inherits ieppa OpenMP behavior unchanged.` (R16 mitigation).
5. Update `src/cp_calib.cpp` body: read inputs from CalibState (replace raw pointer args), construct `CpCalibResult` with `CalibResult base + n_cells + algorithm_requested + algorithm_used + A_norm_estimate + n_power_iter + final_theta + final_tau + final_sigma + fell_back_to_pdhg`. NO per-iter trace storage (production diagnostics = final-iteration scalars only).
6. `src/leafblower.h`: add `RK_ALG_CP = 12,` to `rk_algorithm_t` enum.
7. `src/r_bridge.cpp`:
   - Line 25-37: insert `{"cp", RK_ALG_CP}` after `{"newton_kl", RK_ALG_NEWTON_KL}`.
   - Line ~617 (after newton_kl arm, before catch-all else block at 618): insert new `else if (strcmp(method_str, "cp") == 0) { ... }` arm following the newton_kl pattern, packing CpCalibResult fields into res_n_cells, res_algorithm_used, etc.
   - Line ~709-718: insert `(res_alg_used == static_cast<int>(RK_ALG_CP)) ? "cp"` before `: "iEPPA"` fallback.
   - Result list size: extend `Rf_allocVector(VECSXP, ...)` from 37 to 46 elements; extend names allocation; populate slots 37-45 per Sec 2 Result SEXP-pack table.
8. `R/harvest.R`:
   - Line ~318 (BEFORE line 319): insert `accelerate_explicit <- !missing(accelerate)`.
   - Line ~318 (BEFORE line 322): insert `if (method == "cp" && !accelerate_explicit) accelerate <- TRUE`.
   - Line 319: add `"cp"` to warning predicate whitelist.
   - Line 322: add `"cp"` to `accelerate_bool` predicate whitelist.
   - Line ~488 (after newton_kl status-code-warning arm; mirror lines 479-488 pattern): add CP-specific warning dispatch per Sec 2 R-side warning table.
   - Line 590: add `"cp"` after `"newton_kl"` in match.arg whitelist.
9. `R CMD INSTALL --preclean .` → exit 0; pre-commit isolation gate green (cp_* now allowed in src/leafblower.so).
10. Smoke test: `Rscript -e 'library(leafblower); set.seed(1); df <- data.frame(x=factor(sample(c("a","b","c"), 100, TRUE))); tgt <- list(x=c(a=0.4, b=0.4, c=0.2)); r <- harvest(df, tgt, method="cp", max_weight=5); cat("alg:", attr(r, "algorithm"), "\n")'` → returns "cp" without warnings (accelerate=TRUE default for cp; whitelist now includes cp; no whitelist-warning fires).
| **K-2** | Algorithm 1 obs-level production parity (port spike with src/ conventions) | K-1 | Gemini | ~2h | T1 + T6 PASS; T5 PASS (KL form distinct from greg) |
| **K-3** | Algorithm 2 accelerated variant + fallback when u_max=Inf | K-2 | Opus | ~3h | T1 PASS; T6 PASS (Alg 1 path); default accelerate=TRUE produces stepstone tighter than Alg 1 |
| **K-4** | Cell compression with bounds_mode="cell" + obs-level fallback at M_cell/n > 0.9 + bounds_mode="unit" | K-2 | Opus | ~4h | T3 PASS (cell ≡ obs); T2 PASS (stepstone tighter than ieppa+sraa) |
| **K-5** | Test suite tests/testthat/test-cp.R (T1-T6) | K-1, K-2, K-3, K-4 | Haiku | ~1h | All 6 PASS via `devtools::test()`; full regression FAIL=0 outside Epic-Dβ T2 documented basin |
| **K-6** | NEWS.md additive bullet + harvest.Rd regen + harvest.R docstring | K-5 | Haiku | ~30min | `devtools::document()` clean; `R CMD check` no NOTES related to cp; NEWS.md bullet text matches draft below |

**K-6 NEWS.md draft text** (place under `## New features` of the development version section):
```
* New `method="cp"` (Chambolle-Pock primal-dual; Chambolle & Pock 2011)
  for moderate-skew K>=5 calibration problems. Outperforms `method="ieppa"`
  with `accelerate=TRUE` on stepstone-class fixtures by ~2x (max_err
  5.08e-5 vs 4.39e-4 on stepstone_K9; rate exponent beta=-1.05, R^2=0.96).
  Cell-compressed by default (M_cell/n <= 0.9, bounds_mode="cell"); falls
  back to obs-level otherwise. `accelerate=TRUE` (default for cp) selects
  the O(1/k^2) accelerated PDHG variant; falls back to O(1/k) plain PDHG
  when `max_weight=Inf`. Opt-in only — AUTO routing unchanged in this
  release. NOT recommended for severe-skew K>=5 (target_skew > 5) where
  CP fails to converge (per Epic-J spike investigation: see
  docs/investigations/2026-05-02-ylsy-cp-ipm-spike-result.md).
```

**K-6 harvest.Rd `@param method` insertion** (adjacent to chebyshev grouping):
```
\code{"cp"} (Chambolle-Pock primal-dual splitting; for moderate-skew
K>=5 problems — outperforms ieppa+sraa on stepstone-class fixtures by
~2x. Cell-compressed by default; obs-level fallback when M_cell/n > 0.9
or bounds_mode="unit". Supports accelerate=TRUE for accelerated PDHG.)
| **K-7** | Final code-review-gate (3 adversarial reviewers) + cleanup commit | K-6 | Opus | ~1h | All 3 reviewers (Feasibility, Completeness, Scope & Alignment per `metaswarm:plan-review-gate` convention) PASS. Cleanup: ensure no residual research/cp_calib edits, no orphan .o/.so artefacts staged. |

**Total wall:** ~13–14h sequential.

**Dependency graph:**
```
K-1 ──► K-2 ──► K-3 ──► K-5 ──► K-6 ──► K-7
            │           ▲
            └─► K-4 ────┘
```

K-3 + K-4 both depend on K-2 (Alg 1 baseline correct first). K-5 hard-dep both K-3 + K-4.

**Reversibility**: each WU = one commit; revert single SHA on bug. `research/` untouched.

**Post-execution close protocol:**
- bd close all 7 WU tickets sequentially.
- bd close Epic-K (`leafblower-<TBD>`).
- File Epic-K.2 (AUTO routing sweep + carve-out) — DEFERRED ticket only.
- Optionally close `leafblower-ylsy` (long-standing kk1204) — but ylsy explicitly about kk1204, which Epic-K does not address. Keep ylsy open or close BLOCKED per user.

## 8. Predecessor & Memory References

- **Predecessor**: Epic-J (`leafblower-y2ls`) FAIL at master `484e1e2`. Side-finding documented in Sec 6 of `docs/investigations/2026-05-02-ylsy-cp-ipm-spike-result.md`.
- **Spike CP at**: `research/cp_calib.{hpp,cpp}` master `4e89769` (Algorithm 1, obs-level only).
- **Spec rev 2 (Epic-J)**: `docs/superpowers/specs/2026-05-02-ylsy-cp-ipm-spike-design.md` Sec 3.1 (CP algorithm reference).
- **Cell-table precedent**: `src/cell_table.cpp` + `src/cell_table.hpp`; pattern in `src/ieppa.cpp`.
- **`bd memory ylsy`**: `bd update leafblower-ylsy --notes` on Epic-J close — verdict pointer.
- **TSVD precedent (not used here, but for reference)**: `src/newton_calib.cpp` lines 350-400 Epic-Dβ WL-1.
