# Epic-K: stepstone-CP Productionization — Design

**Date:** 2026-05-02
**Status:** Design — pre-plan
**Predecessor:** Epic-J (`leafblower-y2ls`) FAIL verdict on kk1204 + side-finding CP wins stepstone (parity 0.45 vs `ieppa+sraa` baseline 1.13e-4)
**Investigation report:** `docs/investigations/2026-05-02-ylsy-cp-ipm-spike-result.md` Sec 2 + Sec 6 (Side-finding)

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
└── r_bridge.cpp                      ← add `else if (strcmp(method_str, "cp") == 0)` dispatch arm

R/
└── harvest.R                         ← match.arg += "cp"; @param method docstring; accelerate
                                       whitelist += "cp" (line ~322)

tests/testthat/
└── test-cp.R                         ← NEW (6 tests: T1-T6)

man/harvest.Rd                        ← regenerate via devtools::document()
NEWS.md                               ← additive entry under "## New features"

research/                             ← stays under .Rbuildignore (Epic-J FAIL artefact)
                                       cp_calib.{hpp,cpp} now duplicated in src/; ABANDONED
                                       per FAIL artefact policy
```

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

`CpCalibResult` in `cp_calib.hpp` mirrors `NewtonCalibResult`:
- `CalibResult base` (status, iterations, max_error, weights, best_weights)
- `int n_cells` — `M_cell` if cell-compressed, `n` if obs-level
- `std::string algorithm_used` — `"alg1"` or `"alg2"`
- `double aA_norm_estimate` — power-iter `‖A‖`
- `int n_power_iter` — power-iter convergence count
- `double final_theta` — Alg 2 only, last θ_k
- `double final_tau`, `final_sigma` — Alg 2 final adaptive step sizes

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

**Decision**: if `M_cell * 10 > n * 9` (M_cell/n > 0.9), skip cell-compression and run obs-level. Else run cell-compressed.

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

**Fallback**: if any $u_i = \infty$, $\gamma = 0$ → set `accelerate = FALSE` automatically + log to `verbose ≥ 1`. `algorithm_used = "alg1"` in result struct.

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

**Stability fallback**: if $\theta_k \to 0$ ($\tau_k$ underflows below 1e-300) at any iter, force Algorithm 1 from that iter onward; record `final_theta` in result.

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

`CpCalibResult` populates:
- `n_cells` — int (M_cell or n)
- `algorithm_used` — `"alg1"` or `"alg2"`
- `aA_norm_estimate` — double
- `n_power_iter` — int
- `final_theta`, `final_tau`, `final_sigma` — Alg 2 only

## 4. Test Suite

`tests/testthat/test-cp.R` — 6 tests:

| Test | Purpose |
|---|---|
| **T1** | K=3 small (n=1000): `status=0`, `max_err < 1e-6` (machine-precision sanity). |
| **T2** | stepstone K=9: `cp` `max_err ≤ 1.5 × ieppa+sraa max_err`. The headline-result regression test. |
| **T3** | Cell-compressed CP (`bounds_mode="cell"`) ≡ obs-level CP (`bounds_mode="unit"`) within 1e-10 weight diff on K=2 fixture. |
| **T4** | Bounds-active fallback: tight `max_weight=1.3` on 95/5 target → finite `max_error` (no NaN propagation). |
| **T5** | KL-form vs chi2 (`greg`): cp weights distinct from greg by >1% rel diff; all cp weights `> 0`. |
| **T6** | `accelerate=FALSE` (Algorithm 1) reaches stepstone parity (max_err < 1e-3) — covers Alg 1 code path. |

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
| R3 | `u_max = Inf` violates γ > 0 precondition | L | M | Auto-fallback to Alg 1; verbose log; `algorithm_used="alg1"` |
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

**Discontinuation triggers:**
- R5 (whitelist drift) → halt; verify match.arg before commit.
- R7 (enum collision) → halt; pick next free slot; audit downstream switch tables.
- R12 (algorithm bug) → halt; revert; re-spike before continuing.

## 7. Implementation Phases (Atomic Work Units)

One bd ticket per WU. Sequential per `superpowers:subagent-driven-development`.

| WU | Title | Hard deps | Model | Wall | Decision Gate |
|---|---|---|---|---|---|
| **K-1** | Move + adapt cp_calib to src/, wire enum + r_bridge dispatch + harvest.R match.arg | — | Gemini | ~2h | Build clean; `harvest(small_df, tgt, method="cp")` returns weights without error |
| **K-2** | Algorithm 1 obs-level production parity (port spike with src/ conventions) | K-1 | Gemini | ~2h | T1 + T6 PASS; T5 PASS (KL form distinct from greg) |
| **K-3** | Algorithm 2 accelerated variant + fallback when u_max=Inf | K-2 | Opus | ~3h | T1 PASS; T6 PASS (Alg 1 path); default accelerate=TRUE produces stepstone tighter than Alg 1 |
| **K-4** | Cell compression with bounds_mode="cell" + obs-level fallback at M_cell/n > 0.9 + bounds_mode="unit" | K-2 | Opus | ~4h | T3 PASS (cell ≡ obs); T2 PASS (stepstone tighter than ieppa+sraa) |
| **K-5** | Test suite tests/testthat/test-cp.R (T1-T6) | K-1, K-2, K-3, K-4 | Haiku | ~1h | All 6 PASS via `devtools::test()`; full regression FAIL=0 outside Epic-Dβ T2 documented basin |
| **K-6** | NEWS.md additive bullet + harvest.Rd regen + harvest.R docstring | K-5 | Haiku | ~30min | `devtools::document()` clean; `R CMD check` no NOTES related to cp |
| **K-7** | Final code-review-gate (3 adversarial reviewers) + cleanup commit | K-6 | Opus | ~1h | All 3 reviewers PASS |

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
