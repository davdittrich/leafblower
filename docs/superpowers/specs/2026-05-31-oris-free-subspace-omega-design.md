# ORIS Free-Subspace θ₂ for Box-Constrained Over-Relaxation — Design

**Date:** 2026-05-31
**Status:** Design (rev 3 — addresses design-review-gate iter 2: Architect/Designer/Security)
**Predecessor:** leafblower-mj1p.2 (spectral optimal-ω) closed NO-GO — Lehmann's *global* θ₂ estimator was slower (40 iters) than fixed ω=1.5 (30 iters) on bounded stepstone.
**Research notebook:** NotebookLM `1e3036a1-fcbb-4d05-bc8e-820854f59d8e` (do not delete; 191 sources).
**Verified code map:** `docs/superpowers/derivations/oris_structure_map.md` (byte-exact, the authoritative anchor — line numbers in THIS spec are advisory; the map + symbol names are load-bearing).

---

## 0. Deliverable framing (PRIMARY vs CONDITIONAL)

The user's request is *"derive the Lehmann formula for box-constraint cases."* The derivation is the
**primary deliverable**; SOR is off-by-default (`sor=NULL`), so a shipped feature is not assumed.

- **Phase 1 (PRIMARY, always done): SageMath derivation + documented finding.** Symbolic+numerical
  artifact that (a) **derives** the box-constrained optimal-ω from the KKT/fixed-point structure —
  not merely confirms a borrowed formula — and (b) produces the exact *free-subspace residual
  functional* and *θ₂ estimator formula* that Phase 2 will implement. Expected-value-positive
  regardless of whether ORIS changes.
- **Phase 2 (CONDITIONAL): ORIS implementation** of the Phase-1-derived formula. Attempted only if
  Phase 1 GOs AND the §5 ship gate passes. Else: revert ORIS, `omega_mode_id` default → 1 (fixed),
  commit the derivation + a documented negative.

**Key reframe (resolves Security BLOCKING-1/2 and Architect loop concerns):** the *exact arithmetic*
of the estimator (the θ₂ formula, the residual denominator, the free-cell aggregation) is an
**output of the Phase-1 derivation**, not guessed here. This spec fixes the **guard envelope and
contracts** (§3.4) the derived formula must satisfy; Phase 2 implements the verified formula inside
that envelope. This prevents shipping premature/incoherent arithmetic (the mj1p.2 failure mode).

## Mechanism
Free-subspace spectral estimate of ω: estimate the local convergence rate from error reduction over
the **free (unclamped) coordinate block**, feed Lehmann `ω = 2/(1+√(1−θ₂))`.

## Forbidden
- Global/all-cell residual ratio driving θ₂ (the mj1p.2 bug).
- Changing the fixed point: ω changes only the iteration path.
- Generic LCP solver substitution (Siconos/PATH) — wrong layer.
- SRAA adapt path (deferred; confounds with Anderson — §7).
- `oris_soft` / ALM mode (`alm_active`) — out of scope (§2, §3.1); free-subspace pinning is
  ill-defined when the soft-capacity Newton solution sits strictly inside the box.
- ω ≥ 2, `√` of a negative, NaN/Inf reaching the formula (§3.4 guard envelope).

## Audit
- Phase-1 SageMath artifact derives `ρ(M_II)` and the estimator; it is the Phase-1 GO gate.
- Orchestrator independently re-runs every benchmark number (subagent numbers not trusted).

---

## 1. Problem, root cause, literature (self-contained)

Lehmann–von Renesse–Sambale–Uschmajew (2022, arXiv:2012.12562) give Young's optimal SOR factor
`ω_opt = 2/(1+√(1−ρ(M)²))`, `ρ(M)` = second eigenvalue of the linearized Sinkhorn/IPF Jacobi at the
fixed point (largest is 1 from scaling indeterminacy).

**Settled literature finding** (NotebookLM 191-source synthesis + direct query; see notebook): no
closed-form box-constrained `ω_opt` exists. With per-cell clamping `L_c ≤ X[c] ≤ U_c`, error
propagates only through the **free coordinate block** `I`; the active set is a hard wall. The
active-set-dependent optimum is `ω_opt(I) = 2/(1+√(1−ρ(M_II)²))`, `M_II` = principal submatrix of
`M` on free coordinates. **Extensions checked:** Soma–Uschmajew 2024 (arXiv:2410.14104) and the
"Numerically stable variants of overrelaxation for operator Sinkhorn iteration" paper — *neither
addresses box/inequality constraints*; so this derivation is new, not a re-import. (Phase-1 §4
re-derives `ρ(M_II)` rather than relying on the synthesis assertion — PM BLOCKING-1/2.)

mj1p.2 fed the *global* residual ratio in; pinned cells give zero reduction → global ratio → 1 →
θ₂ → 1 → ω → 1.9 → over-relaxes free coords → oscillation → slower than fixed.

## 2. ORIS code reality (VERIFIED — supersedes rev-1/rev-2 citations)

Source: byte-exact map in `docs/superpowers/derivations/oris_structure_map.md`. Corrects two prior
errors: (a) the water-fill clamp is **SRAA-only**, not flat-path; (b) **there is NO shared
`compute_errRp_ct` helper** — that was a rev-1/rev-2 hallucination; the flat residual is computed
**inline**.

- **Two sibling branches** (not if/else of the iteration body): SRAA `if (sraa_active_lvl)` at
  `oris.cpp:822` (body ~822–995); shared clamp/sync `~997–1046`; flat `if (!sraa_active_lvl)` at
  `oris.cpp:1050` (body 1050–1840). Clean boundary — flat-only change cannot perturb SRAA.
- **Flat capacity clamp (active-set origin):** element-wise `X[c]=std::clamp(X_tilde,L_cell[c],
  U_cell[c])` at `oris.cpp:1362` (linear) / `1444` (log). ALM branch `1433–1442` sets `X=X_alm`
  (Newton, inside box) — **only when `alm_active` (oris_soft), which is OUT OF SCOPE**; for
  `method="oris"` `alm_active==false`, so pinning is well-defined. (Architect BLOCKING-2.)
- **Flat residual (θ₂ driver):** computed **inline** in the convergence sweep `oris.cpp:1504–1549`
  — `S_lin[j]` built inline at `1512–1516`, `errRp_k → per_k_errRp_cache[k]` at `1518–1524`, gated
  on `iter % kErrCheckInterval`. Consumed by the omega-adapt block.
- **Flat omega-adapt block:** `oris.cpp:1588–1642`. `sor_omega[k]` written at `1616` (damp), `1628`
  (spectral), `1631` (fixed), `1634` (recovery). This is the only edit site for the ω logic.
- **θ₂ machinery (mj1p.2):** `kSorSpectralCeiling=1.99` (`:414`), `estimate_theta2` (`:415`, returns
  `(curr/prev)²` — **lag-1 squared**), `omega_from_theta2` (`:425`). Solve-local SOR state:
  `sor_omega` (`:430`), `sor_prev_errRp` (`:431`), `sor_prev_decreasing` (`:432`),
  `sor_theta2_ema` (`:435`) — all `std::vector(st.K)`, no struct fields.
- `L_cell` `:158/161` (level-invariant); `U_cell` `:516–519` (per homotopy level). Both live at the
  clamp and the adapt sites.

**Two distinct loops, stated explicitly (Architect BLOCKING-1):** `is_pinned[c]` is determined at
the **clamp loop** (`1362`/`1444`); the **free residual** is accumulated in the **convergence sweep**
(`1512–1516`, gated on `kErrCheckInterval`). These are different loops separated by the ALM/metrics
blocks. `is_pinned[]` is a solve-local `std::vector<char>(M_cell)` written in the clamp loop and read
in the sweep — no fused single loop is assumed.

## 3. Phase-2 change (flat path, `method="oris"` only)

### 3.1 Active set (new solve-local computation)
Add `std::vector<char> is_pinned(ct.M_cell)`, set in the flat capacity-clamp loop (after `1362`/
`1444`): `is_pinned[c] = (X[c] >= U_cell[c]*(1-kPinTol)) || (X[c] <= L_cell[c]*(1+kPinTol))`.
`kPinTol = 1e-9` (relative; the flat clamp has no bisection tol to inherit — Architect BLOCKING-2).
`alm_active` path excluded by scope (§2).

### 3.2 Free residual (parallel accumulator in the existing sweep — no new pass)
In the convergence sweep `1512–1516`, alongside `S_lin[j] += X[c]`, accumulate
`if (!is_pinned[c]) S_lin_free[j] += X[c];` and a per-margin free-cell count `n_free_k`. The exact
**free-residual functional `errF_k`** (which cells, which denominator) is the **Phase-1 derived
formula** (§4) — Phase 2 implements that. Cost: O(M_cell), folded into the existing sweep.

### 3.3 θ₂ estimator = Phase-1 formula, NOT a reuse of `estimate_theta2`
`estimate_theta2` returns the **lag-1 squared** ratio `(curr/prev)²`. The free-subspace estimator
uses the Phase-1-derived lag-2 (or derived) formula on `errF`. **It must not silently reuse
`estimate_theta2`** (different lag → wrong exponent — Security BLOCKING-1). Phase 1 fixes the single
unambiguous formula and proves it lands θ₂ ∈ [0,1). Add solve-local lag history
(`errF_prev`, `errF_prev2`, `std::vector(st.K)` near `:435`).

### 3.4 Guard envelope (the safety contract — every guard a precondition, ordered)
Per margin `k`, before computing/using θ₂:
1. **Scope:** `alm_active` → not reached (oris_soft excluded).
2. **Free-set gate:** `n_free_k == 0` (fully pinned) → `ω_k = 1`; skip. (Security BLOCKING-2 part.)
3. **Warm-up gate:** `iter < sor_burnin` OR lag history unfilled → `ω_k = 1`. Cold-start:
   `errF_prev/errF_prev2` init to `+inf`; first two informative samples → uninformative → `ω_k = 1`.
   (Security BLOCKING-3.)
4. **Finiteness/mass gate:** `!isfinite(errF_p)` OR `W_total_free < kMinSafeTotalWeight (1e-100)` →
   `ω_k = 1`. (Security BLOCKING-2.)
5. **Converged-denominator gate:** `errF_prev2[k] < kResidFloor (1e-12)` → `ω_k = 1`; no division.
   `kResidFloor` is an **intentional absolute floor** on the relative residual (documented;
   fail-safe direction = under-relaxation). (Security note.)
6. **Ratio + clamp BEFORE sqrt:** compute the Phase-1 θ₂ from the gated ratio; reject `ratio ≥ 1`
   (residual grew) → treat as uninformative → drive ω **down toward 1** (NOT clamp-to-ceiling — the
   mj1p.2 trap); else `θ₂ = clamp(value, 0, 1−1e-9)`. (Security BLOCKING-1.)
7. **EMA:** smooth θ₂. **On every ω=1 fallback (gates 2–5), decay `sor_theta2_ema[k]` toward 0**
   (and require lag re-warm-up after any `n_free_k` transition) — prevents stale-EMA over-relaxation
   after active-set churn. (Security BLOCKING-3.)
8. **Formula:** `ω_k = 2/(1+√(1−θ₂))`, then `min(ω_k, 1.99)`.
9. **Oscillation damp** (`kSorOscillationDamp`, `oris.cpp:1616`) still fires on sign-flip and reads
   the **free-subspace residual `errF`** so damp and estimator agree on "diverging".
10. **Monotone hard-fallback (deterministic, independent of sign-flip):** if `errF_p > errF_prev[k]`
    for `m=3` consecutive adapted steps → force `ω_k = 1` for a cooldown window. Guarantees no
    active-set limit cycle (Security BLOCKING-4).

### 3.5 Degenerate-state → ω (complete; matches §3.4 arithmetic)
| State | ω | Note |
|-------|---|------|
| `alm_active` | n/a | oris_soft out of scope |
| Fully pinned (`n_free_k=0`) | 1 | gate 2 |
| Cold/warm-up/unfilled lag | 1 | gate 3 |
| `!isfinite(errF)` / tiny free mass | 1 | gate 4 |
| `errF_prev2<kResidFloor` (converged) | 1 | gate 5 |
| `ratio ≥ 1` (residual grew) | driven toward 1 | gate 6 (NOT ceiling) |
| 3× consecutive `errF` increase | 1 (cooldown) | gate 10 |
| Normal | `2/(1+√(1−θ₂))`, ≤1.99 | gate 8; **EMA-lagged** — may stay elevated for a few iters after a spike (soft guard) |

### 3.6 Convergence safety (asserted, not just tested)
ω depends on `I` (feedback loop). Safety rests on three asserted invariants, not just "test no
regression": (a) the oscillation damp (gate 9) + monotone hard-fallback (gate 10) bound ω; (b)
best-iterate tracking + max-iter floor guarantee **termination regardless of ω trajectory**
(asserted invariant — Security BLOCKING-4); (c) **baseline snapshot**: Phase-2 step 1 records the
current fixed-ω converged/NOCONV status of all fixtures as the reference; the gate is "no fixture
flips converged→NOCONV". (CTO note.)

### 3.7 No new INPUT ABI; binding parity asserted
`is_pinned`, `S_lin_free`, `n_free_k`, lag buffers = solve-local; `CalibSorCfg` gains no field;
`omega_mode_id=2` reused. Logic in `oris.cpp` core, reached identically by `r_bridge.cpp` and
`c_api.cpp` (no per-binding wiring for the estimator). **Parity test requirement (Designer note):**
the Python parity suite MUST exercise `omega_mode_id=2` explicitly so a c_api wiring miss (the
mj1p.1 failure) is caught before merge.

### 3.8 Observability (OUTPUT ABI — additive; Designer BLOCKING + Architect BLOCKING-4)
Add to the result struct: `sor_omega_mean` (mean realized ω over adapted steps) and **per-gate
fallback counters** — at minimum `sor_n_pinned_fallback` (gate 2) separate from
`sor_n_omega1_other` (gates 3–5,10) so a silent all-fallback is diagnosable (pinned-everywhere vs
warm-up-misconfig vs non-monotone). These are **output-ABI additive** changes: update the result
struct AND the `leafblower.h` ABI-size comment/`EXPECTED_RK_RESULT_BYTES`, wired through both
`r_bridge.cpp` and `c_api.cpp`/pybind11, same commit.

### 3.9 mode-2 default + docs/migration
mj1p.2 shipped `omega_mode_id=2` ("spectral global") as DEFAULT. Decision rule, pre-stated and tied
to **§5 ship gate**: redefined mode 2 ships as default **only if the §5 ship gate passes**; else
default reverts to mode 1 (fixed), free-subspace stays opt-in (`omega_mode_id=2L`). Either way, same
commit updates: roxygen `@param` (`R/harvest.R`), Python binding docstring, `docs/methods/oris.md`
(guarantees-table spectral row), and a NEWS entry noting mode 2's behavioral change.

## 4. Phase-1 SageMath artifact (PRIMARY deliverable + GO gate + estimator source)

`docs/superpowers/derivations/2026-05-31-free-subspace-omega.sage` (+ committed result table):
1. Build a small synthetic box-constrained IPF (3×3 seed, prescribed margins, one cell forced to
   clamp at `U`).
2. **Derive** (symbolically, from the KKT/fixed-point stationarity) the linearized Jacobi `M`,
   extract the free submatrix `M_II`, compute `ρ(M)`, `ρ(M_II)` exactly; **emit the explicit
   free-subspace residual functional and the θ₂ estimator formula** (which cells, denominator, lag)
   that Phase 2 §3.2–3.3 will implement.
3. Verify numerically: free-coordinate error ratio → `ρ(M_II)²`; global ratio → 1 (reproduces the
   bug); `ω_opt(I)` beats fixed ω on iteration count.

**Phase-1 GO** = derivation produced AND all three numerical checks pass. **NO-GO** = stop; commit
artifact + documented negative; do not start Phase 2.

## 5. Phase-2 fixtures + pre-registered ship gate (CTO BLOCKING-1/2/3)

**Early exit (cheap, first):** before building the synthetic fixture, evaluate the bounded win on the
**existing stepstone mw=5** fixture. If free-subspace cannot even tie fixed there, stop (saves
fixture work). mw=5 remains **no-regression-only**, not a ship gate (it was the mj1p.2 NO-GO case).

**Slow-unconstrained fixture (the ship gate; deterministic generator):**
```
set.seed(20260531); n=5000; K=8 binary margins; crossed near-conflicting structure:
  margins 1-4 probs c(0.85,0.15); margins 5-8 same latent flipped c(0.15,0.85);
  targets = opposite skew of the sample; max_weight=1000 (≈unbounded, 0 pinned), min_weight=0.
Plan step: VERIFY this yields ≫50 unconstrained iters; if not, tune skew/K until it does, then
RECORD the achieved ω=1 baseline iter count as the test's pre-registered reference (committed
before the ω change). The SHIP criterion (spectral < fixed on the SAME fixture) is independent of
the absolute count, so tuning cannot manufacture a pass.
```

| Fixture | Purpose | Pass criterion |
|---------|---------|----------------|
| SageMath synthetic (Phase 1) | theory gate | free ratio → ρ(M_II)² |
| stepstone mw=5 (early exit) | bounded no-regression | spectral iters ≤ fixed |
| Slow unconstrained (generator) | **SHIP gate** | spectral iters < fixed iters |
| All converged fixtures (§3.6) | convergence safety | none flip → NOCONV |
| Fixed-point check, all | correctness | weights 1e-8, marginal_kl 1e-9 vs baseline |
| R `devtools::test()` + Python parity (mode-2 exercised) | no regression | 0 fail; rtol=1e-6 |

**SHIP** = Phase-1 GO AND spectral < fixed on the slow-unconstrained fixture AND mw=5 no-regression
AND no converged fixture → NOCONV AND per-iter wall-clock not regressed beyond the iteration win.
**NO-SHIP** = any fail → revert ORIS, default → mode 1, commit derivation + documented second
negative (an acceptable, expected outcome). A **tie** in the bounded production regime (mw 3–5) is
pre-declared **NO-SHIP** (parity does not justify the per-iter cost — CTO BLOCKING-2).

**Alternatives considered (CTO BLOCKING-3):** (a) **ship nothing / documented negative** — lowest
cost; chosen as the NO-SHIP fallback and a legitimate outcome given SOR-off-by-default. (b) Adaptive
PSOR (Wolfe step adaptation) — needs per-step gradients ORIS lacks; deferred (§7). (c) generic LCP
solver — wrong layer, rejected. The phased derivation-first path is the chosen EV-positive route:
the expensive ORIS work is gated behind a cheap decisive symbolic check.

## 6. Software (NotebookLM-selected)
| Tool | Role | Why |
|------|------|-----|
| **SageMath** | symbolic `M_II` derivation + numerical verify | one free tool; SymPy/NumPy-backed (synthesis top pick) |
| NumPy/SciPy | iteration simulation in artifact | already a dep |
| Julia PATHSolver / Complementarity.jl | *optional* ground-truth active set | only if synthetic active set ambiguous |
| Lean/Coq/Isabelle | — | rejected: rate proof is overkill |

## 7. Out of scope (this cycle)
- **SRAA adapt path** (`oris.cpp:857–893`): global `err_rp` + Anderson confound. Follow-up; relates
  to leafblower-e65t.1 (SRAA×over-relaxation).
- **oris_soft / ALM mode** (`alm_active`): pinning ill-defined inside the soft box.
- **Adaptive PSOR**: deferred fallback (needs gradients ORIS lacks).
- `omega_max`/fixed-mode removal — only after a successful spectral ship (user directive).

## 8. Risks
- **Per-iter cost:** §3.2 parallel accumulator is O(M_cell); the `is_pinned` write is in the clamp
  loop, the free accumulation in the (gated) sweep — both already-visited loops. Ship gate includes a
  wall-clock-not-regressed check (§5).
- **Active-set churn:** non-stationary θ₂(I); gates 7 (EMA decay/re-warm) + 10 (monotone fallback) +
  §3.6 termination invariant bound it. Persistent churn that merely fails to beat fixed = clean
  NO-SHIP.
- **Third negative:** may, like mj1p.2, not beat fixed in any shippable regime. §0/§5 make that a
  documented, EV-positive outcome (the derivation stands).

## 9. Commit coherence
- Phase-1: SageMath artifact + result table together (gate evidence).
- Phase-2 (if reached): `oris.cpp` change + new fixture/generator + observability fields +
  result-ABI bump + roxygen/Python docstring/oris.md/NEWS + Python parity test exercising mode-2 —
  all one commit cycle (CLAUDE.md: docs/tests with code).
