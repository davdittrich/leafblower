# src/oris.cpp — Verified Structure Map

File: `/home/dd/Gemini/leafblower/src/oris.cpp` (1860 lines, 102257 bytes, 0 CR — LF-only).
All line numbers + quotes verified against a byte-exact copy `/tmp/oris_raw.cpp`.

---

## 1. SRAA vs flat branch structure

- `const bool sraa_active_lvl = st.accelerate;` declared at **L802**.
- SRAA branch opens: `if (sraa_active_lvl) {` at **L822**.
- Flat (non-SRAA) branch opens: `if (!sraa_active_lvl) {` at **L1050**.
- Flat branch closes: `}  // end if (!sraa_active_lvl)` at **L1840**.

NOT a single if/else. They are two SEPARATE guarded blocks at the same scope:
- SRAA body: **L822 … ~L995** (the `if (sraa_active_lvl)` block; ends before the
  shared "Clamp cell masses" sync at L997+).
- Between the two: shared cell-mass clamp/sync (L997–L1046) + scratch decls
  (`per_k_errRp_cache` L1047, `per_k_errRp_valid` L1048).
- Flat body: **L1050 … L1840** (`if (!sraa_active_lvl) { … }`).

Note: inside the flat body there are still inner `if (sraa_active_lvl)` guards
(L1209, L1385, L1481) — dead under `!sraa_active_lvl` but present for the shared
per-margin sweep code path.

After both: `oris_finalize(st, res, ct, X, X_init, L_cell, U_cell, …)` at **L1853**.

---

## 2. Flat-path omega adaptation

Block: **L1568 … L1641** ("BLOCK 1b — MARGINAL_KL best-iterate" header L1568;
omega update loop body L1594–L1639).

Per-margin errRp recompute / reuse:
```
1594:  double errRp_k;
1595:  if (per_k_errRp_valid) {
1597:      errRp_k = per_k_errRp_cache[k];   // reuse from convergence sweep (773f.7)
1600:  } else { ... recompute via S_lin loop ... }
```
sor_omega[k] updates (verbatim):
```
1612:  bool decreasing = (errRp_k < sor_prev_errRp[k]);
1613:  bool sign_flip  = !decreasing && sor_prev_decreasing[k];
1616:      sor_omega[k] = std::max(omega_min_v, sor_omega[k] * kSorOscillationDamp);
1622:      double theta2_raw_k = estimate_theta2(sor_prev_errRp[k], errRp_k);
1624:      sor_theta2_ema[k] = (sor_theta2_ema[k] < 0.0)
1626:          : kSorEmaAlpha * theta2_raw_k + (1.0 - kSorEmaAlpha) * sor_theta2_ema[k];
1628:      sor_omega[k] = omega_from_theta2(sor_theta2_ema[k], kSorSpectralCeiling);
1631:      sor_omega[k] = omega_max_v;
1634:      sor_omega[k] = std::min(omega_max_v, sor_omega[k] * kSorRecoveryGrowth);
1637:  if (sor_omega[k] < sor_min_omega) sor_min_omega = sor_omega[k];
1638:  sor_prev_decreasing[k] = decreasing;
1639:  sor_prev_errRp[k]      = errRp_k;
```
- Oscillation damp: L1616 `sor_omega[k] * kSorOscillationDamp` (kSorOscillationDamp=0.7, L104).
- Recovery toward omega_max: L1634 `std::min(omega_max_v, sor_omega[k] * kSorRecoveryGrowth)`.
- Spectral/theta2 path: L1622–L1628.

The omega used by the NEXT sweep is read at L569 (`eff_omega = sor_omega[k]`,
linear) and L665 (`eff_omega_log = sor_omega[k]`, log). Adaptation writes here, consume there.

---

## 3. Per-margin residual computation (FLAT path)

errRp computed **INLINE** in oris.cpp, in the convergence sweep at **L1504–L1549**:
```
1504:  double errRp = 0.0;
1512:  std::fill(S_lin.begin(), S_lin.begin() + nj, 0.0);
1516:  if (j >= 0 && j < nj) S_lin[j] += X[c];      // S_lin aggregation INLINE
1518:  double errRp_k = 0.0;                          // 773f.7
1520:  double e = std::fabs(S_lin[j] / W_total - st.targets[k][j]);
1521:  if (e > errRp)   errRp   = e;
1522:  if (e > errRp_k) errRp_k = e;                  // per-margin capture
1524:  per_k_errRp_cache[k] = errRp_k;                // cached for omega block
1549:  per_k_errRp_valid = (use_linear && W_total > 0.0);
```
- formula `errRp_k = max_j |S_lin[j]/W_total - target[k][j]|` is computed inline (L1518–L1522).
- `S_lin[j]` aggregation is built INLINE in the oris hot loop (L1516: `S_lin[j] += X[c]`).
  `S_lin` itself declared at L359: `std::vector<double> S_lin(max_cat, 0.0);` (function-scope scratch).

Two file-LOCAL lambda helpers also exist (used in the per-margin SWEEP, not the
final convergence pass), defined inline in oris.cpp:
- `auto compute_margin_errRp_linear = [&](int k) -> double {` at **L1073** (body L1073–L1091; S_lin built inside, L1075/L1082).
- `auto compute_margin_errRp_log = [&](int k) -> double {` at **L1093**.
Callers (oris.cpp only): L1130, L1146 (linear); L1260, L1280 (log).

NO external/shared helper named `compute_errRp_ct`. The residual is NOT shared
with raking.cpp/sinkhorn.cpp — these are `[&]` lambdas local to `oris_solve`,
cannot be referenced cross-TU. The only SHARED metrics helper is
`lbw::compute_cell_metrics(st, ct, X, W_total, S_lin)` (called L1694, and in SRAA
branch L931/L1694-context) which recomputes full metrics incl. errRp but is the
EXTRA-metrics path, separate from the inline per-margin errRp_k.

---

## 4. Cell mass + bounds in scope at flat omega-adaptation site (L1594–L1639)

In scope and live:
- `X[c]` — cell mass; used in the immediately-preceding errRp loop (L1516 `S_lin[j] += X[c]`)
  and the else-recompute branch (L1604 `if (...) S_lin[j] += X[c]`).
- `L_cell[c]`, `U_cell[c]` — declared at **L158**: `std::vector<double> L_cell(ct.M_cell), U_cell(ct.M_cell);` (function scope).
  - `L_cell[c] = lo * ct.n_per_cell[c];` set once at L161 (level-invariant, min_weight-only).
  - `U_cell[c]` recomputed per homotopy level: L516 comment, L519 `U_cell[c] = hi * ct.n_per_cell[c];`.
  - Last USE before the omega site: clamps at L1352/L1362 (X_tilde_c path), L1434/L1444
    (X_tilde[] path), L1471 (`std::max(X[c]-U_cell[c], L_cell[c]-X[c])`). All
    live and in scope at L1594.

So at the flat omega-adaptation site, `X[c]`, `L_cell[c]`, `U_cell[c]` are ALL in
scope (function-scope vectors), but the omega block itself only reads `X[c]`
(via S_lin); it does not currently read L_cell/U_cell.

---

## 5. theta2 / EMA / spectral (mj1p.2)

Constant: **L414** `static constexpr double kSorSpectralCeiling = 1.99;  // Thibault 2021: convergence strict <2`.

Lambdas (function-scope, defined once):
- **L415** `auto estimate_theta2 = [](double prev, double curr) -> double {` (capture-less; body L415–~423; returns -1 when non-informative per L407).
- **L425** `auto omega_from_theta2 = [&](double theta2, double ceiling) -> double {` (falls back to ceiling on non-informative theta2 per L408/L424).

EMA state: **L435** `std::vector<double> sor_theta2_ema(st.K, -1.0);`
(comment L434: "SRAA uses index [0] (global errRp); flat BCD uses per-margin index [k].")

Spectral mode-2 logic appears in BOTH branches:
- SRAA branch: L872–L884 (`estimate_theta2(sor_prev_errRp[0], curr_errRp)` L875;
  EMA update L877–L879; `omega_from_theta2(sor_theta2_ema[0], kSorSpectralCeiling)` L882; uses index **[0]**).
- Flat branch: L1622–L1628 (`estimate_theta2(sor_prev_errRp[k], errRp_k)` L1622;
  EMA L1624–L1626; `omega_from_theta2(sor_theta2_ema[k], kSorSpectralCeiling)` L1628; uses index **[k]**).

---

## 6. Solve-local SOR state declarations

All function-scope `std::vector` (solve-local, NOT struct fields):
- **L430** `std::vector<double> sor_omega(st.K, omega_init_v);`
- **L431** `std::vector<double> sor_prev_errRp(st.K, std::numeric_limits<double>::infinity());`
- **L432** `std::vector<bool>   sor_prev_decreasing(st.K, false);`
- **L435** `std::vector<double> sor_theta2_ema(st.K, -1.0);`

Confirmed solve-local std::vector sized to `st.K`; none are `st.`/`res.`/struct members.
