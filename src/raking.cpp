// raking.cpp — bounded raking via cyclic IPF with water-filling box projections.
//
// Algorithm: F_eval performs one complete IPF sweep with per-category water-filling
// (autumn single_adjust) that enforces X[c] ∈ [L_c, U_c] within each margin step.
// No Dykstra correction vectors; bounds are enforced inline by water_fill_cat.
//
// Both the flat loop and SRAA-m-accelerated path call F_eval directly.
// Water-filling convergence: cyclic KL projections onto bounded margin constraints
// converge to the bounded KL minimum (Csiszar-Tusnady 1984).
//
// References (classical IPF family):
//   Deming W. E. & Stephan F. F. (1940), "On a Least Squares Adjustment of
//     a Sampled Frequency Table", Ann. Math. Stat. 11, 427-444.
//   Csiszar I. (1975), "I-Divergence Geometry of Probability Distributions
//     and Minimization Problems", Ann. Probab. 3, 146-158.

#include "lbw_config.h"
#include "lbw_math.hpp"
#include "raking.hpp"
#include "cell_table.hpp"
#include "calib_dispatch.hpp"
#include "sraa.hpp"
#include "leafblower.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <limits>
#include <numeric>
#include <vector>

namespace lbw {

// Constrained raking solver: cyclic IPF with per-category water-filling box projections.
// Each margin step applies water_fill_cat() — KL projection onto
// {Σ_c Xv[c]=T_kj, L_c≤Xv[c]≤U_c} (Csiszar-Tusnady 1984; autumn single_adjust).
// No Dykstra correction vectors; stateless F enables SRAA-m Anderson Acceleration.
// inner_max_iter is the single iteration budget; outer_max_iter is unused.
RakingResult raking_solve(CalibState& st) {
    static constexpr double kAbsoluteZeroThreshold = 1e-15;  // bucket_j / free_sum is genuinely zero
    static constexpr double kRelativeZeroFraction  = 1e-15;  // bucket[j] < fraction * W_total
    static constexpr int    kErrCheckInterval     = 10;
    static constexpr int    kMaxNoImprove         = 5;
    // CR-C5b (jy0m): feasibility bar for accepting a metric-plateau as RK_OK.
    // errRp = max_kj |S_p − T| (compute_cell_metrics: fabs(S_p−T)) is the max
    // ABSOLUTE proportion-point margin miss, NOT a relative residual. A calibration
    // is USABLE when every marginal constraint is met to within 1 percentage point —
    // a standard practical acceptance tolerance for survey weighting. This clears the
    // structural floor of feasible problems (stepstone K=9 flat plateaus at
    // errRp≈1.6–2.3e-3; near-identity fixtures at ≈4e-4) with margin, while catching
    // bounds-blocked plateaus (errRp≈0.1–0.5). NOT tol_abs (1e-6): that would force
    // feasible runs to over-iterate and mislabel usable structural-floor results as
    // STALL. LIMITATION (jy0m.1): being absolute, a rare-category infeasibility whose
    // pp miss is < 1pp (e.g. target 0.005 pinned at 0.014) still passes — tracked.
    static constexpr double kRakingFeasTol        = 1e-2;  // 1 percentage point

    RakingResult res;
    res.base.status     = RK_ERR_BUDGET;  // initial; overwritten by criterion/stall; remains if budget exhausted
    res.base.iterations = 0;
    res.base.max_error  = 1.0;

    // Build cell table + initial masses + bounds.
    CellTable ct;
    std::vector<double> X_init;
    double hi_eff;
    std::vector<double> L_cell, U_cell;
    if (lbw::solver_setup_ct_base(st, ct, X_init, hi_eff, L_cell, U_cell, res) != RK_OK)
        return res;

    // Working copy of cell masses.
    std::vector<double> X(X_init);

    // Per-(margin, category) cell index lists for water-filling.
    // cells_per_cat[k][j] = cells where g_per_cell[k][c] == j.
    // Built once: O(M_cell × K). Memory: ~M_cell × K ints.
    auto cells_per_cat = lbw::build_cells_per_cat(ct, st.K, st.cat_counts);

    // Pre-allocated scratch for water-filling inner loop.
    int wf_max_cat = 0;
    for (int k = 0; k < st.K; k++)
        for (int j = 0; j < st.cat_counts[k]; j++)
            wf_max_cat = std::max(wf_max_cat, (int)cells_per_cat[k][j].size());
    std::vector<double>  wf_x_orig(wf_max_cat);
    std::vector<uint8_t> wf_status(wf_max_cat);

    const double hi = hi_eff;

    // No Dykstra correction vectors: water-filling enforces bounds within F_eval.

    bool is_infeasible = false;

    // Pre-computed reciprocals: n_per_cell is constant; avoids division per cell per super-step.
    std::vector<double> inv_n_per_cell(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++)
        inv_n_per_cell[c] = 1.0 / static_cast<double>(ct.n_per_cell[c]);

    int max_cats = lbw::max_cats_count(st.K, st.cat_counts);
    std::vector<double> bucket(max_cats);

    // Descent monitor — tracks solver loss function (weight KL for flat loop, errRp for SRAA-m).
    double min_loss_window = std::numeric_limits<double>::infinity();
    int n_no_improve = 0;

    // pct/l1 tracking at cell level
    std::vector<double> X_prev(X);

    // Convergence rule state
    double prev_metric_for_rule = std::numeric_limits<double>::infinity();

    // G8b: best-iterate tracking via BestIterTracker (replaces ad-hoc vars).
    // best.best_weights stores cell-level X snapshot at the best observed metric.
    BestIterTracker best;

    // Scratch buffers for compute_weight_kl vectorized log path.
    std::vector<double> kl_ratio_scratch(ct.M_cell, 0.0);
    std::vector<double> kl_weight_scratch(ct.M_cell, 0.0);

    // SOR: wire st.sor_cfg into raking's IPF step (same API as oris).
    // Apply only when bounds are active (oscillation risk).
    const bool sor_active  = st.sor_cfg.enabled &&
                              (st.min_weight > 0.0 || hi < kUnboundedSentinel);
    const bool sor_auto    = st.sor_cfg.auto_adapt;
    const double omega_min = st.sor_cfg.omega_min;    // default 0.3
    const double omega_init = st.sor_cfg.omega_init;  // default 1.0
    std::vector<double> sor_omega(st.K, omega_init);
    std::vector<double> sor_prev_errRp(st.K, std::numeric_limits<double>::infinity());

    // Greedy scheduler: per-margin residuals + sort order.
    // errRp_k[k] updated each iter during sweep using bucket[] sums.
    // Survives across iterations — first iter uses uniform priority (1/K init).
    std::vector<double> errRp_k(st.K, 1.0 / st.K);
    std::vector<int> margin_order(st.K);
    std::iota(margin_order.begin(), margin_order.end(), 0);
    const bool use_greedy = (st.scheduler.mode == SchedulerMode::GREEDY);
    // R8: greedy reordering breaks SRAA-m's fixed-point geometry; demote silently.
    bool use_greedy_effective = use_greedy;
    if (st.accelerate && use_greedy_effective) {
        use_greedy_effective = false;
        res.sraa_demoted = true;
        st.log("[raking] greedy scheduler disabled under SRAA-m acceleration; using round_robin");
    }


    // water_fill_cat: KL projection of category j (margin k) onto
    // {Σ_{c∈j} Xv[c] = T_kj, L_cell[c] ≤ Xv[c] ≤ U_cell[c]}.
    // Algorithm: autumn single_adjust() (rake.R:63-93) at cell level.
    // wf_x_orig[] and wf_status[] are pre-allocated scratch (re-used every call).
    auto water_fill_cat = [&](int k, int j, double T_kj, double bucket_j,
                               std::vector<double>& Xv) {
        const auto& cells = cells_per_cat[k][j];
        const int n = static_cast<int>(cells.size());
        if (n == 0) return;
        if (bucket_j < kAbsoluteZeroThreshold) {
            if (T_kj > 0.0) {
                is_infeasible = true;
                // Zero cells to lower bounds — do not leave stale values that
                // survive post-normalization and silently violate the constraint.
                for (int ci = 0; ci < n; ++ci)
                    Xv[cells[ci]] = L_cell[cells[ci]];
            }
            return;
        }

        // Save original weights and mark all cells free
        for (int ci = 0; ci < n; ci++) {
            wf_x_orig[ci] = Xv[cells[ci]];
            wf_status[ci] = 0;
        }

        double clamped_sum = 0.0;
        double free_sum    = bucket_j;  // Σ_{free} X_orig[c]

        // Bound is pass < n (not <= n): at most n clamping passes are possible
        // before the free pool empties (each pass that clamps ≥1 cell shrinks
        // the pool by ≥1).  The only clean-success exit is the !any_clamped
        // commit-and-return inside the loop.  Pass index n is never needed for
        // a real commit; it only re-applies the eps-infeasibility guard
        // (free_sum < kAbsoluteZeroThreshold) which can spuriously flag
        // 0 < free_sum < eps as infeasible — a case the post-loop handler
        // (free_sum <= 0.0) would correctly treat as best-effort.
        for (int pass = 0; pass < n; ++pass) {
            if (free_sum < kAbsoluteZeroThreshold) { is_infeasible = true; break; }
            const double T_free = T_kj - clamped_sum;
            if (T_free <= 0.0) break;
            const double m = T_free / free_sum;

            bool any_clamped = false;
            for (int ci = 0; ci < n; ci++) {
                if (wf_status[ci] != 0) continue;
                const double proposed = wf_x_orig[ci] * m;
                if (proposed > U_cell[cells[ci]]) {
                    wf_status[ci] = 1;
                    clamped_sum += U_cell[cells[ci]];
                    free_sum    -= wf_x_orig[ci];
                    any_clamped  = true;
                } else if (proposed < L_cell[cells[ci]]) {
                    wf_status[ci] = 2;
                    clamped_sum += L_cell[cells[ci]];
                    free_sum    -= wf_x_orig[ci];
                    any_clamped  = true;
                }
            }

            if (!any_clamped) {
                // m applies cleanly — commit final values and return
                for (int ci = 0; ci < n; ci++) {
                    const int c = cells[ci];
                    if      (wf_status[ci] == 0) Xv[c] = wf_x_orig[ci] * m;
                    else if (wf_status[ci] == 1) Xv[c] = U_cell[c];
                    else                          Xv[c] = L_cell[c];
                }
                return;
            }
        }
        // All passes exhausted (infeasible category) — commit best-effort values
        if (free_sum <= 0.0) {
            is_infeasible = true;
            return;  // cannot redistribute — leave cells at bounds
        }
        const double T_final = T_kj - clamped_sum;
        const double m_final = (T_final > 0.0) ? T_final / free_sum : 0.0;
        for (int ci = 0; ci < n; ci++) {
            const int c = cells[ci];
            if      (wf_status[ci] == 0) Xv[c] = wf_x_orig[ci] * m_final;
            else if (wf_status[ci] == 1) Xv[c] = U_cell[c];
            else                          Xv[c] = L_cell[c];
        }
    };

    // last_F_metrics: populated by F_eval via compute_cell_metrics; reused in SRAA-m convergence check.
    // Declared here (outer scope) so the [&] lambda capture includes it.
    lbw::CellMetrics last_F_metrics;

    // 773f.6: controls whether F_eval runs full compute_cell_metrics or fast errRp-only path.
    // Set true before check-interval iterations and always for SRAA (needs full metrics).
    bool f_eval_full_metrics = false;

    // F_eval: one complete bounded IPF iteration using water-filling.
    // Water-filling enforces X[c] ∈ [L[c], U[c]] within each margin step.
    // No correction vectors — F is stateless: enables SRAA-m L2 halving.
    // Mathematical basis: cyclic KL projections onto bounded margin constraints
    // converge to the bounded KL minimum (Csiszar-Tusnady 1984).
    auto F_eval = [&](std::vector<double>& Xv) -> double {
        double W_total = 0.0;
        for (int c = 0; c < ct.M_cell; c++) W_total += Xv[c];

        if (use_greedy_effective)
            std::sort(margin_order.begin(), margin_order.end(),
                      [&](int a, int b){ return errRp_k[a] > errRp_k[b]; });

        for (int ki = 0; ki < st.K; ki++) {
            const int k = use_greedy_effective ? margin_order[ki] : ki;

            // Compute pre-water-fill bucket (needed for errRp_k and SOR)
            std::fill(bucket.begin(), bucket.begin() + st.cat_counts[k], 0.0);
            for (int c = 0; c < ct.M_cell; c++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k]) bucket[g] += Xv[c];
            }

            // errRp_k for Greedy: pre-water-fill residual (how bad margin is NOW).
            // Must be computed BEFORE water-fill — post-fill residual is always near zero.
            if (W_total > 0.0) {
                double ek = 0.0;
                for (int j = 0; j < st.cat_counts[k]; j++) {
                    double e = std::fabs(bucket[j] / W_total - st.targets[k][j]);
                    if (e > ek) ek = e;
                }
                errRp_k[k] = ek;
            }

            // SOR effective omega (under-relaxation toward current bucket)
            const double eff_omega = sor_active ? sor_omega[k] : 1.0;

            // Apply water-filling per category
            for (int j = 0; j < st.cat_counts[k]; j++) {
                double Tkj = st.targets[k][j] * W_total;
                if (bucket[j] < kRelativeZeroFraction * W_total) {
                    if (Tkj > 0.0) {
                        is_infeasible = true;
                        // Clamp cells to lower bounds — stale Xv values would survive
                        // the post-normalization step and silently violate the margin.
                        const auto& cells_kj = cells_per_cat[k][j];
                        for (int ci = 0; ci < (int)cells_kj.size(); ++ci)
                            Xv[cells_kj[ci]] = L_cell[cells_kj[ci]];
                        // 67sk: accumulate post-clamp bucket[j] so SOR ek-loop sees correct value
                        bucket[j] = 0.0;
                        for (int ci2 = 0; ci2 < (int)cells_kj.size(); ci2++)
                            bucket[j] += L_cell[cells_kj[ci2]];
                    } else {
                        // 67sk: zero-target cells clamped to L_cell; update bucket[j] for ek-loop
                        bucket[j] = 0.0;
                        const auto& cells_kj = cells_per_cat[k][j];
                        for (int ci2 = 0; ci2 < (int)cells_kj.size(); ci2++)
                            bucket[j] += L_cell[cells_kj[ci2]];
                    }
                    continue;
                }
                if (eff_omega != 1.0) {
                    // SOR: under-relax target — equivalent to pow(scale, omega) at cell level
                    const double s0 = Tkj / bucket[j];
                    Tkj = bucket[j] * std::pow(s0, eff_omega);
                }
                water_fill_cat(k, j, Tkj, bucket[j], Xv);
                // 67sk: inline bucket[j] update replaces aggregate_to_margin for SOR ek-loop
                bucket[j] = 0.0;
                {
                    const auto& cells_kj_wf = cells_per_cat[k][j];
                    for (int ci2 = 0; ci2 < (int)cells_kj_wf.size(); ci2++)
                        bucket[j] += Xv[cells_kj_wf[ci2]];
                }
            }

            // SOR adaptation: post-water-fill residual for omega adjustment
            if (sor_active && sor_auto && W_total > 0.0) {
                double ek = 0.0;
                for (int j = 0; j < st.cat_counts[k]; j++) {
                    double e = std::fabs(bucket[j] / W_total - st.targets[k][j]);
                    if (e > ek) ek = e;
                }
                if (sor_prev_errRp[k] < ek)
                    sor_omega[k] = std::max(omega_min, sor_omega[k] * 0.7);
                else
                    sor_omega[k] = std::min(1.0, sor_omega[k] * 1.05);
                sor_prev_errRp[k] = ek;
            }
        }

        // Hyperplane: normalize sum to n (no correction vector needed;
        // water-fill preserves category sums, so this is near-no-op for feasible problems)
        double s_hp = 0.0;
        for (int c = 0; c < ct.M_cell; c++) s_hp += Xv[c];
        if (s_hp > 0.0) {
            const double sc_hp = static_cast<double>(st.n) / s_hp;
            for (int c = 0; c < ct.M_cell; c++) Xv[c] *= sc_hp;
        }

        // 773f.6: full metrics only on check iterations; fast errRp-only otherwise.
        if (f_eval_full_metrics) {
            last_F_metrics = lbw::compute_cell_metrics(st, ct, Xv, static_cast<double>(st.n), bucket);
            return last_F_metrics.errRp;
        }
        // Fast path: compute errRp only (single O(K×M_cell) pass, no chi2/kl/grake).
        double errRp_fast = 0.0;
        {
            const double W_fast = static_cast<double>(st.n);
            for (int k = 0; k < st.K; k++) {
                const int nj = st.cat_counts[k];
                std::fill(bucket.begin(), bucket.begin() + nj, 0.0);
                for (int c = 0; c < ct.M_cell; c++) {
                    int g = ct.g_per_cell[k][c];
                    if (g >= 0 && g < nj) bucket[g] += Xv[c];
                }
                for (int j = 0; j < nj; j++) {
                    double e = std::fabs(bucket[j] / W_fast - st.targets[k][j]);
                    if (e > errRp_fast) errRp_fast = e;
                }
            }
        }
        return errRp_fast;
    };

    // SRAA-m state (replaces SRAA-m step-halving while-loop)
    lbw::SRAAState rk_sraa;
    if (st.accelerate) rk_sraa.init(ct.M_cell, lbw::kSRAAm);

    if (st.accelerate) {
        // 773f.6: SRAA always needs full metrics (last_F_metrics reused in convergence check below).
        f_eval_full_metrics = true;
        // SRAA-m replaced by SRAA-m. X_prev_sq (stall detection) removed — SRAA-m-specific.
        // Loop mirrors greenkhorn: sraa_step called once per outer iteration until budget exhausted.
        int f_eval_budget = st.inner_max_iter;
        int f_evals_used = 0;
        int rk_outer_stall_count = 0;
        while (f_evals_used + 1 <= f_eval_budget) {
            rk_sraa.F_cur = X;  // seed F_cur with current X before each sraa_step call
            // B3 fix: save metrics before sraa_step — if AA is attempted then rejected,
            // the second f_eval(X_AA) overwrites last_F_metrics with the rejected iterate.
            // Restore on rejection so convergence check uses metrics from accepted X_k.
            lbw::CellMetrics saved_metrics = last_F_metrics;
            auto r = lbw::sraa_step(F_eval, X, L_cell, U_cell, rk_sraa);
            if (!r.aa_accepted && r.f_evals == 2) last_F_metrics = saved_metrics;
            f_evals_used += r.f_evals;
            // hhmk: use last_F_metrics.errRp (accepted iterate) instead of removed
            // r.err_rp fast-proxy. select_metric() honors cfg.metric.
            const double sraa_metric_val = lbw::select_metric(st.convergence_cfg.metric, last_F_metrics);
            res.base.max_error  = last_F_metrics.errRp;
            res.base.iterations = f_evals_used;

            // CR-C6 (kxna.6): l1 weight change on the SRAA path (last_F_metrics has
            // no l1 — obs-level, not computed by compute_cell_metrics). Mirrors the
            // flat loop (:437-442) reusing the function-scope X_prev; final value is
            // the last accepted step's Σ|ΔX|/n.
            {
                double l1_sraa = 0.0;
                for (int c = 0; c < ct.M_cell; c++) l1_sraa += std::fabs(X[c] - X_prev[c]);
                res.base.l1_weight_change = l1_sraa / static_cast<double>(st.n);
                for (int c = 0; c < ct.M_cell; c++) X_prev[c] = X[c];
            }

            // Best-iterate tracking (uses cfg-specified metric, not max_err proxy)
            if (sraa_metric_val < best.best_metric) {
                best.update(sraa_metric_val,
                            lbw::compute_weight_kl(X, X_init, ct.M_cell, st.n, kl_ratio_scratch.data(), kl_weight_scratch.data()),
                            f_evals_used, X);
            }

            // Outer revert: if quality has degraded significantly, revert to best
            {
                double curr_quality_rk = sraa_metric_val;
                if (curr_quality_rk > best.best_metric * (1.0 + lbw::kSRAAOuterSlack)) {
                    if (++rk_outer_stall_count >= lbw::kSRAAOuterStallWindow) {
                        X = best.best_weights;    // revert; no S_flat rebuild — F_eval handles it
                        rk_sraa.clear();
                        rk_outer_stall_count = 0;
                    }
                } else { rk_outer_stall_count = 0; }
            }

            // Convergence check
            {
                lbw::CellMetrics m_conv = last_F_metrics;
                // m_conv.errRp already populated by F_eval (accepted iterate);
                // no need to override with proxy.
                if (lbw::check_convergence(st.convergence_cfg, m_conv,
                                           prev_metric_for_rule, st.tol_abs)) {
                    // CR-C5b (jy0m): same feasibility gate as the flat loop — a
                    // scale-blind metric plateau is only RK_OK if margins are met
                    // (errRp ≤ 1pp) or the user's explicit absolute tol is met.
                    // A bounds-blocked plateau here does NOT mark_converged; the SRAA
                    // loop then runs to f_evals exhaustion → BUDGET (this path has no
                    // STALL emitter; promoting it to STALL(5) like the flat path is
                    // tracked in jy0m.2). Either way the false RK_OK is removed.
                    // PRE-EXISTING GAP (jy0m.3, NOT introduced here): m_conv=last_F_metrics
                    // and compute_cell_metrics never writes .l1, so with metric="l1_weight"
                    // (or "pct" alias) + absolute_tol>0, user_abs_met reads a phantom 0 and
                    // is spuriously true — bypassing this gate. check_convergence's own
                    // c_abs fires on the same phantom, so main already false-OK'd that
                    // config; this gate is neutral, not a regression. Fixed via l1 in
                    // last_F_metrics under jy0m.3. Flat site is immune (m_conv.l1 real).
                    const bool user_abs_met =
                        st.convergence_cfg.absolute_tol > 0.0 &&
                        lbw::select_metric(st.convergence_cfg.metric, m_conv)
                            < st.convergence_cfg.absolute_tol;
                    if (last_F_metrics.errRp <= kRakingFeasTol || user_abs_met) {
                        lbw::mark_converged(res, st.convergence_cfg, f_evals_used, st.tol_abs);
                        break;
                    }
                }
            }

            if (st.verbose >= 1) {
                char msg[256];
                std::snprintf(msg, 256, "raking[sraa] f_evals=%d errRp=%.2e aa=%d",
                              f_evals_used, last_F_metrics.errRp, (int)r.aa_accepted);
                st.log(msg);
            }
        }
    } else {
        for (int iter = 1; iter <= st.inner_max_iter; iter++) {
            res.base.iterations = iter;

            f_eval_full_metrics = (iter == 1 || iter % kErrCheckInterval == 0 || iter == st.inner_max_iter);
            double errRp = F_eval(X);

            if (f_eval_full_metrics) {
                res.base.max_error = errRp;

                // 6uyk.2: compute wkl once per check iteration; reuse for best-iterate
                // tracking (MAX_ERR and non-MAX_ERR branches) and KL stall detection.
                const double wkl_flat = lbw::compute_weight_kl(X, X_init, ct.M_cell, st.n, kl_ratio_scratch.data(), kl_weight_scratch.data());

                // Best-iterate tracking (MAX_ERR metric)
                if (st.convergence_cfg.metric == lbw::CalibMetric::MAX_ERR) {
                    if (errRp < best.best_metric) {
                        best.update(errRp, wkl_flat, iter, X);
                    }
                }

                // L1 weight change from previous check
                double l1_weight = 0.0;
                for (int c = 0; c < ct.M_cell; c++)
                    l1_weight += std::fabs(X[c] - X_prev[c]);
                l1_weight /= static_cast<double>(st.n);
                for (int c = 0; c < ct.M_cell; c++) X_prev[c] = X[c];

                // Extra metrics needed for: non-MAX_ERR convergence, OR final
                // iter (full diagnostics), OR MAX_ERR with chi2 propagation
                // (leafblower-a605: chi2 was stale on early-converged MAX_ERR
                // when last_F_metrics-write iter ≠ converging iter).
                const lbw::CalibMetric metric = st.convergence_cfg.metric;
                double mean_err = 0.0, kl_max = 0.0, chi2_total = 0.0, grake_norm = 0.0;
                double marginal_kl_sum = 0.0;
                const bool need_extra = (metric == lbw::CalibMetric::MAX_ERR ||
                                         metric == lbw::CalibMetric::MEAN_ERR ||
                                         metric == lbw::CalibMetric::KL ||
                                         metric == lbw::CalibMetric::MARGINAL_KL ||
                                         metric == lbw::CalibMetric::CHI2 ||
                                         metric == lbw::CalibMetric::GRAKE_NORM ||
                                         iter == st.inner_max_iter);
                if (need_extra) {
                    // xc1s.4: reuse F_eval's full-metrics result instead of a second
                    // O(K·M_cell) compute_cell_metrics pass. F_eval already computed
                    // last_F_metrics for this same post-hyperplane X, normalized by
                    // W=st.n — matching errRp's basis (the prior recompute used
                    // W=ΣX≈n, a ~1e-13 inconsistency; max_error is unaffected).
                    const lbw::CellMetrics& m = last_F_metrics;
                    mean_err   = m.mean_err;
                    kl_max     = m.kl;
                    chi2_total = m.chi2;
                    grake_norm = m.grake_norm;
                    marginal_kl_sum = m.marginal_kl;
                    if (metric != lbw::CalibMetric::MAX_ERR) {
                        const double curr_best = lbw::select_metric(metric, m);
                        if (std::isfinite(curr_best) && curr_best < best.best_metric) {
                            best.update(curr_best, wkl_flat, iter, X);
                        }
                    }
                }

                res.base.mean_error       = mean_err;
                res.base.kl               = kl_max;
                res.base.chi2             = chi2_total;
                res.base.grake_norm       = grake_norm;
                res.base.l1_weight_change = l1_weight;

                if (st.verbose >= 1) {
                    char msg[256];
                    std::snprintf(msg, 256, "raking iter %d: errRp=%.2e", iter, errRp);
                    st.log(msg);
                }

                lbw::CellMetrics m_conv;
                m_conv.errRp = errRp; m_conv.mean_err = mean_err;
                m_conv.kl = kl_max; m_conv.chi2 = chi2_total;
                m_conv.grake_norm = grake_norm; m_conv.l1 = l1_weight;
                m_conv.marginal_kl = marginal_kl_sum;
                if (lbw::check_convergence(st.convergence_cfg, m_conv,
                                           prev_metric_for_rule, st.tol_abs)) {
                    // CR-C5b (jy0m): the improvement/pct rule fires on a metric
                    // PLATEAU, which is scale-blind — it cannot tell a feasibly-
                    // converged run (margins met) from a bounds-blocked one whose
                    // metric froze because weights are pinned at their bounds
                    // (X≈X_init from iter 1, residual never falls). Mirroring the
                    // logit eb79.22 fix: accepting the plateau as RK_OK requires
                    // feasibility — the max abs marginal margin miss errRp within the
                    // practical bar kRakingFeasTol (1pp). Feasible IPF plateaus with
                    // margins met well inside 1pp (stepstone ≈2.3e-3, near-identity
                    // ≈4e-4), so this never demotes a genuine convergence; a plateau
                    // with materially unmet margins (errRp≫1pp, e.g. bounds-blocked at
                    // 0.1–0.5) falls through to the STALL detector below (constrained
                    // optimum, status=5), instead of a false RK_OK. user_abs_met: an
                    // EXPLICIT user absolute tolerance that is genuinely met is the
                    // user's own feasibility bar — respect it, do NOT override with 1pp
                    // (else convergence=list(absolute=0.05) on errRp=0.03 wrongly STALLs).
                    const bool user_abs_met =
                        st.convergence_cfg.absolute_tol > 0.0 &&
                        lbw::select_metric(st.convergence_cfg.metric, m_conv)
                            < st.convergence_cfg.absolute_tol;
                    if (errRp <= kRakingFeasTol || user_abs_met) {
                        // Converged — do NOT override with INFEAS here. water_fill_cat may
                        // transiently set is_infeasible when cells temporarily hit U_cell during
                        // convergence. INFEAS only overrides on stall (post-loop check below).
                        lbw::mark_converged(res, st.convergence_cfg, iter, st.tol_abs);
                        break;
                    }
                    // else: metric plateaued but margins unmet — fall through to STALL.
                    // prev_metric_for_rule is intentionally left frozen here: the STALL
                    // detector below tracks its OWN independent min_loss_window/n_no_improve,
                    // so termination does not depend on advancing the improvement-rule
                    // baseline; the rule simply re-fires (≈0 improvement) each check and is
                    // re-gated on errRp until n_no_improve trips.
                }

                // Stall detection: track the active convergence metric (not always wkl).
                // wkl is monotone for KL-configured runs (Csiszar-Tusnady); for other metrics
                // (MAX_ERR, MEAN_ERR, etc.) use select_metric so stall detection matches the
                // user's convergence criterion.
                // CR-C5 (kxna.5): removed the `metric==KL && wkl_flat <= tol_abs`
                // RK_OK shortcut. wkl_flat = Σ X·log(X/X_init)/n measures distance
                // FROM the design weights, not the achieved-margin KL, so when
                // max_weight clamps bind (X ≈ X_init ⇒ wkl_flat ≈ 0) it declared
                // RK_OK on grossly unmet margins and bypassed mark_converged. Real
                // convergence is already caught by check_convergence above; unmet
                // margins now correctly fall through to STALL detection below.
                const double active_metric = lbw::select_metric(st.convergence_cfg.metric, m_conv);
                if (!std::isfinite(min_loss_window)) {
                    min_loss_window = active_metric; n_no_improve = 0;
                } else if (active_metric < min_loss_window * (1.0 - st.convergence_cfg.pct_tol)) {
                    min_loss_window = active_metric; n_no_improve = 0;
                } else {
                    n_no_improve++;
                }

                if (n_no_improve >= kMaxNoImprove) {
                    res.base.status           = RK_ERR_STALL;
                    res.base.stall_kind       = 2;  // kl: metric/loss plateau (plain-IPF path)
                    // leafblower-8eod: stall_kind=2 (kl) because this branch is
                    // !st.accelerate (raking.cpp:392); SRAA path does not emit RK_ERR_STALL.
                    if (st.verbose >= 1) {
                        char msg[256];
                        std::snprintf(msg, 256,
                            "raking: stalled for %d consecutive checks "
                            "(metric=%.2e); aborting at iter %d.",
                            n_no_improve, active_metric, iter);
                        st.log(msg);
                    }
                    break;
                }
            }
        }
    }  // end else flat loop

    // leafblower-236s: propagate chi2 from last full-metrics eval.
    // leafblower-24f7 audit (false-positive): this assignment is always fresh.
    //   Flat loop:   convergence is checked only inside the `if (f_eval_full_metrics)`
    //     block, and on that same iter F_eval already refreshed last_F_metrics
    //     (compute_cell_metrics over the current X). res.base.chi2 was also set just
    //     above from chi2_total (the same compute_cell_metrics result), so this overwrite
    //     is redundant but correct — both use the same X and W≈n.
    //   SRAA loop:   f_eval_full_metrics is always true → F_eval refreshes last_F_metrics
    //     every iter; this line is the ONLY chi2 assignment on the SRAA path.
    // Test guard: test-raking-chi2-freshness.R verifies all scenarios.
    if (best.has_best()) {
        // CR-C6 (kxna.6): assign the full diagnostic set from last_F_metrics so the
        // SRAA path returns the same completeness as the flat path (which sets these
        // at :478-482). Redundant-but-equal for the flat path — same X, W≈n as the
        // convergence-iter compute_cell_metrics (see chi2 audit above). l1_weight_change
        // is set per-path (flat :482; SRAA in-loop, as last_F_metrics has no l1).
        res.base.chi2       = last_F_metrics.chi2;
        res.base.mean_error = last_F_metrics.mean_err;
        res.base.kl         = last_F_metrics.kl;
        res.base.grake_norm = last_F_metrics.grake_norm;
    }

    // Water-filling detects partial infeasibility (some categories can't reach targets
    // within bounds). With Dykstra this was silently masked by bound violations.
    // For stalled iterations: return STALL (status=5) + best weights rather than
    // hard-erroring — caller can use the best achievable calibration.

    // Note: is_infeasible (water_fill_cat capacity exhaustion) is intentionally
    // NOT promoted to RK_ERR_INFEAS here. Upper-bound capacity capping is normal
    // bounded calibration — the solver returns the best achievable constrained
    // optimum with binding capacity constraints, not an infeasibility error.
    // True structural INFEAS (lower bounds exceed target) is caught at validation.

    // Post-loop: normalize sum to n (water-filling already enforces bounds)
    {
        double s_post = 0.0;
        for (int c = 0; c < ct.M_cell; c++) s_post += X[c];
        if (s_post > 0.0) {
            const double sc_post = static_cast<double>(st.n) / s_post;
            for (int c = 0; c < ct.M_cell; c++) X[c] *= sc_post;
        }
    }

    // G8b: expand best.best_weights (cell-level X snapshot) to obs-level.
    res.base.convergence_solver_objective = best.best_objective;
    res.base.convergence_minimized_metric = static_cast<int>(st.convergence_cfg.metric);
    res.base.best_error = best.best_metric;
    res.base.best_iter  = best.best_iter;
    if (best.has_best()) {
        // Normalize cell snapshot, then expand to obs.
        std::vector<double> w_snap = best.best_weights;  // copy — normalize in place
        double s = 0.0;
        for (int c = 0; c < ct.M_cell; c++) s += w_snap[c];
        if (s > 0.0) {
            const double sc = static_cast<double>(st.n) / s;
            for (int c = 0; c < ct.M_cell; c++) w_snap[c] *= sc;
        }
        res.base.best_weights.resize(st.n);
        for (int i = 0; i < st.n; i++) {
            int c = ct.cell_of[i];
            double mult = (X_init[c] > 0.0) ? w_snap[c] / X_init[c] : 1.0;
            res.base.best_weights[i] = st.weights[i] * mult;   // no per-obs clamp
        }
        // CR-D11 (j7x8.11): Σw=n + bounds_mode dispatch via the canonical
        // finalize_weights contract (cell: count-only; unit: water-fill).
        int b_nbv = 0, b_nbc = 0;
        lbw::finalize_weights_buf(res.base.best_weights.data(), st.n, st, ct, b_nbv, b_nbc);
    } else {
        res.base.best_weights.assign(st.n, 0.0);
    }

    // CR-D11 (j7x8.11): obs expansion with NO per-obs clamp; finalize_weights
    // enforces Σw=n and the bounds_mode contract (cell mode leaves per-obs
    // violations as a diagnostic; clamping here distorts marginals — measured).
    lbw::expand_obs(ct, X, X_init, st.n, st.weights);
    lbw::finalize_weights(st, ct, res.n_bounds_violated, res.n_bounds_clamped);

    return res;
}
} // namespace lbw
