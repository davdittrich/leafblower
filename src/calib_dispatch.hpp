#pragma once
// calib_dispatch.hpp — shared CalibMetric / CalibRule helpers (WU-B).
// Eliminates triplicated switch blocks across oris, raking, and lbfgsb_solver.
//
// Semantics:
//   select_metric — map CalibMetric enum → the pre-computed scalar value.
//   apply_rule    — apply a CalibRule stopping condition, updating prev in-place.
//
// apply_rule behaviour (all rules require tol > 0 and finite curr):
//   THRESHOLD  : curr < tol
//   IMPROVEMENT: |curr - prev| / prev < tol
//                Special: curr <= 1e-15 → trivially converged (metric at machine zero).
//                Skip (false) when prev is non-finite or prev <= 1e-15 (first check).
//   PLATEAU    : !(curr < prev * (1 - tol))  — i.e. metric did not improve by ≥tol fraction.
//                Skip (false) when prev is non-finite (first check).
//
// The caller is responsible for passing the correct local variable names for each
// metric; lbfgsb passes pct_change for L1_WEIGHT because it is a batch solver.

#include "leafblower.h"
#include "types.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include "cell_table.hpp"
#include "calib_validate.hpp"
#include <vector>
#include <algorithm>
#include <cassert>
#include <cstring>
#include "lbw_math.hpp"
#include "sinkhorn.hpp"
#include "greg.hpp"
#include "greenkhorn.hpp"
#include "logit_calib.hpp"
#include "oris.hpp"
#include "chebyshev.hpp"
#include "raking.hpp"
#include "newton_calib.hpp"

namespace lbw {

// Select the active metric value from the pre-computed metric set.
// The caller provides all seven metric values in canonical order; only the
// one matching `metric` is used. marginal_kl defaults to 0.0 for callers
// that pre-date MARGINAL_KL=6 and do not compute it.
inline double select_metric(
    CalibMetric metric,
    double max_err,
    double mean_err,
    double kl,
    double chi2,
    double grake_norm,
    double l1_weight,
    double marginal_kl = 0.0) noexcept
{
    switch (metric) {
        case CalibMetric::MAX_ERR:     return max_err;
        case CalibMetric::MEAN_ERR:    return mean_err;
        case CalibMetric::KL:          return kl;
        case CalibMetric::CHI2:        return chi2;
        case CalibMetric::GRAKE_NORM:  return grake_norm;
        case CalibMetric::L1_WEIGHT:   return l1_weight;
        case CalibMetric::MARGINAL_KL: return marginal_kl;
    }
    return max_err;  // unreachable; exhaustive enum silences -Wreturn-type
}

// Apply the stopping rule.
// `prev` is read here and updated to curr ONLY when the rule did NOT fire
// (converged == false).  When convergence fires, prev is left at its
// pre-call value so the caller can distinguish "at-convergence baseline"
// from the halting curr.  Callers that continue iterating despite a fired
// rule (e.g. stop_when=ALL waiting on a second criterion) must update prev
// themselves (see check_convergence below).
//
// Note: lbfgsb_solver does NOT use this helper because it is a batch solver
// with no per-iteration baseline; it retains its own THRESHOLD-only logic.
inline bool apply_rule(
    CalibRule rule,
    double    curr,
    double&   prev,   // in-out: updated to curr on return ONLY when not converged
    double    tol) noexcept
{
    if (!std::isfinite(curr) || tol <= 0.0) return false;

    bool converged = false;
    switch (rule) {
        case CalibRule::THRESHOLD:
            converged = (curr < tol);
            break;

        case CalibRule::IMPROVEMENT: {
            // When curr is at machine-zero the metric cannot improve further;
            // treat relative change as 0 (trivially converged).
            if (curr <= 1e-15) {
                converged = true;
            } else if (std::isfinite(prev) && prev > 1e-15) {
                const double rel = std::fabs(curr - prev) / prev;
                converged = (rel < tol);
            }
            // else: first check (prev==inf) or prev near-zero → skip
            break;
        }

        case CalibRule::PLATEAU:
            // Converge when curr did NOT drop by at least tol fraction vs prev.
            // Skip on first check (prev==inf).
            // B17: when prev=0 and curr>0, !(curr < 0*(1-tol)) = !(curr < 0) fires
            // spuriously — a non-zero metric after a zero baseline is a rebound,
            // not a plateau. Only evaluate when prev>0 or when both prev=0 and
            // curr=0 (metric genuinely flatlined at zero → legitimate plateau).
            if (std::isfinite(prev) && (prev > 0.0 || curr <= 0.0)) {
                converged = !(curr < prev * (1.0 - tol));
            }
            break;
    }

    // Route A: early-return on short-circuit halt — prev retains at-convergence
    // baseline value, NOT the post-halt curr.  This prevents prev from leaking
    // the halting curr into any re-entry path (homotopy restart, SRAA, diagnostics).
    if (converged) return true;

    prev = curr;   // sliding-window baseline for next call (only when not converged)
    return false;
}

struct CellMetrics {
    double errRp       = 0.0;
    double mean_err    = 0.0;
    double kl          = 0.0;   // MAX over margins of per-margin KL (back-compat)
    double chi2        = 0.0;
    double grake_norm  = 0.0;
    double l1          = 0.0;
    double marginal_kl = 0.0;   // SUM over margins of per-margin KL (flat-path marginal_kl)
};

/// Unified best-iterate tracker shared across all solvers.
/// Each solver maintains one instance and calls update() every iteration.
/// Replaces ad-hoc best_metric_seen / best_iter_val / best_objective_seen
/// variables that were duplicated (with inconsistent init semantics) across
/// oris, raking, sinkhorn, greenkhorn, and newton_calib.
struct BestIterTracker {
    double best_metric    = std::numeric_limits<double>::infinity();
    double best_objective = std::numeric_limits<double>::infinity(); // weight-KL at best iter
    int    best_iter      = -1;
    std::vector<double> best_weights;

    /// Returns true if this is a new best (metric strictly improved).
    bool update(double metric, double objective, int iter,
                const std::vector<double>& weights) {
        if (metric < best_metric) {
            best_metric    = metric;
            best_objective = objective;
            best_iter      = iter;
            best_weights   = weights;
            return true;
        }
        return false;
    }

    /// Reset for SRAA path restarts (oris calls this mid-run).
    void reset() {
        best_metric    = std::numeric_limits<double>::infinity();
        best_objective = std::numeric_limits<double>::infinity();
        best_iter      = -1;
        best_weights.clear();
    }

    bool has_best() const { return best_iter >= 0; }
};

/// Convenience overload: select the convergence metric value from a CellMetrics struct.
inline double select_metric(CalibMetric metric, const CellMetrics& m) noexcept {
    return select_metric(metric, m.errRp, m.mean_err, m.kl, m.chi2,
                         m.grake_norm, m.l1, m.marginal_kl);
    // marginal_kl now tracked (Σ_k per-margin KL) — previously omitted, which
    // defaulted to 0.0 and instant-triggered convergence for metric="marginal_kl".
}

// Mathematical objective for NON-KL solvers only.
// KL-minimizing solvers (oris, sinkhorn, raking) use compute_weight_kl inline.
inline double select_solver_objective(int alg_id, const lbw::CellMetrics& m) {
    switch (alg_id) {
    case RK_ALG_GREG:      return m.chi2;
    case RK_ALG_CHEBYSHEV: return m.errRp;
    default:
        // CR-D12 (j7x8.12): CRAN forbids R packages writing to stderr directly.
        // Route through REprintf under the R build; keep stderr for the LBW_NO_R
        // (Python/standalone) build where REprintf is unavailable — matching the
        // established convention in types.hpp.
#ifndef LBW_NO_R
        REprintf("internal: unknown alg_id %d in select_solver_objective\n", alg_id);
#else
        std::fprintf(stderr, "internal: unknown alg_id %d in select_solver_objective\n", alg_id);
#endif
        return std::numeric_limits<double>::quiet_NaN();
    }
}

// Returns true when the convergence criterion fires. Updates prev_metric via apply_rule.
// tol_abs_fallback: value of CalibState::tol_abs, used when neither pct_tol nor absolute_tol set.
//
// prev_metric contract (post-call):
//   - Halting (returns true):  prev_metric retains pre-call value (at-convergence baseline).
//   - Continuing (returns false): prev_metric = curr (sliding-window updated for next check).
//
// For stop_when=ALL: apply_rule fires when the pct criterion alone is met, but the driver
// continues until abs also fires.  In that case apply_rule leaves prev un-updated (Route A);
// we restore prev_metric = curr here so the next check computes improvement/plateau correctly.
inline bool check_convergence(
    const CalibConvergenceCfg& cfg,
    const CellMetrics& m,
    double& prev_metric,
    double tol_abs_fallback) noexcept
{
    const double curr = select_metric(cfg.metric, m);
    const bool have_pct = (cfg.pct_tol > 0.0), have_abs = (cfg.absolute_tol > 0.0);
    const bool c_abs = have_abs && (curr < cfg.absolute_tol);
    const bool c_pct = have_pct && apply_rule(cfg.rule, curr, prev_metric, cfg.pct_tol);
    if (have_pct && have_abs) {
        const bool halting = (cfg.stop_when == CalibStopWhen::ALL)
                             ? (c_pct && c_abs) : (c_pct || c_abs);
        // apply_rule left prev un-updated when c_pct fired (Route A).  If the driver is
        // continuing (not halting), advance the sliding window so the next check has a
        // correct baseline.
        if (!halting && c_pct) prev_metric = curr;
        return halting;
    }
    if (have_pct) return c_pct;
    if (have_abs) return c_abs;
    return (curr < tol_abs_fallback);
}

/// Aggregate cell masses X[] into margin k bucket[0..nj).
/// bucket must be pre-allocated to at least nj elements.
/// Zeroes bucket[0..nj) before accumulating.
inline void aggregate_to_margin(
    const CellTable& ct,
    const std::vector<double>& X,
    int k, int nj,
    double* bucket) noexcept
{
    std::fill(bucket, bucket + nj, 0.0);
    for (int c = 0; c < ct.M_cell; c++) {
        int g = ct.g_per_cell[k][c];
        if (g >= 0 && g < nj) bucket[g] += X[c];
    }
}

// Compute 5 calibration metrics over all K margins from cell-level weight vector X.
// W = sum(X) — passed in to avoid recomputation.
// bucket: pre-allocated scratch of size >= max(st.cat_counts[k]).
inline CellMetrics compute_cell_metrics(
    const CalibState& st, const CellTable& ct,
    const std::vector<double>& X, double W,
    std::vector<double>& bucket) noexcept
{
    constexpr double kMetricEps = 1e-10;
    CellMetrics m;
    double mean_sum = 0.0;
    for (int k = 0; k < st.K; k++) {
        const int nj = st.cat_counts[k];
        assert(static_cast<int>(bucket.size()) >= nj);
        aggregate_to_margin(ct, X, k, nj, bucket.data());
        double max_k = 0.0, kl_k = 0.0;
        for (int j = 0; j < nj; j++) {
            double S_p = bucket[j] / W, T = st.targets[k][j];
            double err = std::fabs(S_p - T);
            if (err > max_k)        max_k = err;
            if (err > m.errRp)      m.errRp = err;
            if (T > 0.0)            kl_k += T * std::log((T + kMetricEps) / (S_p + kMetricEps));
            double obs = bucket[j], pop = T * W;
            // Pearson chi2: (obs - pop)^2 / pop. Skip bins with pop <= kMetricEps:
            // T_j = 0 means the bin is not part of the test (no expected count), and the
            // chi2 contribution is undefined under standard Pearson semantics. Prior
            // kChi2Floor=1.0 Laplace smoothing distorted contributions for rare bins
            // (pop ~ O(1)) by halving them; spec-correct treatment is exclusion.
            if (pop > kMetricEps)   m.chi2 += (obs - pop) * (obs - pop) / pop;
            double nm = std::fabs(obs - pop) / (1.0 + std::fabs(pop));
            if (nm > m.grake_norm)  m.grake_norm = nm;
        }
        mean_sum += max_k;
        if (kl_k > m.kl) m.kl = kl_k;   // MAX over margins (back-compat)
        m.marginal_kl += kl_k;          // SUM over margins (flat-path marginal_kl)
    }
    m.mean_err = (st.K > 0) ? mean_sum / static_cast<double>(st.K) : 0.0;
    return m;
}

// Post-solve obs expansion (NO per-obs clamp): w[i] ← w[i] × X[cell]/X_init[cell].
// Yields Σw = Σ_c X[c] (since Σ_{i∈c} w_init[i] = X_init[c]); finalize_weights
// then enforces Σw=n and the bounds_mode contract (cell: count-only diagnostic;
// unit: per-cell water-fill, which restores Σw=n when redistribution capacity
// suffices). Per the canonical cell contract (see finalize_weights doc),
// clamping per-obs HERE distorts marginals and breaks Σw=n whenever it binds
// (CR-D11 / j7x8.11: measured 13pp margin drift, Σw 3392 vs n 4000).
// Guard: X_init[c] > 1e-10 matches greg's kEps and chebyshev's hardcoded threshold.
// Functionally identical to > 0.0 for all realistic inputs (X_init[c] = sum of initial
// weights for obs in cell c; always >= min_weight when the cell is non-empty).
inline void expand_obs(
    const CellTable& ct,
    const std::vector<double>& X,
    const std::vector<double>& X_init,
    int n, double* weights) noexcept
{
    for (int i = 0; i < n; i++) {
        int c = ct.cell_of[i];
        double mult = (X_init[c] > 1e-10) ? X[c] / X_init[c] : 1.0;
        weights[i] *= mult;
    }
}

/// Sentinel for "unbounded" max_weight upper bound. Chosen as 1e300 (not
/// std::numeric_limits<double>::infinity()) to preserve bit-exact numerical
/// behavior under clamp/multiply: a finite sentinel saturates products and
/// keeps downstream weight × cell-count arithmetic finite. parity_bench
/// rtol=1e-6 vs prior literal-1e300 build is the gate.
inline constexpr double kUnboundedSentinel = 1e300;

/// Minimum total weight for which the Σw=n renormalization (norm = n/total_w) is
/// applied. Below this, total_w is treated as degenerate-zero: weights are left
/// unchanged rather than scaled by an enormous norm (a subnormal total_w ~1e-310
/// would overflow norm to +inf). 1e-100 keeps n/total_w well within double range
/// for any realistic n while rejecting only physically-meaningless near-zero sums.
#ifndef LBW_KMIN_SAFE_TOTAL_WEIGHT_DEFINED
#define LBW_KMIN_SAFE_TOTAL_WEIGHT_DEFINED
inline constexpr double kMinSafeTotalWeight = 1e-100;
#endif

/// Returns st.max_weight if finite, else kUnboundedSentinel.
inline double resolve_hi(const CalibState& st) noexcept {
    return std::isfinite(st.max_weight) ? st.max_weight : kUnboundedSentinel;
}

/// Fills L[c] = lo*n_per_cell[c], U[c] = hi*n_per_cell[c] for c in [0, M_cell).
inline void compute_cell_bounds(
    const CellTable& ct,
    double lo, double hi,
    std::vector<double>& L,
    std::vector<double>& U) noexcept
{
    L.resize(ct.M_cell); U.resize(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) {
        L[c] = lo * static_cast<double>(ct.n_per_cell[c]);
        U[c] = hi * static_cast<double>(ct.n_per_cell[c]);
    }
}

/// Exit-time weight finalization shared by solvers (canonical: oris_finalize).
/// (1) Enforces Σw=n via one pre-bounds scale (sanctioned: applied BEFORE bounds
///     post-processing so unit-mode water-fill sees final-scale weights — NOT a
///     forbidden post-water-fill renormalize).
/// (2) bounds_mode==CELL: counts per-obs violations (diagnostic only, no clamp —
///     cell contract is aggregate X[c] <= U_cell, and clamping per-obs distorts
///     marginals). bounds_mode==UNIT: per-cell water-fill enforcing
///     [min_weight, max_weight] per obs while preserving Σw=n.
/// Degenerate total_w (0 / non-finite): weights left unchanged.
/// Writes the per-obs violation / clamp counts to the out-params (solver result
/// structs differ in whether they surface these — callers pass their own fields
/// or discard).
///
/// _buf variant operates on an explicit weight buffer (w, n) using st only for
/// bounds_mode / min_weight / max_weight / K — used for both the final iterate
/// (st.weights) and the best-iterate snapshot (res.base.best_weights).
inline void finalize_weights_buf(double* w, int n, const CalibState& st,
                                 const CellTable& ct,
                                 int& n_bounds_violated, int& n_bounds_clamped) {
    double total_w = 0.0;
    for (int i = 0; i < n; i++) total_w += w[i];
    if (std::isfinite(total_w) && total_w > kMinSafeTotalWeight) {
        const double norm = static_cast<double>(n) / total_w;
        for (int i = 0; i < n; i++) w[i] *= norm;
    }

    if (st.bounds_mode == RK_BOUNDS_CELL) {
        int violations = 0;
        for (int i = 0; i < n; i++)
            if (w[i] > st.max_weight || w[i] < st.min_weight) ++violations;
        n_bounds_violated = violations;
        n_bounds_clamped  = 0;
        return;
    }

    // Unit mode: per-cell water-fill (mirror of oris_finalize unit branch).
    std::vector<std::vector<int>> cells_of_obs(ct.M_cell);
    for (int i = 0; i < n; i++) cells_of_obs[ct.cell_of[i]].push_back(i);
    const int kWaterFillMaxIter = std::max(50, st.K * 10);
    int total_clamped = 0;
    for (int c = 0; c < ct.M_cell; c++) {
        const auto& idxs = cells_of_obs[c];
        if (idxs.empty()) continue;
        for (int it = 0; it < kWaterFillMaxIter; it++) {
            double excess = 0.0, free_sum = 0.0;
            int n_free = 0; bool any_violation = false;
            for (int i : idxs) {
                if (w[i] > st.max_weight) {
                    excess += w[i] - st.max_weight;
                    w[i] = st.max_weight; any_violation = true; ++total_clamped;
                } else if (w[i] < st.min_weight) {
                    excess -= st.min_weight - w[i];
                    w[i] = st.min_weight; any_violation = true; ++total_clamped;
                } else if (w[i] == st.max_weight || w[i] == st.min_weight) {
                    // Pinned from prior iter (assigned, not computed — FP-equality safe).
                    // Excluded from free_sum so factor distributes excess in full.
                } else {
                    free_sum += w[i]; ++n_free;
                }
            }
            if (!any_violation) break;
            if (n_free == 0 || free_sum <= 0.0) break;  // no room to redistribute
            const double factor = 1.0 + excess / free_sum;
            for (int i : idxs)
                if (w[i] > st.min_weight && w[i] < st.max_weight)
                    w[i] *= factor;
        }
        // CR-D11a (j7x8.15): on kWaterFillMaxIter exhaustion the final
        // redistribution multiply can leave a free weight past a bound with no
        // trailing clamp — the loop's clamp runs at the TOP of the next iter,
        // which never comes. This final clamp-only pass enforces the unit-mode
        // bound contract for every finite weight. No-op when the loop converged
        // (all weights already in bounds); in the rare non-converged case bounds
        // win over an O(excess) Σ drift for this one cell (bounds are the hard
        // promise of bounds_mode="unit"). NaN iterates are deliberately NOT
        // clamped here — a NaN weight is a solver failure that must surface, not
        // be masked into a valid-looking bound.
        for (int i : idxs) {
            if (w[i] > st.max_weight)      { w[i] = st.max_weight; ++total_clamped; }
            else if (w[i] < st.min_weight) { w[i] = st.min_weight; ++total_clamped; }
        }
    }
    n_bounds_violated = 0;
    n_bounds_clamped  = total_clamped;
}

/// Convenience overload operating on st.weights (the final iterate).
inline void finalize_weights(CalibState& st, const CellTable& ct,
                             int& n_bounds_violated, int& n_bounds_clamped) {
    finalize_weights_buf(st.weights, st.n, st, ct, n_bounds_violated, n_bounds_clamped);
}

/// kxna.20 (CR-C7c): re-gate a converged unit-mode result on the POST-finalize
/// returned margin error. In bounds_mode="unit" the per-obs water-fill
/// (finalize_weights_buf) can leave a fully-pinned cell (feasible mass >
/// M_cell*max_weight) off-target AFTER the solver already set RK_OK on the
/// pre-finalize cell iterate, so the RETURNED weights miss margins while status
/// reads 0. Demote RK_OK -> RK_ERR_STALL (best-effort valid weights at a
/// bounds-constrained optimum; NOT hard INFEAS — matches jy0m.1/.2 + harvest.R
/// STALL vocabulary) when the returned error exceeds the solver's tolerance bar.
///
/// `post_finalize_max_error` MUST be the PROPORTION-space margin error of the
/// RETURNED weights (e.g. logit abs/n_d, or compute_cell_metrics(...).errRp on
/// the finalized obs weights aggregated to cells). The threshold is st.tol_abs,
/// also proportion-space — this mirrors each solver's own convergence bar (e.g.
/// logit's max_abs_resid <= tol_abs*n  ==  max_error <= tol_abs). A feasible
/// problem that reaches RK_OK drives the returned error to ~machine-zero (the
/// improvement rule converges to precision), so this never false-demotes a
/// genuinely-calibrated unit run; only water-fill drift on a bounds-infeasible
/// cell produces RK_OK + a large returned error.
///
/// Passed explicitly (not read from base.max_error) because greenkhorn locks
/// base.max_error to a configured-metric scale (CXX.2) — a units mismatch. A NaN
/// error does NOT demote (the NaN surfaces via the weights + reported metric;
/// masking it into a status is wrong, per finalize_weights_buf NaN contract).
/// Cell-mode has no water-fill mass shift, so it is untouched. Never clobbers an
/// already-non-OK status (e.g. an INFEAS promotion that ran first wins).
inline void regate_unit_status(CalibResult& base, const CalibState& st,
                               double post_finalize_max_error) noexcept {
    if (base.status == RK_OK && st.bounds_mode == RK_BOUNDS_UNIT &&
        std::isfinite(post_finalize_max_error) &&
        post_finalize_max_error > st.tol_abs) {
        base.status           = RK_ERR_STALL;
        base.stall_kind       = 2;   // at a bounds-constrained optimum (harvest: "stall_kl")
        base.convergence_iter = -1;  // no valid convergence — matches other non-OK exits
    }
}

/// Builds prefix-sum cat_offset[k] and returns n_cats_total.
/// cat_offset[k] = sum of cat_counts[0..k-1].
inline int build_cat_offset(int K, const int* cat_counts,
                             std::vector<int>& cat_offset) noexcept
{
    cat_offset.resize(K);
    int nct = 0;
    for (int k = 0; k < K; k++) { cat_offset[k] = nct; nct += cat_counts[k]; }
    return nct;
}

/// Returns max(cat_counts[0..K-1]), or 0 when K==0.
inline int max_cats_count(int K, const int* cat_counts) noexcept {
    if (K == 0) return 0;
    return *std::max_element(cat_counts, cat_counts + K);
}

/// Populate standard convergence metadata on a result struct.
/// ResT must have: .base.status (int), .base.convergence_metric (int),
///                 .base.convergence_rule (int), .base.convergence_tol (double),
///                 .base.convergence_iter (int).
template <typename ResT>
inline void mark_converged(ResT& res, const CalibConvergenceCfg& cfg, int iter,
                           double tol_abs_fallback) noexcept {
    res.base.status             = RK_OK;
    res.base.convergence_metric = static_cast<int>(cfg.metric);
    res.base.convergence_rule   = static_cast<int>(cfg.rule);
    // CR-C10b (kxna.22): when NEITHER cfg tol is set, check_convergence halts on the
    // st.tol_abs fallback (`curr < tol_abs_fallback`), so report that same value here
    // instead of a misleading pct_tol=0.0. Other cases unchanged (abs-only -> abs;
    // any pct set -> pct), preserving existing reporting.
    res.base.convergence_tol    = (cfg.absolute_tol > 0.0 && cfg.pct_tol == 0.0) ? cfg.absolute_tol
                                : (cfg.pct_tol > 0.0)                            ? cfg.pct_tol
                                : tol_abs_fallback;
    res.base.convergence_iter   = iter;
}

// Full cell-table setup: build_cell_table + X_init + bounds + cat_offset + validate.
// ResT must have: res.base.status (int) and res.message (char[]).
// Returns RK_OK on success, otherwise sets res.base.status + res.message and returns error code.
template <typename ResT>
inline int solver_setup_ct(
    CalibState&              st,
    CellTable&               ct,
    std::vector<double>&     X_init,
    double&                  hi_eff,
    std::vector<double>&     L_cell,
    std::vector<double>&     U_cell,
    std::vector<int>&        cat_offset,
    int&                     n_cats_total,
    ResT&                    res) noexcept
{
    static_assert(
        std::is_class<typename std::remove_reference<decltype(
            std::declval<ResT>().base)>::type>::value,
        "solver_setup_ct<ResT>: ResT must have a .base member (CalibResultBase). "
        "Valid types: RakingResult, ORISResult, ChebyshevResult, etc."
    );
    if (build_cell_table(st.n, st.K, st.group_ids, st.cat_counts, st.weights, ct) != 0) {
        res.base.status = RK_ERR_BADARG;
        std::strncpy(res.message, "build_cell_table failed", sizeof(res.message) - 1);
        res.message[sizeof(res.message) - 1] = '\0';
        return RK_ERR_BADARG;
    }
    X_init.assign(ct.M_cell, 0.0);
    for (int i = 0; i < st.n; i++) X_init[ct.cell_of[i]] += st.weights[i];
    hi_eff = lbw::resolve_hi(st);
    lbw::compute_cell_bounds(ct, st.min_weight, hi_eff, L_cell, U_cell);
    n_cats_total = lbw::build_cat_offset(st.K, st.cat_counts, cat_offset);
    rk_result_t tmp = {};
    if (calib_validate_preentry(ct, st, &tmp, X_init.data(), n_cats_total) != RK_OK) {
        res.base.status = tmp.status;
        std::strncpy(res.message, tmp.message, sizeof(res.message) - 1);
        res.message[sizeof(res.message) - 1] = '\0';
        return tmp.status;
    }
    return RK_OK;
}

// Partial cell-table setup for raking: build_cell_table + X_init + bounds +
// capacity/negativity pre-entry validation. No cat_offset. ResT only needs
// res.base.status (no message field required — the validate message is dropped
// because RakingResult carries none, matching ORISResult).
// Returns RK_OK on success, otherwise sets res.base.status and returns error code.
template <typename ResT>
inline int solver_setup_ct_base(
    CalibState&              st,
    CellTable&               ct,
    std::vector<double>&     X_init,
    double&                  hi_eff,
    std::vector<double>&     L_cell,
    std::vector<double>&     U_cell,
    ResT&                    res) noexcept
{
    if (build_cell_table(st.n, st.K, st.group_ids, st.cat_counts, st.weights, ct) != 0) {
        res.base.status = RK_ERR_BADARG;
        return RK_ERR_BADARG;
    }
    X_init.assign(ct.M_cell, 0.0);
    for (int i = 0; i < st.n; i++) X_init[ct.cell_of[i]] += st.weights[i];
    hi_eff = lbw::resolve_hi(st);
    lbw::compute_cell_bounds(ct, st.min_weight, hi_eff, L_cell, U_cell);
    // CR-D4 (j7x8.4): run the same capacity/negativity feasibility check the 6
    // full-setup solvers run, so raking rejects negative-weight / structurally-
    // infeasible input up front (RK_ERR_BADARG/INFEAS) instead of letting a negative
    // mass reach log() and burning the full NOCONV budget. n_cats_total=0 skips the
    // cat-count-limit gate — raking IS the high-cardinality fallback, so it does not
    // apply. tmp.message is discarded (RakingResult has no message field).
    rk_result_t tmp = {};
    int vrc = calib_validate_preentry(ct, st, &tmp, X_init.data(), /*n_cats_total=*/0);
    if (vrc != RK_OK) {
        res.base.status = tmp.status;
        return vrc;
    }
    return RK_OK;
}

/// Weight-space KL: Σ_c X[c]*log(X[c]/X_init[c]) / n.
/// ratio_buf and weight_buf are caller-owned scratch of size >= M_cell.
/// Returns 0.0 when result is non-finite.
inline double compute_weight_kl(
    const std::vector<double>& X,
    const std::vector<double>& X_init,
    int M_cell, int n,
    double* ratio_buf,
    double* weight_buf) noexcept
{
    const double inv_n = 1.0 / static_cast<double>(n);
    int valid = 0;
    for (int c = 0; c < M_cell; c++) {
        if (X_init[c] > 0.0 && X[c] > 0.0) {
            ratio_buf[valid]  = X[c] / X_init[c];
            weight_buf[valid] = X[c];
            valid++;
        }
    }
    if (valid == 0) return 0.0;
    lbw::bulk_log(ratio_buf, ratio_buf, valid);
    double wkl = 0.0;
    for (int i = 0; i < valid; i++) wkl += weight_buf[i] * ratio_buf[i] * inv_n;
    return std::isfinite(wkl) ? wkl : 0.0;
}

// ── Shared solver dispatch (SC1, leafblower-rywn) ──────────────────────────
//
// Neutral, internal, ABI-unconstrained result type shared by both FFI
// callers' solver dispatch. A strict superset of rk_result_t's field set
// (see plan 02-01 Task 1's inventory, bd comment on leafblower-rywn):
// c_api.cpp narrows it into rk_result_t (frozen 536B ABI); r_bridge.cpp
// marshals it in full into the SEXP result list (unchanged 49-element
// surface). Every member default-initializes to the same value the
// corresponding res_* local in r_bridge.cpp::C_rk_calibrate used before
// this migration.
struct DispatchResult {
    int             status                       = RK_ERR_NOCONV;
    int             iterations                   = 0;
    double          max_error                    = 1.0;
    rk_algorithm_t  alg_used                     = RK_ALG_ORIS;
    double          mean_error                   = 0.0;
    double          kl                           = 0.0;
    double          chi2                         = 0.0;
    double          l1_weight_change             = 0.0;
    double          grake_norm                   = 0.0;
    int             convergence_metric           = 0;
    int             convergence_rule             = 1;
    double          convergence_tol              = 0.001;
    int             convergence_iter             = -1;
    double          convergence_solver_objective = 0.0;
    int             convergence_minimized_metric = 0;
    double          best_error                   = std::numeric_limits<double>::infinity();
    int             best_iter                    = 0;
    double          metric_first_check           = std::numeric_limits<double>::infinity();
    double          metric_prev_check            = std::numeric_limits<double>::infinity();
    int             prev_check_iter              = -1;
    int             stall_kind                   = 0;
    int             n_bounds_violated            = 0;
    int             n_bounds_clamped             = 0;
    char            solver_message[256]          = {};
    /* SUPERSET-ONLY (no rk_result_t field, Task 1 inventory): newton_kl diagnostics. */
    int             n_projected_dims             = 0;
    double          lm_mu_final                  = 0.0;
    /* SUPERSET-ONLY: oris/oris_soft/raking SRAA scheduler-demotion flag. */
    bool            sraa_demoted                 = false;
    /* ORIS/ORIS_SOFT-only diagnostics (plan 02-05). rk_result_t (leafblower.h)
     * already carries every one of these — pack_oris_result_c wrote them
     * directly from ORISResult before this migration. pack_dispatch_result_c
     * (c_api.cpp, shared by every OTHER migrated solver) deliberately does
     * NOT copy them, so a non-oris solver correctly leaves rk_result_init's
     * memset/1.0 defaults untouched; c_api.cpp's oris/oris_soft branches call
     * a dedicated pack_dispatch_oris_extras_c after pack_dispatch_result_c.
     * Defaults mirror r_bridge.cpp's res_* locals / ORISResult's own defaults. */
    int             n_xcur_writes_per_iter_last  = 0;
    double          min_alpha_seen               = 1.0;
    double          final_alpha                  = 1.0;
    int             homotopy_levels_used         = 0;
    double          homotopy_final_factor        = 1.0;
    int             greedy_sweeps_taken          = 0;
    double          eta_final                    = 0.0;
    double          sor_min_omega                = 1.0;
    int             sor_n_damped                 = 0;
    double          sor_omega_mean               = 1.0;
    int             sor_any_latched              = 0;
    int             sor_n_pinned_fb              = 0;
    int             sor_n_warmup_fb              = 0;
    int             sor_n_conv_fb                = 0;
    int             sor_n_resid_grew             = 0;
    int             sor_n_monotone_cd            = 0;
    /* SUPERSET-ONLY: SRAA-m (Anderson Acceleration) diagnostic. No
     * rk_result_t counterpart — aa_accepted_count is absent from
     * leafblower.h and pack_oris_result_c never surfaced it either, a
     * pre-existing gap flagged in plan 02-01 (bd comment on leafblower-rywn),
     * not fixed by this migration. */
    int             aa_accepted_count            = 0;
    /* ALM diagnostics (oris_soft only; zero elsewhere). */
    double          alm_capacity_mu_final        = 0.0;
    int             alm_n_growth_events          = 0;
    double          alm_max_dual_norm            = 0.0;
    double          alm_sum_drift                = 0.0;
    /* SUPERSET-ONLY: obs-level snapshot, distinct from rk_calibrate's in-place
     * `weights` output param. Heap-backed — callers holding a function-scope
     * DispatchResult across an Rf_error() call site must swap-release this
     * member alongside their other RAII locals (RESEARCH.md Pitfall 4). */
    std::vector<double> best_weights;
};

// Shared solver-selection + result-extraction table. Both c_api.cpp's
// rk_calibrate() and r_bridge.cpp's C_rk_calibrate() dispatch through this
// single function instead of maintaining two independent
// {enum/string -> solver -> result} chains (SC1, leafblower-rywn). Covers
// RK_ALG_SINKHORN (plan 02-01, the tracer slice), RK_ALG_GREG,
// RK_ALG_GREENKHORN, RK_ALG_LOGIT (plan 02-03), RK_ALG_CHEBYSHEV,
// RK_ALG_RAKING (plan 02-04), RK_ALG_ORIS, RK_ALG_ORIS_SOFT (plan 02-05), and
// RK_ALG_NEWTON_KL (plan 02-06) — every named method is now migrated; only
// RK_ALG_AUTO's routing logic (plan 07) still bypasses this table. Every
// other enum value leaves `out` untouched and returns without acting, so
// not-yet-migrated callers keep running their existing branch unchanged.
// Solver-by-solver migration (D-01, 02-CONTEXT.md) adds one case arm per plan.
inline void dispatch_solver(rk_algorithm_t alg, CalibState& st, DispatchResult& out) {
    switch (alg) {
        case RK_ALG_SINKHORN: {
            auto res = lbw::sinkhorn_solve(st);
            out.status                       = res.base.status;
            out.iterations                   = res.base.iterations;
            out.max_error                    = res.base.max_error;
            out.mean_error                   = res.base.mean_error;
            out.kl                           = res.base.kl;
            out.chi2                         = res.base.chi2;
            out.l1_weight_change             = res.base.l1_weight_change;
            out.grake_norm                   = res.base.grake_norm;
            out.convergence_metric           = res.base.convergence_metric;
            out.convergence_rule             = res.base.convergence_rule;
            out.convergence_tol              = res.base.convergence_tol;
            out.convergence_iter             = res.base.convergence_iter;
            out.convergence_solver_objective = res.base.convergence_solver_objective;
            out.convergence_minimized_metric = res.base.convergence_minimized_metric;
            out.best_error                   = res.base.best_error;
            out.best_iter                    = res.base.best_iter;
            out.metric_first_check           = res.base.metric_first_check;
            out.metric_prev_check            = res.base.metric_prev_check;
            out.prev_check_iter              = res.base.prev_check_iter;
            out.stall_kind                   = res.base.stall_kind;
            out.n_bounds_violated            = res.n_bounds_violated;
            out.n_bounds_clamped             = res.n_bounds_clamped;
            std::snprintf(out.solver_message, sizeof(out.solver_message), "%s", res.message);
            out.alg_used                     = RK_ALG_SINKHORN;
            out.best_weights                 = std::move(res.base.best_weights);
            break;
        }
        case RK_ALG_GREG: {
            auto res = lbw::greg_solve(st);
            out.status                       = res.base.status;
            out.iterations                   = res.base.iterations;
            out.max_error                    = res.base.max_error;
            out.mean_error                   = res.base.mean_error;
            out.kl                           = res.base.kl;
            out.chi2                         = res.base.chi2;
            out.l1_weight_change             = res.base.l1_weight_change;
            out.grake_norm                   = res.base.grake_norm;
            out.convergence_metric           = res.base.convergence_metric;
            out.convergence_rule             = res.base.convergence_rule;
            out.convergence_tol              = res.base.convergence_tol;
            out.convergence_iter             = res.base.convergence_iter;
            out.convergence_solver_objective = res.base.convergence_solver_objective;
            out.convergence_minimized_metric = res.base.convergence_minimized_metric;
            out.best_error                   = res.base.best_error;
            out.best_iter                    = res.base.best_iter;
            out.metric_first_check           = res.base.metric_first_check;
            out.metric_prev_check            = res.base.metric_prev_check;
            out.prev_check_iter              = res.base.prev_check_iter;
            out.stall_kind                   = res.base.stall_kind;
            out.n_bounds_violated            = res.n_bounds_violated;
            out.n_bounds_clamped             = res.n_bounds_clamped;
            std::snprintf(out.solver_message, sizeof(out.solver_message), "%s", res.message);
            out.alg_used                     = RK_ALG_GREG;
            out.best_weights                 = std::move(res.base.best_weights);
            break;
        }
        case RK_ALG_GREENKHORN: {
            auto res = lbw::greenkhorn_solve(st);
            out.status                       = res.base.status;
            out.iterations                   = res.base.iterations;
            out.max_error                    = res.base.max_error;
            out.mean_error                   = res.base.mean_error;
            out.kl                           = res.base.kl;
            out.chi2                         = res.base.chi2;
            out.l1_weight_change             = res.base.l1_weight_change;
            out.grake_norm                   = res.base.grake_norm;
            out.convergence_metric           = res.base.convergence_metric;
            out.convergence_rule             = res.base.convergence_rule;
            out.convergence_tol              = res.base.convergence_tol;
            out.convergence_iter             = res.base.convergence_iter;
            out.convergence_solver_objective = res.base.convergence_solver_objective;
            out.convergence_minimized_metric = res.base.convergence_minimized_metric;
            out.best_error                   = res.base.best_error;
            out.best_iter                    = res.base.best_iter;
            out.metric_first_check           = res.base.metric_first_check;
            out.metric_prev_check            = res.base.metric_prev_check;
            out.prev_check_iter              = res.base.prev_check_iter;
            out.stall_kind                   = res.base.stall_kind;
            out.n_bounds_violated            = res.n_bounds_violated;
            out.n_bounds_clamped             = res.n_bounds_clamped;
            std::snprintf(out.solver_message, sizeof(out.solver_message), "%s", res.message);
            out.alg_used                     = RK_ALG_GREENKHORN;
            out.best_weights                 = std::move(res.base.best_weights);
            break;
        }
        case RK_ALG_LOGIT: {
            auto res = lbw::logit_calibrate(st);
            out.status                       = res.base.status;
            out.iterations                   = res.base.iterations;
            out.max_error                    = res.base.max_error;
            out.mean_error                   = res.base.mean_error;
            out.kl                           = res.base.kl;
            out.chi2                         = res.base.chi2;
            out.l1_weight_change             = res.base.l1_weight_change;
            out.grake_norm                   = res.base.grake_norm;
            out.convergence_metric           = res.base.convergence_metric;
            out.convergence_rule             = res.base.convergence_rule;
            out.convergence_tol              = res.base.convergence_tol;
            out.convergence_iter             = res.base.convergence_iter;
            out.convergence_solver_objective = res.base.convergence_solver_objective;
            out.convergence_minimized_metric = res.base.convergence_minimized_metric;
            out.best_error                   = res.base.best_error;
            out.best_iter                    = res.base.best_iter;
            out.metric_first_check           = res.base.metric_first_check;
            out.metric_prev_check            = res.base.metric_prev_check;
            out.prev_check_iter              = res.base.prev_check_iter;
            out.stall_kind                   = res.base.stall_kind;
            out.n_bounds_violated            = res.n_bounds_violated;
            out.n_bounds_clamped             = res.n_bounds_clamped;
            std::snprintf(out.solver_message, sizeof(out.solver_message), "%s", res.message);
            out.alg_used                     = RK_ALG_LOGIT;
            out.best_weights                 = std::move(res.base.best_weights);
            break;
        }
        case RK_ALG_CHEBYSHEV: {
            // Oris warm-start, moved here from both bridges (previously
            // hand-duplicated in r_bridge.cpp and c_api.cpp — xc1s.15/kxna.23)
            // so it runs exactly once per solve. Cold code (once-per-solve,
            // not per-iteration), so living in this header costs nothing.
            std::vector<double> w_warm_obs;
            {
                // SAFETY: weights_copy protects st.weights from oris_solve
                // mutation. st_warm must NOT escape this block (dangling
                // pointer into weights_copy once it is destroyed).
                std::vector<double> weights_copy(st.weights, st.weights + st.n);
                lbw::CalibState st_warm = st;
                st_warm.weights = weights_copy.data();
                st_warm.inner_max_iter = std::max(5, std::min(100, st.inner_max_iter / 10));
                auto oris_res = lbw::oris_solve(st_warm);
                if (!oris_res.base.best_weights.empty() &&
                    static_cast<int>(oris_res.base.best_weights.size()) == st.n &&
                    std::isfinite(oris_res.base.max_error)) {
                    w_warm_obs = std::move(oris_res.base.best_weights);
                }
            }
            auto res = lbw::chebyshev_ipm(st, w_warm_obs);
            out.status                       = res.base.status;
            out.iterations                   = res.base.iterations;
            out.max_error                    = res.base.max_error;
            out.mean_error                   = res.base.mean_error;
            out.kl                           = res.base.kl;
            out.chi2                         = res.base.chi2;
            out.l1_weight_change             = res.base.l1_weight_change;
            out.grake_norm                   = res.base.grake_norm;
            out.convergence_metric           = res.base.convergence_metric;
            out.convergence_rule             = res.base.convergence_rule;
            out.convergence_tol              = res.base.convergence_tol;
            out.convergence_iter             = res.base.convergence_iter;
            out.convergence_solver_objective = res.base.convergence_solver_objective;
            out.convergence_minimized_metric = res.base.convergence_minimized_metric;
            out.best_error                   = res.base.best_error;
            out.best_iter                    = res.base.best_iter;
            out.metric_first_check           = res.base.metric_first_check;
            out.metric_prev_check            = res.base.metric_prev_check;
            out.prev_check_iter              = res.base.prev_check_iter;
            out.stall_kind                   = res.base.stall_kind;
            out.n_bounds_violated            = res.n_bounds_violated;
            out.n_bounds_clamped             = res.n_bounds_clamped;
            std::snprintf(out.solver_message, sizeof(out.solver_message), "%s", res.message);
            out.alg_used                     = RK_ALG_CHEBYSHEV;
            // kxna.15/xl44: the violation guard (inner_max_iter<1, or a
            // solver_setup_ct failure) leaves best_weights empty; fall back to
            // a zero-filled sentinel of length st.n — same fallback the R
            // bridge applied per-callsite before SC1.
            if (!res.base.best_weights.empty())
                out.best_weights = std::move(res.base.best_weights);
            else
                out.best_weights.assign(st.n, 0.0);
            break;
        }
        case RK_ALG_RAKING: {
            auto res = lbw::raking_solve(st);
            out.status                       = res.base.status;
            out.iterations                   = res.base.iterations;
            out.max_error                    = res.base.max_error;
            out.mean_error                   = res.base.mean_error;
            out.kl                           = res.base.kl;
            out.chi2                         = res.base.chi2;
            out.l1_weight_change             = res.base.l1_weight_change;
            out.grake_norm                   = res.base.grake_norm;
            out.convergence_metric           = res.base.convergence_metric;
            out.convergence_rule             = res.base.convergence_rule;
            out.convergence_tol              = res.base.convergence_tol;
            out.convergence_iter             = res.base.convergence_iter;
            out.convergence_solver_objective = res.base.convergence_solver_objective;
            out.convergence_minimized_metric = res.base.convergence_minimized_metric;
            out.best_error                   = res.base.best_error;
            out.best_iter                    = res.base.best_iter;
            out.metric_first_check           = res.base.metric_first_check;
            out.metric_prev_check            = res.base.metric_prev_check;
            out.prev_check_iter              = res.base.prev_check_iter;
            out.stall_kind                   = res.base.stall_kind;
            out.n_bounds_violated            = res.n_bounds_violated;
            out.n_bounds_clamped             = res.n_bounds_clamped;
            out.solver_message[0]            = '\0';  // RakingResult carries no message field
            out.alg_used                     = RK_ALG_RAKING;
            // SUPERSET-ONLY: r_bridge's res_sraa_demoted local is R-only
            // (rk_result_t has no counterpart — see pack_dispatch_result_c).
            out.sraa_demoted                 = res.sraa_demoted;
            out.best_weights                 = std::move(res.base.best_weights);
            break;
        }
        case RK_ALG_ORIS: {
            // Plan 02-05 Task 1: absorbs what r_bridge.cpp's pack_oris_result
            // lambda + c_api.cpp's pack_oris_result_c did today (diffed
            // field-by-field; no divergence found — see 02-05-SUMMARY.md).
            // st.oris_auto_selected stays on the CALLER side (this arm takes
            // st as already-configured state; routing policy is not this
            // arm's job — keeps plan 07's AUTO consolidation a pure move).
            auto res = lbw::oris_solve(st);
            out.status                       = res.base.status;
            out.iterations                   = res.base.iterations;
            out.max_error                    = res.base.max_error;
            out.mean_error                   = res.base.mean_error;
            out.kl                           = res.base.kl;
            out.chi2                         = res.base.chi2;
            out.l1_weight_change             = res.base.l1_weight_change;
            out.grake_norm                   = res.base.grake_norm;
            out.convergence_metric           = res.base.convergence_metric;
            out.convergence_rule             = res.base.convergence_rule;
            out.convergence_tol              = res.base.convergence_tol;
            out.convergence_iter             = res.base.convergence_iter;
            out.convergence_solver_objective = res.base.convergence_solver_objective;
            out.convergence_minimized_metric = res.base.convergence_minimized_metric;
            out.best_error                   = res.base.best_error;
            out.best_iter                    = res.base.best_iter;
            out.metric_first_check           = res.base.metric_first_check;
            out.metric_prev_check            = res.base.metric_prev_check;
            out.prev_check_iter              = res.base.prev_check_iter;
            out.stall_kind                   = res.base.stall_kind;
            out.n_bounds_violated            = res.n_bounds_violated;
            out.n_bounds_clamped             = res.n_bounds_clamped;
            out.solver_message[0]            = '\0';  // ORISResult carries no message field
            out.alg_used                     = RK_ALG_ORIS;
            out.n_xcur_writes_per_iter_last  = res.n_xcur_writes_per_iter_last;
            out.min_alpha_seen               = res.min_alpha_seen;
            out.final_alpha                  = res.final_alpha;
            out.homotopy_levels_used         = res.homotopy_levels_used;
            out.homotopy_final_factor        = res.homotopy_final_factor;
            out.greedy_sweeps_taken          = res.greedy_sweeps_taken;
            out.eta_final                    = res.eta_final;
            out.sor_min_omega                = res.sor_min_omega;
            out.sor_n_damped                 = res.sor_n_damped;
            out.sor_omega_mean               = res.sor_omega_mean;
            out.sor_any_latched              = res.sor_any_latched;
            out.sor_n_pinned_fb              = res.sor_n_pinned_fb;
            out.sor_n_warmup_fb              = res.sor_n_warmup_fb;
            out.sor_n_conv_fb                = res.sor_n_conv_fb;
            out.sor_n_resid_grew             = res.sor_n_resid_grew;
            out.sor_n_monotone_cd            = res.sor_n_monotone_cd;
            out.aa_accepted_count            = res.aa_accepted_count;
            out.sraa_demoted                 = res.sraa_demoted;
            // ORIS never sets alm_* (ALM only runs under oris_soft's
            // use_admm_capacity); leave DispatchResult's 0.0/0 defaults.
            out.best_weights                 = std::move(res.base.best_weights);
            break;
        }
        case RK_ALG_ORIS_SOFT: {
            // Plan 02-05 Task 2: identical field set to RK_ALG_ORIS above
            // (oris_soft is oris_solve() run under ALM capacity penalty),
            // plus the four ALM diagnostics that only this path populates.
            // st.use_admm_capacity / st.oris_auto_selected / capacity_penalty
            // auto-resolution all stay on the CALLER side — the estimate_M_cell
            // memoization question (single vs. duplicate evaluation) is a
            // caller-side concern this arm neither creates nor solves; see
            // 02-05-SUMMARY.md for the two-bridge comparison plan 07 needs.
            auto res = lbw::oris_solve(st);
            out.status                       = res.base.status;
            out.iterations                   = res.base.iterations;
            out.max_error                    = res.base.max_error;
            out.mean_error                   = res.base.mean_error;
            out.kl                           = res.base.kl;
            out.chi2                         = res.base.chi2;
            out.l1_weight_change             = res.base.l1_weight_change;
            out.grake_norm                   = res.base.grake_norm;
            out.convergence_metric           = res.base.convergence_metric;
            out.convergence_rule             = res.base.convergence_rule;
            out.convergence_tol              = res.base.convergence_tol;
            out.convergence_iter             = res.base.convergence_iter;
            out.convergence_solver_objective = res.base.convergence_solver_objective;
            out.convergence_minimized_metric = res.base.convergence_minimized_metric;
            out.best_error                   = res.base.best_error;
            out.best_iter                    = res.base.best_iter;
            out.metric_first_check           = res.base.metric_first_check;
            out.metric_prev_check            = res.base.metric_prev_check;
            out.prev_check_iter              = res.base.prev_check_iter;
            out.stall_kind                   = res.base.stall_kind;
            out.n_bounds_violated            = res.n_bounds_violated;
            out.n_bounds_clamped             = res.n_bounds_clamped;
            out.solver_message[0]            = '\0';  // ORISResult carries no message field
            out.alg_used                     = RK_ALG_ORIS_SOFT;
            out.n_xcur_writes_per_iter_last  = res.n_xcur_writes_per_iter_last;
            out.min_alpha_seen               = res.min_alpha_seen;
            out.final_alpha                  = res.final_alpha;
            out.homotopy_levels_used         = res.homotopy_levels_used;
            out.homotopy_final_factor        = res.homotopy_final_factor;
            out.greedy_sweeps_taken          = res.greedy_sweeps_taken;
            out.eta_final                    = res.eta_final;
            out.sor_min_omega                = res.sor_min_omega;
            out.sor_n_damped                 = res.sor_n_damped;
            out.sor_omega_mean               = res.sor_omega_mean;
            out.sor_any_latched              = res.sor_any_latched;
            out.sor_n_pinned_fb              = res.sor_n_pinned_fb;
            out.sor_n_warmup_fb              = res.sor_n_warmup_fb;
            out.sor_n_conv_fb                = res.sor_n_conv_fb;
            out.sor_n_resid_grew             = res.sor_n_resid_grew;
            out.sor_n_monotone_cd            = res.sor_n_monotone_cd;
            out.aa_accepted_count            = res.aa_accepted_count;
            out.sraa_demoted                 = res.sraa_demoted;
            out.alm_capacity_mu_final        = res.alm_capacity_mu_final;
            out.alm_n_growth_events          = res.alm_n_growth_events;
            out.alm_max_dual_norm            = res.alm_max_dual_norm;
            out.alm_sum_drift                = res.alm_sum_drift;
            out.best_weights                 = std::move(res.base.best_weights);
            break;
        }
        case RK_ALG_NEWTON_KL: {
            auto res = lbw::newton_calibrate(st);
            out.status                       = res.base.status;
            out.iterations                   = res.base.iterations;
            out.max_error                    = res.base.max_error;
            out.mean_error                   = res.base.mean_error;
            out.kl                           = res.base.kl;
            out.chi2                         = res.base.chi2;
            out.l1_weight_change             = res.base.l1_weight_change;
            out.grake_norm                   = res.base.grake_norm;
            out.convergence_metric           = res.base.convergence_metric;
            out.convergence_rule             = res.base.convergence_rule;
            out.convergence_tol              = res.base.convergence_tol;
            out.convergence_iter             = res.base.convergence_iter;
            out.convergence_solver_objective = res.base.convergence_solver_objective;
            out.convergence_minimized_metric = res.base.convergence_minimized_metric;
            out.best_error                   = res.base.best_error;
            out.best_iter                    = res.base.best_iter;
            out.metric_first_check           = res.base.metric_first_check;
            out.metric_prev_check            = res.base.metric_prev_check;
            out.prev_check_iter              = res.base.prev_check_iter;
            out.stall_kind                   = res.base.stall_kind;
            out.n_bounds_violated            = res.n_bounds_violated;
            out.n_bounds_clamped             = res.n_bounds_clamped;
            std::snprintf(out.solver_message, sizeof(out.solver_message), "%s", res.message);
            out.alg_used                     = RK_ALG_NEWTON_KL;
            // SUPERSET-ONLY (no rk_result_t field): newton_kl's own diagnostics.
            out.n_projected_dims             = res.n_projected_dims;
            out.lm_mu_final                  = res.lm_mu_final;
            // Violation guard (newton_calib.cpp: frac_violated > 0.05) leaves
            // res.base.best_weights empty/default; fall back to a zero-filled
            // sentinel of length st.n — the same fallback both bridges applied
            // per-callsite before this migration.
            if (!res.base.best_weights.empty())
                out.best_weights = std::move(res.base.best_weights);
            else
                out.best_weights.assign(st.n, 0.0);
            // newton never sets ORIS-only diagnostics (n_xcur_writes_per_iter_last,
            // min_alpha_seen, ..., alm_sum_drift); DispatchResult's default-
            // constructed values for those members ARE the documented non-ORIS
            // defaults (1.0/0/1.0/... — see struct comment above), so this arm
            // leaves them untouched by design.
            break;
        }
        default:
            break;  // not yet migrated (D-01); caller's existing branch handles it
    }
}

// Lazily-cached lbw::estimate_M_cell (eb79.15). ANY caller that needs
// M_cell/n -- AUTO routing (route_auto below) or oris_soft's
// capacity_penalty auto-resolution (both bridges) -- goes through this, so
// the O(n*K) estimate runs at most once per rk_calibrate()/C_rk_calibrate()
// call regardless of how many of those callers a single call touches.
// m_cell_est_cache: in/out, -1 = uncomputed; whichever caller runs first
// computes and stores it, every later caller in the same call reuses it.
inline int resolve_m_cell_est(int n, int K,
                               const int32_t* const* group_ids,
                               const int* cat_counts,
                               int& m_cell_est_cache) {
    if (m_cell_est_cache < 0) {
        m_cell_est_cache = lbw::estimate_M_cell(n, K, group_ids, cat_counts);
    }
    return m_cell_est_cache;
}

// Result of the AUTO routing DECISION (never solves; see route_auto below).
struct AutoRouteResult {
    rk_algorithm_t algorithm       = RK_ALG_ORIS;
    bool           auto_selected   = true;   // always true when route_auto ran
    bool           force_accelerate = false; // true iff the severe-skew ORIS+SRAA branch fired
};

// AUTO routing (Epic-H WH-g): the single routing DECISION, shared by both
// c_api.cpp's rk_calibrate() and r_bridge.cpp's C_rk_calibrate() (SC1, plan
// 07 -- consolidates what was previously two independently-maintained
// inline copies, r_bridge.cpp:663-788 and c_api.cpp:292-397). This function
// only DECIDES which solver AUTO picks; it never solves. Keeping the
// decision separate from the solve is what lets the auto-fallback path be
// expressed as a second lbw::dispatch_solver() call by the caller, instead
// of nested routing logic.
//
// Thresholds preserved EXACTLY (this migration moves the logic, it does not
// retune it): the exact-integer form of the 0.9 compression comparison
// (M_cell_est*10 >= n*9, PAR.1's `>=` not `>`), K>=5, target_skew > 5.
//
// m_cell_est_cache: in/out, -1 = uncomputed; see resolve_m_cell_est above --
// a caller that already resolved M_cell/n for another purpose (e.g.
// oris_soft's capacity_penalty auto-resolution) passes its cached value in
// and this function reuses it instead of recomputing; a caller with nothing
// cached passes a fresh local initialized to -1.
//
// Takes raw (n, K, group_ids, cat_counts, targets) rather than a CalibState
// so it can be called from c_api.cpp's algorithm-resolution switch, which
// runs BEFORE that caller's CalibState is built (r_bridge.cpp's CalibState
// already exists by the time it calls this and simply passes st.n/st.K/
// st.group_ids/st.cat_counts/st.targets through).
inline AutoRouteResult route_auto(int n, int K,
                                   const int32_t* const* group_ids,
                                   const int* cat_counts,
                                   const double* const* targets,
                                   int& m_cell_est_cache) {
    const int M_cell_est = resolve_m_cell_est(n, K, group_ids, cat_counts, m_cell_est_cache);
    // Exact integer comparison: M_cell_est / n >= 0.9  <=>  M_cell_est*10 >= n*9
    const bool zero_compression =
        (static_cast<int64_t>(M_cell_est) * 10 >= static_cast<int64_t>(n) * 9);

    AutoRouteResult out;
    if (zero_compression && K >= 5) {
        double max_target = 0.0, min_target = 1.0;
        for (int k = 0; k < K; ++k) {
            for (int j = 0; j < cat_counts[k]; ++j) {
                const double t = targets[k][j];
                if (t > max_target) max_target = t;
                if (t > 1e-12 && t < min_target) min_target = t;
            }
        }
        const double target_skew = max_target / std::max(min_target, 1e-12);
        if (target_skew > 5.0) {
            // Epic-H WH-g: kk1204 K=20 evidence -- Newton-KL stalls at
            // gap~=6.24e-2 on severe-skew dual landscape; ORIS+SRAA converges.
            out.algorithm        = RK_ALG_ORIS;
            out.force_accelerate = true;
        } else {
            out.algorithm = RK_ALG_NEWTON_KL;  // moderate skew
        }
    } else if (zero_compression) {
        out.algorithm = RK_ALG_RAKING;  // K<5, zero-compression
    } else {
        out.algorithm = RK_ALG_ORIS;    // compressed regime, any K
    }
    return out;
}

} // namespace lbw
