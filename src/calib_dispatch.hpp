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
        std::fprintf(stderr, "internal: unknown alg_id %d in select_solver_objective\n", alg_id);
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

// Post-solve obs expansion: w[i] ← clamp(w[i] × X[cell]/X_init[cell], lo, hi)
// Guard: X_init[c] > 1e-10 matches greg's kEps and chebyshev's hardcoded threshold.
// Functionally identical to > 0.0 for all realistic inputs (X_init[c] = sum of initial
// weights for obs in cell c; always >= min_weight when the cell is non-empty).
inline void apply_obs_expansion(
    const CellTable& ct,
    const std::vector<double>& X,
    const std::vector<double>& X_init,
    int n, double lo, double hi,
    double* weights) noexcept
{
    for (int i = 0; i < n; i++) {
        int c = ct.cell_of[i];
        double mult = (X_init[c] > 1e-10) ? X[c] / X_init[c] : 1.0;
        weights[i] = std::clamp(weights[i] * mult, lo, hi);
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
    }
    n_bounds_violated = 0;
    n_bounds_clamped  = total_clamped;
}

/// Convenience overload operating on st.weights (the final iterate).
inline void finalize_weights(CalibState& st, const CellTable& ct,
                             int& n_bounds_violated, int& n_bounds_clamped) {
    finalize_weights_buf(st.weights, st.n, st, ct, n_bounds_violated, n_bounds_clamped);
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
inline void mark_converged(ResT& res, const CalibConvergenceCfg& cfg, int iter) noexcept {
    res.base.status             = RK_OK;
    res.base.convergence_metric = static_cast<int>(cfg.metric);
    res.base.convergence_rule   = static_cast<int>(cfg.rule);
    res.base.convergence_tol    = (cfg.absolute_tol > 0.0 && cfg.pct_tol == 0.0)
                                  ? cfg.absolute_tol : cfg.pct_tol;
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

// Partial cell-table setup for raking: build_cell_table + X_init + bounds only.
// No cat_offset, no validate. ResT only needs res.base.status (no message field required).
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

} // namespace lbw
