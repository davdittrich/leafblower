#pragma once
// calib_dispatch.hpp — shared CalibMetric / CalibRule helpers (WU-B).
// Eliminates triplicated switch blocks across ieppa, raking, and lbfgsb_solver.
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
#include <limits>
#include "cell_table.hpp"
#include "calib_validate.hpp"
#include <vector>
#include <algorithm>
#include <cassert>
#include <cstring>

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
// `prev` is read and updated by this function; initialise to +infinity before
// the first iteration.  Returns false on the first call when prev==+inf.
//
// Note: lbfgsb_solver does NOT use this helper because it is a batch solver
// with no per-iteration baseline; it retains its own THRESHOLD-only logic.
inline bool apply_rule(
    CalibRule rule,
    double    curr,
    double&   prev,   // in-out: updated to curr on return when finite
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

    prev = curr;   // always update — sliding-window baseline for next call (intentional side effect)
    return converged;
}

struct CellMetrics {
    double errRp      = 0.0;
    double mean_err   = 0.0;
    double kl         = 0.0;
    double chi2       = 0.0;
    double grake_norm = 0.0;
    double l1         = 0.0;
};

/// Convenience overload: select the convergence metric value from a CellMetrics struct.
/// marginal_kl defaults to 0.0 (not tracked in CellMetrics).
inline double select_metric(CalibMetric metric, const CellMetrics& m) noexcept {
    return select_metric(metric, m.errRp, m.mean_err, m.kl, m.chi2,
                         m.grake_norm, m.l1);
    // marginal_kl omitted — defaults to 0.0 in the 7-arg overload
}

// Mathematical objective for NON-KL solvers only.
// KL-minimizing solvers (ieppa, sinkhorn, raking) use compute_weight_kl inline.
inline double select_solver_objective(int alg_id, const lbw::CellMetrics& m) {
    switch (alg_id) {
    case RK_ALG_GREG:      return m.chi2;
    case RK_ALG_CHEBYSHEV: return m.errRp;
    default:               return m.errRp;
    }
}

// Returns true when the convergence criterion fires. Updates prev_metric via apply_rule.
// tol_abs_fallback: value of CalibState::tol_abs, used when neither pct_tol nor absolute_tol set.
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
    if (have_pct && have_abs)
        return (cfg.stop_when == CalibStopWhen::ALL)
               ? (c_pct && c_abs) : (c_pct || c_abs);
    if (have_pct) return c_pct;
    if (have_abs) return c_abs;
    return (m.errRp < tol_abs_fallback);
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
    constexpr double kChi2Floor = 1.0;
    CellMetrics m;
    double mean_sum = 0.0;
    for (int k = 0; k < st.K; k++) {
        const int nj = st.cat_counts[k];
        assert(static_cast<int>(bucket.size()) >= nj);
        std::fill(bucket.begin(), bucket.begin() + nj, 0.0);
        for (int c = 0; c < ct.M_cell; c++) {
            int g = ct.g_per_cell[k][c];
            if (g >= 0 && g < nj) bucket[g] += X[c];
        }
        double max_k = 0.0, kl_k = 0.0;
        for (int j = 0; j < nj; j++) {
            double S_p = bucket[j] / W, T = st.targets[k][j];
            double err = std::fabs(S_p - T);
            if (err > max_k)        max_k = err;
            if (err > m.errRp)      m.errRp = err;
            if (T > 0.0)            kl_k += T * std::log((T + kMetricEps) / (S_p + kMetricEps));
            double obs = bucket[j], pop = T * W;
            m.chi2 += (obs - pop) * (obs - pop) / (pop + kChi2Floor);
            double nm = std::fabs(obs - pop) / (1.0 + std::fabs(pop));
            if (nm > m.grake_norm)  m.grake_norm = nm;
        }
        mean_sum += max_k;
        if (kl_k > m.kl) m.kl = kl_k;
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

/// Returns st.max_weight if finite, else 1e300.
inline double resolve_hi(const CalibState& st) noexcept {
    return std::isfinite(st.max_weight) ? st.max_weight : 1e300;
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
    if (build_cell_table(st.n, st.K, st.group_ids, st.cat_counts, st.weights, ct) != 0) {
        res.base.status = RK_ERR_BADARG;
        std::strncpy(res.message, "build_cell_table failed", sizeof(res.message) - 1);
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

} // namespace lbw
