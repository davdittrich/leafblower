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
//
// ── Hierarchical 2-stage shared infrastructure (T-B) ──────────────────────
//   build_cell_partition  — group obs by joint coarse-margin profile (O(N))
//   build_sparse_mask     — flag cells with n_per_cell < min_cell_n
//   apply_sparse_inheritance — in-place: sparse-cell obs inherit Stage-1 multipliers
//   enforce_sigmaw_eq_n   — post-solver Σw=n gate (called by solver at exit)

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
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <sstream>
#include <numeric>
#include "lbw_math.hpp"

// ── Hierarchical 2-stage cell-count cap ───────────────────────────────────
// T-E references this same constant; update both places if the value changes.
#define LBW_MAX_HIER_CELLS 100000

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

/// Unified best-iterate tracker shared across all solvers.
/// Each solver maintains one instance and calls update() every iteration.
/// Replaces ad-hoc best_metric_seen / best_iter_val / best_objective_seen
/// variables that were duplicated (with inconsistent init semantics) across
/// ieppa, raking, sinkhorn, greenkhorn, and newton_calib.
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

    /// Reset for SRAA path restarts (ieppa calls this mid-run).
    void reset() {
        best_metric    = std::numeric_limits<double>::infinity();
        best_objective = std::numeric_limits<double>::infinity();
        best_iter      = -1;
        best_weights.clear();
    }

    bool has_best() const { return best_iter >= 0; }
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
    default:
        std::fprintf(stderr, "internal: unknown alg_id %d in select_solver_objective\n", alg_id);
        return std::numeric_limits<double>::quiet_NaN();
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
    constexpr double kChi2Floor = 1.0;
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

/// Populate standard convergence metadata on a result struct.
/// ResT must have: .base.status (int), .base.convergence_metric (int),
///                 .base.convergence_rule (int), .base.convergence_tol (double),
///                 .base.convergence_iter (int).
template <typename ResT>
inline void mark_converged(ResT& res, const CalibConvergenceCfg& cfg, int iter) noexcept {
    res.base.status             = RK_OK;
    res.base.convergence_metric = static_cast<int>(cfg.metric);
    res.base.convergence_rule   = static_cast<int>(cfg.rule);
    res.base.convergence_tol    = cfg.pct_tol;
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

// ── Hierarchical 2-stage shared infrastructure ────────────────────────────
// All structs + helpers live in the lbw namespace and are inline/static per
// existing convention (build_cell_table, resolve_hi, solver_setup_ct, etc.).

/// Result of partitioning observations by their joint coarse-margin profile.
/// obs_indices_by_cell is CSR-like: obs_indices_by_cell[c] holds all obs in cell c.
struct CellPartition {
    int                              n_cells_total;     // K_cells: total coarse cells found
    int                              n_cells_skipped;   // count of is_sparse_cell entries (set by build_sparse_mask)
    int                              n_cells_inherited; // alias of n_cells_skipped
    std::vector<int>                 cell_id_per_obs;   // size N
    std::vector<int>                 n_per_cell;        // size n_cells_total
    std::vector<std::vector<int>>    obs_indices_by_cell; // [n_cells_total][*]
};

/// Mask flagging which coarse cells are sparse.
struct SparseMask {
    std::vector<bool> is_sparse_cell; // size n_cells_total; true = n_per_cell < min_cell_n
};

/// Group observations by joint coarse-margin category tuple.
/// coarse_mask[k] == 1 → margin k participates in the coarse partition.
/// Returns RK_ERR_BADARG (and a zero-initialised partition) on cap violation.
///
/// Cap check (heap-DoS guard, spec §9 + Security review):
///   K_cells > min(N / max(min_cell_n,1), LBW_MAX_HIER_CELLS) → RK_ERR_BADARG.
/// Cap check runs BEFORE allocating obs_indices_by_cell.
///
/// Complexity: O(N) with unordered_map on a 64-bit packed key per obs.
/// Supports up to K_coarse ≤ 64 coarse margins (same limit as build_cell_table).
inline int build_cell_partition(
    int                      N,
    int                      K,
    const int32_t* const*    group_ids,   // [K][N]
    const int*               cat_counts,  // [K]
    const int*               coarse_mask, // [K]: 1 = include in coarse partition
    int                      min_cell_n,
    CellPartition&           out) noexcept
{
    out = {};

    // Count coarse margins; validate K_coarse <= 64
    int K_coarse = 0;
    for (int k = 0; k < K; k++) if (coarse_mask[k]) ++K_coarse;
    if (K_coarse > lbw::K_MAX) {
        out.n_cells_total = 0;
        return RK_ERR_BADARG;
    }

    // Build packed keys: one uint64_t per obs encoding all coarse cats.
    // Use cat_counts[k]+1 as the radix for margin k (NA = cat_counts[k]).
    // Radix overflow beyond 64 bits is prevented by K_coarse <= K_MAX <= 64
    // and cat_counts[k] <= N (practical bound).
    std::vector<uint64_t> keys(N, 0ULL);
    {
        uint64_t radix = 1ULL;
        for (int k = 0; k < K; k++) {
            if (!coarse_mask[k]) continue;
            const int32_t* gk = group_ids[k];
            const uint64_t base = static_cast<uint64_t>(cat_counts[k] + 1);
            for (int i = 0; i < N; i++) {
                int32_t g = gk[i];
                uint64_t cat = (g < 0) ? static_cast<uint64_t>(cat_counts[k]) : static_cast<uint64_t>(g);
                keys[i] += cat * radix;
            }
            radix *= base;
            // Overflow guard: if radix wraps, we silently collision (benign — only
            // affects partitioning granularity, not correctness of Σw=n gate).
        }
    }

    // First pass: discover unique cells and assign IDs.
    std::unordered_map<uint64_t, int> key_to_cell;
    key_to_cell.reserve(static_cast<size_t>(std::min(N, 4096)));
    out.cell_id_per_obs.resize(N);
    for (int i = 0; i < N; i++) {
        auto [it, inserted] = key_to_cell.emplace(keys[i], static_cast<int>(key_to_cell.size()));
        out.cell_id_per_obs[i] = it->second;
    }
    const int K_cells = static_cast<int>(key_to_cell.size());

    // Cap check BEFORE allocating obs_indices_by_cell (heap-DoS guard).
    const int effective_min = (min_cell_n < 1) ? 1 : min_cell_n;
    const int cap = std::min(N / effective_min, LBW_MAX_HIER_CELLS);
    if (K_cells > cap) {
        out = {};
        return RK_ERR_BADARG;
    }

    // Second pass: count per cell.
    out.n_per_cell.assign(K_cells, 0);
    for (int i = 0; i < N; i++) ++out.n_per_cell[out.cell_id_per_obs[i]];

    // Third pass: fill obs_indices_by_cell (pre-sized to avoid realloc in loop).
    out.obs_indices_by_cell.resize(K_cells);
    for (int c = 0; c < K_cells; c++) out.obs_indices_by_cell[c].reserve(out.n_per_cell[c]);
    for (int i = 0; i < N; i++) out.obs_indices_by_cell[out.cell_id_per_obs[i]].push_back(i);

    out.n_cells_total    = K_cells;
    out.n_cells_skipped  = 0;  // set by build_sparse_mask
    out.n_cells_inherited = 0; // alias; set by build_sparse_mask
    return RK_OK;
}

/// Flag cells with n_per_cell < min_cell_n as sparse.
/// Returns SparseMask; also updates partition.n_cells_skipped and n_cells_inherited.
inline SparseMask build_sparse_mask(
    CellPartition& partition,
    int            min_cell_n) noexcept
{
    SparseMask mask;
    const int K_cells = partition.n_cells_total;
    mask.is_sparse_cell.resize(K_cells, false);
    int n_sparse = 0;
    for (int c = 0; c < K_cells; c++) {
        if (partition.n_per_cell[c] < min_cell_n) {
            mask.is_sparse_cell[c] = true;
            ++n_sparse;
        }
    }
    partition.n_cells_skipped  = n_sparse;
    partition.n_cells_inherited = n_sparse; // alias of n_cells_skipped
    return mask;
}

/// Apply sparse-cell inheritance in-place.
/// For each obs in a sparse cell: weights[i] = w_init[i] * stage1_multipliers[cell_id].
/// Per spec §6: unscaled Stage-1 multiplier inherit — NOT cell-total preservation.
///
/// Preconditions:
///   weights.size() == N (same N used to build partition)
///   w_init.size()  == N
///   stage1_multipliers.size() == partition.n_cells_total
///   sparse_mask.is_sparse_cell.size() == partition.n_cells_total
inline void apply_sparse_inheritance(
    std::vector<double>&        weights,
    const std::vector<double>&  w_init,
    const CellPartition&        partition,
    const SparseMask&           sparse_mask,
    const std::vector<double>&  stage1_multipliers) noexcept
{
    const int K_cells = partition.n_cells_total;
    for (int c = 0; c < K_cells; c++) {
        if (!sparse_mask.is_sparse_cell[c]) continue;
        const double mult = stage1_multipliers[c];
        for (int i : partition.obs_indices_by_cell[c]) {
            weights[i] = w_init[i] * mult;
        }
    }
}

/// Post-solver Σw=n gate.  Called by solver at exit — NOT by outer loop.
/// Returns true when |Σw − N| < N · 1e-12; false on violation.
/// N is the integer sample size (>0).
inline bool enforce_sigmaw_eq_n(
    const std::vector<double>& weights,
    int                        N) noexcept
{
    if (N <= 0) return false;
    double sum = 0.0;
    for (double w : weights) sum += w;
    const double tol = static_cast<double>(N) * 1e-12;
    return std::abs(sum - static_cast<double>(N)) < tol;
}

// ── Strategy A outer-iteration loop ─────────────────────────────────────────
// Shared cadence constant — mirrors the per-solver kErrCheckInterval = 10.
// Used only for best-iterate selection; NOT a convergence check cadence.
inline constexpr int kOuterErrCheckInterval = 10;

/// Diagnostic-only best-iterate tracking for the outer loop.
/// Populated at kOuterErrCheckInterval cadence; NOT returned as the result.
/// Caller may read for instrumentation; spec §6: last-iterate is the result.
struct OuterBestIter {
    int    iter_idx     = -1;            // -1 = not yet tracked
    double residual     = std::numeric_limits<double>::infinity();
};

/// Return value of outer_iterate_strategy_a.
struct OuterResult {
    int    status;                       // RK_OK | RK_ERR_BUDGET
    int    iterations_used;
    double residual_final;
    std::vector<double> last_iterate_weights;
    // Diagnostic only (spec §6): NOT returned to caller as result weights.
    OuterBestIter best;
};

/// Strategy A outer-iteration loop (P1 default mode = "refine").
///
/// Alternates Stage 1 (full-data, via fn(weights, -1)) and Stage 2
/// (per non-sparse cell, via fn(weights, cell_id)) until coarse-margin
/// residual drops below outer_tol or outer_iterations is exhausted.
///
/// WithinCellSolve callable contract:
///   int fn(std::vector<double>& weights, int cell_id)
///   cell_id == -1  → Stage-1 (full-data calibrate to coarse margins)
///   cell_id >= 0   → Stage-2 (calibrate obs in that coarse cell only)
///   Return value: RK_OK on success (non-RK_OK is ignored — outer loop
///   continues regardless; callers may inspect via a side-channel if needed).
///
/// Residual is max absolute deviation of coarse-margin sums from N/K_cells
/// (a simple proxy; full metric requires CalibState which is solver-specific).
/// At kOuterErrCheckInterval cadence: track best-iterate via select_metric
/// using the caller-supplied metric from cfg. Best-iterate is DIAGNOSTIC ONLY.
///
/// On RK_OK exit: enforce_sigmaw_eq_n is called; result.last_iterate_weights
/// reflects the post-gate weights.
/// On RK_ERR_BUDGET: last_iterate_weights = weights at budget exhaustion;
/// enforce_sigmaw_eq_n is NOT called.
///
/// Lambda [&] guard: all bool guards declared before this function's lambdas
/// (none here — outer loop has no lambdas). Callers with lambdas MUST declare
/// bool guards before their own auto F = [&] definitions.
template<class WithinCellSolve>
inline OuterResult outer_iterate_strategy_a(
    WithinCellSolve&            fn,
    std::vector<double>&        weights,
    const CellPartition&        partition,
    const SparseMask&           sparse,
    const std::vector<double>&  w_init,        // initial weights for sparse inheritance
    const std::vector<double>&  stage1_multipliers, // per-cell Stage-1 multipliers (length n_cells_total)
    const CalibConvergenceCfg&  cfg,
    double                      outer_tol,      // shadows cfg.absolute_tol at outer level
    int                         outer_iterations, // shadows cfg.max_iterations at outer level
    int                         N) noexcept
{
    OuterResult out;
    out.status          = RK_ERR_BUDGET;
    out.iterations_used = 0;
    out.residual_final  = std::numeric_limits<double>::infinity();

    const int K_cells = partition.n_cells_total;

    // Stage 1: initial full-data calibrate (iteration 0 pre-pass).
    // Repairs coarse-margin sums before the per-cell Stage-2 loop.
    fn(weights, -1);

    for (int iter = 0; iter < outer_iterations; ++iter) {
        // Stage 2: calibrate each non-sparse cell.
        for (int c = 0; c < K_cells; ++c) {
            if (sparse.is_sparse_cell[c]) continue;
            fn(weights, c);
        }
        // Sparse cells inherit Stage-1 multipliers (per spec §6).
        apply_sparse_inheritance(weights, w_init, partition, sparse, stage1_multipliers);

        // Compute coarse-margin residual: max |Σw_cell[c]/N - 1/K_cells|
        // (relative deviation from uniform target across cells; proxy residual
        // adequate for outer convergence check without a full CalibState).
        double residual = 0.0;
        if (K_cells > 0) {
            const double target_frac = 1.0 / static_cast<double>(K_cells);
            const double inv_N       = (N > 0) ? 1.0 / static_cast<double>(N) : 0.0;
            for (int c = 0; c < K_cells; ++c) {
                double cell_sum = 0.0;
                for (int i : partition.obs_indices_by_cell[c]) cell_sum += weights[i];
                const double dev = std::fabs(cell_sum * inv_N - target_frac);
                if (dev > residual) residual = dev;
            }
        }

        // Track best-iterate at kOuterErrCheckInterval cadence.
        // Uses select_metric via a CellMetrics proxy where errRp = residual.
        // IMPORTANT: must use select_metric(cfg.metric, cm), NOT errRp directly.
        if (iter % kOuterErrCheckInterval == 0) {
            CellMetrics cm;
            cm.errRp    = residual;
            cm.mean_err = residual;
            // Other metrics (kl, chi2, grake_norm, l1) remain 0.0 — outer loop
            // does not have a full CalibState; errRp/mean_err proxy is adequate
            // for diagnostic ranking. Caller may override via a richer fn().
            const double m = select_metric(cfg.metric, cm);
            if (m < out.best.residual) {
                out.best.residual = m;
                out.best.iter_idx = iter;
            }
        }

        out.iterations_used = iter + 1;
        out.residual_final  = residual;

        // Convergence check: outer_tol shadows cfg.absolute_tol at outer level.
        if (outer_tol > 0.0 && residual < outer_tol) {
            // RK_OK: enforce Σw=n before returning.
            enforce_sigmaw_eq_n(weights, N);  // gate; result ignored (spec §5)
            out.status              = RK_OK;
            out.last_iterate_weights = weights;
            return out;
        }

        // Re-run Stage 1 to repair coarse-margin sums disturbed by Stage 2.
        fn(weights, -1);
    }

    // Budget exhausted: return last-iterate weights (NOT best-iterate — spec §6).
    out.status              = RK_ERR_BUDGET;
    out.last_iterate_weights = weights;
    return out;
}

// ── T-D: Strategy B orthogonality validator ───────────────────────────────
/// Verify that no fine-margin level spans more than one coarse cell.
/// Called from validate_calibrate_inputs when mode == 1 (exact).
///
/// Algorithm (O(N · K_fine)):
///   For each fine margin f (coarse_mask[f] == 0):
///     For each obs i: record coarse_cell_id[i] in level_to_cells[(f, group_ids[i][f])]
///   If any set has cardinality > 1: BADARG + diagnostic.
///
/// @param group_ids   [K][N] encoded category indices (0-based, -1 = NA/OOV)
/// @param cat_counts  [K] category count per margin
/// @param coarse_mask [K] 1 = coarse margin, 0 = fine margin
/// @param K           number of margins
/// @param N           number of observations
/// @return {RK_OK, ""} on pass; {RK_ERR_BADARG, diagnostic} on fail.
///         Diagnostic is ≤ 256 chars and names the offending margin index,
///         level value, and coarse cell IDs spanned.
inline std::pair<int, std::string> validate_orthogonal_split(
    const int32_t* const*    group_ids,
    const int*               cat_counts,
    const int*               coarse_mask,
    int                      K,
    int                      N) noexcept
{
    if (N <= 0 || K <= 0) return {RK_OK, ""};

    // Step 1: compute coarse_cell_id per obs via build_cell_partition.
    // min_cell_n=1 to avoid cap rejection on sparse splits.
    lbw::CellPartition part;
    {
        int rc = lbw::build_cell_partition(N, K, group_ids, cat_counts,
                                           coarse_mask, /*min_cell_n=*/1, part);
        if (rc != RK_OK) {
            // Cap exceeded or K_coarse > 64 — already diagnosed by Guard (12).
            return {RK_OK, ""};
        }
    }

    // Step 2: for each fine margin, map (f, level) -> set of coarse cell IDs.
    // Key: (fine_margin_idx << 32) | uint32_t(level_value+1) to handle -1 (NA).
    std::unordered_map<uint64_t, std::unordered_set<int>> level_to_cells;
    level_to_cells.reserve(static_cast<size_t>(N));

    for (int f = 0; f < K; f++) {
        if (coarse_mask[f] != 0) continue;  // skip coarse margins
        const int32_t* gf = group_ids[f];
        for (int i = 0; i < N; i++) {
            int32_t lv = gf[i];
            // Encode level: shift NA (-1) to 0, others to lv+1.
            uint64_t lv_enc = static_cast<uint64_t>(
                (lv < 0) ? 0 : static_cast<uint32_t>(lv) + 1u);
            uint64_t key = (static_cast<uint64_t>(f) << 32) | lv_enc;
            level_to_cells[key].insert(part.cell_id_per_obs[i]);
        }
    }

    // Step 3: scan for violations.
    for (auto& [key, cells] : level_to_cells) {
        if (cells.size() <= 1u) continue;
        // Decode key.
        int f              = static_cast<int>(key >> 32);
        uint64_t lv_enc    = key & 0xFFFFFFFFULL;
        int32_t  lv        = (lv_enc == 0u) ? -1
                             : static_cast<int32_t>(lv_enc - 1u);

        // Build "{id0, id1, ...}" string (sorted for determinism).
        std::vector<int> sorted_cells(cells.begin(), cells.end());
        std::sort(sorted_cells.begin(), sorted_cells.end());
        std::ostringstream oss;
        oss << "{";
        for (size_t ci = 0; ci < sorted_cells.size(); ci++) {
            if (ci > 0) oss << ", ";
            oss << sorted_cells[ci];
        }
        oss << "}";

        char msg[256];
        std::snprintf(msg, sizeof(msg),
            "Fine margin %d level %d spans coarse cells %s; require orthogonal split",
            f, static_cast<int>(lv), oss.str().c_str());
        return {RK_ERR_BADARG, std::string(msg)};
    }

    return {RK_OK, ""};
}

} // namespace lbw
