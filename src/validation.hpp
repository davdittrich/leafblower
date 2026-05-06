#pragma once
#include "leafblower.h"
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <climits>
#include <limits>
#include <vector>

// T-B owns canonical site; remove duplicate when T-B merges.
// TODO(T-B): promote to calib_dispatch.hpp as LBW_MAX_HIER_CELLS
#ifndef LBW_MAX_HIER_CELLS
#  define LBW_MAX_HIER_CELLS 100000
#endif

namespace lbw {

// Validate inputs for rk_calibrate. Returns RK_OK or RK_ERR_BADARG.
// On RK_ERR_BADARG, sets result->status and result->message if result != nullptr.
// This is the single source of truth for input validation — called from both
// c_api.cpp (rk_calibrate C ABI path) and r_bridge.cpp (direct C++ path).
inline int validate_calibrate_inputs(int n, int K,
    const double* weights, const int32_t** group_ids,
    const int* cat_counts, const double** targets,
    const rk_params_t* p, rk_result_t* result, rk_algorithm_t alg)
{
    auto err = [&](const char* msg) -> int {
        if (result) {
            result->status = RK_ERR_BADARG;
            snprintf(result->message, 256, "%s", msg);
        }
        return RK_ERR_BADARG;
    };

    if (!weights)   return err("weights pointer is NULL");
    if (!group_ids) return err("group_ids pointer is NULL");
    if (!cat_counts)return err("cat_counts pointer is NULL");
    if (!targets)   return err("targets pointer is NULL");
    if (n <= 0)     return err("n must be > 0");
    if (K <= 0)     return err("K must be > 0");
    if (K > 64)     return err("K exceeds maximum (64); too many margin columns");

    if (p->min_weight >= p->max_weight)
        return err("min_weight must be strictly less than max_weight");

    // cat_counts checks + overflow guard
    size_t total_cats = 0;
    for (int k = 0; k < K; k++) {
        if (cat_counts[k] <= 0)
            return err("cat_counts[k] must be > 0 for all k");
        if (cat_counts[k] > n)
            return err("cat_counts[k] > n: more categories than observations");
        total_cats += (size_t)cat_counts[k];
    }
    if (total_cats > 0 && (size_t)n > SIZE_MAX / 2 / (size_t)total_cats)
        return err("problem too large for platform size_t");

    // Initial weight checks
    double total_w = 0.0;
    for (int i = 0; i < n; i++) {
        if (!std::isfinite(weights[i]))
            return err("NaN or Inf in initial weights[]");
        total_w += weights[i];
    }
    if (total_w < 1e-15)
        return err("total weight is zero or negative");

    // targets checks
    for (int k = 0; k < K; k++) {
        if (!targets[k]) return err("targets[k] is NULL");
        double sum = 0.0;
        for (int j = 0; j < cat_counts[k]; j++) {
            if (!std::isfinite(targets[k][j]))
                return err("NaN or Inf in targets[]");
            if (targets[k][j] < 0.0)
                return err("targets[k][j] < 0");
            sum += targets[k][j];
        }
        if (std::fabs(sum - 1.0) > 1e-6)
            return err("targets[k] does not sum to 1 (within 1e-6)");
    }

    // group_ids range validation — full O(n*K) pass before any weight modification
    for (int k = 0; k < K; k++) {
        if (!group_ids[k]) return err("group_ids[k] is NULL");
        for (int i = 0; i < n; i++) {
            int g = group_ids[k][i];
            if (g < -1)
                return err("group_ids[k][i] < -1: only -1 (NA) is valid");
            if (g >= cat_counts[k])
                return err("group_ids[k][i] >= cat_counts[k]");
        }
    }

    // B13: NA-only category with positive target → structural INFEAS.
    // O(n×K) two-pass bitmap: replaces O(K×ΣC×n) triple-nested scan.
    // Pass 1: mark which (k,j) combinations appear in data — O(n×K)
    std::vector<std::vector<bool>> seen(K);
    for (int k = 0; k < K; k++) seen[k].assign(cat_counts[k], false);
    for (int i = 0; i < n; i++)
        for (int k = 0; k < K; k++) {
            int g = group_ids[k][i];
            if (g >= 0 && g < cat_counts[k]) seen[k][g] = true;
        }
    // Pass 2: check targeted but not-in-data categories — O(K×max_cats)
    for (int k = 0; k < K; k++) {
        for (int j = 0; j < cat_counts[k]; j++) {
            if (!seen[k][j] && targets[k][j] > 1e-12) {
                char msg[256];
                std::snprintf(msg, sizeof(msg),
                    "margin %d, category %d: target=%.6g but 0 observations assigned "
                    "(all NA or missing); structurally infeasible problem",
                    k, j, targets[k][j]);
                if (result) {
                    result->status = RK_ERR_INFEAS;
                    std::snprintf(result->message, 256, "%s", msg);
                }
                return RK_ERR_INFEAS;
            }
        }
    }

    if (!std::isfinite(p->tol_abs) || p->tol_abs <= 0.0)
        return err("tol_abs must be finite and positive");

    /* Enum range guards: prevent silent UB from out-of-range static_cast */
    if (p->metric < 0 || p->metric > 6)
        return err("metric out of range [0,6]: 0=MAX_ERR 1=MEAN_ERR 2=KL 3=CHI2 4=GRAKE_NORM 5=L1_WEIGHT 6=MARGINAL_KL");
    if (p->rule < 0 || p->rule > 2)
        return err("rule out of range [0,2]: 0=THRESHOLD 1=IMPROVEMENT 2=PLATEAU");
    if (p->stop_when < 0 || p->stop_when > 1)
        return err("stop_when out of range [0,1]: 0=ANY 1=ALL");

    // ── T-E: Hierarchical 2-stage input validation ────────────────────────────
    // Guard (1): early-out when hierarchical disabled.
    if (p->hierarchical_enabled != 0) {
        // Helpers
        const int*  mask      = p->hierarchical_coarse_mask;
        const int   min_cn    = p->hierarchical_min_cell_n;
        const int   mode      = p->hierarchical_mode;
        const double outer_tol= p->hierarchical_outer_tol;
        const int   outer_it  = p->hierarchical_outer_iterations;
        char buf[256];

        // NOTE: coarse_mask length == K is checked in r_bridge.cpp (via SEXP
        // LENGTH) before this function is called. The c_api.cpp caller must
        // enforce length independently (it owns the array length contract).

        // Guard (3): coarse_mask must be non-empty (at least one 1-bit).
        {
            int popcount = 0;
            for (int k = 0; k < K; k++) popcount += (mask[k] != 0 ? 1 : 0);
            if (popcount == 0)
                return err("hierarchical: coarse_margins is empty (all mask values are 0)");

            // Guard (4): coarse_mask must be a proper subset (not all K).
            if (popcount == K) {
                return err("hierarchical: coarse_margins == margins; nothing to refine "
                           "-- disable hierarchical or exclude at least one margin");
            }
        }

        // Guard (5): min_cell_n must be >= 1.
        // Guard (6): min_cell_n > N emits a warning (not BADARG) — handled by the
        // caller (r_bridge.cpp) via Rf_warning after this function returns RK_OK.
        if (min_cn < 1) {
            std::snprintf(buf, sizeof(buf),
                "min_cell_n=%d invalid; must satisfy 1 <= min_cell_n <= N=%d",
                min_cn, n);
            return err(buf);
        }

        // Guard (7): mode ∈ {0 = refine, 1 = exact}.
        if (mode < 0 || mode > 1) {
            std::snprintf(buf, sizeof(buf),
                "hierarchical_mode=%d invalid; must be 0 (refine) or 1 (exact)", mode);
            return err(buf);
        }

        // Guard (8): outer_tol must be finite, positive, >= machine_eps * N.
        {
            const double machine_eps = std::numeric_limits<double>::epsilon();
            const double min_tol     = machine_eps * static_cast<double>(n);
            if (!std::isfinite(outer_tol) || outer_tol <= 0.0) {
                std::snprintf(buf, sizeof(buf),
                    "hierarchical_outer_tol=%.6g invalid; must be finite and positive",
                    outer_tol);
                return err(buf);
            }
            if (outer_tol < min_tol) {
                std::snprintf(buf, sizeof(buf),
                    "hierarchical_outer_tol=%.6g < machine_eps*N=%.6g; "
                    "increase outer_tol to avoid infinite outer loops",
                    outer_tol, min_tol);
                return err(buf);
            }
        }

        // Guard (9): outer_iterations in [1, 10000].
        if (outer_it < 1 || outer_it > 10000) {
            std::snprintf(buf, sizeof(buf),
                "hierarchical_outer_iterations=%d invalid; must satisfy "
                "1 <= outer_iterations <= 10000",
                outer_it);
            return err(buf);
        }

        // Guard (10): algorithm × mode constraint.
        // P2 methods (LOGIT=10, GREG=6) cannot use mode=0 (refine) — Newton
        // re-factorization per fine cell is cost-prohibitive.
        if (mode == 0 &&
            (alg == RK_ALG_LOGIT || alg == RK_ALG_GREG)) {
            const char* alg_name = (alg == RK_ALG_LOGIT) ? "logit" : "greg";
            std::snprintf(buf, sizeof(buf),
                "algorithm '%s' cannot use mode='refine' (Newton re-factor cost "
                "prohibitive); pass mode='exact' and ensure orthogonal coarse/fine split",
                alg_name);
            return err(buf);
        }

        // Guard (11): bounds_mode='unit' incompatible with hierarchical.
        // Check fires AT validate entry, not at solver dispatch (per spec §V.4).
        if (p->bounds_mode == RK_BOUNDS_UNIT) {
            return err("bounds_mode='unit' incompatible with hierarchical calibration; "
                       "use bounds_mode='cell' or disable hierarchical");
        }

        // Guard (12): preliminary cell-count cap.
        // cells_estimated = product of cat_counts[k] for coarse dimensions only.
        // Cap = min(N / min_cell_n, LBW_MAX_HIER_CELLS).
        // Skip when min_cn > n: that case degenerates to single-stage (Guard 6 warning
        // is emitted by r_bridge.cpp); integer division N/min_cn would floor to 0.
        if (min_cn <= n) {
            long long cells_est = 1LL;
            bool overflow = false;
            for (int k = 0; k < K; k++) {
                if (mask[k] == 0) continue;
                if (cells_est > (long long)LBW_MAX_HIER_CELLS) { overflow = true; break; }
                cells_est *= (long long)cat_counts[k];
                if (cells_est > (long long)LBW_MAX_HIER_CELLS) { overflow = true; break; }
            }
            // cap = min(N / min_cell_n, LBW_MAX_HIER_CELLS); always >= 1 since min_cn <= n.
            long long cap = (long long)LBW_MAX_HIER_CELLS;
            {
                long long n_over_mincn = (long long)n / (long long)min_cn;
                if (n_over_mincn < cap) cap = n_over_mincn;
            }
            if (overflow || cells_est > cap) {
                long long ce = overflow ? (long long)LBW_MAX_HIER_CELLS + 1LL : cells_est;
                std::snprintf(buf, sizeof(buf),
                    "hierarchical: estimated coarse cells=%lld exceeds cap=%lld "
                    "(N=%d, min_cell_n=%d, LBW_MAX_HIER_CELLS=%d); "
                    "increase min_cell_n or reduce coarse_margins cardinality",
                    (long long)ce, (long long)cap, n, min_cn, LBW_MAX_HIER_CELLS);
                return err(buf);
            }
        }

        // Guard (13): Strategy B orthogonality — delegated to T-D validator.
        // TODO(T-D): call validate_orthogonal_split(mask, K, cat_counts, result) here
        // when mode==1. T-D adds the implementation; this call site is the hook.
    }
    // ── End T-E hierarchical validation ──────────────────────────────────────

    return RK_OK;
}

} // namespace lbw
