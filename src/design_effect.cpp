// Henry & Valliant (2015) calibration design effect + Kish (1965) weighting deff.
// Survey Methodology 41(2), 315-331 (Statistics Canada Catalogue No. 12-001-X).
// Eq 3.5 (zero-correlation approximation): deff_H ≈ deff_K · σ̂²_u / σ̂²_y.

#include "design_effect.hpp"
#include "calib_linalg.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <limits>
#include <vector>

namespace lbw {

namespace {

// Kish (1965) Eq 2.4: n * Σwᵢ² / [Σwᵢ]².
inline double kish_deff(const double* w, int n) noexcept {
    double sum_w = 0.0, sum_w2 = 0.0;
    for (int i = 0; i < n; ++i) { sum_w += w[i]; sum_w2 += w[i] * w[i]; }
    return static_cast<double>(n) * sum_w2 / (sum_w * sum_w);
}

// H&V (2015) Eq 2.7 survey-weighted population unit variance:
//   σ̂²_z = Σwᵢ(zᵢ - z̄_w)² / Σwᵢ,   z̄_w = Σwᵢzᵢ / Σwᵢ.
inline double weighted_var(const double* z, const double* w, int n) noexcept {
    double sum_w = 0.0, sum_wz = 0.0;
    for (int i = 0; i < n; ++i) { sum_w += w[i]; sum_wz += w[i] * z[i]; }
    const double zbar = sum_wz / sum_w;
    double num = 0.0;
    for (int i = 0; i < n; ++i) {
        const double d = z[i] - zbar;
        num += w[i] * d * d;
    }
    return num / sum_w;
}

}  // anonymous namespace

int design_effect_compute(
    const double* weights,
    const double* outcome,
    const int*    data_codes,
    const int*    cat_counts,
    int n, int K,
    rk_design_effect_result_t* out
) noexcept {
    out->deff_K = std::numeric_limits<double>::quiet_NaN();
    out->deff_H = std::numeric_limits<double>::quiet_NaN();
    out->rank_def = 0;
    out->message[0] = '\0';

    // ---- Validation ----
    if (n < 1) {
        std::snprintf(out->message, sizeof(out->message), "n must be >= 1");
        return RK_ERR_BADARG;
    }
    double sum_w = 0.0;
    for (int i = 0; i < n; ++i) {
        if (!std::isfinite(weights[i])) {
            std::snprintf(out->message, sizeof(out->message), "weights[%d] not finite", i);
            return RK_ERR_BADARG;
        }
        sum_w += weights[i];
    }
    if (sum_w <= 0.0) {
        std::snprintf(out->message, sizeof(out->message),
            "sum(weights) must be positive (got %g)", sum_w);
        return RK_ERR_BADARG;
    }

    // ---- Kish deff (always) ----
    out->deff_K = kish_deff(weights, n);

    // ---- 1-arg path: outcome NULL → Kish only ----
    if (outcome == nullptr) return RK_OK;

    // ---- Outcome validation ----
    for (int i = 0; i < n; ++i) {
        if (!std::isfinite(outcome[i])) {
            std::snprintf(out->message, sizeof(out->message), "outcome[%d] not finite", i);
            return RK_ERR_BADARG;
        }
    }

    // ---- K=0 short-circuit → deff_H = deff_K ----
    if (K == 0) { out->deff_H = out->deff_K; return RK_OK; }

    // ---- K-margin validation ----
    if (data_codes == nullptr || cat_counts == nullptr) {
        std::snprintf(out->message, sizeof(out->message),
            "K=%d requires non-null data_codes + cat_counts", K);
        return RK_ERR_BADARG;
    }
    // j7x8.22: accumulate the residual dof in size_t — each term is <= n <= INT_MAX,
    // but the running SUM Σ(cat_counts[k]-1) overflows int on many-margins input (e.g.
    // K=1e5, n=1e5 → ~1e10). Signed int overflow is UB; a wrapped-negative p would skip
    // the p>0 guard and sign-extend into a huge alloc. Reject p > INT_MAX before the
    // int narrowing below (all downstream uses index/dim with int p).
    std::size_t p_acc = 1;  // leafblower-xfz4: seed with the constant column (see below)
    for (int k = 0; k < K; ++k) {
        if (cat_counts[k] < 2) {
            std::snprintf(out->message, sizeof(out->message),
                "cat_counts[%d]=%d must be >= 2", k, cat_counts[k]);
            return RK_ERR_BADARG;
        }
        // CR-D14: cap cat_counts at n (mirror validation.cpp:47).
        if (cat_counts[k] > n) {
            std::snprintf(out->message, sizeof(out->message),
                "cat_counts[%d]=%d > n=%d: more categories than observations", k, cat_counts[k], n);
            return RK_ERR_BADARG;
        }
        p_acc += static_cast<std::size_t>(cat_counts[k] - 1);
    }
    if (p_acc > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        std::snprintf(out->message, sizeof(out->message),
            "total residual dof p=%zu exceeds INT_MAX: too many margins for platform", p_acc);
        return RK_ERR_BADARG;
    }
    const int p = static_cast<int>(p_acc);
    // CR-D14 (j7x8.14): guard n*K (data_codes index) and n*p (design matrix) against
    // size_t overflow before any indexing, mirroring validation.cpp's overflow guard.
    if ((K > 0 && static_cast<std::size_t>(n) > SIZE_MAX / static_cast<std::size_t>(K)) ||
        (p > 0 && static_cast<std::size_t>(n) > SIZE_MAX / static_cast<std::size_t>(p))) {
        std::snprintf(out->message, sizeof(out->message),
            "n=%d * (K=%d, p=%d) overflows size_t: problem too large for platform", n, K, p);
        return RK_ERR_BADARG;
    }
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < K; ++k) {
            const int c = data_codes[static_cast<std::size_t>(i) * K + k];
            if (c < 0 || c >= cat_counts[k]) {
                std::snprintf(out->message, sizeof(out->message),
                    "data_codes[%d,%d]=%d out of range [0,%d)", i, k, c, cat_counts[k]);
                return RK_ERR_BADARG;
            }
        }
    }

    // ---- Build X (n × p, column-major): constant column + drop-first-level dummies ----
    // Column 0 is the constant vector 1 (all n rows); margin k's dummy for category
    // c >= 1 lands at column col_off + c - 1, with col_off starting at 1 and advancing by
    // cat_counts[k] - 1 per margin. The calibration constraints always include the
    // population total, so the constant belongs to the calibration column space: without
    // it, u_i is not a true GREG residual (the reference cell's residual equals its raw
    // outcome), which can push σ̂²_u above σ̂²_y and violate the deff_H <= deff_K
    // invariant that Eq 3.5 guarantees once 1 ∈ col(X). Henry & Valliant (2015) Eq 3.5;
    // leafblower-xfz4.
    std::vector<double> X(static_cast<std::size_t>(n) * p, 0.0);
    {
        for (int i = 0; i < n; ++i) X[i] = 1.0;  // column 0: constant
        int col_off = 1;
        for (int k = 0; k < K; ++k) {
            for (int i = 0; i < n; ++i) {
                const int c = data_codes[static_cast<std::size_t>(i) * K + k];
                if (c >= 1) X[static_cast<std::size_t>(col_off + c - 1) * n + i] = 1.0;
            }
            col_off += cat_counts[k] - 1;
        }
    }

    // ---- X^T W X (p × p, column-major, LOWER triangle for uplo='L') ----
    // calib_linalg.hpp cholesky_factor_inplace uses dpotrf(uplo='L') → fill lower.
    // Lower triangle: row l >= col j → store at XtWX[j * p + l].
    std::vector<double> XtWX(static_cast<std::size_t>(p) * p, 0.0);
    std::vector<double> XtWy(p, 0.0);
    for (int j = 0; j < p; ++j) {
        const double* X_j = &X[static_cast<std::size_t>(j) * n];
        for (int i = 0; i < n; ++i) {
            const double wXji = weights[i] * X_j[i];
            // lower triangle: l >= j
            for (int l = j; l < p; ++l) {
                XtWX[static_cast<std::size_t>(j) * p + l] += wXji * X[static_cast<std::size_t>(l) * n + i];
            }
            XtWy[j] += wXji * outcome[i];
        }
    }

    // ---- Cholesky factor + solve: β̂ in XtWy ----
    constexpr double kEpsPerturb = 1e-12;
    const int factor_status = cholesky_factor_inplace(
        XtWX.data(), static_cast<std::size_t>(p), kEpsPerturb);
    if (factor_status != RK_OK) {
        out->rank_def = 1;
        out->deff_H = out->deff_K;
        std::snprintf(out->message, sizeof(out->message),
            "calibration margins rank-deficient; deff_H=deff_K");
        return RK_OK;
    }
    cholesky_solve(XtWX.data(), static_cast<std::size_t>(p), XtWy.data());

    // ---- GREG residuals û = y - X β̂ ----
    std::vector<double> u(n);
    for (int i = 0; i < n; ++i) {
        double xi_beta = 0.0;
        for (int j = 0; j < p; ++j)
            xi_beta += X[static_cast<std::size_t>(j) * n + i] * XtWy[j];
        u[i] = outcome[i] - xi_beta;
    }

    // ---- σ̂²_y, σ̂²_u → deff_H (H&V 2015 Eq 3.5) ----
    const double sigma2_y = weighted_var(outcome, weights, n);
    if (sigma2_y < std::numeric_limits<double>::epsilon()) {
        out->deff_H = out->deff_K;
        return RK_OK;
    }
    out->deff_H = out->deff_K * weighted_var(u.data(), weights, n) / sigma2_y;
    return RK_OK;
}

}  // namespace lbw
