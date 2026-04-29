#pragma once
// src/sraa.hpp — Safeguarded Regularized Anderson Acceleration (Type II AA)
// Used by greenkhorn.cpp and raking.cpp to replace SQUAREM CBB acceleration.
#include "calib_linalg.hpp"
#include "lbw_config.h"
#include <R.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace lbw {

static constexpr int    kSRAAMaxM         = 10;    // stack array bound; init() enforces window <= this
static constexpr int    kSRAAm            = 5;     // default window: ~176 MB at stepstone scale
static constexpr int    kSRAAMinCount     = 2;     // min DX/DR pairs before AA fires
static constexpr double kSRAAdeltaReg     = 1e-10; // relative Tikhonov on Gram matrix
static constexpr double kSRAARestartGamma = 2.0;   // restart when ||R_k||^2 > 4 x prev_norm
static constexpr double kSRAAOuterSlack       = 0.10; // 10% above best_errRp → stall
static constexpr int    kSRAAOuterStallWindow = 5;    // 5 stall iters → revert+restart

struct SRAAStepResult {
    bool   aa_accepted;
    int    f_evals;     // 1 (plain) or 2 (AA attempted)
    double err_result;  // max_error of accepted step
};

struct SRAAState {
    int M = 0, m = 0, head = 0, count = 0;
    bool has_prev = false;
    int aa_accepted_count = 0; // cumulative; NOT reset by clear()

    // Slot-contiguous: buf[slot*M + cell]. Each slot is contiguous M doubles -> SIMD dot products.
    std::vector<double> dX_buf;   // m x M: ΔX differences
    std::vector<double> dR_buf;   // m x M: ΔR differences
    std::vector<double> R_prev;   // M
    std::vector<double> X_prev;   // M
    std::vector<double> F_cur;    // M: F(X_k)
    std::vector<double> scratch;  // M: R_k temp, then X_AA, then F(X_AA)

    double gram[kSRAAMaxM * kSRAAMaxM] = {};
    double rhs[kSRAAMaxM] = {};
    double gamma_[kSRAAMaxM] = {}; // gamma_ avoids <cmath> gamma conflict

    double prev_resid_norm = 0.0;  // read only when has_prev=true

    void init(int M_cell, int window) {
        if (window > kSRAAMaxM)
            Rf_error("SRAA window %d exceeds kSRAAMaxM=%d", window, kSRAAMaxM);
        M = M_cell; m = window;
        try {
            dX_buf.assign((size_t)m * M, 0.0);
            dR_buf.assign((size_t)m * M, 0.0);
            R_prev.assign(M, 0.0); X_prev.assign(M, 0.0);
            F_cur.assign(M, 0.0);  scratch.assign(M, 0.0);
        } catch (std::bad_alloc&) {
            Rf_error("SRAA: out of memory allocating %.0f MB (m=%d, M=%d)",
                     (2.0*m + 4.0) * M * 8.0 / 1e6, m, M);
        }
        clear();
    }

    // Resets history. aa_accepted_count NOT reset (cumulative).
    // prev_resid_norm reset to 0 (safe sentinel; never read when has_prev=false).
    void clear() { head = 0; count = 0; has_prev = false; prev_resid_norm = 0.0; }
};

// sraa_step: one AA super-step.
// f_eval: (std::vector<double>&) -> double
//   Modifies Xv in-place (K steps), returns max_errRp.
//   INVARIANT: must receive state.F_cur or state.scratch, NEVER the outer X.
// X is updated in-place via std::swap (O(1)).
template<typename FEval>
SRAAStepResult sraa_step(
    FEval& f_eval,
    std::vector<double>& X,
    const std::vector<double>& L_cell,
    const std::vector<double>& U_cell,
    SRAAState& state)
{
    const int M = state.M;

    // --- Step 1: F(X_k) -> F_cur; compute R_k into scratch; compute norm ---
    double err_plain = f_eval(state.F_cur);
    double norm_k = 0.0;
    for (int c = 0; c < M; c++) {
        double rk = state.F_cur[c] - X[c];
        state.scratch[c] = rk;  // temporarily holds R_k
        norm_k += rk * rk;
    }

    // --- Step 2: Restart check (guarded by has_prev — never fires on first call) ---
    if (state.has_prev &&
        norm_k > kSRAARestartGamma * kSRAARestartGamma * state.prev_resid_norm) {
        state.clear();
        std::swap(X, state.F_cur);
        return {false, 1, err_plain};
    }

    // --- Step 3: Append DX, DR to circular buffer (only when has_prev) ---
    if (state.has_prev) {
        double* dX_s = state.dX_buf.data() + state.head * M;
        double* dR_s = state.dR_buf.data() + state.head * M;
        for (int c = 0; c < M; c++) {
            dX_s[c] = X[c]              - state.X_prev[c];
            dR_s[c] = state.scratch[c]  - state.R_prev[c];  // R_k - R_{k-1}
        }
        state.head  = (state.head + 1) % state.m;
        state.count = std::min(state.count + 1, state.m);
    }

    // --- Step 4: Update prev state ---
    for (int c = 0; c < M; c++) {
        state.X_prev[c] = X[c];
        state.R_prev[c] = state.scratch[c];  // R_k
    }
    state.prev_resid_norm = norm_k;
    state.has_prev = true;

    // --- Step 5: Not enough history -> plain ---
    if (state.count < kSRAAMinCount) {
        std::swap(X, state.F_cur);
        return {false, 1, err_plain};
    }

    // --- Step 6: Gram matrix G[i][j] = DR[i].DR[j]; RHS rhs[i] = DR[i].R_k ---
    const int n = state.count;
    double max_diag = 0.0;
    for (int i = 0; i < n; i++) {
        const double* dRi = state.dR_buf.data() + i * M;
        for (int j = i; j < n; j++) {
            const double* dRj = state.dR_buf.data() + j * M;
            double dot = 0.0;
            for (int c = 0; c < M; c++) dot += dRi[c] * dRj[c];
            state.gram[i * kSRAAMaxM + j] = state.gram[j * kSRAAMaxM + i] = dot;
        }
        max_diag = std::max(max_diag, state.gram[i * kSRAAMaxM + i]);
    }
    for (int i = 0; i < n; i++) {
        const double* dRi = state.dR_buf.data() + i * M;
        double dot = 0.0;
        for (int c = 0; c < M; c++) dot += dRi[c] * state.scratch[c];
        state.rhs[i] = dot;
    }

    // --- Step 7: Regularize + LDLT solve (n x n submatrix) ---
    double eps = kSRAAdeltaReg * (max_diag > 0.0 ? max_diag : 1.0);
    double G[kSRAAMaxM * kSRAAMaxM];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            G[i * n + j] = state.gram[i * kSRAAMaxM + j];
    for (int i = 0; i < n; i++) G[i * n + i] += eps;

    if (lbw::ldlt_factor_inplace(G, (size_t)n, 0.0) != RK_OK) {
        state.clear();
        std::swap(X, state.F_cur);
        return {false, 1, err_plain};
    }
    for (int i = 0; i < n; i++) state.gamma_[i] = state.rhs[i];
    lbw::ldlt_solve(G, (size_t)n, state.gamma_);

    // --- Step 8: Extrapolate + clamp into scratch ---
    // scratch currently holds R_k; overwrite in-place with X_AA
    for (int c = 0; c < M; c++) {
        double Rk_c = state.scratch[c];  // read R_k before overwrite
        double corr = 0.0;
        for (int i = 0; i < n; i++)
            corr += state.gamma_[i] * (state.dX_buf[i * M + c] + state.dR_buf[i * M + c]);
        state.scratch[c] = std::clamp(X[c] + Rk_c - corr, L_cell[c], U_cell[c]);
    }

    // --- Step 9: F(X_AA); scratch -> F(X_AA) ---
    double err_AA = f_eval(state.scratch);

    // NaN/inf guard (includes Wv=0 case from degenerate AA input — see greenkhorn.cpp)
    if (!std::isfinite(err_AA)) {
        state.clear();
        std::swap(X, state.F_cur);
        return {false, 1, err_plain};
    }

    // --- Step 10: Safeguard ---
    // Accept AA if it matches or beats the plain F-step quality.
    // Note: both err_AA and err_plain are from f_eval (fixed K-step sort).
    // On K>=3 overlapping-margin problems, SRAA may still converge to a
    // KL-optimal (vs max_err-optimal) fixed point — this is a known limitation.
    // Full fix requires outer quality tracking from the adaptive-sort greenkhorn loop.
    if (err_AA <= err_plain) {
        std::swap(X, state.scratch);    // O(1); X = F(X_AA)
        state.aa_accepted_count++;
        return {true, 2, err_AA};
    } else {
        state.clear();                  // bad step; reset history
        std::swap(X, state.F_cur);     // X = F(X_k)
        return {false, 2, err_plain};
    }
}

} // namespace lbw
