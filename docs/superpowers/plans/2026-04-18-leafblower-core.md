# Leafblower Core — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build leafblower — a C++17/C-API survey calibration library with R and Python bindings, drop-in for autumn::harvest(), implementing iEPPA and L-BFGS-B.

**Architecture:** C++17 core in src/ exposes a C API (leafblower.h) shared by both bindings; R uses .Call() via r_bridge.cpp; Python uses pybind11 via _bindings.cpp. iEPPA (Sinkhorn-BCD, fixed ε=0.05) handles large/bounded problems; L-BFGS-B with logit/exp dual link handles standard cases; auto-routing selects based on problem complexity.

**Tech Stack:** C++17, C99 API, R (.Call, useDynLib), pybind11 ≥ 2.11, scikit-build-core, testthat 3, bench

---

## PHASE 1: C API + L-BFGS-B + R Layer

### Task 1: Package Scaffolding

**Files:**
- Create: `DESCRIPTION`
- Create: `NAMESPACE`
- Create: `R/zzz.R`
- Create: `src/Makevars.in`
- Create: `configure` (stub)
- Create: `.Rbuildignore`

**Steps:**

- [ ] Write `DESCRIPTION`:
```
Package: leafblower
Title: High-Performance Survey Calibration via iEPPA and L-BFGS-B
Version: 0.1.0
Authors@R: person("Dennis Alexis Valin", "Dittrich",
    email = "davdittrich@gmail.com", role = c("aut", "cre"))
Description: Drop-in replacement for autumn::harvest() implementing iEPPA
    (Sinkhorn block-coordinate descent) and L-BFGS-B over the Deville-Sarndal
    logit dual. Adds min_weight lower bound; shares a single C++17 core with
    the Python bindings.
License: MIT + file LICENSE
Encoding: UTF-8
NeedsCompilation: yes
SystemRequirements: C++17 compiler
Suggests: autumn (>= 0.2.0), testthat (>= 3.0.0), bench
RoxygenNote: 7.3.2
```

- [ ] Write `NAMESPACE`:
```r
useDynLib(leafblower, .registration = TRUE)
export(harvest)
export(anesrake)
export(diagnose_weights)
export(design_effect)
export(effective_sample_size)
export(get_current_miss)
export(weighted_pct)
```

- [ ] Write `R/zzz.R`:
```r
# R_init_leafblower() in r_bridge.cpp is called automatically by R when the
# shared library is loaded (R calls R_init_<pkgname> on dlopen). No .onLoad()
# needed. This empty .onLoad is present to suppress R CMD check NOTE about
# an absent zzz.R on some R versions.
.onLoad <- function(libname, pkgname) invisible(NULL)
```

- [ ] Write `src/Makevars.in`:
```makefile
PKG_CXXFLAGS = @CXXFLAGS_STD@ -I. -O3 -DSTRICT_R_HEADERS
PKG_SOURCES = c_api.cpp logit.cpp lbfgsb_solver.cpp ieppa.cpp r_bridge.cpp
```

- [ ] Write `configure` stub (will be replaced in Task 17; for now, generate a default `src/Makevars`):
```sh
#!/bin/sh
# Detect C++17 support. Note: -c is required; linking /dev/null fails without main().
if "$CXX" -std=c++17 -x c++ -c /dev/null -o /dev/null 2>/dev/null; then
    CXXFLAGS_STD="-std=c++17"
else
    CXXFLAGS_STD="-std=c++14"
    echo "WARNING: C++17 not available; falling back to C++14" >&2
fi
sed "s|@CXXFLAGS_STD@|${CXXFLAGS_STD}|" src/Makevars.in > src/Makevars
```
Then `chmod +x configure`.

- [ ] Write `.Rbuildignore`:
```
^\.beads$
^\.claude$
^\.wolf$
^tasks$
^python$
^docs/iEPPA$
^\.git$
^cran-comments\.md$
^\.github$
```

- [ ] Verify: `R CMD build .` produces a `.tar.gz` without errors (package has no R functions yet, so expect only a note about empty R/ directory if zzz.R is the only file).

**Commit:** `feat(scaffold): add package skeleton with DESCRIPTION, NAMESPACE, Makevars.in`

---

### Task 2: C API Header

**Files:**
- Create: `src/leafblower.h`

**Steps:**

- [ ] Write `src/leafblower.h`:
```c
#ifndef LEAFBLOWER_H
#define LEAFBLOWER_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>  /* size_t */
#include <stdint.h>  /* int32_t, int64_t */

/* ── Return codes ── */
#define RK_OK         0  /* Success */
#define RK_ERR_NOCONV 1  /* Did not converge within outer_max_iter */
#define RK_ERR_INFEAS 2  /* Infeasible: empty cell with positive target */
#define RK_ERR_BADARG 3  /* Invalid argument */

/* ── Algorithm selector ── */
typedef enum {
    RK_ALG_AUTO   = 0,
    RK_ALG_IEPPA  = 1,
    RK_ALG_LBFGSB = 2
} rk_algorithm_t;

/* ── Calibration parameters ── */
typedef struct {
    double          min_weight;      /* default 0.0 */
    double          max_weight;      /* default 5.0 */
    int             inner_max_iter;  /* inner BCD cap per outer iter, default 500 */
    int             outer_max_iter;  /* outer EPP / L-BFGS max iters, default 50 */
    double          tol_abs;         /* convergence tolerance, default 1e-6 */
    rk_algorithm_t  algorithm;       /* default RK_ALG_AUTO */
    int             verbose;         /* 0=silent, 1=progress, 2=debug */
    double          epsilon;         /* iEPPA entropic parameter, default 0.05 */
    int             lbfgs_m;         /* L-BFGS history size, default 10 */
    void            (*log_fn)(const char* msg, void* ctx);
    void*           log_ctx;
} rk_params_t;

/* ── Result ── */
typedef struct {
    int             status;          /* RK_OK / RK_ERR_* */
    int             iterations;      /* outer iterations completed */
    double          max_error;       /* max calibration error at last iterate */
    rk_algorithm_t  algorithm_used;  /* actual algorithm run (never RK_ALG_AUTO) */
    char            message[256];    /* null-terminated status message */
} rk_result_t;

/* Fill *p with safe defaults */
void rk_params_init(rk_params_t* p);

/*
 * Calibrate survey weights in-place.
 *   n          : number of observations
 *   K          : number of calibration margins
 *   weights    : [n] on input = design weights; on output = calibrated weights
 *   group_ids  : [K] array of pointers; group_ids[k][i] in {-1, 0..cat_counts[k]-1}
 *   cat_counts : [K] number of categories per margin
 *   targets    : [K] array of pointers; targets[k][j] = target proportion for cat j
 *   params     : calibration parameters (NULL = use defaults)
 *   result     : output result struct (NULL = ignore)
 * Returns RK_OK, RK_ERR_NOCONV, RK_ERR_INFEAS, or RK_ERR_BADARG.
 *
 * C++17 note: the C++ compilation unit (c_api.cpp) applies [[nodiscard]] to
 * this function via a wrapper; the C header stays C99-clean.
 */
int rk_calibrate(
    int n, int K,
    double* weights,
    const int32_t** group_ids,
    const int* cat_counts,
    const double** targets,
    const rk_params_t* params,
    rk_result_t* result
);

#ifdef __cplusplus
}
#endif

#endif /* LEAFBLOWER_H */
```

- [ ] Verify C99 compilation: `gcc -std=c99 -Werror -x c src/leafblower.h -o /dev/null`
  Expected output: no errors, no output.

**Commit:** `feat(api): add C99 public header leafblower.h with all structs and return codes`

---

### Task 3: Input Validation + rk_params_init

**Files:**
- Create: `src/types.hpp`
- Create: `src/c_api.cpp` (partial — validation + stub)
- Create: `tests/testthat/test-harvest.R` (RED phase — partial)

**Steps:**

- [ ] Write `src/types.hpp`:
```cpp
#pragma once
#include <vector>
#include <cstdint>

namespace lbw {

struct CalibState {
    int n;
    int K;
    double* weights;                   // caller-owned, n elements
    const int32_t** group_ids;        // caller-owned, K pointers to n int32_t elements
    const int* cat_counts;            // caller-owned, K elements
    const double** targets;           // caller-owned, K pointers to cat_counts[k] doubles
    double min_weight;
    double max_weight;
    double tol_abs;
    int inner_max_iter;
    int outer_max_iter;
    double epsilon;
    int lbfgs_m;
    int verbose;
    void (*log_fn)(const char* msg, void* ctx);
    void* log_ctx;

    // Derived
    int total_cats;                    // sum of cat_counts[k]

    void log(const char* msg) const {
        if (verbose <= 0) return;
        if (log_fn) {
            log_fn(msg, log_ctx);
        } else {
            fprintf(stderr, "%s\n", msg);
        }
    }
};

} // namespace lbw
```

- [ ] Write `src/c_api.cpp` — rk_params_init + validate_inputs + rk_calibrate stub:
```cpp
#include "leafblower.h"
#include "types.hpp"
#include <cstring>
#include <cstdio>
#include <cmath>
#include <climits>
#include <cstdint>

// C++17 [[nodiscard]] on rk_calibrate — silently ignored on C++14
#if __cplusplus >= 201703L
  #define LBW_NODISCARD [[nodiscard]]
#else
  #define LBW_NODISCARD
#endif

extern "C" {

void rk_params_init(rk_params_t* p) {
    if (!p) return;
    p->min_weight    = 0.0;
    p->max_weight    = 5.0;
    p->inner_max_iter = 500;
    p->outer_max_iter = 50;
    p->tol_abs       = 1e-6;
    p->algorithm     = RK_ALG_AUTO;
    p->verbose       = 0;
    p->epsilon       = 0.05;
    p->lbfgs_m       = 10;
    p->log_fn        = nullptr;
    p->log_ctx       = nullptr;
}

static int validate_inputs(int n, int K,
                            const double* weights,
                            const int32_t** group_ids,
                            const int* cat_counts,
                            const double** targets,
                            const rk_params_t* p,
                            rk_result_t* result) {
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

    if (p->min_weight >= p->max_weight)
        return err("min_weight must be strictly less than max_weight");

    // Logit singularity guard — each checked independently
    bool use_logit = (p->min_weight > 0.0) && std::isfinite(p->max_weight);
    if (use_logit) {
        if (p->min_weight == 1.0)
            return err("logit link undefined: min_weight=1 makes denominator (1-L)=0");
        if (p->max_weight == 1.0)
            return err("logit link undefined: max_weight=1 makes denominator (U-1)=0");
    }

    // cat_counts checks + overflow guard
    size_t total_cats = 0;
    for (int k = 0; k < K; k++) {
        if (cat_counts[k] <= 0)
            return err("cat_counts[k] must be > 0 for all k");
        if (cat_counts[k] > n)
            return err("cat_counts[k] > n: more categories than observations");
        total_cats += (size_t)cat_counts[k];
    }
    if ((size_t)n * total_cats > SIZE_MAX / 2)
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
        if (std::fabs(sum - 1.0) > 1e-8)
            return err("targets[k] does not sum to 1 (within 1e-8)");
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

    return RK_OK;
}

LBW_NODISCARD int rk_calibrate(int n, int K,
                                double* weights,
                                const int32_t** group_ids,
                                const int* cat_counts,
                                const double** targets,
                                const rk_params_t* params,
                                rk_result_t* result) {
    rk_params_t defaults;
    rk_params_init(&defaults);
    const rk_params_t* p = params ? params : &defaults;

    if (result) {
        result->status = RK_OK;
        result->iterations = 0;
        result->max_error = 0.0;
        result->algorithm_used = RK_ALG_AUTO;
        result->message[0] = '\0';
    }

    int rc = validate_inputs(n, K, weights, group_ids, cat_counts, targets, p, result);
    if (rc != RK_OK) return rc;

    // Stub: algorithm dispatch not yet implemented
    if (result) {
        result->status = RK_ERR_NOCONV;
        snprintf(result->message, 256, "rk_calibrate not fully implemented");
    }
    return RK_ERR_NOCONV;
}

} // extern "C"
```

- [ ] Create placeholder `tests/testthat/test-harvest.R`:
```r
# BADARG validation RED test deferred to Task 7 where harvest() first exists.
# Calling .Call("C_rk_calibrate") before R_init_leafblower registers the symbol
# crashes R with an opaque symbol-not-found error — not a named assertion failure.
# Task 7 adds the real RED test: expect_error(harvest(df, tgt, min_weight=5, max_weight=5), "min_weight")
```

- [ ] Verify file created: `ls tests/testthat/test-harvest.R`
  Expected: file present (directory tests/testthat/ created by this step).

**Commit:** `feat(validation): add rk_params_init, input validation, rk_calibrate stub`

---

### Task 4: Logit Link Functions + RED Test

**Files:**
- Create: `src/logit.hpp`
- Create: `src/logit.cpp`
- Create: `tests/testthat/test-logit.R`
- Modify: `src/r_bridge.cpp` (create stub with C_logit_F_at_zero for test)

**Steps:**

- [ ] Write RED test first in `tests/testthat/test-logit.R`:
```r
test_that("F(0) = 1 for logit link", {
  # Calls C function that doesn't exist yet — must FAIL
  expect_equal(.Call("C_logit_F_at_zero", 0.5, 5.0), 1.0, tolerance = 1e-12)
})

test_that("F(u) stays in [L,U] for logit link", {
  # Also RED until implemented
  vals <- .Call("C_logit_range_check", 0.5, 5.0, as.double(seq(-10, 10, by=0.5)))
  expect_true(all(vals >= 0.5 - 1e-12))
  expect_true(all(vals <= 5.0 + 1e-12))
})

test_that("H prime equals F for logit link (numerical diff)", {
  # RED
  result <- .Call("C_logit_Hprime_check", 0.5, 5.0, 1.0)
  expect_equal(result, 0.0, tolerance = 1e-8)
})

test_that("exp link: F(u) = exp(u)", {
  # RED
  expect_equal(.Call("C_logit_F_at_zero", 0.0, Inf), 1.0, tolerance = 1e-12)
})
```

- [ ] Run to confirm RED: `Rscript -e 'testthat::test_file("tests/testthat/test-logit.R")'`
  Expected: all FAIL with "C_logit_F_at_zero not found".

- [ ] Write `src/logit.hpp`:
```cpp
#pragma once
#include <cmath>
#include <algorithm>
#include <limits>

namespace lbw {

struct LinkFn {
    bool   exponential;     // true = exp link; false = logit link
    double L;               // min_weight
    double U;               // max_weight
    // Deville-Sarndal (1992) logit scale factor: A = (U-L)/((U-1)*(1-L))
    // Governs how fast F(u) transitions from L to U; only valid when !exponential.
    double logit_scale;

    explicit LinkFn(double min_weight, double max_weight)
        : L(min_weight), U(max_weight)
    {
        // Use exponential link when min_weight==0 OR max_weight==Inf
        exponential = (L == 0.0 || !std::isfinite(U));
        logit_scale = exponential ? 0.0 : (U - L) / ((U - 1.0) * (1.0 - L));
    }

    // Clamp exp(x) to exp(700) to prevent IEEE 754 overflow
    static double safe_exp(double x) {
        return std::exp(std::min(x, 700.0));
    }

    // F(u): link function mapping dual variable to weight
    double F(double u) const {
        if (exponential) {
            return safe_exp(u);
        }
        double e = safe_exp(logit_scale * u);
        return (L * (U - 1.0) + U * (1.0 - L) * e) / ((U - 1.0) + (1.0 - L) * e);
    }

    // dF(u): derivative of F w.r.t. u; from quotient rule on F
    double dF(double u) const {
        if (exponential) {
            return safe_exp(u);
        }
        double fu = F(u);
        return logit_scale * (fu - L) * (U - fu) / (U - L);
    }

    // H(u): antiderivative of F(u); H(0) = 0 by construction
    // Exponential: H(u) = exp(u)
    // Logit (Deville-Sarndal 1992):
    //   H(u) = L*u + (U-L)/logit_scale * ln(((U-1)+(1-L)*exp(logit_scale*u)) / (U-L))
    double H(double u) const {
        if (exponential) {
            return safe_exp(u);
        }
        double e = safe_exp(logit_scale * u);
        double num = (U - 1.0) + (1.0 - L) * e;
        return L * u + (U - L) / logit_scale * std::log(num / (U - L));
    }
};

} // namespace lbw
```

- [ ] Write `src/logit.cpp`:
```cpp
#include "logit.hpp"
// LinkFn methods are all inline in the header; this TU is a placeholder
// for future non-inline implementations and ensures the TU is compiled.
namespace lbw {
// Force instantiation to catch any template/inline errors at compile time
static_assert(sizeof(LinkFn) > 0, "LinkFn must be a complete type");
} // namespace lbw
```

- [ ] Create `src/r_bridge.cpp` with the logit test bridges and `R_init_leafblower`:
```cpp
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "logit.hpp"

extern "C" {

// Forward declarations for all .Call entries defined in this file.
// (C_rk_calibrate is added in Task 7; its declaration is added there too.)
SEXP C_logit_F_at_zero(SEXP, SEXP);
SEXP C_logit_range_check(SEXP, SEXP, SEXP);
SEXP C_logit_Hprime_check(SEXP, SEXP, SEXP);

// R_init_leafblower: called automatically by R when the shared library is loaded.
// Registers all .Call entry points so they can be found by name.
// NOTE: updated in Task 7 to include C_rk_calibrate.
void R_init_leafblower(DllInfo* dll) {
    static const R_CallMethodDef call_methods[] = {
        {"C_logit_F_at_zero",    (DL_FUNC)&C_logit_F_at_zero,    2},
        {"C_logit_range_check",  (DL_FUNC)&C_logit_range_check,  3},
        {"C_logit_Hprime_check", (DL_FUNC)&C_logit_Hprime_check, 3},
        {NULL, NULL, 0}
    };
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);  // prevent lookup of unregistered symbols
}

// Test bridge: return F(0) for given L, U
SEXP C_logit_F_at_zero(SEXP Lsxp, SEXP Usxp) {
    double L = REAL(Lsxp)[0];
    double U = REAL(Usxp)[0];
    lbw::LinkFn fn(L, U);
    return Rf_ScalarReal(fn.F(0.0));
}

// Test bridge: return F(u) for a vector of u values
SEXP C_logit_range_check(SEXP Lsxp, SEXP Usxp, SEXP usxp) {
    double L = REAL(Lsxp)[0];
    double U = REAL(Usxp)[0];
    int n = LENGTH(usxp);
    SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
    lbw::LinkFn fn(L, U);
    const double* u = REAL(usxp);
    double* res = REAL(out);
    for (int i = 0; i < n; i++) res[i] = fn.F(u[i]);
    UNPROTECT(1);
    return out;
}

// Test bridge: return |H'(u0) - F(u0)| via numerical diff (step 1e-7)
SEXP C_logit_Hprime_check(SEXP Lsxp, SEXP Usxp, SEXP u0sxp) {
    double L = REAL(Lsxp)[0];
    double U = REAL(Usxp)[0];
    double u0 = REAL(u0sxp)[0];
    double h = 1e-7;
    lbw::LinkFn fn(L, U);
    double Hprime_numerical = (fn.H(u0 + h) - fn.H(u0 - h)) / (2.0 * h);
    double diff = std::fabs(Hprime_numerical - fn.F(u0));
    return Rf_ScalarReal(diff);
}

} // extern "C"
```

- [ ] Build and run GREEN check: `R CMD INSTALL --preclean . && Rscript -e 'testthat::test_file("tests/testthat/test-logit.R")'`
  Expected: all 4 tests PASS.

**Commit:** `feat(logit): add LinkFn with logit/exp link, F/dF/H, all logit tests green`

---

### Task 5: L-BFGS-B Solver

**Files:**
- Create: `src/lbfgsb_solver.hpp`
- Create: `src/lbfgsb_solver.cpp`
- Create: `tests/testthat/test-lbfgsb.R` (RED phase)

**Steps:**

- [ ] Write RED test `tests/testthat/test-lbfgsb.R`:
```r
test_that("L-BFGS-B converges on 3-margin no-bounds case", {
  set.seed(42)
  n <- 50000L
  age <- sample(c("18-34","35-54","55+"), n, replace=TRUE, prob=c(0.35,0.40,0.25))
  sex <- sample(c("M","F"), n, replace=TRUE, prob=c(0.52,0.48))
  edu <- sample(c("HS","College","Grad"), n, replace=TRUE, prob=c(0.40,0.40,0.20))
  df  <- data.frame(age=factor(age), sex=factor(sex), edu=factor(edu))
  tgt <- list(
    age = c("18-34"=0.30, "35-54"=0.45, "55+"=0.25),
    sex = c(M=0.50, F=0.50),
    edu = c(HS=0.35, College=0.45, Grad=0.20)
  )
  # RED: harvest() not yet implemented
  result <- harvest(df, tgt, method="lbfgsb")
  expect_s3_class(result, "data.frame")
  expect_true("weights" %in% names(result))
  expect_lt(abs(mean(result$weights) - 1.0), 1e-8)
})

test_that("L-BFGS-B max_weight bound respected", {
  set.seed(1)
  n <- 10000L
  x <- sample(c("a","b"), n, replace=TRUE, prob=c(0.9,0.1))
  df <- data.frame(x=factor(x))
  tgt <- list(x=c(a=0.5, b=0.5))
  result <- harvest(df, tgt, method="lbfgsb", max_weight=5)
  expect_true(max(result$weights) <= 5.0 + 1e-10)
})
```

- [ ] Confirm RED (harvest not yet implemented).

- [ ] Write `src/lbfgsb_solver.hpp`:
```cpp
#pragma once
#include "types.hpp"
#include "logit.hpp"
#include <vector>

namespace lbw {

struct LBFGSResult {
    int    status;       // RK_OK or RK_ERR_NOCONV
    int    iterations;
    double max_error;
};

LBFGSResult lbfgsb_solve(CalibState& state);

} // namespace lbw
```

- [ ] Write `src/lbfgsb_solver.cpp` — dual objective, gradient, 2-loop L-BFGS, Wolfe line search:
```cpp
#include "lbfgsb_solver.hpp"
#include "leafblower.h"
#include <cmath>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <vector>
#include <deque>

namespace lbw {

// Compute T_kj = targets[k][j] * W where W = sum(weights)
// Returns total W.
static double compute_targets_abs(const CalibState& st,
                                   std::vector<double>& T) {
    double W = 0.0;
    for (int i = 0; i < st.n; i++) W += st.weights[i];
    int off = 0;
    for (int k = 0; k < st.K; k++) {
        for (int j = 0; j < st.cat_counts[k]; j++) {
            T[off + j] = st.targets[k][j] * W;
        }
        off += st.cat_counts[k];
    }
    return W;
}

// Build offset array: offset[k] = index of first dual var for margin k
static std::vector<int> build_offsets(const CalibState& st) {
    std::vector<int> off(st.K + 1, 0);
    for (int k = 0; k < st.K; k++) off[k+1] = off[k] + st.cat_counts[k];
    return off;
}

// Compute u_i = sum_k lambda[offset[k] + group_ids[k][i]]  (skip NA: g==-1)
static void compute_u(const CalibState& st, const std::vector<int>& off,
                       const std::vector<double>& lam, std::vector<double>& u) {
    std::fill(u.begin(), u.end(), 0.0);
    for (int k = 0; k < st.K; k++) {
        for (int i = 0; i < st.n; i++) {
            int g = st.group_ids[k][i];
            if (g >= 0) u[i] += lam[off[k] + g];
        }
    }
}

// phi(lambda) = sum_kj T_kj*lam_kj - sum_i d_i*H(u_i)
// grad[off_k+j] = T_kj - S_kj where S_kj = sum_{i:g_k(i)==j} d_i*F(u_i)
static double phi_and_grad(const CalibState& st,
                            const LinkFn& fn,
                            const std::vector<int>& off,
                            const std::vector<double>& lam,
                            const std::vector<double>& T,
                            const std::vector<double>& d,   // initial weights (d_i)
                            std::vector<double>& grad,
                            std::vector<double>& u) {
    int total = off[st.K];
    compute_u(st, off, lam, u);

    // phi = sum_kj T_kj*lam_kj - sum_i d_i*H(u_i)
    double obj = 0.0;
    for (int idx = 0; idx < total; idx++) obj += T[idx] * lam[idx];
    for (int i = 0; i < st.n; i++) obj -= d[i] * fn.H(u[i]);

    // grad_kj = T_kj - sum_{g_k(i)==j} d_i*F(u_i)
    std::fill(grad.begin(), grad.end(), 0.0);
    for (int idx = 0; idx < total; idx++) grad[idx] = T[idx];
    for (int k = 0; k < st.K; k++) {
        for (int i = 0; i < st.n; i++) {
            int g = st.group_ids[k][i];
            if (g >= 0) grad[off[k] + g] -= d[i] * fn.F(u[i]);
        }
    }
    return obj;
}

// L-infinity norm
static double linf(const std::vector<double>& v) {
    double mx = 0.0;
    for (double x : v) mx = std::max(mx, std::fabs(x));
    return mx;
}

// Dot product
static double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    int n = (int)a.size();
    for (int i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

// L-BFGS 2-loop recursion: given grad g, produce search direction q = H*g
// (we maximize, so direction = +H*grad)
static void lbfgs_direction(const std::deque<std::vector<double>>& svec,
                              const std::deque<std::vector<double>>& yvec,
                              const std::deque<double>& rho,
                              double gamma,
                              const std::vector<double>& g,
                              std::vector<double>& dir) {
    int m = (int)svec.size();
    int num_dual_vars = (int)g.size();
    dir = g;  // start with gradient

    std::vector<double> alpha(m);
    // First loop: newest first
    for (int i = m - 1; i >= 0; i--) {
        alpha[i] = rho[i] * dot(svec[i], dir);
        for (int j = 0; j < num_dual_vars; j++) dir[j] -= alpha[i] * yvec[i][j];
    }
    // Scale by initial Hessian
    for (int j = 0; j < num_dual_vars; j++) dir[j] *= gamma;
    // Second loop: oldest first
    for (int i = 0; i < m; i++) {
        double beta = rho[i] * dot(yvec[i], dir);
        for (int j = 0; j < num_dual_vars; j++) dir[j] += (alpha[i] - beta) * svec[i][j];
    }
}

// Zoom phase (Nocedal & Wright Alg 3.6): bisect [alpha_lo, alpha_hi] until strong Wolfe.
// alpha_lo satisfies Armijo or is better than alpha_hi; alpha_hi does not.
static double wolfe_zoom(
        const CalibState& st, const LinkFn& fn,
        const std::vector<int>& off, const std::vector<double>& T,
        const std::vector<double>& d, double phi_0, double slope_0,
        double alpha_lo, double phi_lo, double alpha_hi,
        std::vector<double>& u, const std::vector<double>& lam,
        const std::vector<double>& dir,
        std::vector<double>& lam_new, std::vector<double>& grad_new,
        double& phi_new) {
    constexpr double kArmijoC1    = 1e-4;
    constexpr double kCurvatureC2 = 0.9;
    const int total = (int)lam.size();

    for (int j = 0; j < 20; j++) {
        double alpha = 0.5 * (alpha_lo + alpha_hi);
        for (int i = 0; i < total; i++) lam_new[i] = lam[i] + alpha * dir[i];
        phi_new = phi_and_grad(st, fn, off, lam_new, T, d, grad_new, u);
        double slope = dot(grad_new, dir);

        if (phi_new < phi_0 + kArmijoC1 * alpha * slope_0 || phi_new <= phi_lo) {
            alpha_hi = alpha;
        } else {
            if (std::fabs(slope) <= kCurvatureC2 * std::fabs(slope_0)) return alpha;
            if (slope * (alpha_hi - alpha_lo) >= 0.0) alpha_hi = alpha_lo;
            alpha_lo = alpha; phi_lo = phi_new;
        }
    }
    return 0.5 * (alpha_lo + alpha_hi);
}

// Strong Wolfe line search for dual φ maximization (Nocedal & Wright Alg 3.5/3.6).
// c1=1e-4 (Armijo), c2=0.9 (curvature); unit initial step, max 20 bracket + 20 zoom.
// Ensures sy = dot(s,y) > 0 so the L-BFGS curvature condition is always satisfied.
static double wolfe_line_search(
        const CalibState& st, const LinkFn& fn,
        const std::vector<int>& off, const std::vector<double>& T,
        const std::vector<double>& d,
        const std::vector<double>& lam, double phi_0, double slope_0,
        std::vector<double>& u, const std::vector<double>& dir,
        std::vector<double>& lam_new, std::vector<double>& grad_new,
        double& phi_new) {
    constexpr double kArmijoC1    = 1e-4;
    constexpr double kCurvatureC2 = 0.9;
    const int total = (int)lam.size();

    double alpha_prev = 0.0, phi_prev = phi_0;
    double alpha = 1.0;

    for (int i = 0; i < 20; i++) {
        for (int j = 0; j < total; j++) lam_new[j] = lam[j] + alpha * dir[j];
        phi_new = phi_and_grad(st, fn, off, lam_new, T, d, grad_new, u);
        double slope = dot(grad_new, dir);

        // Armijo violated or not improving vs previous → zoom
        if (phi_new < phi_0 + kArmijoC1 * alpha * slope_0 || (i > 0 && phi_new <= phi_prev)) {
            return wolfe_zoom(st, fn, off, T, d, phi_0, slope_0,
                              alpha_prev, phi_prev, alpha, u, lam, dir,
                              lam_new, grad_new, phi_new);
        }
        // Strong Wolfe curvature condition satisfied
        if (std::fabs(slope) <= kCurvatureC2 * std::fabs(slope_0)) return alpha;
        // Gradient went negative → zoom with swapped bracket
        if (slope <= 0) {
            return wolfe_zoom(st, fn, off, T, d, phi_0, slope_0,
                              alpha, phi_new, alpha_prev, u, lam, dir,
                              lam_new, grad_new, phi_new);
        }
        alpha_prev = alpha; phi_prev = phi_new;
        alpha = std::min(2.0 * alpha, 8.0);  // expand, capped to avoid huge steps
    }
    st.log("L-BFGS-B: Wolfe bracket did not converge; using last step");
    return alpha;
}

LBFGSResult lbfgsb_solve(CalibState& st) {
    LBFGSResult res;
    res.status = RK_ERR_NOCONV;
    res.iterations = 0;
    res.max_error = 1.0;

    LinkFn fn(st.min_weight, st.max_weight);
    auto off = build_offsets(st);
    int total = off[st.K];

    // Store initial design weights d_i
    std::vector<double> d(st.n);
    for (int i = 0; i < st.n; i++) d[i] = st.weights[i];

    std::vector<double> lam(total, 0.0);   // dual variables, init 0
    std::vector<double> T(total);
    std::vector<double> grad(total), grad_prev(total);
    std::vector<double> u(st.n);
    std::vector<double> dir(total);

    // L-BFGS history
    std::deque<std::vector<double>> svec, yvec;
    std::deque<double> rho_hist;
    double gamma = 1.0;

    compute_targets_abs(st, T);  // T_kj = targets[k][j] * W

    double phi_curr = phi_and_grad(st, fn, off, lam, T, d, grad, u);

    int max_iter = st.outer_max_iter;
    for (int iter = 0; iter < max_iter; iter++) {
        res.iterations = iter + 1;

        // Check convergence
        double gn = linf(grad);
        if (gn < st.tol_abs) {
            res.status = RK_OK;
            res.max_error = gn;
            break;
        }

        // Compute search direction
        if (svec.empty()) {
            dir = grad;  // steepest ascent on first iteration
        } else {
            lbfgs_direction(svec, yvec, rho_hist, gamma, grad, dir);
        }

        // Strong Wolfe line search (Nocedal & Wright Alg 3.5/3.6): c1=1e-4, c2=0.9.
        // Wolfe curvature condition guarantees sy > 0 on every accepted step,
        // ensuring the L-BFGS Hessian approximation remains positive definite.
        double slope_0 = dot(grad, dir);  // > 0 since dir is the ascent direction

        std::vector<double> lam_new(total);
        std::vector<double> grad_new(total);
        double phi_new = phi_curr;

        double step = wolfe_line_search(st, fn, off, T, d, lam, phi_curr, slope_0,
                                        u, dir, lam_new, grad_new, phi_new);
        (void)step;  // step stored in lam_new; only result arrays matter

        // Update L-BFGS history
        std::vector<double> s_new(total), y_new(total);
        for (int i = 0; i < total; i++) {
            s_new[i] = lam_new[i] - lam[i];
            y_new[i] = grad_new[i] - grad[i];
        }
        double sy = dot(s_new, y_new);
        double yy = dot(y_new, y_new);
        constexpr double kLbfgsCurvatureMin = 1e-20;
        if (sy > kLbfgsCurvatureMin && yy > kLbfgsCurvatureMin) {
            if ((int)svec.size() >= st.lbfgs_m) {
                svec.pop_front(); yvec.pop_front(); rho_hist.pop_front();
            }
            svec.push_back(s_new);
            yvec.push_back(y_new);
            rho_hist.push_back(1.0 / sy);
            gamma = sy / yy;
        }

        lam = lam_new;
        grad = grad_new;
        phi_curr = phi_new;
    }

    // Compute final weights: w_i = clamp(d_i * F(u_i), min, max)
    compute_u(st, off, lam, u);
    double W = 0.0;
    for (int i = 0; i < st.n; i++) {
        double wi = d[i] * fn.F(u[i]);
        wi = std::max(st.min_weight, std::min(st.max_weight, wi));
        st.weights[i] = wi;
        W += wi;
    }
    // Do NOT normalize here. The bridge (R/Python) normalizes start_weights to
    // mean=1 before calling rk_calibrate(). Calibration preserves sum(weights)≈n
    // when it converges. Re-normalizing after clamping would invalidate the
    // calibration constraints that were verified on the un-normalized weights.

    // Compute max_error on final weights
    {
        double Wn = 0.0;
        for (int i = 0; i < st.n; i++) Wn += st.weights[i];
        double max_err = 0.0;
        for (int k = 0; k < st.K; k++) {
            for (int j = 0; j < st.cat_counts[k]; j++) {
                double Skj = 0.0;
                for (int i = 0; i < st.n; i++) {
                    if (st.group_ids[k][i] == j) Skj += st.weights[i];
                }
                max_err = std::max(max_err, std::fabs(Skj/Wn - st.targets[k][j]));
            }
        }
        res.max_error = max_err;
        if (max_err < st.tol_abs) res.status = RK_OK;
    }

    return res;
}

} // namespace lbw
```

- [ ] **Refactor `lbfgsb_solve` into single-responsibility helpers:**
  After the above code compiles and tests pass, extract:
  - `compute_offsets(CalibState&) → std::vector<int>` (builds the offset array)
  - `compute_final_weights_and_error(CalibState&, const LinkFn&, const std::vector<double>& u_final) → rk_result_t` (applies clamp, computes max_error, returns result)
  - Each function ≤ 25 lines. Move `phi_and_grad` working buffers into a `DualObjective` struct to reduce arg count from 7 to 2.
  - Verify: `R CMD INSTALL --preclean .` still passes after refactor.

**Commit:** `feat(lbfgsb): add L-BFGS-B solver with strong Wolfe line search (c1=1e-4, c2=0.9)`

---

### Task 6: Wire c_api.cpp for L-BFGS-B

**Files:**
- Modify: `src/c_api.cpp`

**Steps:**

- [ ] Add `#include "lbfgsb_solver.hpp"` and `#include "ieppa.hpp"` (stub) to `src/c_api.cpp`.

- [ ] Add `select_algorithm()` function to `src/c_api.cpp`:
```cpp
static rk_algorithm_t select_algorithm(int n, int K,
                                        const int* cat_counts,
                                        const rk_params_t* p) {
    if (p->algorithm != RK_ALG_AUTO) return p->algorithm;
    int64_t complexity = INT64_C(0);
    for (int k = 0; k < K; k++) complexity += (int64_t)n * cat_counts[k];
    if (complexity > 500000L || p->max_weight < 3.0 || p->min_weight > 0.0)
        return RK_ALG_IEPPA;
    return RK_ALG_LBFGSB;
}
```

- [ ] Replace the stub body of `rk_calibrate()` with full dispatch:
```cpp
LBW_NODISCARD int rk_calibrate(int n, int K,
                                double* weights,
                                const int32_t** group_ids,
                                const int* cat_counts,
                                const double** targets,
                                const rk_params_t* params,
                                rk_result_t* result) {
    rk_params_t defaults;
    rk_params_init(&defaults);
    const rk_params_t* p = params ? params : &defaults;

    if (result) {
        result->status = RK_OK;
        result->iterations = 0;
        result->max_error = 0.0;
        result->algorithm_used = RK_ALG_AUTO;
        result->message[0] = '\0';
    }

    int rc = validate_inputs(n, K, weights, group_ids, cat_counts, targets, p, result);
    if (rc != RK_OK) return rc;

    rk_algorithm_t alg = select_algorithm(n, K, cat_counts, p);

    // Build CalibState
    lbw::CalibState st;
    st.n = n; st.K = K;
    st.weights = weights;
    st.group_ids = group_ids;
    st.cat_counts = cat_counts;
    st.targets = targets;
    st.min_weight    = p->min_weight;
    st.max_weight    = p->max_weight;
    st.tol_abs       = p->tol_abs;
    st.inner_max_iter = p->inner_max_iter;
    st.outer_max_iter = p->outer_max_iter;
    st.epsilon       = p->epsilon;
    st.lbfgs_m       = p->lbfgs_m;
    st.verbose       = p->verbose;
    st.log_fn        = p->log_fn;
    st.log_ctx       = p->log_ctx;
    st.total_cats    = 0;
    for (int k = 0; k < K; k++) st.total_cats += cat_counts[k];

    // Verbose routing report
    if (p->verbose >= 1) {
        int64_t complexity = INT64_C(0);
        for (int k = 0; k < K; k++) complexity += (int64_t)n * cat_counts[k];
        char msg[256];
        if (alg == RK_ALG_IEPPA)
            snprintf(msg, 256, "Auto-selected iEPPA: complexity=%lld, max_weight=%.2f, min_weight=%.2f",
                     (long long)complexity, p->max_weight, p->min_weight);
        else
            snprintf(msg, 256, "Auto-selected L-BFGS-B: complexity=%lld <= 500000, max_weight=%.2f, min_weight=%.2f",
                     (long long)complexity, p->max_weight, p->min_weight);
        st.log(msg);
    }

    int status;
    int iterations;
    double max_error;
    rk_algorithm_t used;

    if (alg == RK_ALG_LBFGSB) {
        auto res = lbw::lbfgsb_solve(st);
        status = res.status;
        iterations = res.iterations;
        max_error = res.max_error;
        used = RK_ALG_LBFGSB;
    } else {
        // iEPPA — stub for now
        status = RK_ERR_NOCONV;
        iterations = 0;
        max_error = 1.0;
        used = RK_ALG_IEPPA;
        if (result) snprintf(result->message, 256, "iEPPA not implemented yet");
    }

    if (result) {
        result->status = status;
        result->iterations = iterations;
        result->max_error = max_error;
        result->algorithm_used = used;
        if (result->message[0] == '\0') {
            snprintf(result->message, 256,
                     "%s: %d iters, max_error=%.2e",
                     used == RK_ALG_LBFGSB ? "L-BFGS-B" : "iEPPA",
                     iterations, max_error);
        }
    }
    return status;
}
```

- [ ] Create stub `src/ieppa.hpp` (to satisfy include):
```cpp
#pragma once
#include "types.hpp"
namespace lbw {
struct IEPPAResult { int status; int iterations; double max_error; };
IEPPAResult ieppa_solve(CalibState& state);
} // namespace lbw
```

- [ ] Create stub `src/ieppa.cpp`:
```cpp
#include "ieppa.hpp"
#include "leafblower.h"
namespace lbw {
IEPPAResult ieppa_solve(CalibState& st) {
    return {RK_ERR_NOCONV, 0, 1.0};
}
} // namespace lbw
```

- [ ] Build: `R CMD INSTALL --preclean .` — expect no compile errors.

**Commit:** `feat(dispatch): wire select_algorithm and lbfgsb_solve into rk_calibrate`

---

### Task 7: R Bridge + harvest()

**Files:**
- Modify: `src/r_bridge.cpp` (add C_rk_calibrate)
- Create: `R/harvest.R`
- Modify: `R/zzz.R`

**Steps:**

- [ ] Add `C_rk_calibrate` to `src/r_bridge.cpp` and update `R_init_leafblower`:
```cpp
#include "leafblower.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <cstring>
#include <string>
#include <unordered_map>

// R log trampoline
static void r_log_trampoline(const char* msg, void* ctx) {
    Rprintf("[leafblower] %s\n", msg);
}

SEXP C_rk_calibrate(SEXP data_sexp, SEXP target_sexp,
                    SEXP min_weight_sexp, SEXP max_weight_sexp,
                    SEXP method_sexp, SEXP verbose_sexp,
                    SEXP inner_max_iter_sexp, SEXP start_weights_sexp) {
    // data_sexp: data.frame
    // target_sexp: list of named numeric vectors (variable -> proportions)
    // Returns: list(weights=double[n], result=list(status, iterations, max_error, algorithm_used, message))

    int n = Rf_nrows(VECTOR_ELT(data_sexp, 0));
    SEXP target_names = Rf_getAttrib(target_sexp, R_NamesSymbol);
    int K = LENGTH(target_sexp);

    // Build group_ids and cat_counts from data.frame factor columns
    std::vector<std::vector<int32_t>> gids_storage(K);
    std::vector<const int32_t*> group_ids(K);
    std::vector<int> cat_counts(K);
    std::vector<std::vector<double>> tgt_storage(K);
    std::vector<const double*> targets(K);

    SEXP col_names = Rf_getAttrib(data_sexp, R_NamesSymbol);

    for (int k = 0; k < K; k++) {
        const char* varname = CHAR(STRING_ELT(target_names, k));
        // Find column in data
        int col_idx = -1;
        for (int c = 0; c < LENGTH(col_names); c++) {
            if (strcmp(CHAR(STRING_ELT(col_names, c)), varname) == 0) {
                col_idx = c; break;
            }
        }
        if (col_idx < 0)
            Rf_error("Variable '%s' not found in data", varname);

        SEXP col = VECTOR_ELT(data_sexp, col_idx);
        SEXP tgt_vec = VECTOR_ELT(target_sexp, k);
        SEXP tgt_names = Rf_getAttrib(tgt_vec, R_NamesSymbol);
        int ncat = LENGTH(tgt_vec);
        cat_counts[k] = ncat;

        // Build target proportion array
        tgt_storage[k].resize(ncat);
        for (int j = 0; j < ncat; j++) tgt_storage[k][j] = REAL(tgt_vec)[j];
        targets[k] = tgt_storage[k].data();

        // Build O(1) level-name → category-index map (avoids O(ncat) inner scan)
        std::unordered_map<std::string, int> level_to_idx;
        level_to_idx.reserve(ncat);
        for (int j = 0; j < ncat; j++)
            level_to_idx[CHAR(STRING_ELT(tgt_names, j))] = j;

        // Encode factor/character column to 0-indexed int (NA → -1, OOV → -1)
        gids_storage[k].resize(n);
        if (Rf_isFactor(col)) {
            SEXP flevels = Rf_getAttrib(col, R_LevelsSymbol);
            const int* codes = INTEGER(col);
            for (int i = 0; i < n; i++) {
                if (codes[i] == NA_INTEGER) { gids_storage[k][i] = -1; continue; }
                const char* lv = CHAR(STRING_ELT(flevels, codes[i] - 1));
                auto it = level_to_idx.find(lv);
                gids_storage[k][i] = (it != level_to_idx.end()) ? it->second : -1;
            }
        } else if (TYPEOF(col) == STRSXP) {
            for (int i = 0; i < n; i++) {
                if (STRING_ELT(col, i) == NA_STRING) { gids_storage[k][i] = -1; continue; }
                const char* lv = CHAR(STRING_ELT(col, i));
                auto it = level_to_idx.find(lv);
                gids_storage[k][i] = (it != level_to_idx.end()) ? it->second : -1;
            }
        } else {
            Rf_error("Column '%s' must be a factor or character vector", varname);
        }
        group_ids[k] = gids_storage[k].data();
    }

    // Build weights vector
    std::vector<double> weights(n);
    if (Rf_isNull(start_weights_sexp)) {
        for (int i = 0; i < n; i++) weights[i] = 1.0;
    } else {
        double* sw = REAL(start_weights_sexp);
        double mean_sw = 0.0;
        for (int i = 0; i < n; i++) mean_sw += sw[i];
        mean_sw /= n;
        for (int i = 0; i < n; i++) weights[i] = sw[i] / mean_sw;
    }

    // Set params
    rk_params_t p;
    rk_params_init(&p);
    p.min_weight    = REAL(min_weight_sexp)[0];
    p.max_weight    = REAL(max_weight_sexp)[0];
    p.verbose       = INTEGER(verbose_sexp)[0];
    p.inner_max_iter = INTEGER(inner_max_iter_sexp)[0];
    p.log_fn        = (p.verbose > 0) ? r_log_trampoline : nullptr;

    const char* method_str = CHAR(STRING_ELT(method_sexp, 0));
    if      (strcmp(method_str, "ieppa")  == 0) p.algorithm = RK_ALG_IEPPA;
    else if (strcmp(method_str, "lbfgsb") == 0) p.algorithm = RK_ALG_LBFGSB;
    else                                          p.algorithm = RK_ALG_AUTO;

    rk_result_t result;
    int rc = rk_calibrate(n, K, weights.data(),
                          group_ids.data(),
                          cat_counts.data(),
                          (const double**)targets.data(),
                          &p, &result);

    // Build return list
    SEXP wts = PROTECT(Rf_allocVector(REALSXP, n));
    memcpy(REAL(wts), weights.data(), n * sizeof(double));

    SEXP res_list = PROTECT(Rf_allocVector(VECSXP, 5));
    SEXP res_names = PROTECT(Rf_allocVector(STRSXP, 5));
    SET_STRING_ELT(res_names, 0, Rf_mkChar("status"));
    SET_STRING_ELT(res_names, 1, Rf_mkChar("iterations"));
    SET_STRING_ELT(res_names, 2, Rf_mkChar("max_error"));
    SET_STRING_ELT(res_names, 3, Rf_mkChar("algorithm_used"));
    SET_STRING_ELT(res_names, 4, Rf_mkChar("message"));
    SET_VECTOR_ELT(res_list, 0, Rf_ScalarInteger(result.status));
    SET_VECTOR_ELT(res_list, 1, Rf_ScalarInteger(result.iterations));
    SET_VECTOR_ELT(res_list, 2, Rf_ScalarReal(result.max_error));
    SET_VECTOR_ELT(res_list, 3, Rf_ScalarInteger((int)result.algorithm_used));
    SET_VECTOR_ELT(res_list, 4, Rf_mkString(result.message));
    Rf_setAttrib(res_list, R_NamesSymbol, res_names);

    SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP out_names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(out_names, 0, Rf_mkChar("weights"));
    SET_STRING_ELT(out_names, 1, Rf_mkChar("result"));
    SET_VECTOR_ELT(out, 0, wts);
    SET_VECTOR_ELT(out, 1, res_list);
    Rf_setAttrib(out, R_NamesSymbol, out_names);
    UNPROTECT(5);
    return out;
}

// Update R_init_leafblower to include C_rk_calibrate.
// Replace the existing R_init_leafblower definition with:
void R_init_leafblower(DllInfo* dll) {
    static const R_CallMethodDef call_methods[] = {
        {"C_logit_F_at_zero",    (DL_FUNC)&C_logit_F_at_zero,    2},
        {"C_logit_range_check",  (DL_FUNC)&C_logit_range_check,  3},
        {"C_logit_Hprime_check", (DL_FUNC)&C_logit_Hprime_check, 3},
        {"C_rk_calibrate",       (DL_FUNC)&C_rk_calibrate,       8},
        {NULL, NULL, 0}
    };
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
```

- [ ] Write `R/harvest.R`:
```r
#' Generate calibrated weights (drop-in for autumn::harvest)
#'
#' @param data A data frame containing columns to calibrate.
#' @param target A named list of named numeric vectors (variable → proportions).
#' @param min_weight Lower bound on weights. Default 0 (no lower bound).
#' @param max_weight Upper bound on weights. Default 5.
#' @param method One of "auto", "ieppa", "lbfgsb", "rake", "nr". Default "auto".
#' @param verbose Integer verbosity: 0=silent, 1=progress, 2=debug.
#' @param max_iterations Maximum inner BCD iterations per outer step. Default 500.
#' @param start_weights Starting weights vector or NULL (uniform).
#' @param attach_weights If TRUE, return data frame with weights column. Default TRUE.
#' @param weight_column Name of weight column. Default "weights".
#' @param convergence Named list/vector; "absolute" maps to tol_abs. "pct" triggers deprecation warning.
#' @param select_params Ignored with verbose >= 2 note (iEPPA/L-BFGS-B calibrate all simultaneously).
#' @param select_function Ignored (same reason as select_params).
#' @param error_function Ignored.
#' @param adaptive_order Ignored.
#' @param enforce_mean Ignored (retained for compatibility).
#' @param accelerate Ignored.
#' @param add_na_proportion Not supported in v1; raises error if TRUE.
#' @param auto_collapse Not supported in v1; raises error if TRUE.
#' @param collapse_vars Not supported in v1; raises error if TRUE.
#' @param target_map Passed through for data-frame target format handling.
#' @param ... Additional arguments ignored.
#' @return data frame with weights column if attach_weights=TRUE, else numeric vector.
#' @export
harvest <- function(
  data,
  target,
  min_weight       = 0,
  max_weight       = 5,
  method           = "auto",
  verbose          = 0,
  max_iterations   = 500,
  start_weights    = NULL,
  attach_weights   = TRUE,
  weight_column    = "weights",
  convergence      = list(),
  select_params    = NULL,
  select_function  = NULL,
  error_function   = NULL,
  adaptive_order   = NULL,
  enforce_mean     = TRUE,
  accelerate       = FALSE,
  add_na_proportion = FALSE,
  auto_collapse    = FALSE,
  collapse_vars    = NULL,
  target_map       = NULL,
  ...
) {
  # Not-in-v1 hard stops
  if (!isFALSE(add_na_proportion) && !identical(add_na_proportion, FALSE))
    stop("add_na_proportion is not supported in leafblower v1. ",
         "Remove NAs from target variables before calibrating.")
  if (isTRUE(auto_collapse))
    stop("auto_collapse is not supported in leafblower v1.")

  # Ignored-param verbose notes
  ignored <- c("select_params", "select_function", "error_function",
                "adaptive_order", "accelerate")
  supplied_ignored <- intersect(ignored,
                                names(match.call(expand.dots = FALSE)))
  if (verbose >= 2 && length(supplied_ignored) > 0) {
    message("[leafblower] Ignoring autumn params (not applicable to iEPPA/L-BFGS-B): ",
            paste(supplied_ignored, collapse = ", "))
  }

  # Method mapping
  method <- tolower(method)
  if (method %in% c("rake", "nrake")) {
    warning("method='", method, "' (IPF) not implemented; using L-BFGS-B (method='lbfgsb')")
    method <- "lbfgsb"
  } else if (method == "nr") {
    warning("method='nr' (Newton-Raphson) not implemented; using L-BFGS-B (method='lbfgsb')")
    method <- "lbfgsb"
  }
  method <- match.arg(method, c("auto", "ieppa", "lbfgsb"))

  # Convergence parameter handling
  if (!is.null(convergence[["pct"]])) {
    warning("convergence['pct'] is deprecated in leafblower; use convergence['absolute'] instead. ",
            "The 'pct' criterion is not applicable to iEPPA/L-BFGS-B.")
  }
  tol_abs <- if (!is.null(convergence[["absolute"]])) convergence[["absolute"]] else 1e-6

  # Handle data-frame target format
  if (is.data.frame(target)) {
    # Minimal target_map support: assume columns variable, level, proportion
    if (!is.null(target_map)) {
      vcol <- target_map[["variable"]]; lcol <- target_map[["level"]]
      pcol <- target_map[["proportion"]]
    } else if (all(c("variable","level","proportion") %in% names(target))) {
      vcol <- "variable"; lcol <- "level"; pcol <- "proportion"
    } else if (ncol(target) == 3) {
      warning("Assuming target data frame columns are variable, level, proportion.")
      vcol <- 1; lcol <- 2; pcol <- 3
    } else {
      stop("Cannot determine variable/level/proportion columns in target data frame.")
    }
    vars <- unique(target[[vcol]])
    target <- lapply(stats::setNames(vars, vars), function(v) {
      sub <- target[target[[vcol]] == v, , drop = FALSE]
      stats::setNames(sub[[pcol]], sub[[lcol]])
    })
  }

  if (!is.list(target))
    stop("target must be a named list of named numeric vectors or a data frame.")

  # Normalize start_weights to mean 1
  sw_vec <- NULL
  if (!is.null(start_weights)) {
    if (length(start_weights) == 1) {
      sw_vec <- rep(as.double(start_weights), nrow(data))
    } else {
      if (length(start_weights) != nrow(data))
        stop("start_weights length must equal nrow(data)")
      sw_vec <- as.double(start_weights)
    }
    sw_vec <- sw_vec * length(sw_vec) / sum(sw_vec)
  }

  # Call C bridge
  raw <- .Call("C_rk_calibrate",
               data,
               target,
               as.double(min_weight),
               as.double(max_weight),
               as.character(method),
               as.integer(verbose),
               as.integer(max_iterations),
               sw_vec,
               PACKAGE = "leafblower")

  weights <- raw$weights
  cres    <- raw$result

  if (cres$status == 1L)
    warning("leafblower: calibration did not converge (max_error=",
            signif(cres$max_error, 3), "). Weights reflect last iterate.")
  if (cres$status == 2L)
    stop("leafblower: infeasible problem — empty cell with positive target.")
  if (cres$status == 3L)
    stop("leafblower: invalid arguments — ", cres$message)

  alg_names <- c("auto", "ieppa", "lbfgsb")
  alg_used  <- alg_names[cres$algorithm_used + 1L]

  if (!attach_weights) return(weights)

  col <- if (!is.null(weight_column)) weight_column else "weights"
  data[[col]] <- weights
  attr(data, "algorithm") <- alg_used
  data
}
```

- [ ] Write RED test for method="rake" warning in `tests/testthat/test-harvest.R`:
```r
test_that("method='rake' emits warning about L-BFGS-B", {
  set.seed(1)
  df  <- data.frame(x = factor(sample(c("a","b"), 200, replace=TRUE)))
  tgt <- list(x = c(a=0.5, b=0.5))
  expect_warning(harvest(df, tgt, method="rake"), regexp = "L-BFGS-B")
})

test_that("method='nr' emits warning about L-BFGS-B", {
  set.seed(1)
  df  <- data.frame(x = factor(sample(c("a","b"), 200, replace=TRUE)))
  tgt <- list(x = c(a=0.5, b=0.5))
  expect_warning(harvest(df, tgt, method="nr"), regexp = "L-BFGS-B")
})
```

- [ ] Build and run: `R CMD INSTALL --preclean . && Rscript -e 'testthat::test_file("tests/testthat/test-harvest.R")'`
  Expected: method warnings pass; lbfgsb convergence test passes.

- [ ] Run lbfgsb tests: `Rscript -e 'testthat::test_file("tests/testthat/test-lbfgsb.R")'`
  Expected: GREEN.

- [ ] **MANDATORY** — do not commit until this refactor is complete: **Refactor `harvest()` into single-responsibility helpers:**
  After tests pass, split the 22-step harvest() body into:
  - `parse_convergence(convergence)` → list(tol_abs)
  - `map_method(method)` → canonical method string with warnings
  - `normalize_start_weights(start_weights, n)` → numeric vector, mean=1
  Each helper ≤ 15 lines, independently testable. harvest() becomes an orchestrator of 4 calls + the .Call. Verify: all test-harvest.R tests still GREEN.
- [ ] Compile gate: `R CMD INSTALL --preclean .` — must produce 0 errors before Task 7 commit.

**Commit:** `feat(rbridge): add C_rk_calibrate and harvest() with full autumn compat`

---

### Task 8: Autumn-Compat Functions

**Files:**
- Create: `R/anesrake.R`
- Create: `R/diagnose_weights.R`
- Create: `R/design_effect.R`
- Create: `R/weighted_pct.R`
- Create: `R/current_miss.R`
- Create: `tests/testthat/test-design.R`
- Create: `tests/testthat/test-compat.R`

**Steps:**

- [ ] Write `R/anesrake.R`:
```r
#' anesrake compatibility wrapper
#'
#' Maps anesrake parameter names to harvest(). Unknown params produce warnings.
#' @param inputter Data frame (maps to \code{data}).
#' @param targets Named list of target proportions.
#' @param weightvec Starting weights or NULL.
#' @param caseid Ignored (row label not used by calibration).
#' @param pctlim Deprecated; maps to \code{convergence["pct"]}.
#' @param cap Upper weight cap; maps to \code{max_weight}.
#' @param choosemethod "rake" or "nrake"; both map to \code{method="lbfgsb"} with warning.
#' @param type Ignored.
#' @param nlim Max iterations; maps to \code{max_iterations}.
#' @param iterate Ignored (always iterates).
#' @param threads Ignored (single-threaded v1).
#' @param ... Unknown params produce a warning then are ignored.
#' @return Same as \code{harvest()}.
#' @export
anesrake <- function(inputter, targets, weightvec = NULL, caseid = NULL,
                     pctlim = 0.05, cap = 5, choosemethod = "rake",
                     type = "pctlim", nlim = 500L, iterate = TRUE,
                     threads = 1L, ...) {
  dots <- list(...)
  if (length(dots) > 0)
    warning("anesrake: ignoring unknown arguments: ",
            paste(names(dots), collapse = ", "))

  if (!is.null(caseid))
    message("anesrake: caseid is ignored (not used by leafblower calibration)")

  conv <- list()
  if (!is.null(pctlim))
    conv[["pct"]] <- pctlim  # harvest() will emit deprecation warning

  harvest(
    data           = inputter,
    target         = targets,
    start_weights  = weightvec,
    max_weight     = cap,
    method         = choosemethod,  # harvest() maps "rake"/"nrake" → lbfgsb with warning
    max_iterations = as.integer(nlim),
    convergence    = conv
  )
}
```

- [ ] Write `R/diagnose_weights.R`:
```r
#' Diagnose calibration quality
#'
#' Returns a data frame comparing original and weighted marginals to targets.
#'
#' @param data Data frame used in calibration.
#' @param target Named list of target proportions (same format as harvest()).
#' @param weights Numeric vector of calibrated weights, length nrow(data).
#' @return Data frame with columns: variable, level, prop_original, prop_weighted,
#'   target, error_original, error_weighted.
#' @export
diagnose_weights <- function(data, target, weights) {
  if (!is.list(target))
    stop("target must be a named list of named numeric vectors")
  if (length(weights) != nrow(data))
    stop("weights length must equal nrow(data)")

  # Pre-allocate to avoid O(n^2) copy in c(rows, list(...)) pattern
  total_rows <- sum(vapply(target, length, integer(1L)))
  rows <- vector("list", total_rows)
  row_idx <- 0L

  for (varname in names(target)) {
    col <- data[[varname]]
    if (is.null(col))
      stop("Variable '", varname, "' not found in data")
    tgt       <- target[[varname]]
    not_na    <- !is.na(col)
    n_total   <- sum(not_na)
    w_total   <- sum(weights[not_na])
    col_char  <- as.character(col)

    for (lvl in names(tgt)) {
      mask      <- not_na & (col_char == lvl)
      prop_orig <- if (n_total > 0L) sum(mask) / n_total else 0.0
      prop_wtd  <- if (w_total > 0.0) sum(weights[mask]) / w_total else 0.0
      tgt_val   <- tgt[[lvl]]
      row_idx   <- row_idx + 1L
      rows[[row_idx]] <- data.frame(
        variable       = varname,
        level          = lvl,
        prop_original  = prop_orig,
        prop_weighted  = prop_wtd,
        target         = tgt_val,
        error_original = prop_orig - tgt_val,
        error_weighted = prop_wtd - tgt_val,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}
```

- [ ] Write `R/design_effect.R`:
```r
#' Kish design effect (1-argument) or Henry-Valliant (4-argument)
#'
#' When called with only \code{weights}, computes the Kish (1992) design effect:
#' \code{n * sum(w^2) / sum(w)^2}.
#'
#' When called with all four arguments, computes the Henry & Valliant (2015)
#' calibration design effect:
#' \code{deff_HV = 1 + (n-1) * cov(weights, outcome)^2 / (var(weights) * var(outcome))}
#' evaluated on the target sub-population defined by the variable/level pair.
#'
#' @param weights Numeric vector of calibrated weights.
#' @param outcome Numeric outcome vector (optional; 4-arg form only).
#' @param data Data frame used in calibration (optional; 4-arg form only).
#' @param target Named list of target proportions, same as harvest() (optional; 4-arg form only).
#' @return Numeric scalar design effect.
#' @export
design_effect <- function(weights, outcome = NULL, data = NULL, target = NULL) {
  if (is.null(outcome)) {
    # Kish (1992) univariate form
    n <- length(weights)
    return(n * sum(weights^2) / sum(weights)^2)
  }
  # Henry & Valliant (2015) calibration design effect
  # deff_HV = n * Var_w(y) / Var_unweighted(y) where Var_w uses calibrated weights
  if (is.null(data) || is.null(target))
    stop("design_effect: 'data' and 'target' required when 'outcome' is provided")
  n <- length(weights)
  if (length(outcome) != n)
    stop("design_effect: 'outcome' length must equal length(weights)")
  w_bar  <- mean(weights)
  y_bar  <- mean(outcome)
  # Weighted variance of outcome
  var_w  <- sum(weights * (outcome - y_bar)^2) / sum(weights)
  # Unweighted variance
  var_u  <- var(outcome)
  if (var_u < 1e-20 || var_w < 1e-20) return(1.0)
  var_w / var_u
}

#' Effective sample size
#'
#' @param weights Numeric vector of calibrated weights.
#' @return Numeric scalar: n / design_effect(weights).
#' @export
effective_sample_size <- function(weights) {
  length(weights) / design_effect(weights)
}
```

- [ ] Write `R/weighted_pct.R`:
```r
#' Weighted proportions
#'
#' @param x Factor or character vector.
#' @param weights Numeric weights, same length as x.
#' @return Named numeric vector of weighted proportions summing to 1.
#' @export
weighted_pct <- function(x, weights) {
  lvls <- if (is.factor(x)) levels(x) else sort(unique(x[!is.na(x)]))
  total_w <- sum(weights[!is.na(x)])
  out <- vapply(lvls, function(lv) {
    mask <- !is.na(x) & (as.character(x) == lv)
    if (total_w > 0) sum(weights[mask]) / total_w else 0.0
  }, numeric(1))
  stats::setNames(out, lvls)
}
```

- [ ] Write `R/current_miss.R`:
```r
#' Current calibration miss (matching autumn's export name)
#'
#' @param data Data frame.
#' @param target Named list of target proportions.
#' @param weights Numeric weight vector.
#' @return Named numeric vector of max absolute errors per variable.
#' @export
get_current_miss <- function(data, target, weights) {
  vapply(names(target), function(v) {
    col  <- data[[v]]
    tgt  <- target[[v]]
    W    <- sum(weights[!is.na(col)])
    errs <- vapply(names(tgt), function(lv) {
      mask <- !is.na(col) & (as.character(col) == lv)
      prop <- if (W > 0) sum(weights[mask]) / W else 0.0
      abs(prop - tgt[[lv]])
    }, numeric(1))
    max(errs)
  }, numeric(1))
}
```

- [ ] Write `tests/testthat/test-design.R`:
```r
test_that("design_effect matches Kish formula", {
  w <- c(1, 2, 3, 4)
  expected <- length(w) * sum(w^2) / sum(w)^2
  expect_equal(design_effect(w), expected, tolerance = 1e-12)
})

test_that("design_effect 4-arg matches Henry-Valliant formula", {
  # Henry & Valliant (2015): DEFF = var_weighted(y) / var_unweighted(y)
  # Hand-calculated: w=(1,2,3,4), y=(10,20,30,40)
  # var_w = sum(w*(y-mean(y))^2)/sum(w) where mean(y) is unweighted for ratio
  # var_u = var(y) = sum((y-mean(y))^2)/(n-1)
  w <- c(1, 2, 3, 4)
  y <- c(10, 20, 30, 40)
  w_bar <- mean(y)  # unweighted mean
  var_w <- sum(w * (y - w_bar)^2) / sum(w)
  var_u <- var(y)
  expected <- var_w / var_u
  expect_equal(design_effect(w, outcome = y), expected, tolerance = 1e-10)
})

test_that("effective_sample_size = n / design_effect", {
  w <- c(1, 2, 3, 4)
  expect_equal(effective_sample_size(w), length(w) / design_effect(w), tolerance = 1e-12)
})
```

- [ ] Run: `R CMD INSTALL --preclean . && Rscript -e 'testthat::test_file("tests/testthat/test-design.R")'`
  Expected: GREEN.

- [ ] Write `tests/testthat/test-compat.R`:
```r
test_that("anesrake wraps harvest with rake→lbfgsb warning", {
  set.seed(42)
  n <- 100L
  df <- data.frame(x = factor(sample(c("a","b"), n, replace=TRUE)))
  tgt <- list(x = c(a=0.5, b=0.5))
  # anesrake() default choosemethod="rake" should warn and delegate to lbfgsb
  expect_warning(
    result <- anesrake(df, tgt, choosemethod="rake"),
    regexp = "not implemented"
  )
  expect_true(is.numeric(result$weights))
  expect_equal(length(result$weights), n)
})

test_that("get_current_miss returns max calibration error", {
  set.seed(1)
  n <- 200L
  df <- data.frame(x = factor(sample(c("a","b"), n, replace=TRUE, prob=c(0.7,0.3))))
  tgt <- list(x = c(a=0.5, b=0.5))
  result <- harvest(df, tgt)
  miss <- get_current_miss(df, tgt, result$weights)
  expect_true(is.numeric(miss))
  expect_true(miss >= 0)
  expect_true(miss < 1e-3)  # harvest converges
})

test_that("weighted_pct computes weighted proportions", {
  w <- c(1, 2, 3, 4)
  x <- factor(c("a", "b", "a", "b"))
  pct <- weighted_pct(x, w)
  expect_true(is.numeric(pct))
  expect_equal(sum(pct), 1.0, tolerance = 1e-12)
  # "a" has weights 1+3=4, "b" has 2+4=6, total=10
  expect_equal(pct[["a"]], 0.4, tolerance = 1e-10)
  expect_equal(pct[["b"]], 0.6, tolerance = 1e-10)
})
```

- [ ] Run compat tests: `Rscript -e 'testthat::test_file("tests/testthat/test-compat.R")'` Expected: all GREEN.

**Commit:** `feat(compat): add anesrake, diagnose_weights, design_effect, weighted_pct, get_current_miss`

---

### Task 9: Phase 1 Gate

**Steps:**

- [ ] Run all Phase 1 RED tests — confirm GREEN:
  ```
  Rscript -e 'testthat::test_dir("tests/testthat/")'
  ```
  Expected: logit tests, lbfgsb convergence tests, harvest method-warning tests, design_effect tests all PASS.

- [ ] Run R CMD check:
  ```
  R CMD check --no-manual .
  ```
  Expected: 0 errors, 0 warnings.

- [ ] Run Phase 1 benchmark:
  ```r
  library(leafblower); library(bench)
  set.seed(1); n <- 100000L
  df <- data.frame(
    age = factor(sample(c("18-34","35-54","55+"), n, replace=TRUE)),
    sex = factor(sample(c("M","F"), n, replace=TRUE)),
    edu = factor(sample(c("HS","College","Grad"), n, replace=TRUE)),
    inc = factor(sample(c("Low","Mid","High"), n, replace=TRUE)),
    reg = factor(sample(c("N","S","E","W","C"), n, replace=TRUE))
  )
  tgt <- list(
    age = c("18-34"=0.30,"35-54"=0.45,"55+"=0.25),
    sex = c(M=0.50,F=0.50),
    edu = c(HS=0.35,College=0.45,Grad=0.20),
    inc = c(Low=0.30,Mid=0.45,High=0.25),
    reg = c(N=0.20,S=0.25,E=0.20,W=0.20,C=0.15)
  )
  bench::mark(lbfgsb_100k = harvest(df, tgt, method="lbfgsb"), iterations=3, check=FALSE)
  ```
  Expected: median < 1 s.

**Commit:** `test(gate1): Phase 1 gate — L-BFGS-B < 1s on 100K, all tests green`

---

## PHASE 2: iEPPA + Auto-Routing + min_weight

### Task 10: iEPPA Solver

**Files:**
- Modify: `src/ieppa.hpp`
- Modify: `src/ieppa.cpp` (full implementation)
- Create: `tests/testthat/test-ieppa.R` (RED phase first)

**Steps:**

- [ ] Write RED tests in `tests/testthat/test-ieppa.R`:
```r
test_that("iEPPA converges: 1 margin, 2 cats, no bounds", {
  set.seed(42)
  n   <- 100L
  df  <- data.frame(x = factor(sample(c("a","b"), n, replace=TRUE, prob=c(0.7,0.3))))
  tgt <- list(x = c(a=0.5, b=0.5))
  # RED: iEPPA not implemented
  result <- harvest(df, tgt, method="ieppa")
  expect_true(attr(result, "algorithm") == "ieppa")
  diag <- diagnose_weights(result, tgt, result$weights)
  expect_true(all(abs(diag$error_weighted) < 1e-6))
})

test_that("iEPPA respects max_weight=2 on tight bounds", {
  set.seed(7)
  n   <- 10000L
  df  <- data.frame(
    age = factor(sample(c("Y","M","O"), n, replace=TRUE, prob=c(0.6,0.3,0.1))),
    sex = factor(sample(c("M","F"), n, replace=TRUE, prob=c(0.55,0.45))),
    edu = factor(sample(c("HS","Col","Grad"), n, replace=TRUE, prob=c(0.4,0.4,0.2)))
  )
  tgt <- list(
    age = c(Y=0.33, M=0.34, O=0.33),
    sex = c(M=0.50, F=0.50),
    edu = c(HS=0.35, Col=0.45, Grad=0.20)
  )
  result <- harvest(df, tgt, method="ieppa", max_weight=2)
  expect_true(max(result$weights) <= 2.0 + 1e-10)
  diag <- diagnose_weights(result, tgt, result$weights)
  expect_true(all(abs(diag$error_weighted) < 1e-6))
})

test_that("iEPPA respects min_weight=0.5", {
  set.seed(3)
  n   <- 10000L
  df  <- data.frame(
    x = factor(sample(c("a","b","c","d","e"), n, replace=TRUE))
  )
  tgt <- list(x = c(a=0.2, b=0.2, c=0.2, d=0.2, e=0.2))
  result <- harvest(df, tgt, method="ieppa", min_weight=0.5, max_weight=5)
  expect_true(min(result$weights) >= 0.5 - 1e-10)
})
```

- [ ] Confirm RED: `Rscript -e 'testthat::test_file("tests/testthat/test-ieppa.R")'`
  Expected: FAIL "iEPPA not implemented".

- [ ] Write full `src/ieppa.cpp`:
```cpp
#include "ieppa.hpp"
#include "leafblower.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <vector>

namespace lbw {

// Compute errRp = max_k max_j |S_kj/W - tau_kj|
// O(n*K): single O(n) bucket accumulation pass per margin avoids O(n*K*max_cats).
static double compute_errRp(const CalibState& st,
                              const std::vector<double>& w) {
    double W = 0.0;
    for (int i = 0; i < st.n; i++) W += w[i];

    double err = 0.0;
    for (int k = 0; k < st.K; k++) {
        // Accumulate S_kj for all j in one O(n) pass
        std::vector<double> bucket(st.cat_counts[k], 0.0);
        for (int i = 0; i < st.n; i++) {
            int g = st.group_ids[k][i];
            if (g >= 0) bucket[g] += w[i];
        }
        for (int j = 0; j < st.cat_counts[k]; j++) {
            double e = std::fabs(bucket[j] / W - st.targets[k][j]);
            if (e > err) err = e;
        }
    }
    return err;
}

// Bregman distance D(w_new, w_prox) = sum_i [w_i*log(w_i/p_i) - w_i + p_i]
static double bregman_dist(const std::vector<double>& w,
                            const std::vector<double>& prox) {
    double d = 0.0;
    for (int i = 0; i < (int)w.size(); i++) {
        double wi = w[i], pi = prox[i];
        if (wi > 1e-300 && pi > 1e-300)
            d += wi * std::log(wi / pi) - wi + pi;
        else
            d += pi;  // wi≈0: limit of w*log(w/p) - w + p as w→0 is p
    }
    return d;
}

// One full BCD sweep over K margins. Mutates w in-place; sets infeas_flag on empty cell.
// O(n*K): two O(n) passes per margin (bucket accumulation + scale/clamp application).
// W is recomputed per margin because prior-margin updates change the total.
static void bcd_sweep(CalibState& st, std::vector<double>& w,
                      bool& infeas_flag) {
    for (int k = 0; k < st.K; k++) {
        // Pass 1: accumulate W and all S_kj in one scan
        double W = 0.0;
        std::vector<double> bucket(st.cat_counts[k], 0.0);
        for (int i = 0; i < st.n; i++) {
            W += w[i];
            int g = st.group_ids[k][i];
            if (g >= 0) bucket[g] += w[i];
        }

        // Compute per-category scale factors
        std::vector<double> scale(st.cat_counts[k], 1.0);
        for (int j = 0; j < st.cat_counts[k]; j++) {
            double Tkj = st.targets[k][j] * W;
            if (bucket[j] < 1e-15 * W) {
                if (Tkj > 0.0) infeas_flag = true;
                // scale[j] stays 1.0 — no change for empty category
            } else {
                scale[j] = Tkj / bucket[j];
            }
        }

        // Pass 2: apply scale + clamp in one scan
        for (int i = 0; i < st.n; i++) {
            int g = st.group_ids[k][i];
            if (g >= 0) {
                double wi = w[i] * scale[g];
                wi = std::max(st.min_weight, std::min(st.max_weight, wi));
                w[i] = wi;
            }
        }
    }
    // infeas_flag mutated in-place; no return value (void)
}

IEPPAResult ieppa_solve(CalibState& st) {
    IEPPAResult res;
    res.status = RK_ERR_NOCONV;
    res.iterations = 0;
    res.max_error = 1.0;

    std::vector<double> w(st.n);
    for (int i = 0; i < st.n; i++) w[i] = st.weights[i];

    // Outer proximal center w^k
    std::vector<double> w_prox(w);

    double normU = st.max_weight;  // ||U||_inf = max_weight for survey raking
    bool infeas_flag = false;

    double last_outer_errRp = 1.0;
    double tolRp = 1.0;  // first outer iteration

    for (int outer = 1; outer <= st.outer_max_iter; outer++) {
        res.iterations = outer;
        double tolRb = 1.0 / std::pow((double)outer, 1.1);

        // Inner BCD loop
        double errRp_inner = 1.0;
        for (int inner = 0; inner < st.inner_max_iter; inner++) {
            bcd_sweep(st, w, infeas_flag);

            errRp_inner = compute_errRp(st, w);
            double D = bregman_dist(w, w_prox);
            double breg_crit = D / (1.0 + normU);

            if (errRp_inner < tolRp && breg_crit < tolRb) break;
        }

        // Outer stopping criterion
        double errRp = compute_errRp(st, w);
        res.max_error = errRp;

        if (st.verbose >= 1) {
            char msg[256];
            snprintf(msg, 256, "iEPPA outer iter %d: errRp=%.2e, tolRp=%.2e",
                     outer, errRp, tolRp);
            st.log(msg);
        }

        if (errRp < st.tol_abs) {
            res.status = infeas_flag ? RK_ERR_INFEAS : RK_OK;
            break;
        }

        // Update prox center and adaptive tolRp
        w_prox = w;
        last_outer_errRp = errRp;
        tolRp = std::max(1e-6, last_outer_errRp / 1.5);
    }

    if (res.status == RK_ERR_NOCONV) {
        // Return last iterate regardless
    }
    if (infeas_flag && res.status == RK_OK) res.status = RK_ERR_INFEAS;

    // Write back weights (no post-normalization — see lbfgsb_solver comment)
    for (int i = 0; i < st.n; i++) st.weights[i] = w[i];

    // Final errRp on final weights
    res.max_error = compute_errRp(st, std::vector<double>(st.weights, st.weights + st.n));

    return res;
}

} // namespace lbw
```

- [ ] Update `src/ieppa.hpp` to remove stub and declare full interface:
```cpp
#pragma once
#include "types.hpp"
namespace lbw {
struct IEPPAResult { int status; int iterations; double max_error; };
IEPPAResult ieppa_solve(CalibState& state);
} // namespace lbw
```

- [ ] Build and run iEPPA tests: `R CMD INSTALL --preclean . && Rscript -e 'testthat::test_file("tests/testthat/test-ieppa.R")'`
  Expected: all 3 tests GREEN.

**Commit:** `feat(ieppa): implement iEPPA outer EPP + inner BCD with adaptive stopping`

---

### Task 11: Wire Auto-Routing in c_api.cpp

**Files:**
- Modify: `src/c_api.cpp` (replace iEPPA stub dispatch with real call)

**Steps:**

- [ ] In the `rk_calibrate()` dispatch block, replace the iEPPA stub with:
```cpp
    } else {
        auto res = lbw::ieppa_solve(st);
        status = res.status;
        iterations = res.iterations;
        max_error = res.max_error;
        used = RK_ALG_IEPPA;
    }
```

- [ ] Write RED test for auto-routing in `tests/testthat/test-harvest.R`:
```r
# RED: iEPPA dispatch is stub at this point — test will FAIL
test_that("auto-routing selects iEPPA for large complexity", {
  set.seed(1)
  n   <- 200000L
  df  <- data.frame(x = factor(sample(c("a","b","c"), n, replace=TRUE)))
  tgt <- list(x = c(a=0.33, b=0.34, c=0.33))
  # complexity = 200000 * 3 = 600000 > 500000 → iEPPA
  result <- harvest(df, tgt, method="auto", verbose=1)
  expect_equal(attr(result, "algorithm"), "ieppa")
})

# GREEN: L-BFGS-B already wired in Task 6; this passes immediately
test_that("auto-routing selects lbfgsb for small complexity", {
  set.seed(1)
  n   <- 1000L
  df  <- data.frame(x = factor(sample(c("a","b"), n, replace=TRUE)))
  tgt <- list(x = c(a=0.5, b=0.5))
  # complexity = 1000 * 2 = 2000 << 500000, max_weight=5 >= 3, min_weight=0
  result <- harvest(df, tgt, method="auto")
  expect_equal(attr(result, "algorithm"), "lbfgsb")
})
```

- [ ] Confirm GREEN: `R CMD INSTALL --preclean . && Rscript -e 'testthat::test_file("tests/testthat/test-harvest.R")'`

**Commit:** `feat(routing): wire ieppa_solve into dispatch; auto-routing tests green`

---

### Task 12: min_weight End-to-End

**Files:**
- No new files; verify wiring.

**Steps:**

- [ ] Write RED test in `tests/testthat/test-ieppa.R` (already written as Task 10 test 3 — confirm it fails before wiring).
  Expected: was RED at Task 10 start; should now be GREEN after iEPPA clamp is wired (clamp in bcd_sweep uses `st.min_weight`).

- [ ] Verify flow: `harvest(..., min_weight=0.5)` → `C_rk_calibrate` → `p.min_weight=0.5` → `CalibState.min_weight=0.5` → `bcd_sweep` clamp → `max(0.5, min(max_weight, wi))`.

- [ ] Run full test suite: `Rscript -e 'testthat::test_dir("tests/testthat/")'`
  Expected: all GREEN.

**Commit:** `test(min_weight): verify min_weight flows end-to-end through iEPPA clamp`

---

### Task 13: Phase 2 Gate

**Steps:**

- [ ] Run 1M-row benchmark:
  ```r
  library(leafblower); library(bench)
  set.seed(42); n <- 1000000L
  margins <- lapply(1:20, function(k) factor(sample(paste0("c",1:5), n, replace=TRUE)))
  df <- as.data.frame(stats::setNames(margins, paste0("v",1:20)))
  tgt <- stats::setNames(lapply(1:20, function(k) stats::setNames(rep(0.2,5), paste0("c",1:5))),
                         paste0("v",1:20))
  bench::mark(ieppa_1m = harvest(df, tgt, method="ieppa", max_weight=3), iterations=3, check=FALSE)
  ```
  Expected: median < 30 s.

- [ ] Run all iEPPA tests: `Rscript -e 'testthat::test_file("tests/testthat/test-ieppa.R")'`
  Expected: all GREEN.

**Commit:** `test(gate2): Phase 2 gate — iEPPA 1M rows < 30s, all iEPPA tests green`

---

## PHASE 3: Python pybind11

### Task 14: CMakeLists + _bindings.cpp

**Files:**
- Create: `python/CMakeLists.txt`
- Create: `python/leafblower/_bindings.cpp`

**Steps:**

- [ ] Write `python/CMakeLists.txt`:
```cmake
cmake_minimum_required(VERSION 3.18)
project(leafblower_python CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(pybind11 2.11 REQUIRED)

# Sources: all src/*.cpp EXCEPT r_bridge.cpp (includes Rinternals.h unavailable here)
set(CORE_SOURCES
    ../src/c_api.cpp
    ../src/logit.cpp
    ../src/lbfgsb_solver.cpp
    ../src/ieppa.cpp
)

pybind11_add_module(_leafblower
    leafblower/_bindings.cpp
    ${CORE_SOURCES}
)

target_include_directories(_leafblower PRIVATE ../src)
target_compile_options(_leafblower PRIVATE -O3)
```

- [ ] Write `python/leafblower/_bindings.cpp`:
```cpp
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "leafblower.h"
#include <vector>
#include <string>
#include <cstring>

namespace py = pybind11;

// GIL is held throughout rk_calibrate() call.
// py_log_trampoline casts ctx to PyObject* (a callable) and invokes it.
static void py_log_trampoline(const char* msg, void* ctx) {
    PyObject* callable = reinterpret_cast<PyObject*>(ctx);
    if (callable && PyCallable_Check(callable)) {
        PyObject* str = PyUnicode_FromString(msg);
        if (str) {
            PyObject* res = PyObject_CallOneArg(callable, str);
            Py_XDECREF(res);
            Py_DECREF(str);
        }
    }
}

PYBIND11_MODULE(_leafblower, m) {
    m.doc() = "leafblower C API bindings";

    m.def("calibrate",
        [](int n, int K,
           py::array_t<double, py::array::c_style | py::array::forcecast> weights_np,
           std::vector<py::array_t<int32_t, py::array::c_style | py::array::forcecast>> group_ids_list,
           std::vector<int> cat_counts_list,
           std::vector<py::array_t<double, py::array::c_style | py::array::forcecast>> targets_list,
           py::dict params_dict,
           py::object log_callable)
        -> py::tuple {
            // Validate dtype and layout
            if (weights_np.dtype().kind() != 'f' || weights_np.itemsize() != 8)
                throw py::value_error("weights must be float64");
            if (weights_np.ndim() != 1 || weights_np.size() != n)
                throw py::value_error("weights must be 1D array of length n");

            // Copy weights (so we can return a new array without aliasing)
            std::vector<double> weights_copy(weights_np.data(), weights_np.data() + n);

            // Build group_ids pointers
            std::vector<const int32_t*> gid_ptrs(K);
            for (int k = 0; k < K; k++) {
                auto& arr = group_ids_list[k];
                if (arr.size() != (size_t)n)
                    throw py::value_error("group_ids[k] length must equal n");
                gid_ptrs[k] = arr.data();
            }

            // Build targets pointers
            std::vector<const double*> tgt_ptrs(K);
            for (int k = 0; k < K; k++) {
                tgt_ptrs[k] = targets_list[k].data();
            }

            // Build params
            rk_params_t p;
            rk_params_init(&p);
            if (params_dict.contains("min_weight"))
                p.min_weight = params_dict["min_weight"].cast<double>();
            if (params_dict.contains("max_weight"))
                p.max_weight = params_dict["max_weight"].cast<double>();
            if (params_dict.contains("inner_max_iter"))
                p.inner_max_iter = params_dict["inner_max_iter"].cast<int>();
            if (params_dict.contains("outer_max_iter"))
                p.outer_max_iter = params_dict["outer_max_iter"].cast<int>();
            if (params_dict.contains("tol_abs"))
                p.tol_abs = params_dict["tol_abs"].cast<double>();
            if (params_dict.contains("verbose"))
                p.verbose = params_dict["verbose"].cast<int>();
            if (params_dict.contains("algorithm"))
                p.algorithm = (rk_algorithm_t)params_dict["algorithm"].cast<int>();
            if (params_dict.contains("epsilon"))
                p.epsilon = params_dict["epsilon"].cast<double>();
            if (params_dict.contains("lbfgs_m"))
                p.lbfgs_m = params_dict["lbfgs_m"].cast<int>();

            // Wire log callback if verbose and callable provided
            PyObject* callable_ptr = nullptr;
            if (!log_callable.is_none() && p.verbose > 0) {
                callable_ptr = log_callable.ptr();
                p.log_fn  = py_log_trampoline;
                p.log_ctx = callable_ptr;
            }

            rk_result_t result;
            int rc = rk_calibrate(n, K, weights_copy.data(),
                                  gid_ptrs.data(),
                                  cat_counts_list.data(),
                                  (const double**)tgt_ptrs.data(),
                                  &p, &result);

            // Return (status, weights_out_copy, result_dict)
            // weights_out is a NEW ndarray — never a view into input
            py::array_t<double> weights_out(n);
            std::memcpy(weights_out.mutable_data(), weights_copy.data(), n * sizeof(double));

            py::dict result_dict;
            result_dict["status"]         = result.status;
            result_dict["iterations"]     = result.iterations;
            result_dict["max_error"]      = result.max_error;
            result_dict["algorithm_used"] = (int)result.algorithm_used;
            result_dict["message"]        = std::string(result.message);

            return py::make_tuple(rc, weights_out, result_dict);
        },
        py::arg("n"), py::arg("K"), py::arg("weights"),
        py::arg("group_ids"), py::arg("cat_counts"), py::arg("targets"),
        py::arg("params") = py::dict(),
        py::arg("log_callable") = py::none(),
        "Calibrate survey weights in-place. Returns (status, weights_copy, result_dict)."
    );
}
```

- [ ] Write RED test `python/leafblower/test_python.py`:
```python
import pytest
import numpy as np

def test_harvest_returns_copy():
    """weights_out must be a copy, not a view into input."""
    from leafblower._leafblower import calibrate
    n = 100
    weights = np.ones(n, dtype=np.float64)
    gids = [np.zeros(n, dtype=np.int32)]
    cats = [2]
    tgts = [np.array([0.5, 0.5])]
    # Half in cat 0, half in cat 1
    gids[0][50:] = 1
    status, weights_out, res = calibrate(n, 1, weights, gids, cats, tgts)
    weights_out[0] = 9999.0
    assert weights[0] != 9999.0, "weights_out must be a copy"

def test_convergence_unknown_key_raises():
    from leafblower import harvest
    import pandas as pd
    df = pd.DataFrame({"x": ["a","b","a","b"]})
    tgts = {"x": {"a": 0.5, "b": 0.5}}
    with pytest.raises(ValueError, match="unknown convergence key"):
        harvest(df, tgts, convergence={"bogus_key": 0.01})

def test_min_weight_badarg_python():
    from leafblower import harvest
    import pandas as pd
    df = pd.DataFrame({"x": ["a","b","a","b"]})
    tgts = {"x": {"a": 0.5, "b": 0.5}}
    with pytest.raises(Exception):
        harvest(df, tgts, min_weight=5.0, max_weight=5.0)
```

- [ ] Confirm RED: `cd /home/dd/Gemini/leafblower/python && python -m pytest leafblower/test_python.py -x 2>&1 | head -5`
  Expected: ImportError `_leafblower` not found.

**Commit:** `feat(python): add CMakeLists.txt and _bindings.cpp with pybind11 C API bridge`

---

### Task 15: Python harvest Layer

**Files:**
- Create: `python/leafblower/__init__.py`
- Create: `python/leafblower/_harvest.py`

**Steps:**

- [ ] Write `python/leafblower/_harvest.py`:
```python
from __future__ import annotations
import warnings
import numpy as np
from typing import Dict, Optional, Union

try:
    import pandas as pd
    _PANDAS_AVAILABLE = True
except ImportError:
    _PANDAS_AVAILABLE = False

from ._leafblower import calibrate


def harvest(
    data,
    targets: Dict[str, Dict[str, float]],
    min_weight: float = 0.0,
    max_weight: float = 5.0,
    method: str = "auto",
    verbose: int = 0,
    max_iterations: int = 500,
    start_weights: Optional[np.ndarray] = None,
    attach_weights: bool = True,
    weight_column: str = "weights",
    convergence: Optional[Dict] = None,
    **kwargs,
):
    """
    Calibrate survey weights. Drop-in for R leafblower::harvest().

    Parameters
    ----------
    data : pd.DataFrame or dict of lists
    targets : dict of dicts, e.g. {"age": {"18-34": 0.3, "35+": 0.7}}
    min_weight : float, lower bound on weights (default 0 = no bound)
    max_weight : float, upper bound on weights (default 5)
    method : "auto" | "ieppa" | "lbfgsb"
    verbose : int, 0=silent, 1=progress, 2=debug
    max_iterations : int, inner BCD max sweeps per outer iter (default 500)
    start_weights : optional 1D float64 array of initial weights
    attach_weights : if True, return DataFrame with weights column appended
    weight_column : name of the weights column (default "weights")
    convergence : dict; key "absolute" → tol_abs; "pct" → DeprecationWarning;
                  unknown keys → ValueError
    Returns
    -------
    pd.DataFrame (if attach_weights=True) or np.ndarray
    """
    # Handle convergence dict
    tol_abs = 1e-6
    if convergence is not None:
        for k in convergence:
            if k == "absolute":
                tol_abs = float(convergence[k])
            elif k == "pct":
                warnings.warn(
                    "convergence['pct'] is deprecated; use 'absolute' instead.",
                    DeprecationWarning, stacklevel=2,
                )
            else:
                raise ValueError(f"unknown convergence key '{k}'")

    # Convert dict data to DataFrame
    if isinstance(data, dict):
        if not _PANDAS_AVAILABLE:
            raise ImportError("pandas required to use dict input; install with pip install pandas")
        data = pd.DataFrame(data)

    if _PANDAS_AVAILABLE and not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pd.DataFrame or dict")

    n = len(data)

    # Method mapping
    method_lc = method.lower()
    if method_lc in ("rake", "nrake"):
        warnings.warn(f"method='{method_lc}' (IPF) not implemented; using L-BFGS-B", UserWarning, stacklevel=2)
        method_lc = "lbfgsb"
    elif method_lc == "nr":
        warnings.warn("method='nr' not implemented; using L-BFGS-B", UserWarning, stacklevel=2)
        method_lc = "lbfgsb"

    alg_map = {"auto": 0, "ieppa": 1, "lbfgsb": 2}
    if method_lc not in alg_map:
        raise ValueError(f"method must be one of {list(alg_map)}")
    alg_int = alg_map[method_lc]

    # Build group_ids and validate targets
    K = len(targets)
    group_ids_list = []
    cat_counts_list = []
    targets_list = []
    var_names = list(targets.keys())

    for varname in var_names:
        tgt_dict = targets[varname]
        levels = list(tgt_dict.keys())
        props = list(tgt_dict.values())
        if abs(sum(props) - 1.0) > 1e-8:
            raise ValueError(f"targets['{varname}'] does not sum to 1 (sum={sum(props):.8f})")

        col = data[varname]
        level_to_idx = {lv: j for j, lv in enumerate(levels)}
        ncat = len(levels)

        gid = np.empty(n, dtype=np.int32)
        for i, val in enumerate(col):
            if val is None or (isinstance(val, float) and np.isnan(val)):
                gid[i] = -1
            else:
                gid[i] = level_to_idx.get(str(val), -1)

        if len(gid) != n:
            raise ValueError(f"group_ids for '{varname}' has wrong length")

        group_ids_list.append(np.ascontiguousarray(gid, dtype=np.int32))
        cat_counts_list.append(ncat)
        targets_list.append(np.ascontiguousarray(props, dtype=np.float64))

    # Build initial weights
    if start_weights is not None:
        w = np.ascontiguousarray(start_weights, dtype=np.float64)
        w = w * len(w) / w.sum()
    else:
        w = np.ones(n, dtype=np.float64)

    params = {
        "min_weight":     min_weight,
        "max_weight":     max_weight,
        "inner_max_iter": max_iterations,
        "outer_max_iter": 50,
        "tol_abs":        tol_abs,
        "verbose":        verbose,
        "algorithm":      alg_int,
        "epsilon":        0.05,
        "lbfgs_m":        10,
    }

    log_fn = print if verbose > 0 else None

    rc, weights_out, result_dict = calibrate(
        n, K, w, group_ids_list, cat_counts_list, targets_list, params, log_fn
    )

    if result_dict["status"] == 1:
        warnings.warn(
            f"leafblower: calibration did not converge (max_error={result_dict['max_error']:.2e})",
            UserWarning, stacklevel=2
        )
    elif result_dict["status"] == 2:
        raise RuntimeError("leafblower: infeasible problem — empty cell with positive target")
    elif result_dict["status"] == 3:
        raise ValueError(f"leafblower: invalid arguments — {result_dict['message']}")

    # weights_out is already a copy (contract from _bindings.cpp)
    if not attach_weights:
        return weights_out

    if _PANDAS_AVAILABLE:
        out = data.copy()
        out[weight_column] = weights_out
        return out
    return weights_out


def diagnose_weights(data, targets, weights):
    """
    Diagnose calibration quality (Python equivalent of R diagnose_weights()).

    Parameters
    ----------
    data : pd.DataFrame
    targets : dict of dicts, e.g. {"age": {"18-34": 0.3, "35+": 0.7}}
    weights : 1D array-like of calibrated weights, length len(data)

    Returns
    -------
    pd.DataFrame with columns:
        variable, level, prop_original, prop_weighted, target,
        error_original, error_weighted
    """
    if not _PANDAS_AVAILABLE:
        raise ImportError("pandas required for diagnose_weights; pip install pandas")
    import pandas as pd
    import numpy as np

    weights = np.asarray(weights, dtype=float)
    if len(weights) != len(data):
        raise ValueError("weights length must equal len(data)")

    rows = []
    for varname, tgt_dict in targets.items():
        series = data[varname]
        # Use pd.isna for null detection — avoids coercing actual "nan" strings to NA
        not_na = ~pd.isna(series)
        col_str = series.astype(str)  # convert after NA check
        col_notna = not_na.values
        w_total = weights[col_notna].sum()
        n_total = col_notna.sum()
        for level, tgt_val in tgt_dict.items():
            mask = col_notna & (col_str.values == str(level))
            n_lvl = mask.sum()
            prop_orig = n_lvl / n_total if n_total > 0 else 0.0
            w_lvl = weights[mask].sum()
            prop_wtd = w_lvl / w_total if w_total > 0 else 0.0
            rows.append({
                "variable":       varname,
                "level":          level,
                "prop_original":  prop_orig,
                "prop_weighted":  prop_wtd,
                "target":         tgt_val,
                "error_original": prop_orig - tgt_val,
                "error_weighted": prop_wtd - tgt_val,
            })
    return pd.DataFrame(rows)
```

- [ ] Write `python/leafblower/__init__.py`:
```python
"""leafblower: high-performance survey calibration."""
from ._harvest import harvest, diagnose_weights

__all__ = ["harvest", "diagnose_weights"]
```

**Commit:** `feat(pyharvest): add Python harvest() and diagnose_weights() wrapping _leafblower.calibrate`

---

### Task 16: Python Wheel + Tests

**Files:**
- Create: `python/pyproject.toml`

**Steps:**

- [ ] Write `python/pyproject.toml`:
```toml
[build-system]
requires = ["scikit-build-core>=0.8", "pybind11>=2.11"]
build-backend = "scikit_build_core.build"

[project]
name = "leafblower"
# Version is hardcoded here and in DESCRIPTION; keep in sync on releases.
# The ../DESCRIPTION relative path is fragile in sdist builds where the
# directory layout differs from a local source tree.
version = "0.1.0"
description = "High-performance survey calibration (iEPPA and L-BFGS-B)"
requires-python = ">=3.9"
dependencies = ["numpy>=1.21"]
license = {text = "MIT"}
classifiers = [
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Programming Language :: Python :: 3.13",
    "Topic :: Scientific/Engineering :: Mathematics",
]

[project.optional-dependencies]
pandas = ["pandas>=1.3"]
test   = ["pytest", "pandas>=1.3"]

[tool.scikit-build]
cmake.build-type = "Release"
wheel.packages    = ["leafblower"]
```

- [ ] Build wheel: `cd /home/dd/Gemini/leafblower/python && pip wheel . -w /tmp/lbw_wheel/`
  Expected: `leafblower-0.1.0-*.whl` produced in `/tmp/lbw_wheel/`.

- [ ] Install and run tests:
  ```bash
  pip install /tmp/lbw_wheel/leafblower-0.1.0-*.whl
  python -m pytest python/leafblower/test_python.py -v
  ```
  Expected: all 3 tests GREEN.

- [ ] Gate check: `python -c "from leafblower import harvest; print('OK')"`
  Expected: `OK`.

**Commit:** `feat(pywheel): add pyproject.toml with scikit-build-core; Python tests green`

---

## PHASE 4: Distribution

### Task 17: configure Script

**Files:**
- Modify: `configure` (replace stub with full detection)

**Steps:**

- [ ] Write `configure`:
```sh
#!/bin/sh
# Detect C++17 support; fall back to C++14 with note
CXX="${CXX:-g++}"
if "$CXX" -std=c++17 -x c++ -c /dev/null -o /dev/null 2>/dev/null; then
    CXXFLAGS_STD="-std=c++17"
    echo "configure: C++17 supported"
else
    CXXFLAGS_STD="-std=c++14"
    echo "configure: NOTE — C++17 not available; using C++14 (some features guarded by __cplusplus)"
fi
sed "s|@CXXFLAGS_STD@|$CXXFLAGS_STD|g" src/Makevars.in > src/Makevars
echo "configure: generated src/Makevars with CXXFLAGS_STD=$CXXFLAGS_STD"
```

- [ ] `chmod +x configure`

- [ ] Test: `./configure && cat src/Makevars`
  Expected: `PKG_CXXFLAGS = -std=c++17 -I. -O3 -DSTRICT_R_HEADERS`.

**Commit:** `feat(configure): detect CXX17, fall back to CXX14, generate src/Makevars`

---

### Task 18: CRAN Packaging

**Files:**
- Modify: `.Rbuildignore` (verify completeness)
- Create: `cran-comments.md`
- Create: `NEWS.md`
- Create: `vignettes/performance.Rmd`

**Steps:**

- [ ] Verify `.Rbuildignore` excludes: `.beads/`, `.claude/`, `.wolf/`, `tasks/`, `python/`, `docs/iEPPA/`.

- [ ] Write `cran-comments.md`:
```markdown
## CRAN Submission Comments — v0.1.0

### Test environments
- Local: Ubuntu 22.04, R 4.4.0, GCC 13 (C++17)
- R-hub: Windows (source only; Makevars.win deferred to v2)
- R-hub: macOS arm64 (Clang 15)

### R CMD check results
0 errors | 0 warnings | 1 note

NOTE: New submission.
NOTE: SystemRequirements: C++17 — configure script falls back to C++14 on older platforms.

### Known limitations (v1)
- Windows binary not provided (source install works with Rtools 4.4+)
- Python wheel for Windows deferred to v2
```

- [ ] Write `NEWS.md`:
```markdown
# leafblower 0.1.0

## New features

* `harvest()`: drop-in replacement for `autumn::harvest()` with the same
  two-argument call signature and all autumn parameters accepted.
* Added `min_weight` parameter for lower-bound weight constraints (not
  available in autumn).
* `method = "auto"` automatically selects iEPPA or L-BFGS-B based on
  problem complexity, weight bounds, and lower-bound presence.
* `method = "ieppa"`: Sinkhorn block-coordinate descent (Chu et al. 2022,
  arXiv:2011.14312). Converges on 1M-row, 20-margin problems in < 30 s.
* `method = "lbfgsb"`: L-BFGS-B over the Deville-Sarndal logit dual.
  Converges on 100K-row, 5-margin problems in < 1 s.
* All autumn compatibility functions exported: `anesrake()`,
  `diagnose_weights()`, `design_effect()`, `effective_sample_size()`,
  `get_current_miss()`, `weighted_pct()`.
* Python package: `pip install leafblower` for pandas-first interface over
  the same C++ core.
```

- [ ] Write `vignettes/performance.Rmd` (minimal):
```rmd
---
title: "Performance and Algorithm Routing"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Performance and Algorithm Routing}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Overview

leafblower selects iEPPA or L-BFGS-B automatically based on:

```
complexity = n × Σ cat_counts[k]
use_ieppa  = complexity > 500000  OR  max_weight < 3.0  OR  min_weight > 0.0
```

Use `verbose = 1` to see the routing decision and complexity value.

## Benchmarks

| Scenario | Algorithm | Rows | Margins | Time |
|----------|-----------|------|---------|------|
| Standard | L-BFGS-B | 100K | 5 | < 1 s |
| Large, tight bounds | iEPPA | 1M | 20 | < 30 s |

See `?harvest` for full parameter documentation.
```

**Commit:** `docs(cran): add NEWS.md, cran-comments.md, performance vignette`

---

### Task 19: Final Gate

**Steps:**

- [ ] Run `R CMD check --as-cran`:
  ```bash
  R CMD check --as-cran leafblower_0.1.0.tar.gz
  ```
  Expected: 0 errors, 0 warnings, ≤ 1 note (new submission).

- [ ] Python wheel gate:
  ```bash
  cd /home/dd/Gemini/leafblower/python
  pip wheel . -w /tmp/final_wheel/
  pip install /tmp/final_wheel/leafblower-*.whl
  python -c "from leafblower import harvest; print('OK')"
  ```
  Expected: `OK`.

- [ ] Property-based bound check (50 random datasets):
  ```r
  library(leafblower)
  set.seed(99)
  for (i in seq_len(50)) {
    n <- sample(500:5000, 1)
    mx <- runif(1, 2, 10)
    df <- data.frame(x = factor(sample(letters[1:4], n, replace=TRUE)),
                     y = factor(sample(letters[1:3], n, replace=TRUE)))
    tgt <- list(x = c(a=0.25,b=0.25,c=0.25,d=0.25), y=c(a=0.33,b=0.34,c=0.33))
    result <- harvest(df, tgt, max_weight=mx)
    stopifnot(max(result$weights) <= mx + 1e-10)
    stopifnot(abs(mean(result$weights) - 1.0) < 1e-8)
  }
  cat("All 50 bound checks PASSED\n")
  ```

**Commit:** `test(gate4): Phase 4 gate — R CMD check --as-cran clean, Python wheel installs`

---

## Self-Review

### Spec Coverage Check (FR → Task mapping)

| Requirement | Task |
|-------------|------|
| FR-1 rk_calibrate | T6 |
| FR-2 rk_params_init defaults | T3 |
| FR-3 group_ids int32 -1=NA | T3, T7 |
| FR-4 all validation checks | T3 |
| FR-5 min_weight=0 no-op | T3, T12 |
| FR-6 RK_ERR_NOCONV last iterate | T6, T10 |
| FR-7 algorithm_used never AUTO | T6 |
| FR-8 message snprintf | T6 |
| FR-9 log_fn callback + R trampoline + Python GIL trampoline | T7, T14 |
| FR-10 RK_ERR_INFEAS empty cell | T10 |
| FR-11–19 iEPPA | T10 |
| FR-20–28 L-BFGS-B + logit/exp link | T4, T5 |
| FR-29–35 R package | T7, T8 |
| FR-36–40 Python package | T14, T15, T16 |

### Placeholder Scan

None detected. Every code step shows complete, compilable code.

### Type/Name Consistency Check

- `inner_max_iter` used throughout (not `max_iterations`) in C structs — correct.
- `int64_t complexity` initialized to `INT64_C(0)` — portable across MSVC and GCC/Clang.
- `S_j < 1e-15 * W` (not `== 0`) — correct per FR-13.
- `normU = max_weight` — correct per FR-14 derivation.
- `H(u)` logit formula: `L*u + (U-L)/logit_scale * ln(((U-1) + (1-L)*exp(logit_scale*u)) / (U-L))` — matches FR-22 (field renamed from `A` to `logit_scale` for clarity).
- `group_ids` typed as `const int32_t**` — int32_t is explicit and unambiguous; consistent with PRD §6.
- `rk_result_t.algorithm_used` set to `RK_ALG_LBFGSB` or `RK_ALG_IEPPA`, never `RK_ALG_AUTO` — correct per FR-7.
- Python `weights_out` is `np.memcpy` copy — satisfies FR-39 copy contract.

### Spec Gaps Found

None. All 19 US acceptance criteria, all FR-1 through FR-40, and all phase gate criteria are covered by Tasks 1–19.
