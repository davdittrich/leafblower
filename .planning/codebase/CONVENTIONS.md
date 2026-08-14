# Coding Conventions

**Analysis Date:** 2026-08-15

## Naming Patterns

**Files:**
- C++ source/header: snake_case (`raking.cpp`, `oris_trajectory.cpp`, `types.hpp`)
- R files: snake_case (`harvest.R`, `design_effect.R`)
- Python modules: snake_case with leading underscore for internal modules (`_harvest.py`, `_design_effect.py`)

**Functions:**
- C++: snake_case (`raking_solve`, `build_cells_per_cat`, `compute_cell_metrics`)
  - Static/file-local: still snake_case, documented with inline comments
  - Exported C API: snake_case prefixed with namespace (`calibrate`, `harvest`)
- R: snake_case (`harvest`, `diagnose_weights`, `design_effect`)
- Python: snake_case with leading underscore for internal (`_parse_convergence`, `_validate_pos_scalar`)

**Variables:**
- C++: snake_case (`X_init`, `L_cell`, `max_cats`, `kl_ratio_scratch`)
- R: snake_case (df, tgt, result, weights)
- Python: snake_case (weights_out, target_mass, sparse_cats)

**Types and Structures:**
- C++: PascalCase for structs/enums (`CalibResult`, `CalibMetric`, `HomotopyConfigLbw`, `RakingResult`)
- Enum values: SCREAMING_SNAKE_CASE (`MAX_ERR`, `RK_OK`, `RK_ERR_NOCONV`)
- Constants: SCREAMING_SNAKE_CASE with `k` prefix (`kAbsoluteZeroThreshold`, `kErrCheckInterval`, `kRakingFeasTol`, `kMinSafeTotalWeight`)
  - Static/constexpr constants in function scope: k-prefixed
  - Global configuration: lowercase with `k` prefix

**Type Hints:**
- Python: use type hints consistently
  - Function signatures: `def _parse_convergence(conv: Optional[Dict]) -> Tuple[...]`
  - Variables: `Optional[float]`, `Dict[str, list]`, `frozenset`

## Code Style

**Formatting:**
- No enforced formatter (clang-format absent)
- Manual style adherence based on existing code:
  - C++: 4-space indentation, `namespace lbw {` opening braces on same line
  - R: 2-space indentation (roxygen2 generated)
  - Python: 4-space indentation (PEP 8 standard)

**Linting:**
- No linter configuration files (`.eslintrc`, `.pylintrc` absent)
- Manual code review for:
  - C++: namespace correctness, const-correctness, RAII patterns
  - Python: type annotation completeness, docstring coverage

## Import Organization

**C++:**
1. `lbw_config.h` (configuration)
2. `lbw_math.hpp` (math utilities)
3. Standard library headers (`<cmath>`, `<vector>`, `<algorithm>`)
4. Project headers (`leafblower.h`, `types.hpp`, `calib_dispatch.hpp`)

Example from `raking.cpp`:
```cpp
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
```

**Python:**
1. Future imports (`from __future__ import annotations`)
2. Standard library (math, warnings, os, sys)
3. Third-party (numpy, pandas)
4. Local imports (`from ._leafblower import calibrate`)

Example from `_harvest.py`:
```python
from __future__ import annotations
import math
import warnings
import numpy as np
import pandas as pd
from typing import Dict, Optional

from ._leafblower import calibrate
```

**R:**
- Uses roxygen2 for documentation
- Package dependencies listed in DESCRIPTION
- No explicit import organization at file level

## Error Handling

**C++:**
- Return status codes from functions: `RK_OK`, `RK_ERR_NOCONV`, `RK_ERR_BUDGET`, `RK_ERR_INFEAS`
- CalibResult struct contains status field for all solver results
- Example from `raking_solve()`:
  ```cpp
  RakingResult res;
  res.base.status     = RK_ERR_BUDGET;  // initial; overwritten by criterion/stall
  if (lbw::solver_setup_ct_base(st, ct, X_init, hi_eff, L_cell, U_cell, res) != RK_OK)
      return res;
  ```
- Guard clauses for bounds validation early in functions

**Python:**
- Raise `ValueError` with descriptive messages for invalid inputs
- Mirroring R's validation: type, range, and domain checks
- Example from `_validate_pos_scalar()`:
  ```python
  if not math.isfinite(val) or val <= 0 or (check_upper and val > 1e15):
      raise ValueError(
          f"{name} must be {null_clause}a positive finite scalar; got: {val!r}"
      )
  ```
- Use `warnings.warn()` for soft failures (e.g., below recommended range)

**R:**
- R errors thrown via `stop()` for fatal validation failures
- Attributes attached to result object for diagnostics: `attr(result, "result")`
- Example:
  ```r
  if (status != 0L) {
    stop("Calibration failed: ", error_message)
  }
  ```

## Logging

**Framework:** console output via stdio

**C++:**
- `CalibState::log()` method for verbose output
- Uses `log_fn` callback if set, else prints to R/C stderr via `Rprintf` or `fprintf`
- Verbose levels: 0 (silent), 1 (progress), 2+ (debug)
- Example:
  ```cpp
  void log(const char* msg) const {
      if (verbose <= 0) return;
      if (log_fn) {
          log_fn(msg, log_ctx);
      }
  }
  ```

**R:**
- `message()` and `warning()` for user-facing output
- Verbose flags passed to C++ via `verbose` parameter
- Example: `tryCatch(capture.output(...), type = "output")`

**Python:**
- `warnings.warn()` for deprecations and soft errors
- Test output via `print()` or pytest's `-s` flag
- Example from tests:
  ```python
  warnings.warn(f"... below recommended range; ...", UserWarning, stacklevel=2)
  ```

## Comments

**When to Comment:**
- Algorithm description at file start (`raking.cpp` lines 1-15: References to Deming, Csiszar)
- Complex mathematical invariants (e.g., "X_tilde[c] = X_init[c] * exp(sum_k lf[k][g_k(c)])")
- Known workarounds or deliberate simplifications marked with `ponytail:` prefix
- Bug/ticket references (e.g., "leafblower-8eod", "CR-C5b", "eb79.1")
- State machine transitions and status codes

**JSDoc/Roxygen2:**
- R functions: comprehensive roxygen2 documentation (`#' @param`, `#' @details`, `#' @return`)
- Python: docstrings for public functions with parameter descriptions
- C++: no formal documentation format; comments in code
- Example from `harvest.R`:
  ```r
  #' Generate calibrated weights (drop-in for autumn::harvest)
  #'
  #' @param data A data frame containing columns to calibrate.
  #' @param target A named list of named numeric vectors (variable -> proportions).
  ```

## Function Design

**Size:** 
- Functions stay under ~300 lines; longer solvers (`oris.cpp` 2246L) split cold logic (`oris_trajectory.cpp`, `oris_finalize.cpp`)
- Inner hot loops (per-iteration) co-located with caller for inlining without LTO

**Parameters:**
- C++: pass large structures by reference (`CalibState& st`), const-correct
- R: use named lists for configuration (`convergence = list(...)`)
- Python: use dicts for configuration, type hints on all params

**Return Values:**
- C++: return struct by value (`RakingResult`, `OrisSoftResult`); base fields at `.base.*`
- R: return vector or list with attributes (`attr(result, "algorithm")`, `attr(result, "result")`)
- Python: return dict or numpy array with metadata keys

## Module Design

**Exports:**
- C++: exported from `leafblower.h` (C API only)
- R: roxygen2 `@export` for public functions
- Python: bare functions in `__init__.py`, internal functions with `_` prefix

**Barrel Files:**
- Python: `__init__.py` re-exports public API
- R: namespace imports via roxygen2 (no explicit barrel files)

**Shared Logic Location:**
- C++: `calib_dispatch.hpp` for solver-shared helpers (finalize_weights, solver_setup_ct_base)
- C++: `cell_table.hpp` for cell-aggregation logic
- Python: `_harvest.py` for parameter parsing and validation

**Single Responsibility:**
- One solver per file (`raking.cpp`, `oris.cpp`, `newton_calib.cpp`)
- One algorithm per handler (not monolithic dispatch)
- Internal helpers marked `static` or in anonymous namespaces

---

*Convention analysis: 2026-08-15*
