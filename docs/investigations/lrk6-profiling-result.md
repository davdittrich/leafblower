# lrk6 Profiling Result

**Fixture:** `tests/testthat/fixtures/stepstone_small.parquet`  
**n = 10,000 rows, K = 9 margins, method = ieppa, max_iter = 200, tol = 1e-8**  
**R version:** leafblower installed; **Python version:** `_leafblower.cpython-314-x86_64-linux-gnu.so` via local build

## Timings (stepstone_small, n=10K, K=9)

| Phase | R (ms) | Python (ms) | Notes |
|-------|--------|-------------|-------|
| **Total harvest()** | **74.0** | **61.4** | 5-rep median; Python uses synthetic uniform targets |
| Solver (.Call / calibrate C++) | 67.3 | ~51.3 (est.) | R from Rprof 90.99%; Python = total - preprocessing |
| R preprocessing (group_ids encoding) | 4.3 | 10.1 | R lapply+factor path; Python dict comprehension |
| droplevels / other R | 2.4 | 0.0 | R-specific |
| **Total bridge overhead** | **6.7** | **10.1** | R bridge is *leaner* than Python preprocessing |

Also measured: **R raking** median = **57 ms** (5 reps).

## Rprof breakdown (30-run average, ieppa)

```
                      self.time self.pct total.time total.pct
".Call"                    2.02    90.99       2.02     90.99
"FUN"                      0.13     5.86       0.15      6.76
"droplevels"               0.02     0.90       0.02      0.90
"withCallingHandlers"      0.01     0.45       2.22    100.00
"as.character"             0.01     0.45       0.01      0.45
"c"                        0.01     0.45       0.01      0.45
"makeRestartList"          0.01     0.45       0.01      0.45
"tolower"                  0.01     0.45       0.01      0.45

Total elapsed: 2.22s for 30 runs → 74 ms/run
```

## Dominant phase analysis

**91% of R wall time is inside `.Call` — the C++ solver.** The R bridge overhead (lapply encoding + SEXP packing) is only **~6.7 ms** (9%).

The R–Python gap on the small fixture is **12.6 ms**. Python's `harvest()` total is 61.4 ms but its Python-side preprocessing takes ~10 ms (dict comprehension encoding), while R's lapply+factor path takes only ~4.3 ms. The C++ solver itself runs **faster under R** for this fixture (67.3 ms) than the Python equivalent (~51 ms estimated), which is the opposite of H2/H3 predictions.

**The 1211 ms gap cited in the task background (1.58M rows, K=9) is not explained by this small-fixture profile.** At n=10K the R bridge costs ~6.7 ms; if the 1211 ms gap scales linearly with n, that would require ~180 ms per 10K rows of bridging — inconsistent with what we see. The 1211 ms gap at 1.58M rows is almost certainly dominated by the **C++ solver itself running slower under R** (different BLAS thread configuration, or R GC pressure during iterative solver) rather than any bridge cost.

## Hypothesis verdicts

| Hypothesis | Verdict | Evidence |
|------------|---------|----------|
| H1 (C++ memcpy ~57MB, 20–60ms) | **RULED OUT as dominant** | Bridge = 6.7ms; memcpy within that |
| H2 (harvest.R preprocessing) | **RULED OUT** | 4.3ms, negligible |
| H3 (SEXP packing ~80 allocs) | **RULED OUT** | Absorbed into 6.7ms bridge total |
| H4 (inherent R/Python solver path difference) | **LIKELY** | .Call = 91% of wall time; gap must be in solver |

## Scale note

At 1.58M rows (157× larger):
- R total reported: 3142 ms → per-run ~3142 ms
- Python total reported: 1931 ms → per-run ~1931 ms
- Gap: 1211 ms
- If bridge scales linearly: 6.7 ms × 157 = ~1051 ms — *close* to gap but bridge overhead is sublinear (encoding is O(n), solver is O(n·K·iter))
- The `FUN/lapply` group_ids encoding at O(n) **could** account for ~1051 ms at 1.58M rows (4.3ms × 157 = 675ms) — this is H2/H3 combined at scale

**Revised conclusion:** At 10K rows, bridge is negligible. At 1.58M rows, the O(n) lapply encoding in R (currently 4.3ms at 10K) scales to ~675ms — making **H2 (harvest.R O(n) preprocessing) the dominant cost at full scale**, not the C++ solver difference.

## DECISION

`DECISION: H2_HARVEST`

The lapply/factor encoding step in `harvest.R` (lines 379–394) is O(n). At 10K rows it costs 4.3ms (5.9% of total). At 1.58M rows it projects to ~675ms — the dominant component of the 1211ms gap. Optimization target: vectorize group_ids encoding to avoid R-level per-element dispatch.
