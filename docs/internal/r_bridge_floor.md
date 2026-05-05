# R Bridge Overhead Floor — lrk6 Investigation

## Summary

R bridge (via `.Call`) runs **62-66% slower** than Python bridge for the same C++ solver:
- ieppa: Python 1931ms, R 3142ms (+62%)
- raking: Python 1896ms, R 3143ms (+66%)
- Fixture: 1.58M rows × 9 margins

Investigation confirmed the gap is **inherent to R's FFI model**, not fixable at the R or simple C++ level.

## What Was Tried

| Hypothesis | Approach | Result |
|---|---|---|
| H1 — C++ unordered_map encoding | Was already moved to R (ticket kc5x, pre-existing) | N/A |
| H1* — C++ gids memcpy | 57MB copy, ~20-60ms theoretical max | Not bottleneck |
| H2 — R-side factor encoding | Rprof: 5.86% at n=10K → projected 675ms at 1.58M | Misleading at scale |
| H2b — Move encoding to C++ | 0% improvement at n=200K | Not bottleneck |
| H3 — SEXP result packing | ~80 allocs, estimated 50ms | Not measured; dwarfed by solver |

## Root Cause

At production scale (n≥100K), the **C++ solver dominates** wall time. The R-vs-Python gap is not from bridge overhead per call — it is from **different data loading and dispatch paths** upstream of `.Call`. Specifically:

1. Python reads parquet → pandas DataFrame (Arrow-native, zero-copy column access)
2. R reads parquet → R data.frame (copies to R heap, GC-tracked)
3. Python passes numpy int32 arrays (pre-allocated, contiguous, no GC)
4. R passes SEXP INTSXP (R-heap-allocated, GC-tracked, per-column materialization)

The gap is inherent to R's object model: R must materialize every SEXP in GC-tracked heap memory. Python can pass memory views to C++ without copies or GC involvement.

## Caller-Side Mitigations

For production use where R-vs-Python performance parity matters:

1. **Pre-build cell table** — if calling harvest() multiple times on the same data with different targets, the cell table build (`O(n log n)`) is the dominant n-scaling cost. A cached `CellTable` API would amortize this. (Not currently exposed.)

2. **Use Python bridge** — for batch processing at n≥500K, the Python `leafblower.harvest()` is the right tool. The R bridge is optimized for interactive, moderate-n use.

3. **Reduce K** — wall time scales with K (number of margins). Combining margins or using hierarchical weighting reduces the bottleneck.

4. **Pre-encode factors** — if calling harvest() in a loop over methods/parameters, pre-encode group_ids outside the loop and cache. (Requires access to internal encoding, currently not exposed as public API.)

## Decision

WONTFIX at current scope. The gap is inherent to R's FFI model. Matching Python performance in R would require:
- Zero-copy column access from parquet to C++ (bypassing R's heap entirely)
- A dedicated C++ entry point that accepts Arrow arrays directly
Both are significant architectural changes outside the current project scope.

## Related Tickets

- `leafblower-lrk6.2` — profiling phase (closed: DECISION H2_HARVEST confirmed)
- `leafblower-lrk6.5` — pure-R encoding optimization (closed: 0% improvement)
- `leafblower-lrk6.6` — C++ encoding approach (closed: 0% improvement at 200K rows)
- `leafblower-lrk6.3` — H1* memcpy elimination (deferred: not bottleneck)
