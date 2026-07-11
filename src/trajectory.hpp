// STUDY-BRANCH-ONLY-DO-NOT-MERGE
#ifndef LBW_TRAJECTORY_HPP
#define LBW_TRAJECTORY_HPP

// Generic (non-ORIS) cold-path trajectory diagnostics, mirroring the
// oris_trajectory.cpp facility (parse_trajectory_iters/write_trajectory_csv)
// but under distinct names (traj_*) so this can live alongside it without a
// symbol collision in the `lbw` namespace. Every leafblower solver can use
// this to emit per-iteration convergence CSVs for RQ3 curves.
//
// Driven by env vars:
//   LBW_TRAJECTORY_AT   — comma-separated iters to sample, e.g. "1,2,5,10".
//   LBW_TRAJECTORY_OUT  — CSV output path. No env => no parsing, no file I/O
//                          (zero cost on the hot path; parsed/written once
//                          per solve, never per-iteration).
//   LBW_TRAJECTORY_KIND — overrides the metric column name passed by the
//                          caller (fallback "metric" if neither is set).

#include <deque>
#include <utility>
#include <vector>

namespace lbw {

// Parses LBW_TRAJECTORY_AT into a sorted-unique vector of target iterations.
// Returns {} when the env var is unset/empty.
std::vector<int> traj_parse_iters();

// Writes samples to LBW_TRAJECTORY_OUT as "iter,<metric_col>\n" rows.
// No-op when samples is empty or the env var is unset. metric_col is the
// caller-supplied column name (e.g. "errRp", "dual_gap"); LBW_TRAJECTORY_KIND,
// if set, overrides it.
void traj_write_csv(const std::vector<std::pair<int, double>>& samples,
                     const char* metric_col);

// Shared probe-queue-advance logic used by every solver's convergence-check
// block: pops all queue entries <= iter, pushing exactly one (iter, metric)
// sample per queue-front reached (matches oris's "nearest check >= target"
// semantics). No-op (queue untouched) when the queue is empty or iter hasn't
// reached the next target — the common case when LBW_TRAJECTORY_AT is unset.
inline void traj_record(std::deque<int>& queue, int iter, double metric,
                         std::vector<std::pair<int, double>>& samples) {
    if (queue.empty() || iter < queue.front()) return;
    bool first = true;
    while (!queue.empty() && iter >= queue.front()) {
        if (first) { samples.emplace_back(iter, metric); first = false; }
        queue.pop_front();
    }
}

}  // namespace lbw

#endif  // LBW_TRAJECTORY_HPP
