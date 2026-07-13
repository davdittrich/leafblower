// STUDY-BRANCH-ONLY-DO-NOT-MERGE
// Generic per-solver cold-path trajectory diagnostics (leafblower-2ouc.42).
// Copied from oris_trajectory.cpp's parse_trajectory_iters/write_trajectory_csv
// logic under distinct names (traj_*) so every solver can emit RQ3 convergence
// curves without disturbing ORIS's own copy. Never on the hot path — parsed/
// written once per solve.

#include "trajectory.hpp"

#include <algorithm>
#include <cstdlib>
#include <exception>
#include <cstdio>
#include <string>

namespace lbw {

std::vector<int> traj_parse_iters() {
    const char* s = std::getenv("LBW_TRAJECTORY_AT");
    if (!s || !*s) return {};
    std::vector<int> out;
    std::string buf;
    for (const char* p = s;; ++p) {
        if (*p == ',' || *p == '\0') {
            if (!buf.empty()) {
                try {
                    out.push_back(std::stoi(buf));
                } catch (const std::exception&) {
                    // skip malformed token
                }
                buf.clear();
            }
            if (*p == '\0') break;
        } else buf.push_back(*p);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

void traj_write_csv(const std::vector<lbw::TrajSample>& samples,
                     const char* metric_col) {
    if (samples.empty()) return;
    const char* path = std::getenv("LBW_TRAJECTORY_OUT");
    if (!path || !*path) return;
    // STUDY-BRANCH-ONLY-DO-NOT-MERGE: C stdio, NOT std::ofstream. The Python
    // _leafblower.so statically links libstdc++, so iostream/locale globals in
    // the module collide with the host process libstdc++ and SIGSEGV on the
    // first operator<<. fopen/fprintf touch no C++ locale/iostream state.
    std::FILE* f = std::fopen(path, "w");
    if (!f) return;
    const char* kind_env = std::getenv("LBW_TRAJECTORY_KIND");
    const char* col = (kind_env && *kind_env) ? kind_env
                      : (metric_col && *metric_col) ? metric_col : "metric";
    std::fprintf(f, "iter,%s,marginal_kl\n", col);
    for (const auto& s : samples)
        std::fprintf(f, "%d,%.17g,%.17g\n", s.iter, s.natural, s.marginal_kl);
    std::fclose(f);
}

}  // namespace lbw
