// Cold trajectory-diagnostics I/O for ORIS, split out of oris.cpp (uu8r.1).
// Driven by env vars LBW_TRAJECTORY_AT (iters to sample) / LBW_TRAJECTORY_OUT
// (csv path). Never on the hot path — called once per solve. Externalized from
// the former file-static definitions so oris_finalize.cpp (uu8r.2) can link to
// write_trajectory_csv.

#include "oris_internal.hpp"

#include <algorithm>
#include <cstdlib>
#include <exception>
#include <cstdio>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

namespace lbw {

std::vector<int> parse_trajectory_iters() {
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


void write_trajectory_csv(  // STUDY-BRANCH-ONLY-DO-NOT-MERGE
    const std::vector<lbw::TrajSample>& samples)
{
    if (samples.empty()) return;
    const char* path = std::getenv("LBW_TRAJECTORY_OUT");
    if (!path || !*path) return;
    // STUDY-BRANCH-ONLY-DO-NOT-MERGE: C stdio to avoid the static-libstdc++
    // iostream-global SIGSEGV when this cold CSV writer runs from the Python
    // _leafblower.so. Exit-only, read-only on weights -> solve byte-identical.
    std::FILE* f = std::fopen(path, "w");
    if (!f) return;
    std::fprintf(f, "iter,errRp,marginal_kl\n");
    for (const auto& s : samples)
        std::fprintf(f, "%d,%.17g,%.17g\n", s.iter, s.natural, s.marginal_kl);
    std::fclose(f);
}

}  // namespace lbw
