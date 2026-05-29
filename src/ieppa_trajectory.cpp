// Cold trajectory-diagnostics I/O for iEPPA, split out of ieppa.cpp (uu8r.1).
// Driven by env vars LBW_TRAJECTORY_AT (iters to sample) / LBW_TRAJECTORY_OUT
// (csv path). Never on the hot path — called once per solve. Externalized from
// the former file-static definitions so ieppa_finalize.cpp (uu8r.2) can link to
// write_trajectory_csv.

#include "ieppa_internal.hpp"

#include <algorithm>
#include <cstdlib>
#include <exception>
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


void write_trajectory_csv(
    const std::vector<std::pair<int,double>>& samples)
{
    if (samples.empty()) return;
    const char* path = std::getenv("LBW_TRAJECTORY_OUT");
    if (!path || !*path) return;
    std::ofstream f(path);
    if (!f.is_open()) return;
    f << "iter,errRp\n";
    for (const auto& p : samples)
        f << p.first << "," << p.second << "\n";
}

}  // namespace lbw
