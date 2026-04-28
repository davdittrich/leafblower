#include "cell_table.hpp"
#include <algorithm>
#include <numeric>
#include <unordered_set>

namespace lbw {

namespace {
inline int bits_needed(int n_vals) {
    int bits = 0;
    int vv = n_vals - 1;
    while (vv > 0) { bits++; vv >>= 1; }
    return (bits < 1) ? 1 : bits;
}
}

// Pack (g_1, ..., g_K) into a 64-bit key when feasible.
// Layout: low-order bits = margin 0, higher bits = later margins.
// NA (-1) encoded as cat_counts[k].
static bool pack_key_fits(int K, const int* cat_counts) {
    if (K > 8) return false;
    uint64_t total_bits = 0;
    for (int k = 0; k < K; k++) {
        total_bits += bits_needed(cat_counts[k] + 1);
        if (total_bits > 64) return false;
    }
    return true;
}

static uint64_t pack_key_compute(int K,
                                  const int32_t* const* group_ids, int i,
                                  const int* cat_counts,
                                  const int* bit_widths) {
    uint64_t key = 0;
    int shift = 0;
    for (int k = 0; k < K; k++) {
        int g = group_ids[k][i];
        int encoded = (g == -1) ? cat_counts[k] : g;
        key |= (uint64_t)encoded << shift;
        shift += bit_widths[k];
    }
    return key;
}

int build_cell_table(int n, int K,
                     const int32_t* const* group_ids,
                     const int* cat_counts,
                     const double* weights,
                     CellTable& out) {
    if (K > K_MAX) return -1;
    if (K <= 0 || n <= 0) return -1;

    // Build keys for each observation.
    const bool use_packed = pack_key_fits(K, cat_counts);

    // Indices sorted by key.
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);

    if (use_packed) {
        std::vector<int> bit_widths(K);
        for (int k = 0; k < K; k++) bit_widths[k] = bits_needed(cat_counts[k] + 1);
        std::vector<uint64_t> keys(n);
        for (int i = 0; i < n; i++)
            keys[i] = pack_key_compute(K, group_ids, i, cat_counts, bit_widths.data());
        std::sort(idx.begin(), idx.end(),
                  [&](int a, int b) { return keys[a] < keys[b]; });
        // Scan to identify cells.
        out.cell_of.resize(n);
        out.n_per_cell.clear();
        out.g_per_cell.assign(K, std::vector<int>());
        int current_cell = -1;
        uint64_t current_key = 0;
        for (int r = 0; r < n; r++) {
            int i = idx[r];
            uint64_t k = keys[i];
            if (current_cell == -1 || k != current_key) {
                current_cell++;
                current_key = k;
                out.n_per_cell.push_back(0);
                for (int m = 0; m < K; m++) {
                    int g = group_ids[m][i];
                    int encoded = (g == -1) ? cat_counts[m] : g;
                    out.g_per_cell[m].push_back(encoded);
                }
            }
            out.cell_of[i] = current_cell;
            out.n_per_cell[current_cell]++;
        }
        out.M_cell = current_cell + 1;
    } else {
        // General path: sort by tuple directly.
        auto tuple_less = [&](int a, int b) {
            for (int k = 0; k < K; k++) {
                int ga = group_ids[k][a];
                int gb = group_ids[k][b];
                int ea = (ga == -1) ? cat_counts[k] : ga;
                int eb = (gb == -1) ? cat_counts[k] : gb;
                if (ea != eb) return ea < eb;
            }
            return false;
        };
        std::sort(idx.begin(), idx.end(), tuple_less);
        out.cell_of.resize(n);
        out.n_per_cell.clear();
        out.g_per_cell.assign(K, std::vector<int>());
        int current_cell = -1;
        auto tuples_equal = [&](int a, int b) {
            for (int k = 0; k < K; k++) {
                int ga = group_ids[k][a];
                int gb = group_ids[k][b];
                int ea = (ga == -1) ? cat_counts[k] : ga;
                int eb = (gb == -1) ? cat_counts[k] : gb;
                if (ea != eb) return false;
            }
            return true;
        };
        int prev_i = -1;
        for (int r = 0; r < n; r++) {
            int i = idx[r];
            if (current_cell == -1 || !tuples_equal(i, prev_i)) {
                current_cell++;
                out.n_per_cell.push_back(0);
                for (int m = 0; m < K; m++) {
                    int g = group_ids[m][i];
                    int encoded = (g == -1) ? cat_counts[m] : g;
                    out.g_per_cell[m].push_back(encoded);
                }
            }
            out.cell_of[i] = current_cell;
            out.n_per_cell[current_cell]++;
            prev_i = i;
        }
        out.M_cell = current_cell + 1;
    }

    out.W_input = 0.0;
    for (int i = 0; i < n; i++) out.W_input += weights[i];
    return 0;
}

int estimate_M_cell(int n, int K,
                    const int32_t* const* group_ids,
                    const int* cat_counts) {
    if (K <= 0 || n <= 0) return 0;
    if (K > 8) {
        int64_t prod = 1;
        for (int k = 0; k < K; k++) {
            prod *= static_cast<int64_t>(cat_counts[k]);
            if (prod >= static_cast<int64_t>(n)) return n;
        }
        return static_cast<int>(prod);
    }
    if (!pack_key_fits(K, cat_counts)) return n;

    std::vector<int> bit_widths(K);
    for (int k = 0; k < K; k++) bit_widths[k] = bits_needed(cat_counts[k] + 1);

    std::unordered_set<uint64_t> seen;
    seen.reserve(static_cast<size_t>(std::min(n, 1 << 20)));  // caps allocation; count remains exact
    for (int i = 0; i < n; i++)
        seen.insert(pack_key_compute(K, group_ids, i, cat_counts, bit_widths.data()));
    return static_cast<int>(seen.size());
}

} // namespace lbw
