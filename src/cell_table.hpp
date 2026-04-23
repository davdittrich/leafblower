#pragma once
#include <cstdint>
#include <vector>

namespace lbw {

struct CellTable {
    int M_cell;                                  // number of unique cells
    std::vector<int> cell_of;                    // size n: obs i -> cell index
    std::vector<int> n_per_cell;                 // size M_cell: count per cell
    std::vector<std::vector<int>> g_per_cell;    // [K][M_cell]: margin-k cat per cell
    double W_input;                              // sum of input weights
};

// Build cell table from group_ids. Returns 0 on success, -1 if K > 64.
// Caller owns all input memory; function allocates output vectors.
// NA (group_ids[k][i] = -1) is treated as category `cat_counts[k]` internally
// (distinct bucket per margin).
int build_cell_table(int n, int K,
                     const int32_t* const* group_ids,
                     const int* cat_counts,
                     const double* weights,
                     CellTable& out);

// Maximum K supported (prevents unbounded key allocation).
inline constexpr int K_MAX = 64;

} // namespace lbw
