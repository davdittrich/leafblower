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
    double capacity_mu_auto = 0.0;  // auto-computed ALM default: M_cell/n (set by build_cell_table)
};

// Build cell table from group_ids. Returns 0 on success, -1 on invalid dims (K > 64, or K <= 0 / n <= 0).
// Caller owns all input memory; function allocates output vectors.
// NA (group_ids[k][i] = -1) is treated as category `cat_counts[k]` internally
// (distinct bucket per margin).
int build_cell_table(int n, int K,
                     const int32_t* const* group_ids,
                     const int* cat_counts,
                     const double* weights,
                     CellTable& out);

// Quick O(n) cell count estimate without building the full CellTable.
// Returns the capped distinct-cell product, clamped to n (so it returns n once the
// product of cat_counts reaches n — e.g. the K > 8 fast-exit path). Used for AUTO routing.
int estimate_M_cell(int n, int K,
                    const int32_t* const* group_ids,
                    const int* cat_counts);

// Maximum K supported (prevents unbounded key allocation).
inline constexpr int K_MAX = 64;

/// Build cells_per_cat[k][j] = sorted list of cell indices c where g_per_cell[k][c] == j.
/// K and cat_counts must match the CellTable that was built from the same inputs.
std::vector<std::vector<std::vector<int>>>
build_cells_per_cat(const CellTable& ct, int K, const int* cat_counts);

} // namespace lbw
