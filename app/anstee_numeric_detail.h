#pragma once

#include "ortools/graph/min_cost_flow.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>

// Internal numeric operations shared with bounded-memory regression tests.
namespace anstee_numeric_detail {

constexpr std::size_t index_limit = static_cast<std::size_t>(
    std::min(std::numeric_limits<int>::max(),
             std::numeric_limits<int32_t>::max()));

inline int checked_partition_size(std::size_t n, int n_u) {
    if (n_u < 0 || static_cast<std::size_t>(n_u) > n)
        throw std::invalid_argument("n_u out of range");
    // The Hitchcock network has 2*n+2 nodes. SimpleMinCostFlow adds two
    // more nodes when checking feasibility; its node counts are int32_t.
    if (n > (index_limit - 4) / 2)
        throw std::length_error("graph too large for Anstee flow node indices");
    return static_cast<int>(n);
}

// n has already passed checked_partition_size. Reserve two arcs per vertex,
// one bypass arc and the solver's two supply/demand arcs before edge arcs.
inline std::size_t max_input_edges(int n) {
    return (index_limit - 3) / 2 - static_cast<std::size_t>(n);
}

inline int64_t symmetrized_flow(
    const operations_research::SimpleMinCostFlow& flow, int arc_fwd, int arc_rev)
{
    // Each flow can reach INT_MAX; the addition itself must be 64-bit.
    const int64_t a_fwd = arc_fwd >= 0 ? flow.Flow(arc_fwd) : 0;
    const int64_t a_rev = arc_rev >= 0 ? flow.Flow(arc_rev) : 0;
    return a_fwd + a_rev;
}

inline int checked_multiplicity(int64_t x2) {
    if (x2 < 0 || x2 / 2 > std::numeric_limits<int>::max())
        throw std::overflow_error("Anstee multiplicity outside int range");
    if (x2 % 2 != 0)
        throw std::logic_error("Anstee multiplicity is not integral");
    return static_cast<int>(x2 / 2);
}

} // namespace anstee_numeric_detail
