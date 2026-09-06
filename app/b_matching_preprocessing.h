#pragma once

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

// Preprocessing for SIMPLE bipartite b-matching: return min(b[v], degree(v)).
// U = [0,n_u), V = [n_u,adj.size()). Read each U-side entry once and count
// both endpoints; reverse rows are optional and parallel entries are distinct
// edges. The input graph and capacities are not modified. Isolates get zero.
// Generate the graph first, then share this vector among all matching solvers.
// Do not use this degree bound for simple=false (repeated use of an edge).
inline std::vector<int> clamp_bipartite_capacities(
    const std::vector<std::vector<std::pair<int, int>>>& adj,
    int n_u, const std::vector<int>& capacity)
{
    if (adj.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::length_error("graph exceeds integer vertex indexing");
    if (n_u < 0 || static_cast<std::size_t>(n_u) > adj.size())
        throw std::invalid_argument("n_u out of range");
    if (capacity.size() != adj.size())
        throw std::invalid_argument("capacity size must equal vertex count");
    if (std::any_of(capacity.begin(), capacity.end(), [](int b) { return b < 0; }))
        throw std::invalid_argument("capacities must be nonnegative");

    std::vector<int> effective(adj.size(), 0);
    for (int u = 0; u < n_u; ++u) {
        for (const auto& entry : adj[u]) {
            const int v = entry.first;
            if (v < n_u || static_cast<std::size_t>(v) >= adj.size())
                throw std::out_of_range("edge endpoint outside V");
            // Saturating at the supplied capacity avoids overflowing a degree
            // counter, while computing exactly min(capacity[v], degree(v)).
            if (effective[u] < capacity[u]) ++effective[u];
            if (effective[v] < capacity[v]) ++effective[v];
        }
    }
    return effective;
}

inline std::vector<int> clamp_bipartite_capacities(
    const std::vector<std::vector<std::pair<int, int>>>& adj,
    int n_u, int capacity)
{
    if (capacity < 0) throw std::invalid_argument("capacity must be nonnegative");
    return clamp_bipartite_capacities(adj, n_u, std::vector<int>(adj.size(), capacity));
}
