#include "greedy_b_matching.h"

#include <algorithm>
#include <stdexcept>

namespace {

std::vector<MatchingEdge> collect_unique_edges(
    const std::vector<std::vector<std::pair<int, int>>>& adj
) {
    std::vector<MatchingEdge> edges;

    for (int u = 0; u < static_cast<int>(adj.size()); ++u) {
        for (const auto& [v, weight] : adj[u]) {
            if (v < 0 || v >= static_cast<int>(adj.size())) {
                throw std::out_of_range("edge endpoint is outside the graph");
            }
            if (u < v) {
                edges.push_back({u, v, weight});
            }
        }
    }

    std::sort(edges.begin(), edges.end(), [](const MatchingEdge& left, const MatchingEdge& right) {
        if (left.weight != right.weight) {
            return left.weight > right.weight;
        }
        if (left.u != right.u) {
            return left.u < right.u;
        }
        return left.v < right.v;
    });

    return edges;
}

} // namespace

MatchingResult greedy_weighted_b_matching(
    const std::vector<std::vector<std::pair<int, int>>>& adj,
    int capacity
) {
    if (capacity < 0) {
        throw std::invalid_argument("capacity must be non-negative");
    }

    return greedy_weighted_b_matching(adj, std::vector<int>(adj.size(), capacity));
}

MatchingResult greedy_weighted_b_matching(
    const std::vector<std::vector<std::pair<int, int>>>& adj,
    const std::vector<int>& capacity
) {
    if (capacity.size() != adj.size()) {
        throw std::invalid_argument("capacity vector size must match graph vertex count");
    }
    if (std::any_of(capacity.begin(), capacity.end(), [](int value) { return value < 0; })) {
        throw std::invalid_argument("capacities must be non-negative");
    }

    const auto candidateEdges = collect_unique_edges(adj);
    MatchingResult result;
    result.degree.assign(adj.size(), 0);

    for (const auto& edge : candidateEdges) {
        if (result.degree[edge.u] >= capacity[edge.u] ||
            result.degree[edge.v] >= capacity[edge.v]) {
            continue;
        }

        result.edges.push_back(edge);
        result.totalWeight += edge.weight;
        ++result.degree[edge.u];
        ++result.degree[edge.v];
    }

    return result;
}
