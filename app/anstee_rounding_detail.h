#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace anstee_detail {

// Optional structural counters for regression tests, independent of timing.
// Each event costs O(1), apart from amortized O(1) adjacency/list insertion.
struct RoundingWork {
    std::size_t vertex_checks = 0;
    std::size_t edge_checks = 0;
    std::size_t adjacency_probes = 0;
    std::size_t traversed_edges = 0;
    std::size_t open_trails = 0;
    std::size_t closed_trails = 0;
};

// Stage 2b on a feasible, optimal half-integral bipartite b-matching.
// Edge must have vertex fields u and v; parallel edges retain distinct indices.
// x2 stores twice each edge value in int64_t without narrowing.
// Integral capacities/bounds give each odd endpoint half a unit of slack and
// each fractional edge room for either rounding sign. Alternating on a maximal
// open trail or an even closed trail is therefore feasible in both directions.
// Optimality implies zero objective change for either sign; use -1,+1,... .
template<class Edge>
void round_half_integral(int n, const std::vector<Edge>& edges,
                         std::vector<int64_t>& x2,
                         RoundingWork* work = nullptr) {
    assert(n >= 0 && x2.size() == edges.size());
    if (work) *work = {};

    // Every fractional edge contributes two adjacency entries. Cursors never
    // move backwards; an active bit retires both entries in constant time.
    std::vector<std::vector<std::pair<int, std::size_t>>> adjacency(n);
    std::vector<std::size_t> degree(n, 0), cursor(n, 0);
    std::vector<unsigned char> active(edges.size(), 0);
    for (std::size_t e = 0; e < edges.size(); ++e) {
        if (work) ++work->edge_checks;
        assert(x2[e] >= 0);
        if (x2[e] % 2 == 0) continue;
        const auto& edge = edges[e];
        adjacency[edge.u].emplace_back(edge.v, e);
        adjacency[edge.v].emplace_back(edge.u, e);
        ++degree[edge.u];
        ++degree[edge.v];
        active[e] = 1;
    }

    std::vector<int> odd_vertices;
    odd_vertices.reserve(n);
    for (int v = 0; v < n; ++v) {
        if (work) ++work->vertex_checks;
        if (degree[v] % 2) odd_vertices.push_back(v);
    }

    // Consume a maximal trail and apply its alternating changes as we go.
    // No edge is revisited, and no temporary trail needs to be copied/stored.
    auto round_trail = [&](int start, bool closed) {
        int current = start;
        std::size_t length = 0;
        for (;;) {
            auto& next = cursor[current];
            const auto& incident = adjacency[current];
            while (next < incident.size()) {
                if (work) ++work->adjacency_probes;
                if (active[incident[next].second]) break;
                ++next;
            }
            if (next == incident.size()) break;
            const auto [neighbor, e] = incident[next++];
            active[e] = 0;
            --degree[current];
            --degree[neighbor];
            x2[e] += (length % 2 == 0) ? -1 : 1;
            ++length;
            if (work) ++work->traversed_edges;
            current = neighbor;
        }
        assert(length > 0 && degree[current] == 0);
        assert(closed ? (current == start && length % 2 == 0)
                      : (current != start));
        if (work) {
            if (closed) ++work->closed_trails;
            else ++work->open_trails;
        }
    };

    // A completed open trail makes its two odd endpoints even; internal
    // vertices lose even degree. Thus no new odd vertices need enqueuing.
    // An entry whose vertex was already used as an endpoint is simply skipped.
    for (int v : odd_vertices) {
        if (work) ++work->vertex_checks;
        if (degree[v] % 2) round_trail(v, false);
    }

    // Remaining degrees are even. A maximal trail closes at its start and
    // exhausts that vertex. A single vertex sweep finds all remaining trails,
    // including cycles attached to vertices inside earlier open/closed trails.
    for (int v = 0; v < n; ++v) {
        if (work) ++work->vertex_checks;
        if (degree[v] != 0) round_trail(v, true);
    }

    // For h fractional edges: <= 3*n vertex checks, m edge checks,
    // <= 2*h adjacency probes and exactly h traversals. Allocations and
    // initialization cost O(n+m), giving O(n+m) time and auxiliary space.
    // This bound covers rounding only, not the flow solver or expansion of
    // every selected multiplicity unit into the public result's edge list.
}

} // namespace anstee_detail
