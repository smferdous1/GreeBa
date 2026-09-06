#include "anstee_b_matching.h"
#include "anstee_numeric_detail.h"
#include "anstee_rounding_detail.h"

#include "ortools/graph/min_cost_flow.h"

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <vector>

// ---------------------------------------------------------------------------
// Anstee bipartite b-matching (per-vertex capacity vector).
// ---------------------------------------------------------------------------
MatchingResult anstee_bipartite_b_matching(
    const std::vector<std::vector<std::pair<int, int>>>& adj,
    int n_u,
    const std::vector<int>& b,
    bool simple)
{
    const int n = anstee_numeric_detail::checked_partition_size(adj.size(), n_u);
    if (b.size() != adj.size())
        throw std::invalid_argument("b.size() must equal adj.size()");
    if (std::any_of(b.begin(), b.end(), [](int value) { return value < 0; }))
        throw std::invalid_argument("capacities must be non-negative");

    // Read each U-side entry once; reverse adjacency entries are ignored.
    struct RawEdge { int u, v, w; };
    std::vector<RawEdge> edges;
    const std::size_t max_edges = anstee_numeric_detail::max_input_edges(n);
    for (int u = 0; u < n_u; ++u) {
        if (adj[u].size() > max_edges - edges.size())
            throw std::length_error("graph too large for Anstee flow arc indices");
        for (auto& [v, w] : adj[u]) {
            if (v < n_u || v >= n)
                throw std::out_of_range("edge references vertex outside V");
            if (w < 0)
                throw std::invalid_argument("negative edge weights not supported");
            edges.push_back({u, v, w});
        }
    }
    const int m = static_cast<int>(edges.size());

    // ---- Stage 1: Hitchcock network + OR-Tools SimpleMinCostFlow ----
    // Node layout: 0=s, 1=t, 2+i=R(i), 2+n+i=S(i)  for i in 0..n-1.
    const int s = 0, t = 1;
    auto R = [&](int i) { return 2 + i; };
    auto S = [&](int i) { return 2 + n + i; };

    operations_research::SimpleMinCostFlow smcf;

    // Two arcs per edge: R(u)->S(v) with cost -w and R(v)->S(u) with cost -w.
    // Costs are negated so that minimum cost maximises matching weight.
    std::vector<int> arc_fwd(m, -1), arc_rev(m, -1);
    for (int e = 0; e < m; ++e) {
        const int cap = simple ? 1 : std::min(b[edges[e].u], b[edges[e].v]);
        if (cap == 0) continue;
        arc_fwd[e] = smcf.AddArcWithCapacityAndUnitCost(R(edges[e].u), S(edges[e].v), cap, -edges[e].w);
        arc_rev[e] = smcf.AddArcWithCapacityAndUnitCost(R(edges[e].v), S(edges[e].u), cap, -edges[e].w);
    }

    // Supply/demand arcs: s->R(x) and S(x)->t, each with capacity b[x].
    int64_t total_supply = 0;
    for (int x = 0; x < n; ++x) {
        if (b[x] > 0) {
            smcf.AddArcWithCapacityAndUnitCost(s, R(x), b[x], 0);
            smcf.AddArcWithCapacityAndUnitCost(S(x), t, b[x], 0);
            total_supply += b[x];
        }
    }
    // Route unused capacity directly to the sink at zero cost. Without this
    // bypass, SolveMaxFlowWithMinCost prioritizes cardinality, which can force
    // two light edges in place of one heavier edge. Fix the total flow and
    // let minimum cost choose how much passes through actual graph edges.
    smcf.AddArcWithCapacityAndUnitCost(s, t, total_supply, 0);
    smcf.SetNodeSupply(s, total_supply);
    smcf.SetNodeSupply(t, -total_supply);

    if (smcf.Solve() != operations_research::SimpleMinCostFlow::OPTIMAL)
        throw std::runtime_error("OR-Tools SimpleMinCostFlow did not reach OPTIMAL");

    // ---- Stage 2a: symmetrize to x2[e] = 2*x_{uv} (always an integer) ----
    std::vector<int64_t> x2(m, 0);
    for (int e = 0; e < m; ++e) {
        x2[e] = anstee_numeric_detail::symmetrized_flow(smcf, arc_fwd[e], arc_rev[e]);
    }

    // ---- Stage 2b: resolve half-integral edges via alternating trails ----
    anstee_detail::round_half_integral(n, edges, x2);

    // ---- Build result ----
    MatchingResult result;
    result.degree.assign(n, 0);
    for (int e = 0; e < m; ++e) {
        const int xval = anstee_numeric_detail::checked_multiplicity(x2[e]);
        if (xval > 0) {
            const auto& [u, v, w] = edges[e];
            for (int k = 0; k < xval; ++k)
                result.edges.push_back({u, v, w});
            result.totalWeight += static_cast<long long>(xval) * w;
            result.degree[u] += xval;
            result.degree[v] += xval;
        }
    }
    return result;
}

// ---------------------------------------------------------------------------
// Convenience overload: uniform capacity.
// ---------------------------------------------------------------------------
MatchingResult anstee_bipartite_b_matching(
    const std::vector<std::vector<std::pair<int, int>>>& adj,
    int n_u,
    int capacity,
    bool simple)
{
    if (capacity < 0)
        throw std::invalid_argument("capacity must be non-negative");
    anstee_numeric_detail::checked_partition_size(adj.size(), n_u);
    return anstee_bipartite_b_matching(
        adj, n_u, std::vector<int>(adj.size(), capacity), simple);
}
