#include "anstee_b_matching.h"

#include "ortools/graph/min_cost_flow.h"

#include <cassert>
#include <stdexcept>
#include <unordered_set>
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
    const int n = static_cast<int>(adj.size());
    const int n_v = n - n_u;
    if (n_u < 0 || n_v < 0)
        throw std::invalid_argument("n_u out of range");
    if (static_cast<int>(b.size()) != n)
        throw std::invalid_argument("b.size() must equal adj.size()");

    // Collect edges (u in U, v in V) without duplicates.
    struct RawEdge { int u, v, w; };
    std::vector<RawEdge> edges;
    for (int u = 0; u < n_u; ++u) {
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
    // (Costs are negated so that min-cost max-flow maximises weight.)
    std::vector<int> arc_fwd(m, -1), arc_rev(m, -1);
    for (int e = 0; e < m; ++e) {
        const int cap = simple ? 1 : std::min(b[edges[e].u], b[edges[e].v]);
        if (cap == 0) continue;
        arc_fwd[e] = smcf.AddArcWithCapacityAndUnitCost(R(edges[e].u), S(edges[e].v), cap, -edges[e].w);
        arc_rev[e] = smcf.AddArcWithCapacityAndUnitCost(R(edges[e].v), S(edges[e].u), cap, -edges[e].w);
    }

    // Supply/demand arcs: s->R(x) and S(x)->t, each with capacity b[x].
    int total_supply = 0;
    for (int x = 0; x < n; ++x) {
        if (b[x] > 0) {
            smcf.AddArcWithCapacityAndUnitCost(s, R(x), b[x], 0);
            smcf.AddArcWithCapacityAndUnitCost(S(x), t, b[x], 0);
            total_supply += b[x];
        }
    }
    smcf.SetNodeSupply(s, total_supply);
    smcf.SetNodeSupply(t, -total_supply);

    if (smcf.SolveMaxFlowWithMinCost() != operations_research::SimpleMinCostFlow::OPTIMAL)
        throw std::runtime_error("OR-Tools SimpleMinCostFlow did not reach OPTIMAL");

    // ---- Stage 2a: symmetrize to x2[e] = 2*x_{uv} (always an integer) ----
    std::vector<int> x2(m, 0);
    for (int e = 0; e < m; ++e) {
        int a_fwd = (arc_fwd[e] >= 0) ? static_cast<int>(smcf.Flow(arc_fwd[e])) : 0;
        int a_rev = (arc_rev[e] >= 0) ? static_cast<int>(smcf.Flow(arc_rev[e])) : 0;
        x2[e] = a_fwd + a_rev;
    }

    // ---- Stage 2b: resolve half-integral edges via alternating trails ----
    // H = {edges with x2[e] odd}.  We work with x2 scaled by 2 throughout:
    // alternation subtracts 1 from even-indexed trail positions and adds 1 to
    // odd-indexed ones.  After resolution every x2[e] is even (= 2 * x_{uv}).

    // Build adjacency for H.
    std::vector<std::vector<std::pair<int, int>>> H_adj(n); // vertex -> [(nbr, edge_idx)]
    std::unordered_set<int> H_rem;
    for (int e = 0; e < m; ++e) {
        if (x2[e] % 2 == 1) {
            H_adj[edges[e].u].emplace_back(edges[e].v, e);
            H_adj[edges[e].v].emplace_back(edges[e].u, e);
            H_rem.insert(e);
        }
    }

    auto deg_in_H = [&](int v) {
        int d = 0;
        for (auto& [nbr, eidx] : H_adj[v])
            if (H_rem.count(eidx)) ++d;
        return d;
    };

    // Walk a maximal trail from `start`, removing traversed edges from H_rem.
    auto find_trail = [&](int start) {
        std::vector<int> trail;
        int cur = start;
        for (;;) {
            int found_e = -1, found_nbr = -1;
            for (auto& [nbr, eidx] : H_adj[cur]) {
                if (H_rem.count(eidx)) { found_e = eidx; found_nbr = nbr; break; }
            }
            if (found_e < 0) break;
            trail.push_back(found_e);
            H_rem.erase(found_e);
            cur = found_nbr;
        }
        return trail;
    };

    // Anstee eq. (7): 0-indexed position i even -> x2 -= 1, odd -> x2 += 1.
    auto apply_alt = [&](const std::vector<int>& trail) {
        for (int i = 0; i < static_cast<int>(trail.size()); ++i) {
            if (i % 2 == 0) x2[trail[i]] -= 1;
            else             x2[trail[i]] += 1;
        }
    };

    // Step 1: eliminate odd-degree vertices in H via maximal trails.
    // Each iteration finds one odd-degree vertex and walks until stuck.
    // By Eulerian trail theory, the walk must end at another odd-degree vertex,
    // reducing the count of odd-degree vertices by 2 per iteration.
    for (;;) {
        int odd_v = -1;
        for (int v = 0; v < n && odd_v < 0; ++v)
            if (deg_in_H(v) & 1) odd_v = v;
        if (odd_v < 0) break;
        auto trail = find_trail(odd_v);
        if (trail.empty()) break;
        apply_alt(trail);
    }

    // Step 2: H now has all even degrees; decompose into closed trails.
    // Bipartite H guarantees all closed trails have even length.
    while (!H_rem.empty()) {
        int e0 = *H_rem.begin();
        auto trail = find_trail(edges[e0].u);
        assert(trail.size() % 2 == 0); // even-length closed trail in bipartite H
        apply_alt(trail);
    }

    // ---- Build result ----
    MatchingResult result;
    result.degree.assign(n, 0);
    for (int e = 0; e < m; ++e) {
        assert(x2[e] >= 0 && x2[e] % 2 == 0);
        const int xval = x2[e] / 2;
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
    return anstee_bipartite_b_matching(
        adj, n_u, std::vector<int>(adj.size(), capacity), simple);
}
