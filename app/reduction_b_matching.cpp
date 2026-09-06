#include "reduction_b_matching.h"

#include "ortools/graph/min_cost_flow.h"
#include <algorithm>
#include <limits>
#include <stdexcept>

BMatchingReduction reduce_b_matching(
    const std::vector<MatchingEdge>& edges, const std::vector<int>& capacity)
{
    const auto max_nodes = static_cast<size_t>(std::numeric_limits<int>::max() - 2);
    if (capacity.size() > max_nodes || edges.size() > max_nodes / 2)
        throw std::length_error("graph exceeds integer vertex indexing");
    size_t nodes = 2 * edges.size();
    for (int b : capacity) {
        if (b < 0) throw std::invalid_argument("capacities must be nonnegative");
        if (static_cast<size_t>(b) > max_nodes - nodes)
            throw std::length_error("reduced graph exceeds integer vertex indexing");
        nodes += b;
    }
    for (const auto& e : edges) {
        if (e.u < 0 || e.v < 0 || static_cast<size_t>(e.u) >= capacity.size() ||
            static_cast<size_t>(e.v) >= capacity.size())
            throw std::out_of_range("edge endpoint outside graph");
        if (e.u == e.v) throw std::invalid_argument("self-loops are not supported");
        if (e.weight < 0) throw std::invalid_argument("weights must be nonnegative");
    }

    BMatchingReduction r;
    r.adj.resize(nodes);
    r.copies.resize(capacity.size());
    r.original_edges = edges;
    int next = 0;
    for (size_t v = 0; v < capacity.size(); ++v)
        for (int j = 0; j < capacity[v]; ++j) r.copies[v].push_back(next++);
    auto add_edge = [&](int u, int v, int w) {
        r.adj[u].emplace_back(v, w);
        r.adj[v].emplace_back(u, w);
    };
    for (const auto& e : edges) {
        const int pu = next++;
        const int pv = next++;
        r.gadgets.emplace_back(pu, pv);
        add_edge(pu, pv, e.weight);
        for (int u : r.copies[e.u]) add_edge(u, pu, e.weight);
        for (int v : r.copies[e.v]) add_edge(v, pv, e.weight);
        r.baseline_weight += e.weight;
    }
    return r;
}

MatchingResult recover_b_matching(
    const BMatchingReduction& r, const std::vector<int>& mate)
{
    if (mate.size() != r.adj.size())
        throw std::invalid_argument("mate size must equal reduced vertex count");
    for (size_t v = 0; v < mate.size(); ++v) {
        const int u = mate[v];
        if (u == -1) continue;
        if (u < 0 || static_cast<size_t>(u) >= mate.size() ||
            mate[u] != static_cast<int>(v))
            throw std::invalid_argument("mate must describe a symmetric 1-matching");
        if (std::none_of(r.adj[v].begin(), r.adj[v].end(),
                        [&](const auto& entry) { return entry.first == u; }))
            throw std::invalid_argument("matched edge is absent from reduced graph");
    }
    MatchingResult result;
    result.degree.assign(r.copies.size(), 0);
    for (size_t i = 0; i < r.gadgets.size(); ++i) {
        const auto [pu, pv] = r.gadgets[i];
        if (mate[pu] == -1 || mate[pv] == -1 || mate[pu] == pv) continue;
        const auto& e = r.original_edges[i];
        result.edges.push_back(e);
        result.totalWeight += e.weight;
        ++result.degree[e.u];
        ++result.degree[e.v];
    }
    return result;
}

MatchingResult reduction_bipartite_b_matching(
    const std::vector<std::vector<std::pair<int, int>>>& adj,
    int n_u, const std::vector<int>& capacity)
{
    if (n_u < 0 || static_cast<size_t>(n_u) > adj.size())
        throw std::invalid_argument("n_u out of range");
    if (capacity.size() != adj.size())
        throw std::invalid_argument("capacity size must equal vertex count");
    std::vector<MatchingEdge> edges;
    for (int u = 0; u < n_u; ++u) {
        for (const auto& [v, w] : adj[u]) {
            if (v < n_u || static_cast<size_t>(v) >= adj.size())
                throw std::out_of_range("edge endpoint outside V");
            edges.push_back({u, v, w});
        }
    }
    const auto r = reduce_b_matching(edges, capacity);
    const int n = static_cast<int>(r.adj.size());
    std::vector<bool> left(n, false);
    for (int u = 0; u < n_u; ++u)
        for (int copy : r.copies[u]) left[copy] = true;
    // U copies -- p_e,u -- p_e,v -- V copies.
    for (const auto& [pu, pv] : r.gadgets) left[pv] = true;

    operations_research::SimpleMinCostFlow flow;
    const int source = n, sink = n + 1;
    int supply = 0;
    for (int v = 0; v < n; ++v) {
        if (left[v]) {
            flow.AddArcWithCapacityAndUnitCost(source, v, 1, 0);
            ++supply;
        } else {
            flow.AddArcWithCapacityAndUnitCost(v, sink, 1, 0);
        }
    }
    struct Arc { int id, u, v; };
    std::vector<Arc> arcs;
    for (int u = 0; u < n; ++u) {
        if (!left[u]) continue;
        for (const auto& [v, w] : r.adj[u])
            arcs.push_back({flow.AddArcWithCapacityAndUnitCost(u, v, 1,
                            -static_cast<long long>(w)), u, v});
    }
    // Unused capacity bypasses the graph. This optimizes weight without
    // forcing maximum cardinality, which can sacrifice heavier edges.
    flow.AddArcWithCapacityAndUnitCost(source, sink, supply, 0);
    flow.SetNodeSupply(source, supply);
    flow.SetNodeSupply(sink, -supply);
    if (flow.Solve() != operations_research::SimpleMinCostFlow::OPTIMAL)
        throw std::runtime_error("reduced 1-matching solver did not reach OPTIMAL");
    std::vector<int> mate(n, -1);
    for (const auto& a : arcs) {
        if (flow.Flow(a.id)) {
            mate[a.u] = a.v;
            mate[a.v] = a.u;
        }
    }
    return recover_b_matching(r, mate);
}

MatchingResult reduction_bipartite_b_matching(
    const std::vector<std::vector<std::pair<int, int>>>& adj, int n_u, int capacity)
{
    if (capacity < 0) throw std::invalid_argument("capacity must be nonnegative");
    return reduction_bipartite_b_matching(adj, n_u,
                                         std::vector<int>(adj.size(), capacity));
}
