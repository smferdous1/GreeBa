#include "b_matching_preprocessing.h"
#include "anstee_b_matching.h"
#include "greedy_b_matching.h"
#include "milp_b_matching.h"
#include "reduction_b_matching.h"
#include "GreeBa/genGraph.h"

#include <climits>
#include <iostream>
#include <random>
#include <tuple>

namespace {
using Adj = std::vector<std::vector<std::pair<int, int>>>;
void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
template<class F> void rejects(F call) {
    bool rejected = false;
    try { call(); } catch (const std::logic_error&) { rejected = true; }
    check(rejected, "invalid preprocessing input accepted");
}

void validate(const MatchingResult& result, const std::vector<MatchingEdge>& edges,
              const std::vector<int>& capacity) {
    auto remaining = edges;
    std::vector<int> degree(capacity.size(), 0);
    long long weight = 0;
    for (const auto& e : result.edges) {
        const auto it = std::find_if(remaining.begin(), remaining.end(), [&](const auto& x) {
            return std::tie(e.u, e.v, e.weight) == std::tie(x.u, x.v, x.weight);
        });
        check(it != remaining.end(), "selected an absent edge or reused an input edge");
        remaining.erase(it);
        check(++degree[e.u] <= capacity[e.u] && ++degree[e.v] <= capacity[e.v],
              "matching exceeds original capacity");
        weight += e.weight;
    }
    check(degree == result.degree && weight == result.totalWeight, "incorrect matching metadata");
}

void verify(const std::vector<MatchingEdge>& edges, int n_u,
            const std::vector<int>& original) {
    Adj adj(original.size());
    std::vector<int> degree(original.size(), 0);
    for (const auto& e : edges) {
        adj[e.u].push_back({e.v, e.weight});
        ++degree[e.u]; ++degree[e.v];
    }
    const auto input_before = adj;
    const auto capacity_before = original;
    const auto effective = clamp_bipartite_capacities(adj, n_u, original);
    for (std::size_t v = 0; v < original.size(); ++v)
        check(effective[v] == std::min(original[v], degree[v]), "incorrect capacity projection");
    check(clamp_bipartite_capacities(adj, n_u, effective) == effective, "not idempotent");
    check(adj == input_before && original == capacity_before, "preprocessing mutated input");
    auto symmetric = adj;
    for (const auto& e : edges) symmetric[e.v].push_back({e.u, e.weight});
    check(clamp_bipartite_capacities(symmetric, n_u, original) == effective,
          "reverse adjacency changed clamping");

    // Independently enumerate every simple edge subset. Clamping must preserve
    // the entire feasible set, not just the optimum found by another solver.
    long long optimum = 0;
    for (unsigned mask = 0; mask < (1u << edges.size()); ++mask) {
        std::vector<int> used(original.size(), 0);
        long long weight = 0;
        for (std::size_t e = 0; e < edges.size(); ++e) {
            if (!(mask & (1u << e))) continue;
            ++used[edges[e].u]; ++used[edges[e].v];
            weight += edges[e].weight;
        }
        bool before = true, after = true;
        for (std::size_t v = 0; v < used.size(); ++v) {
            before = before && used[v] <= original[v];
            after = after && used[v] <= effective[v];
        }
        check(before == after, "clamping changed the feasible set");
        if (before) optimum = std::max(optimum, weight);
    }
    const auto greedy = greedy_weighted_b_matching(adj, effective);
    check(greedy.totalWeight == greedy_weighted_b_matching(adj, original).totalWeight,
          "clamping changed greedy weight");
    validate(greedy, edges, original);
    for (const auto& graph : {adj, symmetric}) {
        for (const auto& result : {anstee_bipartite_b_matching(graph, n_u, effective),
                                  reduction_bipartite_b_matching(graph, n_u, effective),
                                  milp_bipartite_b_matching(graph, n_u, effective)}) {
            check(result.totalWeight == optimum, "preprocessed solver differs from original optimum");
            validate(result, edges, original);
        }
    }
}
} // namespace

int main() {
    try {
        verify({}, 0, {});
        verify({}, 2, {INT_MAX, 16, 0, INT_MAX});
        verify({{0, 1, 7}}, 1, {INT_MAX, INT_MAX});
        verify({{0, 1, 9}, {0, 1, 9}, {0, 1, 3}}, 1, {INT_MAX, 2});
        verify({{0, 2, 100}, {0, 3, 1}, {1, 2, 1}}, 2, {1, INT_MAX, 1, INT_MAX});
        verify({{0, 2, 0}, {0, 3, 7}, {1, 3, 4}}, 2, {0, INT_MAX, INT_MAX, 16});
        std::mt19937 rng(20260906);
        const int capacities[] = {0, 1, 2, 16, INT_MAX};
        for (int trial = 0; trial < 100; ++trial) {
            const int n_u = 1 + rng() % 3, n_v = 1 + rng() % 3;
            std::vector<int> b(n_u + n_v);
            for (int& value : b) value = capacities[rng() % 5];
            std::vector<MatchingEdge> edges;
            const int m = rng() % 9;
            for (int e = 0; e < m; ++e)
                edges.push_back({static_cast<int>(rng() % n_u),
                                 n_u + static_cast<int>(rng() % n_v),
                                 static_cast<int>(rng() % 20)});
            verify(edges, n_u, b);
        }

        // A concrete size regression proves the reduction actually receives
        // fewer copies, with the graph generator's b left at its original 16.
        const auto [adj, n_u, n_v] = make_greeba_bpt(1, 16, 100);
        const auto before = adj;
        const auto effective = clamp_bipartite_capacities(adj, n_u, 16);
        std::vector<MatchingEdge> edges;
        for (int u = 0; u < n_u; ++u)
            for (auto [v, w] : adj[u]) edges.push_back({u, v, w});
        const auto reduced = reduce_b_matching(edges, effective);
        check(adj == before && adj.size() == 544 && edges.size() == 768,
              "clamping changed the generated graph");
        check(reduced.adj.size() == 2560, "reduced vertex count did not shrink");
        std::size_t reduced_edges_twice = 0;
        for (const auto& row : reduced.adj) reduced_edges_twice += row.size();
        check(reduced_edges_twice == 2 * 17664, "reduced edge count did not shrink");

        const Adj one_edge{{{1, 7}}, {}};
        check(clamp_bipartite_capacities(one_edge, 1, INT_MAX) == std::vector<int>({1, 1}),
              "scalar overload failed");
        // Raw solver APIs still support edge multiplicity: callers do not
        // apply this simple-graph preprocessing when simple=false.
        check(anstee_bipartite_b_matching(one_edge, 1, 3, false).totalWeight == 21 &&
              milp_bipartite_b_matching(one_edge, 1, 3, false).totalWeight == 21,
              "edge-reuse behavior changed");
        rejects([&] { clamp_bipartite_capacities(one_edge, INT_MIN, 1); });
        rejects([&] { clamp_bipartite_capacities(one_edge, 3, 1); });
        rejects([&] { clamp_bipartite_capacities(one_edge, 1, -1); });
        rejects([&] { clamp_bipartite_capacities(one_edge, 1, std::vector<int>{1}); });
        rejects([&] { clamp_bipartite_capacities(one_edge, 1, std::vector<int>{0, -1}); });
        for (int v : {-1, 0, 2, INT_MAX})
            rejects([&] { clamp_bipartite_capacities(Adj{{{v, 7}}, {}}, 1, 0); });
        std::cout << "Capacity preprocessing passed: 106 feasible-set/solver comparisons, "
                     "optional reverse rows, parallel edges, isolates, INT_MAX, and graph-size regression.\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
