#include "reduction_b_matching.h"
#include "anstee_b_matching.h"
#include "milp_b_matching.h"
#include "GreeBa/genGraph.h"
#include "b_matching_preprocessing.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <tuple>

namespace {
void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

template<class F> void rejects(F f) {
    bool rejected = false;
    try { f(); } catch (const std::logic_error&) { rejected = true; }
    check(rejected, "invalid input was accepted");
}

MatchingResult brute_force(const std::vector<MatchingEdge>& edges,
                           const std::vector<int>& b) {
    MatchingResult best;
    best.degree.assign(b.size(), 0);
    for (unsigned mask = 0; mask < (1u << edges.size()); ++mask) {
        MatchingResult current;
        current.degree.assign(b.size(), 0);
        bool feasible = true;
        for (size_t i = 0; i < edges.size(); ++i) {
            if (!(mask & (1u << i))) continue;
            const auto& e = edges[i];
            if (++current.degree[e.u] > b[e.u] || ++current.degree[e.v] > b[e.v]) {
                feasible = false;
                break;
            }
            current.edges.push_back(e);
            current.totalWeight += e.weight;
        }
        if (feasible && current.totalWeight > best.totalWeight) best = current;
    }
    return best;
}

void validate_result(const MatchingResult& actual,
                     const std::vector<MatchingEdge>& edges,
                     const std::vector<int>& b) {
    auto remaining = edges;
    std::vector<int> degree(b.size(), 0);
    long long weight = 0;
    for (const auto& e : actual.edges) {
        const auto it = std::find_if(remaining.begin(), remaining.end(),
            [&](const auto& x) { return std::tie(x.u, x.v, x.weight) ==
                                       std::tie(e.u, e.v, e.weight); });
        check(it != remaining.end(), "edge selected more than once or absent");
        remaining.erase(it);
        ++degree[e.u];
        ++degree[e.v];
        weight += e.weight;
    }
    check(degree == actual.degree && weight == actual.totalWeight,
          "incorrect result metadata");
    for (size_t v = 0; v < b.size(); ++v)
        check(degree[v] <= b[v], "capacity exceeded");
}

void verify(const std::vector<MatchingEdge>& edges, int n_u,
            const std::vector<int>& b) {
    std::vector<std::vector<std::pair<int, int>>> adj(b.size());
    for (const auto& e : edges) adj[e.u].emplace_back(e.v, e.weight);
    const auto reduced = reduction_bipartite_b_matching(adj, n_u, b);
    const auto anstee = anstee_bipartite_b_matching(adj, n_u, b);
    const auto milp = milp_bipartite_b_matching(adj, n_u, b);
    const auto optimum = brute_force(edges, b).totalWeight;
    check(reduced.totalWeight == optimum, "reduction disagrees with exhaustive optimum");
    check(anstee.totalWeight == optimum, "Anstee disagrees with exhaustive optimum");
    check(milp.totalWeight == optimum, "MILP disagrees with exhaustive optimum");
    validate_result(reduced, edges, b);
    validate_result(anstee, edges, b);
    validate_result(milp, edges, b);
}

void verify_generated(int n, int b, int M) {
    const auto [adj, n_u, n_v] = make_greeba_bpt(n, b, M);
    const auto effective = clamp_bipartite_capacities(adj, n_u, b);
    const auto anstee = anstee_bipartite_b_matching(adj, n_u, effective);
    const auto reduced = reduction_bipartite_b_matching(adj, n_u, effective);
    const auto milp = milp_bipartite_b_matching(adj, n_u, effective);
    if (anstee.totalWeight != reduced.totalWeight || anstee.totalWeight != milp.totalWeight) {
        std::cerr << "Generated graph n=" << n << " b=" << b << " M=" << M
                  << ": Anstee=" << anstee.totalWeight
                  << " reduction=" << reduced.totalWeight
                  << " MILP=" << milp.totalWeight << '\n';
        throw std::runtime_error("generated graph weights differ");
    }
    std::vector<MatchingEdge> edges;
    for (int u = 0; u < n_u; ++u)
        for (const auto& [v, w] : adj[u]) edges.push_back({u, v, w});
    const std::vector<int> capacity(adj.size(), b);
    validate_result(anstee, edges, capacity);
    validate_result(reduced, edges, capacity);
    validate_result(milp, edges, capacity);
}

void verify_multiplicity(const std::vector<MatchingEdge>& edges, int n_u,
                         const std::vector<int>& b) {
    std::vector<std::vector<std::pair<int, int>>> adj(b.size());
    std::vector<MatchingEdge> expanded;
    for (const auto& e : edges) {
        adj[e.u].emplace_back(e.v, e.weight);
        for (int i = 0; i < std::min(b[e.u], b[e.v]); ++i) expanded.push_back(e);
    }
    const auto actual = anstee_bipartite_b_matching(adj, n_u, b, false);
    const auto milp = milp_bipartite_b_matching(adj, n_u, b, false);
    check(actual.totalWeight == brute_force(expanded, b).totalWeight,
          "Anstee multiplicity mode disagrees with exhaustive optimum");
    validate_result(actual, expanded, b);
    check(milp.totalWeight == actual.totalWeight, "MILP multiplicity mode disagrees with optimum");
    validate_result(milp, expanded, b);
}
} // namespace

int main() {
    try {
        verify_generated(4, 1, 1); // formerly Anstee=0, reduction=1
        verify_generated(8, 1, 2); // formerly Anstee=4, reduction=5
        verify_generated(24, 2, 100);
        for (int b : {1, 2, 3, 4})
            for (int M : {1, 2, 3, 100})
                for (int blocks : {1, 2, 3, 5})
                    for (int offset : {-1, 0, 1})
                        verify_generated(blocks * (2*b*b + 2*b) + offset, b, M);

        // Heavy edge beats the two-edge maximum-cardinality matching.
        verify({{0, 2, 100}, {0, 3, 1}, {1, 2, 1}}, 2, {1, 1, 1, 1});
        verify({}, 0, {});
        verify({}, 1, {2, 0});
        verify({{0, 1, 9}}, 1, {0, 3});
        verify({{0, 1, 9}}, 1, {3, 3}); // one input edge stays single
        verify({{0, 1, 9}, {0, 1, 7}}, 1, {2, 2}); // parallel edges
        verify({{0, 1, 0}}, 1, {1, 1});
        verify({{0, 2, 2000000000}, {1, 3, 2000000000}}, 2, {1, 1, 1, 1});
        std::mt19937 rng(42);
        for (int trial = 0; trial < 300; ++trial) {
            std::vector<int> b(6);
            for (int& cap : b) cap = rng() % 4;
            std::vector<MatchingEdge> edges;
            for (int u = 0; u < 3; ++u)
                for (int v = 3; v < 6; ++v)
                    if (rng() % 3) edges.push_back({u, v, static_cast<int>(rng() % 20)});
            verify(edges, 3, b);
        }

        verify_multiplicity({{0, 1, 9}}, 1, {3, 3});
        // The total supply can exceed INT_MAX even for a tiny simple graph.
        check(anstee_bipartite_b_matching({{{1, 7}}, {}}, 1,
                  std::vector<int>(2, std::numeric_limits<int>::max())).totalWeight == 7,
              "large aggregate capacity failed");
        verify_multiplicity({{0, 2, 100}, {0, 3, 1}, {1, 2, 1}}, 2, {2, 2, 2, 2});
        for (int trial = 0; trial < 100; ++trial) {
            std::vector<int> b(4);
            for (int& cap : b) cap = rng() % 3;
            std::vector<MatchingEdge> edges;
            for (int u = 0; u < 2; ++u)
                for (int v = 2; v < 4; ++v)
                    if (rng() % 3) edges.push_back({u, v, static_cast<int>(rng() % 20)});
            verify_multiplicity(edges, 2, b);
        }

        const auto r = reduce_b_matching({{0, 1, 7}}, {2, 1});
        check(r.adj.size() == 5 && r.baseline_weight == 7, "wrong gadget size");
        const auto [pu, pv] = r.gadgets[0];
        std::vector<int> mate(r.adj.size(), -1);
        check(recover_b_matching(r, mate).edges.empty(), "empty matching recovery");
        mate[pu] = pv; mate[pv] = pu;
        check(recover_b_matching(r, mate).edges.empty(), "middle edge recovery");
        mate.assign(r.adj.size(), -1);
        mate[pu] = r.copies[0][0]; mate[r.copies[0][0]] = pu;
        check(recover_b_matching(r, mate).edges.empty(), "single outer edge recovery");
        mate[pv] = r.copies[1][0]; mate[r.copies[1][0]] = pv;
        check(recover_b_matching(r, mate).totalWeight == 7, "two outer edge recovery");

        // General graph construction/recovery: triangle, solved exhaustively.
        const auto triangle = reduce_b_matching({{0, 1, 8}, {1, 2, 5}, {0, 2, 4}}, {1, 1, 1});
        std::vector<MatchingEdge> expanded;
        for (int u = 0; u < static_cast<int>(triangle.adj.size()); ++u)
            for (const auto& [v, w] : triangle.adj[u])
                if (u < v) expanded.push_back({u, v, w});
        const auto optimum = brute_force(expanded, std::vector<int>(triangle.adj.size(), 1));
        mate.assign(triangle.adj.size(), -1);
        for (const auto& e : optimum.edges) { mate[e.u] = e.v; mate[e.v] = e.u; }
        const auto recovered = recover_b_matching(triangle, mate);
        check(recovered.totalWeight == 8 &&
              optimum.totalWeight == triangle.baseline_weight + recovered.totalWeight,
              "general graph optimum or baseline identity failed");

        rejects([] { reduce_b_matching({}, {-1}); });
        rejects([] { milp_bipartite_b_matching({{}}, 1, std::vector<int>{-1}); });
        rejects([] { milp_bipartite_b_matching({{}}, 1, std::vector<int>{}); });
        rejects([] { milp_bipartite_b_matching({}, -1, 1); });
        rejects([] { milp_bipartite_b_matching({}, 1, 1); });
        rejects([] { milp_bipartite_b_matching({}, 0, -1); });
        rejects([] { milp_bipartite_b_matching({{{0, 1}}}, 1, 1); });
        rejects([] { milp_bipartite_b_matching({{{2, 1}}, {}}, 1, 1); });
        rejects([] { milp_bipartite_b_matching({{{1, -1}}, {}}, 1, 1); });
        rejects([] { anstee_bipartite_b_matching({{}}, 1, std::vector<int>{-1}); });
        rejects([] { reduce_b_matching({{0, 1, -1}}, {1, 1}); });
        rejects([] { reduce_b_matching({{0, 0, 1}}, {1}); });
        rejects([] { reduce_b_matching({{0, 2, 1}}, {1, 1}); });
        rejects([] { reduction_bipartite_b_matching({{}}, 2, 1); });
        rejects([] { reduction_bipartite_b_matching({{}}, 1, std::vector<int>{}); });
        rejects([] { reduction_bipartite_b_matching({{{0, 1}}}, 1, 1); });
        rejects([] { reduction_bipartite_b_matching({}, 0, -1); });
        rejects([&] { recover_b_matching(r, {}); });
        rejects([&] { recover_b_matching(r, {3, -1, -1, -1, -1}); });
        rejects([&] { recover_b_matching(r, {1, 0, -1, -1, -1}); });
        std::cout << "Anstee, reduction, and MILP agree: 195 generated-graph comparisons "
                     "and 300 random exhaustive comparisons.\n"
                     "Anstee and MILP also pass 100 random exhaustive multiplicity comparisons.\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
