#include "GreeBa/genGraph.h"
#include "greedy_b_matching.h"
#include "anstee_b_matching.h"
#include "reduction_b_matching.h"
#include "milp_b_matching.h"

#include <algorithm>
#include <iostream>

int main() {
    int n = 24;
    int b = 2;
    int M = 100;
    auto [adj, n_u, n_v] = make_greeba_bpt(n, b, M);
    std::cout << "Number of vertices: " << adj.size() << std::endl;
    std::cout << "n_u: " << n_u << ", n_v: " << n_v << std::endl;
    std::cout << "Adjacency list for first 5 vertices:" << std::endl;
    for (size_t i = 0; i < std::min<size_t>(adj.size(), 200); ++i) {
        std::cout << "Vertex " << i << ":";
        for (auto& p : adj[i]) {
            std::cout << " (" << p.first << "," << p.second << ")";
        }
        std::cout << std::endl;
    }

    const auto greedy = greedy_weighted_b_matching(adj, b);
    std::cout << "\nGreedy weighted b-matching:" << std::endl;
    std::cout << "  matched edges: " << greedy.edges.size() << std::endl;
    std::cout << "  total weight:  " << greedy.totalWeight << std::endl;
    std::cout << "  first matched edges:" << std::endl;
    for (size_t i = 0; i < std::min<size_t>(greedy.edges.size(), 10); ++i) {
        const auto& edge = greedy.edges[i];
        std::cout << "    " << edge.u << " -- " << edge.v
                  << " (w=" << edge.weight << ")" << std::endl;
    }

    const auto anstee = anstee_bipartite_b_matching(adj, n_u, b, /*simple=*/true);
    std::cout << "\nAnstee exact b-matching:" << std::endl;
    std::cout << "  matched edges: " << anstee.edges.size() << std::endl;
    std::cout << "  total weight:  " << anstee.totalWeight << std::endl;
    std::cout << "  first matched edges:" << std::endl;
    for (size_t i = 0; i < std::min<size_t>(anstee.edges.size(), 10); ++i) {
        const auto& edge = anstee.edges[i];
        std::cout << "    " << edge.u << " -- " << edge.v
                  << " (w=" << edge.weight << ")" << std::endl;
    }

    std::cout << "\nAnstee weight >= Greedy weight: "
              << (anstee.totalWeight >= greedy.totalWeight ? "YES" : "NO") << std::endl;

    const auto reduced = reduction_bipartite_b_matching(adj, n_u, b);
    std::cout << "\nReduction to exact 1-matching:\n"
              << "  matched edges: " << reduced.edges.size() << '\n'
              << "  total weight:  " << reduced.totalWeight << '\n';
    for (const auto& edge : reduced.edges)
        std::cout << "    " << edge.u << " -- " << edge.v
                  << " (w=" << edge.weight << ")\n";

    const auto milp = milp_bipartite_b_matching(adj, n_u, b);
    std::cout << "\nMILP b-matching (SCIP):\n"
              << "  matched edges: " << milp.edges.size() << '\n'
              << "  total weight:  " << milp.totalWeight << '\n';
    for (const auto& edge : milp.edges)
        std::cout << "    " << edge.u << " -- " << edge.v
                  << " (w=" << edge.weight << ")\n";

    const bool same_weight = anstee.totalWeight == reduced.totalWeight &&
                             anstee.totalWeight == milp.totalWeight;
    std::cout << "\nAnstee weight == Reduction weight == MILP weight: "
              << (same_weight ? "YES" : "NO") << '\n';
    return same_weight ? 0 : 1;
}
