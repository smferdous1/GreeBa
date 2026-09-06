#pragma once

#include "greedy_b_matching.h"
#include <utility>
#include <vector>

// Appendix B.1-B.2 of Ferdous, Algorithms for Degree-Constrained Subgraphs
// and Applications. Every vertex in adj has matching capacity ONE.
struct BMatchingReduction {
    std::vector<std::vector<std::pair<int, int>>> adj;
    std::vector<std::vector<int>> copies;
    std::vector<std::pair<int, int>> gadgets; // (p_e,u, p_e,v), in edge order
    std::vector<MatchingEdge> original_edges;
    long long baseline_weight = 0; // sum of original edge weights
};

// Each input edge is supplied once. Supports general loopless graphs,
// including parallel edges as distinct selectable edges. Weights and
// capacities must be nonnegative; vertex IDs index capacity.
BMatchingReduction reduce_b_matching(
    const std::vector<MatchingEdge>& edges, const std::vector<int>& capacity);

// mate[v] is -1 or the matched neighbor of v in the reduced graph.
// Validates the 1-matching and recovers edges with TWO matched outer edges.
// Single outer edges can be replaced by their middle edges at equal weight;
// this normalization is implicit and does not modify mate.
// Optimality requires a maximum-WEIGHT reduced matching, not just maximality.
MatchingResult recover_b_matching(
    const BMatchingReduction& reduction, const std::vector<int>& mate);

// Exact simple b-matching through the reduction and a bipartite 1-matching
// solver. U = [0,n_u), V = [n_u,adj.size()). Only U rows are read; reverse
// adjacency entries are optional. Parallel entries are distinct input edges.
MatchingResult reduction_bipartite_b_matching(
    const std::vector<std::vector<std::pair<int, int>>>& adj,
    int n_u, const std::vector<int>& capacity);

MatchingResult reduction_bipartite_b_matching(
    const std::vector<std::vector<std::pair<int, int>>>& adj,
    int n_u, int capacity);
