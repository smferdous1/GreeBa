#pragma once

#include "greedy_b_matching.h"
#include <vector>

// Maximum-weight b-matching on a bipartite graph via Anstee's algorithm
// (Information Processing Letters 24 (1987) 153-157).
//
// adj    : adjacency list over all vertices 0..n-1.
//          U = {0..n_u-1}, V = {n_u..n-1}.
//          Only edges from U to V are processed; both directions may be
//          present in adj (the function reads only u < n_u rows).
// n_u    : size of the left partition.
// b      : per-vertex capacity vector (length == adj.size()).
// simple : if true, each edge may be selected at most once (x_e in {0,1}).
//          if false, x_e in {0..min(b_u, b_v)}.
// Weights and capacities must be nonnegative. Unused capacity is allowed;
// the objective maximizes total weight without a cardinality requirement.
//
// Returns a MatchingResult whose edges list contains one MatchingEdge per
// unit of multiplicity (e.g. if x_{uv}=2, two identical entries appear).

MatchingResult anstee_bipartite_b_matching(
    const std::vector<std::vector<std::pair<int, int>>>& adj,
    int n_u,
    const std::vector<int>& b,
    bool simple = true
);

MatchingResult anstee_bipartite_b_matching(
    const std::vector<std::vector<std::pair<int, int>>>& adj,
    int n_u,
    int capacity,
    bool simple = true
);
