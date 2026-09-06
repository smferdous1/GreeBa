#pragma once

#include "greedy_b_matching.h"
#include <utility>
#include <vector>

// Independent maximum-weight bipartite b-matching MILP, solved by SCIP
// through OR-Tools. U = [0,n_u), V = [n_u,adj.size()). Only U rows are read;
// each entry is a distinct edge, so reverse adjacency entries are optional.
// Weights and capacities must be nonnegative integers.
//
// Maximize sum_e w_e*x_e, subject to sum_{e incident to v} x_e <= b_v.
// simple=true: x_e is binary. Otherwise x_e is integer in [0,min(b_u,b_v)].
// Returns only a solver-certified optimum; solver failures throw runtime_error.
// MatchingResult contains one edge entry per unit of selected multiplicity.
MatchingResult milp_bipartite_b_matching(
    const std::vector<std::vector<std::pair<int, int>>>& adj,
    int n_u, const std::vector<int>& capacity, bool simple = true);

MatchingResult milp_bipartite_b_matching(
    const std::vector<std::vector<std::pair<int, int>>>& adj,
    int n_u, int capacity, bool simple = true);
