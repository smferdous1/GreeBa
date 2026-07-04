#pragma once

#include <vector>

struct MatchingEdge {
    int u;
    int v;
    int weight;
};

struct MatchingResult {
    std::vector<MatchingEdge> edges;
    long long totalWeight = 0;
    std::vector<int> degree;
};

MatchingResult greedy_weighted_b_matching(
    const std::vector<std::vector<std::pair<int, int>>>& adj,
    int capacity
);

MatchingResult greedy_weighted_b_matching(
    const std::vector<std::vector<std::pair<int, int>>>& adj,
    const std::vector<int>& capacity
);
