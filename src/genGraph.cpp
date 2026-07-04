#include <vector>
#include <tuple>
#include <algorithm>
#include<iostream>
#include <cmath>
#include "GreeBa/genGraph.h"

std::tuple<std::vector<std::vector<std::pair<int, int>>>, int, int> make_greeba_bpt(int n, int b, int M) {
    int denom = 2 * b * b + 2 * b;
    n = std::ceil(n / (double)denom) * denom;
    int n_cpy = n / denom;
    std::vector<std::tuple<int, int, int>> edgelist;
    int n_u = 0;
    int n_v = 0;
    int t_nv = 0;
    std::vector<int> outerU, outerV;
    for (int i = 0; i < n_cpy; ++i) {
        for (int j = 0; j < b; ++j) {
            for (int k = 0; k <b; ++k) {
                edgelist.push_back({j+n_u, k+n_v, M});
            }
        }
        n_u = n_u + b;
        n_v = n_v + b;
        t_nv = n_v; //temporary n_v for correct index in the second block
        if(i % 2 ==1) outerV.clear();
        for (int j = 0; j < b; ++j) {
            for (int k = 0; k < b; ++k) {
                edgelist.push_back({n_u - b + j, n_v, M - 1});
                if (i % 2 == 1) outerV.push_back(n_v);
                ++n_v;
            }
        }
        if (i%2==0) outerU.clear();
        for (int j = 0; j < b; ++j) {
            for (int k = 0; k < b; ++k) {
                edgelist.push_back({n_u, t_nv- b + j, M - 1});
                if (i % 2 == 0) outerU.push_back(n_u);
                ++n_u;
            }
        }
        if (i > 0) {
            for (size_t l = 0; l < outerU.size(); ++l) {
                edgelist.push_back({outerU[l], outerV[l], 1});
            }
        }
        
    }
    std::vector<std::vector<std::pair<int, int>>> adj(n_u+n_v);
    for (auto& e : edgelist) {
        auto [u, v, w] = e;
        adj[u].push_back({n_v+v, w});
        adj[n_v+v].push_back({u, w});
    }
    return std::make_tuple(adj, n_u, n_v);
}