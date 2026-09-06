#include "anstee_rounding_detail.h"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

namespace {
struct Edge { int u, v, w; };
using Values = std::vector<int64_t>;
std::size_t checked_states = 0;
std::size_t open_trails = 0;
std::size_t closed_trails = 0;

void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

anstee_detail::RoundingWork verify(int n_u, const std::vector<Edge>& edges,
                                 const std::vector<int>& b, Values x2,
                                 bool simple) {
    const auto before = x2;
    check(x2.size() == edges.size(), "fixture size mismatch");
    Values degree(b.size(), 0);
    int64_t weight_before = 0;
    std::size_t fractional = 0;
    for (std::size_t e = 0; e < edges.size(); ++e) {
        const auto& edge = edges[e];
        check(edge.u >= 0 && edge.u < n_u && edge.v >= n_u &&
              edge.v < static_cast<int>(b.size()), "fixture not bipartite");
        const int64_t cap = simple ? 1 : std::min(b[edge.u], b[edge.v]);
        check(x2[e] >= 0 && x2[e] <= 2 * cap, "infeasible fixture edge");
        degree[edge.u] += x2[e];
        degree[edge.v] += x2[e];
        weight_before += x2[e] * edge.w;
        fractional += x2[e] % 2;
    }
    for (std::size_t v = 0; v < b.size(); ++v)
        check(degree[v] <= 2LL * b[v], "infeasible fixture degree");

    // This is the very same helper called by anstee_b_matching.cpp, with
    // instrumentation at its actual scans/traversals, not a model algorithm.
    anstee_detail::RoundingWork work;
    anstee_detail::round_half_integral(static_cast<int>(b.size()), edges, x2, &work);

    std::fill(degree.begin(), degree.end(), 0);
    int64_t weight_after = 0;
    for (std::size_t e = 0; e < edges.size(); ++e) {
        const auto& edge = edges[e];
        const int64_t cap = simple ? 1 : std::min(b[edge.u], b[edge.v]);
        check(x2[e] >= 0 && x2[e] % 2 == 0 && x2[e] <= 2 * cap,
              "rounded edge not integral or out of bounds");
        if (before[e] % 2 == 0)
            check(x2[e] == before[e], "integral edge changed");
        else
            check(x2[e] == before[e] - 1 || x2[e] == before[e] + 1,
                  "fractional edge not rounded to an adjacent integer");
        degree[edge.u] += x2[e];
        degree[edge.v] += x2[e];
        weight_after += x2[e] * edge.w;
    }
    for (std::size_t v = 0; v < b.size(); ++v)
        check(degree[v] <= 2LL * b[v], "rounded degree exceeds capacity");
    check(weight_before == weight_after, "rounding changed objective weight");

    // Each of these events does constant work (amortized for vector pushes).
    // The initial classification sees m edges, vertex sweeps/list see <=3n,
    // cursors see each of the 2h entries once, and each odd edge is retired once.
    check(work.edge_checks == edges.size(), "edge classification rescanned");
    check(work.vertex_checks <= 3 * b.size(), "vertex scans not linear");
    check(work.adjacency_probes <= 2 * fractional, "adjacency rescanned");
    check(work.traversed_edges == fractional, "edge revisited or left fractional");
    check(work.open_trails + work.closed_trails <= fractional,
          "empty trail or too many trails");
    ++checked_states;
    open_trails += work.open_trails;
    closed_trails += work.closed_trails;
    return work;
}

void fixed_cases() {
    for (bool simple : {true, false}) {
        verify(0, {}, {}, {}, simple);
        verify(2, {}, {0, 3, 1}, {}, simple);
        verify(1, {{0, 1, 7}}, {1, 1}, {2}, simple);

        // Opposite optimal choices in a two-edge star give x=1/2 on each edge.
        auto work = verify(1, {{0, 1, 7}, {0, 2, 7}}, {1, 1, 1}, {1, 1}, simple);
        check(work.open_trails == 1 && work.closed_trails == 0, "star branch");

        // Two perfect matchings of equal weight 12 on an even cycle.
        work = verify(2, {{0, 2, 3}, {0, 3, 7}, {1, 2, 5}, {1, 3, 9}},
                      {1, 1, 1, 1}, {1, 1, 1, 1}, simple);
        check(work.closed_trails == 1 && work.open_trails == 0, "cycle branch");

        // Parallel input entries are separate edges and form a length-2 cycle.
        work = verify(1, {{0, 1, 4}, {0, 1, 4}}, {1, 1}, {1, 1}, simple);
        check(work.closed_trails == 1 && work.open_trails == 0, "parallel cycle");

        // Even open path: two disjoint optimum matchings of weight 2.
        work = verify(2, {{0, 2, 1}, {0, 3, 1}, {1, 3, 1}, {1, 4, 1}},
                      {1, 1, 1, 1, 1}, {1, 1, 1, 1}, simple);
        check(work.open_trails == 1, "even open path");
        // Odd open path: the middle weight 3 ties outer weights 1+2.
        work = verify(2, {{0, 2, 1}, {1, 2, 3}, {1, 3, 2}},
                      {1, 1, 1, 1}, {1, 1, 1}, simple);
        check(work.open_trails == 1, "odd open path");
        verify(1, {{0, 1, 0}}, {1, 1}, {1}, simple);

        // The first open trail leaves a cycle attached to its internal vertex.
        work = verify(2, {{0, 2, 1}, {0, 3, 1}, {0, 4, 1}, {1, 4, 1},
                          {1, 5, 1}, {0, 5, 1}},
                      {2, 1, 1, 1, 1, 1}, Values(6, 1), simple);
        check(work.open_trails == 1 && work.closed_trails == 1,
              "residual cycle after open trail");
        // Likewise, a maximal closed trail can leave another cycle at an
        // internal vertex. The final vertex sweep must find it.
        work = verify(3, {{0, 3, 1}, {1, 3, 1}, {1, 4, 1}, {0, 4, 1},
                          {1, 5, 1}, {2, 5, 1}, {2, 6, 1}, {1, 6, 1}},
                      {1, 2, 1, 1, 1, 1, 1}, Values(8, 1), simple);
        check(work.closed_trails == 2, "residual cycle after closed trail");
    }

    // A multiplicity optimum with values above INT_MAX in the scaled state.
    // The assignments (q-2,q+1) and (q-1,q) tie and saturate the center.
    // This only calls rounding; it never expands billions of output entries.
    const int q = 1 << 30;
    verify(1, {{0, 1, 5}, {0, 2, 5}},
           {std::numeric_limits<int>::max(), q - 1, q + 1},
           {2LL * q - 3, 2LL * q + 1}, false);
    verify(1, {{0, 1, 5}}, {q, q}, {2LL * q}, false);

    // Mixed integral and fractional entries: only the latter are touched.
    verify(1, {{0, 1, 5}, {0, 2, 5}, {0, 3, 5}}, {4, 2, 2, 2},
           {2, 3, 3}, false);
}

// Independent integer oracle: enumerate edge multiplicities with capacity
// pruning and retain all optima. It uses no flow or half-integral rounding.
struct Optima {
    int64_t weight = -1;
    std::vector<Values> assignments;
};

void enumerate(const std::vector<Edge>& edges, bool simple, std::size_t e,
               std::vector<int>& remaining, Values& values,
               int64_t weight, Optima& optima) {
    if (e == edges.size()) {
        if (weight > optima.weight) {
            optima.weight = weight;
            optima.assignments.clear();
        }
        if (weight == optima.weight) optima.assignments.push_back(values);
        return;
    }
    const auto& edge = edges[e];
    int upper = std::min(remaining[edge.u], remaining[edge.v]);
    if (simple) upper = std::min(upper, 1);
    for (int value = 0; value <= upper; ++value) {
        remaining[edge.u] -= value;
        remaining[edge.v] -= value;
        values[e] = value;
        enumerate(edges, simple, e + 1, remaining, values,
                  weight + static_cast<int64_t>(value) * edge.w, optima);
        remaining[edge.u] += value;
        remaining[edge.v] += value;
    }
}

void verify_optimal_midpoints(int n_u, const std::vector<Edge>& edges,
                             const std::vector<int>& b, bool simple) {
    Optima optima;
    auto remaining = b;
    Values values(edges.size());
    enumerate(edges, simple, 0, remaining, values, 0, optima);
    // Each midpoint is a valid optimal symmetrization: use one optimum for
    // the forward orientation and the other for the reverse orientation.
    for (std::size_t i = 0; i < optima.assignments.size(); ++i) {
        for (std::size_t j = i; j < optima.assignments.size(); ++j) {
            Values x2(edges.size());
            int64_t weight = 0;
            for (std::size_t e = 0; e < edges.size(); ++e) {
                x2[e] = optima.assignments[i][e] + optima.assignments[j][e];
                weight += x2[e] * edges[e].w;
            }
            check(weight == 2 * optima.weight, "oracle midpoint weight");
            verify(n_u, edges, b, x2, simple);
        }
    }
}

void exhaustive_midpoints() {
    // All K(2,2) edge states absent/weight 0/weight 1/weight 3, capacities
    // in {0,1,2}, both modes, and every unordered pair of integral optima.
    for (int code = 0; code < 256; ++code) {
        std::vector<Edge> edges;
        int edge_code = code;
        for (int u = 0; u < 2; ++u) {
            for (int v = 2; v < 4; ++v) {
                int state = edge_code % 4;
                edge_code /= 4;
                if (state) edges.push_back({u, v, state == 3 ? 3 : state - 1});
            }
        }
        for (int cap_code = 0; cap_code < 81; ++cap_code) {
            std::vector<int> b(4);
            int digits = cap_code;
            for (int& cap : b) { cap = digits % 3; digits /= 3; }
            for (bool simple : {true, false})
                verify_optimal_midpoints(2, edges, b, simple);
        }
    }

    // Additional graphs include parallel entries and differently ordered
    // adjacency at revisited vertices, exercising nontrivial maximal trails.
    std::mt19937 rng(20260906);
    for (int trial = 0; trial < 500; ++trial) {
        const int n_u = 1 + rng() % 3, n_v = 1 + rng() % 3;
        std::vector<int> b(n_u + n_v);
        for (int& cap : b) cap = rng() % 3;
        std::vector<Edge> edges;
        const int m = rng() % 8;
        for (int e = 0; e < m; ++e)
            edges.push_back({static_cast<int>(rng() % n_u),
                             n_u + static_cast<int>(rng() % n_v),
                             static_cast<int>(rng() % 5)});
        for (bool simple : {true, false})
            verify_optimal_midpoints(n_u, edges, b, simple);
    }
}

void print_work(const char* family, int size, std::size_t n, std::size_t m,
                const anstee_detail::RoundingWork& work) {
    std::cout << family << '=' << size << " n=" << n << " m=" << m
              << " vertex_checks=" << work.vertex_checks
              << " edge_checks=" << work.edge_checks
              << " adjacency_probes=" << work.adjacency_probes
              << " traversed_edges=" << work.traversed_edges
              << " open=" << work.open_trails << " closed=" << work.closed_trails
              << '\n';
}

void scaling() {
    for (int k : {100, 200, 400, 800, 1600, 3200, 6400, 12800}) {
        // Reproduce the original sparse quadratic fixture on the real helper.
        std::vector<Edge> stars;
        for (int i = 0; i < k; ++i) {
            stars.push_back({i, k + 2 * i, 1});
            stars.push_back({i, k + 2 * i + 1, 1});
        }
        auto work = verify(k, stars, std::vector<int>(3 * k, 1), Values(2 * k, 1), true);
        check(work.open_trails == static_cast<std::size_t>(k) &&
              work.closed_trails == 0, "scaling stars trail count");
        print_work("stars", k, 3 * k, 2 * k, work);

        // A high-degree vertex revisited once per parallel-edge cycle exposes
        // any restart of its adjacency scan. Each leaf contributes weight w.
        std::vector<Edge> cycles;
        std::vector<int> b(k + 1, 1);
        b[0] = k;
        for (int v = 1; v <= k; ++v) {
            cycles.push_back({0, v, 1 + v % 7});
            cycles.push_back({0, v, 1 + v % 7});
        }
        for (bool simple : {true, false}) {
            work = verify(1, cycles, b, Values(2 * k, 1), simple);
            check(work.open_trails == 0 && work.closed_trails == 1,
                  "scaling parallel cycles trail count");
            if (simple) print_work("parallel_cycles", k, k + 1, 2 * k, work);
        }
    }
}
} // namespace

int main() {
    try {
        fixed_cases();
        exhaustive_midpoints();
        scaling();
        check(open_trails > 0 && closed_trails > 0, "missing trail branch coverage");
        std::cout << checked_states << " half-integral states passed; open trails="
                  << open_trails << " closed trails=" << closed_trails << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Rounding regression failed: " << error.what() << '\n';
        return 1;
    }
}
