#include "anstee_b_matching.h"
#include "anstee_numeric_detail.h"

#include <climits>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
using Adj = std::vector<std::vector<std::pair<int, int>>>;
namespace numeric = anstee_numeric_detail;

void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

template<class Exception, class F> void rejects(F f) {
    try { f(); } catch (const Exception&) { return; }
    throw std::runtime_error("expected exception was not thrown");
}

// Solve the real two-vertex Hitchcock network, then use the same numeric
// operation as production Stage 2a. Never build an expanded MatchingResult:
// capacities near INT_MAX must still require only six nodes and seven arcs.
void check_symmetrization(int capacity, bool simple) {
    operations_research::SimpleMinCostFlow flow;
    const int cap = simple ? 1 : capacity;
    const int fwd = cap ? flow.AddArcWithCapacityAndUnitCost(2, 5, cap, -1) : -1;
    const int rev = cap ? flow.AddArcWithCapacityAndUnitCost(3, 4, cap, -1) : -1;
    for (int v = 0; v < 2; ++v) {
        if (capacity == 0) continue;
        flow.AddArcWithCapacityAndUnitCost(0, 2 + v, capacity, 0);
        flow.AddArcWithCapacityAndUnitCost(4 + v, 1, capacity, 0);
    }
    const int64_t supply = int64_t{2} * capacity;
    flow.AddArcWithCapacityAndUnitCost(0, 1, supply, 0);
    flow.SetNodeSupply(0, supply);
    flow.SetNodeSupply(1, -supply);
    check(flow.Solve() == operations_research::SimpleMinCostFlow::OPTIMAL,
          "bounded Hitchcock flow did not solve");
    const int expected = capacity == 0 ? 0 : cap;
    if (fwd >= 0) check(flow.Flow(fwd) == expected, "incorrect forward flow");
    if (rev >= 0) check(flow.Flow(rev) == expected, "incorrect reverse flow");
    const std::vector<int64_t> x2{numeric::symmetrized_flow(flow, fwd, rev)};
    check(x2[0] == int64_t{2} * expected, "symmetrized flow overflowed");
    check(numeric::checked_multiplicity(x2[0]) == expected,
          "checked multiplicity conversion failed");
    check(numeric::symmetrized_flow(flow, -1, rev) == expected,
          "missing forward arc mishandled");
    check(numeric::symmetrized_flow(flow, fwd, -1) == expected,
          "missing reverse arc mishandled");
    check(numeric::symmetrized_flow(flow, -1, -1) == 0,
          "missing arcs mishandled");
    if (!simple && capacity == (1 << 30))
        std::cout << "1073741824 + 1073741824 = " << x2[0]
                  << " (bounded-memory production symmetrization)\n";
}

void check_empty(const MatchingResult& result, std::size_t n) {
    check(result.edges.empty() && result.totalWeight == 0 &&
              result.degree == std::vector<int>(n, 0),
          "incorrect empty result");
}

void check_public_api(bool simple) {
    for (int n : {0, 1, 4}) {
        const Adj adj(n);
        for (int n_u = 0; n_u <= n; ++n_u) {
            for (int cap : {0, 3, INT_MAX}) {
                check_empty(anstee_bipartite_b_matching(adj, n_u, cap, simple), n);
                check_empty(anstee_bipartite_b_matching(
                    adj, n_u, std::vector<int>(n, cap), simple), n);
            }
        }
        for (int invalid : {INT_MIN, -1, n + 1, INT_MAX}) {
            rejects<std::invalid_argument>([&] {
                anstee_bipartite_b_matching(adj, invalid, 0, simple);
            });
            rejects<std::invalid_argument>([&] {
                anstee_bipartite_b_matching(adj, invalid, std::vector<int>(n), simple);
            });
        }
    }

    const Adj one_edge{{{1, INT_MAX}}, {{0, INT_MAX}}};
    for (int cap : {0, 1, 7}) {
        const auto result = anstee_bipartite_b_matching(one_edge, 1, cap, simple);
        const int expected = simple && cap > 0 ? 1 : cap;
        check(result.edges.size() == static_cast<std::size_t>(expected) &&
                  result.degree == std::vector<int>(2, expected) &&
                  result.totalWeight == static_cast<long long>(expected) * INT_MAX,
              "incorrect result for zero/unit/ordinary capacity or large weight");
    }
    for (const auto& capacities : {std::vector<int>{INT_MAX, 0},
                                  std::vector<int>{0, INT_MAX},
                                  std::vector<int>{INT_MAX, 3},
                                  std::vector<int>{3, INT_MAX}}) {
        const auto result = anstee_bipartite_b_matching(one_edge, 1, capacities, simple);
        const int bound = std::min(capacities[0], capacities[1]);
        const int expected = simple && bound > 0 ? 1 : bound;
        check(result.edges.size() == static_cast<std::size_t>(expected) &&
                  result.degree == std::vector<int>(2, expected) &&
                  result.totalWeight == static_cast<long long>(expected) * INT_MAX,
              "large-capacity bottleneck failed");
    }
    if (simple) {
        const auto result = anstee_bipartite_b_matching(one_edge, 1, INT_MAX, true);
        check(result.edges.size() == 1 && result.totalWeight == INT_MAX &&
                  result.degree == std::vector<int>(2, 1),
              "large uniform simple capacity failed");
    }
    const auto parallel = anstee_bipartite_b_matching(
        Adj{{{1, 6}, {1, 4}}, {}}, 1, 5, simple);
    check(parallel.edges.size() == (simple ? 2u : 5u) &&
              parallel.totalWeight == (simple ? 10 : 30) &&
              parallel.degree == std::vector<int>(2, simple ? 2 : 5),
          "parallel-edge multiplicity semantics changed");

    rejects<std::invalid_argument>([&] {
        anstee_bipartite_b_matching(one_edge, 1, std::vector<int>{1}, simple);
    });
    rejects<std::invalid_argument>([&] {
        anstee_bipartite_b_matching(one_edge, 1, std::vector<int>{-1, 1}, simple);
    });
    rejects<std::invalid_argument>([&] {
        anstee_bipartite_b_matching(one_edge, 1, -1, simple);
    });
    rejects<std::invalid_argument>([&] {
        anstee_bipartite_b_matching(Adj{{{1, -1}}, {}}, 1, 1, simple);
    });
    for (int invalid : {INT_MIN, -1, 0, 2, INT_MAX}) {
        rejects<std::out_of_range>([&] {
            anstee_bipartite_b_matching(Adj{{{invalid, 1}}, {}}, 1, 1, simple);
        });
    }
}

void check_numeric_limits() {
    // Exercise the production size guard without allocating huge graphs.
    const std::size_t largest_n = (numeric::index_limit - 4) / 2;
    check(numeric::checked_partition_size(largest_n, 0) == static_cast<int>(largest_n),
          "supported node-count boundary rejected");
    check(numeric::checked_partition_size(largest_n, static_cast<int>(largest_n)) ==
              static_cast<int>(largest_n), "supported partition boundary rejected");
    for (std::size_t n : {largest_n + 1, numeric::index_limit,
                          std::numeric_limits<std::size_t>::max()}) {
        rejects<std::length_error>([&] { numeric::checked_partition_size(n, 0); });
    }
    const int n = static_cast<int>(largest_n);
    check(int64_t{2} + n + (n - 1) + 1 + 2 <=
              static_cast<int64_t>(numeric::index_limit),
          "node indices or solver node count exceed supported range");
    for (int vertices : {0, 2, n}) {
        const auto m = numeric::max_input_edges(vertices);
        const uint64_t arcs = uint64_t{2} * m + uint64_t{2} * vertices + 3;
        check(arcs <= numeric::index_limit && arcs + 2 > numeric::index_limit,
              "arc-count bound does not cover network and solver arcs");
    }
    check(numeric::checked_multiplicity(0) == 0, "zero narrowing failed");
    check(numeric::checked_multiplicity(int64_t{2} * INT_MAX) == INT_MAX,
          "maximum int multiplicity rejected");
    for (int64_t x2 : {int64_t{-2}, int64_t{2} * INT_MAX + 2,
                       std::numeric_limits<int64_t>::max()}) {
        rejects<std::overflow_error>([&] { numeric::checked_multiplicity(x2); });
    }
    rejects<std::logic_error>([] { numeric::checked_multiplicity(1); });
}
} // namespace

int main() {
    try {
        check_numeric_limits();
        for (bool simple : {true, false}) {
            for (int capacity : {0, 1, 7, (1 << 30) - 1, 1 << 30, INT_MAX})
                check_symmetrization(capacity, simple);
            check_public_api(simple);
        }
        std::cout << "Anstee numeric and partition regressions passed in both modes.\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
