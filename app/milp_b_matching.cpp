#include "milp_b_matching.h"

#include "ortools/linear_solver/linear_solver.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

MatchingResult milp_bipartite_b_matching(
    const std::vector<std::vector<std::pair<int, int>>>& adj,
    int n_u, const std::vector<int>& capacity, bool simple)
{
    if (adj.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
        throw std::length_error("graph exceeds integer vertex indexing");
    const int n = static_cast<int>(adj.size());
    if (n_u < 0 || n_u > n) throw std::invalid_argument("n_u out of range");
    if (capacity.size() != adj.size())
        throw std::invalid_argument("capacity size must equal vertex count");
    if (std::any_of(capacity.begin(), capacity.end(), [](int b) { return b < 0; }))
        throw std::invalid_argument("capacities must be nonnegative");
    std::vector<MatchingEdge> edges;
    for (int u = 0; u < n_u; ++u) {
        for (const auto& [v, w] : adj[u]) {
            if (v < n_u || v >= n) throw std::out_of_range("edge endpoint outside V");
            if (w < 0) throw std::invalid_argument("weights must be nonnegative");
            edges.push_back({u, v, w});
        }
    }

    using operations_research::MPSolver;
    std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("SCIP"));
    if (!solver) throw std::runtime_error("OR-Tools was built without the SCIP MILP backend");

    std::vector<operations_research::MPConstraint*> degree_constraints;
    for (int v = 0; v < n; ++v)
        degree_constraints.push_back(solver->MakeRowConstraint(
            0, capacity[v], "degree_" + std::to_string(v)));

    auto* objective = solver->MutableObjective();
    objective->SetMaximization();
    std::vector<operations_research::MPVariable*> variables;
    for (const auto& e : edges) {
        const std::string name = "x_" + std::to_string(variables.size());
        auto* x = simple ? solver->MakeBoolVar(name) : solver->MakeIntVar(
            0, std::min(capacity[e.u], capacity[e.v]), name);
        variables.push_back(x);
        degree_constraints[e.u]->SetCoefficient(x, 1);
        degree_constraints[e.v]->SetCoefficient(x, 1);
        objective->SetCoefficient(x, e.weight);
    }
    // Do not accept a merely feasible incumbent or stop at a nonzero MIP gap.
    operations_research::MPSolverParameters parameters;
    parameters.SetDoubleParam(operations_research::MPSolverParameters::RELATIVE_MIP_GAP, 0);
    const auto status = solver->Solve(parameters);
    if (status != MPSolver::OPTIMAL)
        throw std::runtime_error("SCIP MILP did not reach OPTIMAL (status " +
                                 std::to_string(static_cast<int>(status)) + ")");
    if (!solver->VerifySolution(1e-6, true))
        throw std::runtime_error("SCIP MILP solution verification failed");

    MatchingResult result;
    result.degree.assign(n, 0);
    for (size_t i = 0; i < edges.size(); ++i) {
        const auto& e = edges[i];
        const double value = variables[i]->solution_value();
        const long long count = std::llround(value);
        const int upper = simple ? 1 : std::min(capacity[e.u], capacity[e.v]);
        if (std::abs(value - count) > 1e-6 || count < 0 || count > upper ||
            count > capacity[e.u] - result.degree[e.u] ||
            count > capacity[e.v] - result.degree[e.v])
            throw std::runtime_error("MILP solution is not an integral feasible b-matching");
        result.degree[e.u] += static_cast<int>(count);
        result.degree[e.v] += static_cast<int>(count);
        result.totalWeight += count * e.weight;
        for (long long k = 0; k < count; ++k) result.edges.push_back(e);
    }
    return result;
}

MatchingResult milp_bipartite_b_matching(
    const std::vector<std::vector<std::pair<int, int>>>& adj,
    int n_u, int capacity, bool simple)
{
    if (capacity < 0) throw std::invalid_argument("capacity must be nonnegative");
    return milp_bipartite_b_matching(adj, n_u,
                                    std::vector<int>(adj.size(), capacity), simple);
}
