# GreeBa

A C++17 library implementing greedy and exact maximum-weight **b-matching** algorithms on bipartite graphs.

## Algorithms

- **Greedy b-matching** — fast heuristic that sorts edges by weight and greedily selects them while respecting per-vertex degree capacities.
- **Anstee b-matching** — exact maximum-weight b-matching on bipartite graphs via Anstee's algorithm (Information Processing Letters 24, 1987), backed by OR-Tools min-cost flow. A zero-cost bypass allows unused capacity so that matching weight is optimized without forcing maximum cardinality.
- **Reduction to 1-matching** — exact maximum-weight b-matching on bipartite graphs using the vertex-copy and edge-gadget construction in Appendix B of Ferdous, *Algorithms for Degree-Constrained Subgraphs and Applications* (pp. 156–160).
- **MILP b-matching** — an independent integer programming formulation solved by SCIP through OR-Tools, used to verify the two exact combinatorial implementations.

### Bipartite b-matching MILP

For a bipartite graph `G = (U ∪ V, E)`, let `x_e` indicate whether edge `e` is selected. The model in [`milp_b_matching.cpp`](app/milp_b_matching.cpp) is:

```text
maximize    sum_{e in E} w_e x_e
subject to  sum_{e incident to v} x_e <= b(v)    for every v in U ∪ V
            x_e in {0, 1}                      for every e in E
```

Each edge can be selected at most once, and vertices may leave capacity unused. The objective is total weight. With `simple=false`, each variable instead has integer bounds `0 <= x_e <= min(b(u), b(v))`, allowing repeated selection of an edge.

```cpp
#include "milp_b_matching.h"

auto milp = milp_bipartite_b_matching(adj, n_u, b);
// Per-vertex capacities are also supported:
auto milp_vector = milp_bipartite_b_matching(adj, n_u, capacities);
```

The MILP reads only U-side adjacency entries, using one variable per entry; reverse entries are optional and parallel entries represent distinct edges. Weights and capacities must be nonnegative integers. It returns `MatchingResult`, requests zero relative MIP gap, requires SCIP's `OPTIMAL` status, and checks the solution's integrality and capacity constraints before returning. Unlike the flow solvers, SCIP uses floating-point numerical tolerances. The comparison tests recompute selected-edge weights as integers and require exact equality of all three totals; optimal edge sets may differ when there are ties.

### Reduction from b-matching to 1-matching

[`reduction_b_matching.cpp`](app/reduction_b_matching.cpp) implements the construction and recovery from Appendix B.1–B.2:

1. Create `b(v)` copies of each original vertex `v`. Every reduced vertex has capacity **1**.
2. For each original edge `e = (u, v)`, create two vertices `p_e,u` and `p_e,v`, join them with a middle edge, and connect each to every copy of its corresponding endpoint. All `b(u) + b(v) + 1` gadget edges have weight `w(e)`.
3. Compute a maximum-weight 1-matching on the reduced graph.
4. Recover `e` precisely when both of its gadget vertices are matched to endpoint copies. A single matched outer edge can be replaced by the middle edge at equal weight and contributes no original edge.

The reduced graph has `sum_v b(v) + 2|E|` vertices and `|E| + sum_v degree(v)b(v)` edges. These counts determine construction time and space. The copy capacity of 1 and the `2|E|` vertex term correct apparent typographical errors in B.1. With nonnegative weights, unmatched gadgets can be filled with middle edges and single outer edges normalized, yielding `weight(M') = sum_e w(e) + weight(M_b)`; thus an optimal reduced matching recovers an optimal b-matching.

```cpp
#include "reduction_b_matching.h"

// U vertices precede V vertices in adj. A scalar capacity also works.
std::vector<int> capacities(adj.size(), b);
auto result = reduction_bipartite_b_matching(adj, n_u, capacities);
```

The bipartite solver uses OR-Tools min-cost flow with unit capacities and a zero-cost bypass for unused capacity, so it optimizes weight without requiring maximum cardinality. It returns the existing `MatchingResult` format. Inputs require nonnegative integer weights and capacities; zero capacities are supported. Only U-side adjacency rows are read, with each entry representing a distinct edge selectable at most once.

For general loopless graphs, `reduce_b_matching(edges, capacities)` exposes the reduced adjacency list, copy IDs, gadget IDs, and baseline weight. Supply each original edge once. Solve that graph with a general maximum-weight 1-matching solver, then pass its symmetric mate vector (`-1` for unmatched vertices) to `recover_b_matching(reduction, mate)`. The bundled exact solver supports bipartite graphs; construction and recovery also support non-bipartite graphs. This implements the b-matching reduction only.

### GreeBa bipartite graph generation

[`make_greeba_bpt(n, b, M)`](src/genGraph.cpp) constructs a deterministic weighted bipartite graph for comparing the matching algorithms. `n` is the requested total vertex count, `b` controls the block size and is also used as the matching capacity in the example, and `M` sets the core edge weight. Use nonnegative `n` and positive `b`; the function does not validate its inputs.

1. Round `n` up to a multiple of `2b² + 2b`, giving `q = ceil(n / (2b² + 2b))` blocks.
2. In each block, create a complete bipartite core with `b` vertices on each side and `b²` edges of weight `M`.
3. Attach `b` distinct new outer vertices to each core vertex on the opposite side, using edges of weight `M - 1`. This adds `b²` outer vertices to each partition and `2b²` edges per block.
4. Connect each consecutive pair of blocks with `b²` edges of weight `1`, pairing their outer vertices by index. With zero-based block indices, these links join the U-side outer vertices of the even block to the V-side outer vertices of the odd block.
5. Convert the edge list into an undirected adjacency list, storing every edge in both directions.

The return value is `(adj, n_u, n_v)`, where each adjacency entry is a `(neighbor, weight)` pair and `n_u = n_v = q(b² + b)`. U vertices occupy indices `[0, n_u)` and V vertices occupy `[n_u, n_u + n_v)`. For `q > 0`, the graph has `q(2b² + 2b)` vertices and `(4q - 1)b²` undirected edges. Construction takes `O(qb²)` time and space.

For example, `make_greeba_bpt(24, 2, 100)` builds two blocks with 12 vertices per partition and 28 edges, with weights 100, 99, and 1.

## Requirements

- CMake ≥ 3.10
- C++17 compiler
- [OR-Tools](https://developers.google.com/optimization) with the SCIP MILP backend (installed via Homebrew: `brew install or-tools`)

## Install OR-Tools

### macOS with Homebrew

1. Install Apple's Command Line Tools if they are not already installed:

   ```bash
   xcode-select --install
   ```

2. Install [Homebrew](https://brew.sh/) if needed, then install CMake and the OR-Tools C++ library:

   ```bash
   brew install cmake or-tools
   ```

   The [Homebrew OR-Tools package](https://formulae.brew.sh/formula/or-tools) includes SCIP as a dependency, which this project's MILP solver requires.

3. From the GreeBa project directory, configure and build using your Homebrew installation paths:

   ```bash
   cmake -S . -B build \
     -DCMAKE_PREFIX_PATH="$(brew --prefix)" \
     -DCMAKE_MODULE_PATH="$(brew --prefix or-tools)/lib/cmake/ortools/modules"
   cmake --build build
   ctest --test-dir build --output-on-failure
   ./build/app/test_app
   ```

   These paths let CMake locate OR-Tools and its dependency modules without relying on the Apple Silicon and version-specific paths currently in `CMakeLists.txt`. The tests exercise both min-cost flow and the SCIP MILP backend. The default example should report weight **1588** for Anstee, reduction, and MILP, followed by an equality check of **YES**.

### Other platforms or custom installations

Follow the official [OR-Tools C++ installation guide](https://developers.google.com/optimization/install/cpp) for binary distributions or source builds. Install the C++ headers, libraries, and CMake package, with SCIP support enabled. Set `CMAKE_PREFIX_PATH` to your installation prefix and, if needed, `CMAKE_MODULE_PATH` to its `lib/cmake/ortools/modules` directory when configuring GreeBa. If CMake cannot find `ortoolsConfig.cmake`, set `ortools_DIR` to the directory containing that file.

## Build

```bash
cmake -S . -B build
cmake --build build
```

## Run the example

```bash
./build/app/test_app
```

The example generates a deterministic GreeBa bipartite graph, runs greedy, Anstee, reduction, and MILP matching, and prints the matched edges and total weight for each. It checks that Anstee, the reduction, and the MILP return equal weights and exits with an error if they differ.

Run the matching tests, including 195 generated-graph comparisons between Anstee, the reduction, and the MILP (`b` from 1 to 4 and `M` in `{1, 2, 3, 100}`), exhaustive comparisons for all three solvers on 300 small random graphs, and 100 random exhaustive checks of Anstee and MILP edge-multiplicity modes:

```bash
ctest --test-dir build --output-on-failure
```

## Project structure

```
.
├── app/                  # Executable and algorithm implementations
│   ├── greedy_b_matching.cpp/h
│   ├── anstee_b_matching.cpp/h
│   ├── reduction_b_matching.cpp/h
│   ├── milp_b_matching.cpp/h
│   ├── reduction_tests.cpp
│   └── test.cpp
├── include/
│   ├── GreeBa/           # Graph generation utilities
│   └── LiteGraph/        # Lightweight graph data structures
└── src/                  # Library source
    └── genGraph.cpp      # GreeBa bipartite graph construction
```
