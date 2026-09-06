# GreeBa

A C++17 library implementing greedy and exact maximum-weight **b-matching** algorithms on bipartite graphs.

## Algorithms

- **Greedy b-matching** — fast heuristic that sorts edges by weight and greedily selects them while respecting per-vertex degree capacities.
- **Anstee b-matching** — exact maximum-weight b-matching on bipartite graphs via Anstee's algorithm (Information Processing Letters 24, 1987), backed by OR-Tools min-cost flow. A zero-cost bypass allows unused capacity so that matching weight is optimized without forcing maximum cardinality.
- **Reduction to 1-matching** — exact maximum-weight b-matching on bipartite graphs using the vertex-copy and edge-gadget construction in Appendix B of Ferdous, *Algorithms for Degree-Constrained Subgraphs and Applications* (pp. 156–160).
- **MILP b-matching** — an independent integer programming formulation solved by SCIP through OR-Tools, used to verify the two exact combinatorial implementations.

### Shared capacity preprocessing

Before running **simple** b-matching, use [`clamp_bipartite_capacities`](app/b_matching_preprocessing.h)
to replace each capacity with `min(b(v), degree(v))`. Generate the graph first,
then compute this vector once and pass it to **greedy, Anstee, reduction, and MILP**:

```cpp
#include "GreeBa/genGraph.h"
#include "b_matching_preprocessing.h"
#include "greedy_b_matching.h"
#include "anstee_b_matching.h"
#include "reduction_b_matching.h"
#include "milp_b_matching.h"

auto [adj, n_u, n_v] = make_greeba_bpt(n, b, M);
const auto capacity = clamp_bipartite_capacities(adj, n_u, b);
// A vector of per-vertex capacities can replace the scalar b above.
auto greedy = greedy_weighted_b_matching(adj, capacity);
auto anstee = anstee_bipartite_b_matching(adj, n_u, capacity, true);
auto reduction = reduction_bipartite_b_matching(adj, n_u, capacity);
auto milp = milp_bipartite_b_matching(adj, n_u, capacity, true);
```

This is a separate preprocessing step; the solver APIs use the capacities they
receive. It preserves every feasible simple matching and its weight because
each incident input edge can be selected at most once. Isolates get capacity
zero. The helper reads U-side entries once, counting both endpoints, so reverse
rows are optional and do not double-count degrees. Parallel input entries count
as distinct edges. Degree counters saturate at the requested capacity to avoid
integer overflow. Time is `O(n+m)` and extra space is `O(n)`.

The generator's `b`, graph topology, weights, and original capacity vector are
not changed. The reduction size formulas below use the **effective capacities
passed to the solver**. Do not apply this degree bound to `simple=false`, which
allows repeated selection of an edge. CTest's `capacity_preprocessing` checks
feasible-set preservation and all four algorithms, including isolates,
parallel edges, optional reverse rows, and `INT_MAX` capacities.

### Anstee rounding and numeric limits

Anstee reads each U-side adjacency entry as a distinct edge, including parallel entries; reverse adjacency rows are optional. A zero-cost source-to-sink bypass permits unused capacity.

The two directed flow values and their symmetrized sum `x2 = 2*x` use `int64_t`. Stage 2 rounds fractional edges using residual degree counters, a list of initially odd vertices, per-vertex adjacency cursors, and per-edge active flags. Each adjacency entry is inspected at most once and each fractional edge is traversed once, so rounding takes `O(n+m)` time and auxiliary space. This bound covers Stage 2 only: OR-Tools uses its own cost-scaling flow algorithm, and this implementation does not claim the paper's full strongly-polynomial bound.

Weights and per-vertex capacities must be nonnegative `int` values. Invalid partitions are rejected before arithmetic or narrowing. Oversized graphs are rejected before constructing node or arc indices; the conservative bounds reserve room for the auxiliary nodes/arcs used by OR-Tools. Output multiplicities are checked before conversion back to `int`. The returned total weight must fit in `long long`; aggregate result-weight overflow is not currently checked. OR-Tools also has internal numeric limits, and failure to reach `OPTIMAL` is reported as an exception.

`MatchingResult::edges` contains one entry per selected unit of multiplicity. If `K` units are selected, result construction requires `O(K)` time and storage; a large capacity can therefore be impractical even though symmetrization itself is safe. The numeric regression solves the real large-capacity flow and calls the production symmetrization helper without expanding the result.

CTest includes `anstee_numeric` and `anstee_rounding` alongside the existing matching suite. The rounding tests exercise open and closed trails, parallel entries, leftover cycles, 64-bit scaled values, averages of independently enumerated optimal matchings, and operation-count bounds on growing fixtures.

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
#include "b_matching_preprocessing.h"

auto effective = clamp_bipartite_capacities(adj, n_u, b);
auto milp = milp_bipartite_b_matching(adj, n_u, effective);
// Per-vertex capacities are also supported:
auto effective_vector = clamp_bipartite_capacities(adj, n_u, capacities);
auto milp_vector = milp_bipartite_b_matching(adj, n_u, effective_vector);
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
#include "b_matching_preprocessing.h"

// U vertices precede V vertices in adj. Preprocess after graph generation.
const auto capacities = clamp_bipartite_capacities(adj, n_u, b);
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
- [OR-Tools](https://developers.google.com/optimization) C++ headers, libraries, and CMake package, built with the SCIP MILP backend
- `pkg-config` (provided by `pkgconf` on Homebrew) when required by the installed OR-Tools dependency packages

`pip install ortools` installs the Python package; it does not provide the C++ headers and CMake package needed to build this project. GreeBa's CMake configuration finds an existing C++ installation; it does not download or install OR-Tools.

## Install OR-Tools

### macOS with Homebrew

1. Install Apple's Command Line Tools if they are not already installed:

   ```bash
   xcode-select --install
   ```

2. Install [Homebrew](https://brew.sh/) if needed, then install CMake and the OR-Tools C++ library:

   ```bash
   brew install cmake pkgconf or-tools
   ```

   The [Homebrew OR-Tools package](https://formulae.brew.sh/formula/or-tools) includes SCIP as a dependency, which this project's MILP solver requires.

3. From the GreeBa project directory, configure and build using your Homebrew installation paths:

   ```bash
   cmake -S . -B build -DCMAKE_PREFIX_PATH="$(brew --prefix)"
   cmake --build build
   ctest --test-dir build --output-on-failure
   ./build/app/test_app
   ```

   The prefix is supplied at configure time; no Homebrew or version-specific paths are embedded in `CMakeLists.txt`. OR-Tools supplies its own dependency module paths. The tests exercise both min-cost flow and the SCIP MILP backend. The default example should report weight **1588** for Anstee, reduction, and MILP, followed by an equality check of **YES**.

### Linux clusters

Use an OR-Tools C++ installation built for the cluster, with SCIP support and its dependencies available. If the cluster provides environment modules, load its compiler, CMake, and OR-Tools modules first; module names depend on the site. Keep that environment loaded when running the executable, including inside a batch job.

#### If the cluster does not provide OR-Tools

Install OR-Tools under your own account once, then point GreeBa at that installation. No administrator access or Homebrew is needed. The recipe below pins OR-Tools to `v9.15` and follows its [Linux source-build instructions](https://developers.google.com/optimization/install/cpp/source_linux) and [versioned CMake requirements](https://github.com/google/or-tools/blob/v9.15/cmake/README.md#requirement).

First load the cluster's C/C++ compiler and CMake modules. This OR-Tools release requires **CMake 3.24 or newer**; its build documentation calls for **GCC 10 or newer**, or an equivalent supported compiler. You also need Git and a build tool such as Make. The clone and configuration steps need internet access to fetch OR-Tools and its dependencies. Use source, build, and installation paths without spaces, with enough storage for a dependency build.

Run the following on the cluster in Bash, using an allocation suitable for compilation. These commands preserve your current working directory:

```bash
greeba_ortools_src="$HOME/src/or-tools-9.15"
greeba_ortools_prefix="$HOME/.local/or-tools-9.15"
mkdir -p "$HOME/src"

# Clone once; on a later build, reuse this source directory.
git clone --depth 1 --branch v9.15 \
  https://github.com/google/or-tools.git "$greeba_ortools_src"

cmake -S "$greeba_ortools_src" -B "$greeba_ortools_src/build" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$greeba_ortools_prefix" \
  -DBUILD_DEPS=ON \
  -DINSTALL_BUILD_DEPS=ON \
  -DUSE_SCIP=ON \
  -DBUILD_SAMPLES=OFF \
  -DBUILD_EXAMPLES=OFF \
  -DBUILD_TESTING=OFF
cmake --build "$greeba_ortools_src/build" --parallel 4
cmake --install "$greeba_ortools_src/build"
```

`BUILD_DEPS=ON` downloads and builds the required dependency libraries; `USE_SCIP=ON` includes the MILP solver used by GreeBa. `INSTALL_BUILD_DEPS=ON` installs those libraries alongside OR-Tools. Samples, examples, and upstream tests are disabled to reduce setup work; GreeBa's own tests remain enabled below. The installation prefix can instead be a persistent project directory accessible from the compute nodes. Adjust `--parallel 4` to the CPU and memory allocation, and use the same compiler environment for both builds.

After installation, run this from the GreeBa source directory:

```bash
cmake -S . -B build/cluster -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$HOME/.local/or-tools-9.15"
cmake --build build/cluster --parallel 4
ctest --test-dir build/cluster --output-on-failure
./build/cluster/app/test_app
```

Use your chosen prefix if it differs from the example. The matching tests exercise both OR-Tools min-cost flow and SCIP. This source-install recipe has been checked against OR-Tools v9.15's build configuration, but has not been executed on your cluster.

If the cluster has no internet access, cloning only the OR-Tools repository elsewhere is insufficient: configuration downloads dependency sources too. A matching [Linux C++ binary distribution](https://developers.google.com/optimization/install/cpp/binary_linux) can be downloaded elsewhere and transferred if its architecture and system-library requirements match the cluster. Otherwise, stage all dependency sources or use a compatible Linux build environment with network access before transferring the installation. Do not transfer the macOS libraries.

#### If OR-Tools is already installed

If the modules make OR-Tools discoverable by CMake, configure without additional paths:

```bash
cmake -S . -B build/cluster -DCMAKE_BUILD_TYPE=Release
```

For an installation in a custom location, supply its prefix instead (replace the example path):

```bash
cmake -S . -B build/cluster -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="/path/to/or-tools"
```

The prefix is the installation directory containing `include` and `lib` or `lib64`. If dependencies are installed separately, pass a quoted semicolon-separated list of their prefixes, for example `-DCMAKE_PREFIX_PATH="/path/to/or-tools;/path/to/dependencies"`. Alternatively, set `-Dortools_DIR="/path/to/directory/containing/ortoolsConfig.cmake"`; dependency prefixes may still be needed. CMake prints the selected OR-Tools package directory during configuration.

Then build, test, and run in a suitable cluster allocation:

```bash
cmake --build build/cluster --parallel 4
ctest --test-dir build/cluster --output-on-failure
./build/cluster/app/test_app
```

Choose a build parallelism that fits the allocated CPUs. Transfer the source and configure a new build on the cluster: the macOS executables and existing `CMakeCache.txt` are specific to the local machine. Changing compiler or dependency installations also calls for a fresh build directory.

Some SCIP packages export a CMake target named `libscip`, while OR-Tools expects `SCIP::libscip`. GreeBa creates a forwarding target only when the former exists and the latter is missing. This avoids OR-Tools' compatibility warning and preserves the installed library's include paths and link requirements; it does not modify or reinstall SCIP or OR-Tools.

### Other platforms or custom installations

Follow the official [OR-Tools C++ installation guide](https://developers.google.com/optimization/install/cpp) for binary distributions or source builds. Install the C++ headers, libraries, and CMake package, with SCIP support enabled. Set `CMAKE_PREFIX_PATH` to your installation and dependency prefixes. If CMake cannot find `ortoolsConfig.cmake`, set `ortools_DIR` to the directory containing that file. The dependency versions and compiler requirements are determined by the OR-Tools installation you use.

## Build

```bash
cmake -S . -B build
cmake --build build
```

## Run the example

```bash
./build/app/test_app
```

The example generates a deterministic GreeBa bipartite graph, clamps capacities once, passes the same vector to greedy, Anstee, reduction, and MILP matching, and prints the matched edges and total weight for each. It checks that Anstee, the reduction, and the MILP return equal weights and exits with an error if they differ.

Run the matching tests, including 195 generated-graph comparisons between Anstee, the reduction, and the MILP (`b` from 1 to 4 and `M` in `{1, 2, 3, 100}`), exhaustive comparisons for all three solvers on 300 small random graphs, and 100 random exhaustive checks of Anstee and MILP edge-multiplicity modes:

```bash
ctest --test-dir build --output-on-failure
```

## Project structure

```
.
├── app/                  # Executable and algorithm implementations
│   ├── greedy_b_matching.cpp/h
│   ├── b_matching_preprocessing.h
│   ├── preprocessing_tests.cpp
│   ├── anstee_b_matching.cpp/h
│   ├── anstee_numeric_detail.h
│   ├── anstee_numeric_tests.cpp
│   ├── anstee_rounding_detail.h
│   ├── anstee_rounding_tests.cpp
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
