# GreeBa

A C++17 library implementing greedy and exact maximum-weight **b-matching** algorithms on bipartite graphs.

## Algorithms

- **Greedy b-matching** — fast heuristic that sorts edges by weight and greedily selects them while respecting per-vertex degree capacities.
- **Anstee b-matching** — exact maximum-weight b-matching on bipartite graphs via Anstee's algorithm (Information Processing Letters 24, 1987), backed by OR-Tools for LP/min-cost flow solving.

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
- [OR-Tools](https://developers.google.com/optimization) (installed via Homebrew: `brew install or-tools`)

## Build

```bash
cmake -S . -B build
cmake --build build
```

## Run the example

```bash
./build/app/test_app
```

The example generates a deterministic GreeBa bipartite graph, runs both matching algorithms, and prints the matched edges and total weight for each.

## Project structure

```
.
├── app/                  # Executable and algorithm implementations
│   ├── greedy_b_matching.cpp/h
│   ├── anstee_b_matching.cpp/h
│   └── test.cpp
├── include/
│   ├── GreeBa/           # Graph generation utilities
│   └── LiteGraph/        # Lightweight graph data structures
└── src/                  # Library source
    └── genGraph.cpp      # GreeBa bipartite graph construction
```
