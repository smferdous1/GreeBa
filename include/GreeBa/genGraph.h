#ifndef GENGRAPH_H
#define GENGRAPH_H

#include <vector>
#include <utility>
#include <tuple>

std::tuple<std::vector<std::vector<std::pair<int, int>>>, int, int> make_greeba_bpt(int n, int b, int M);

#endif // GENGRAPH_H