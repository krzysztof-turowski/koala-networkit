#pragma once

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

using planar_embedding_t = std::unordered_map<NetworKit::node, std::vector<NetworKit::node>>;
using node_subsets_t = std::vector<std::vector<NetworKit::node>>;

namespace Koala {

node_subsets_t findSuitableRDivision(NetworKit::Graph& graph, int r, int c);

} /* namespace Koala */
