#pragma once

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

using planar_embedding_t = std::unordered_map<NetworKit::node, std::vector<NetworKit::node>>;
using nodeSubsets_t = std::vector<std::vector<NetworKit::node>>;

namespace Koala {

    nodeSubsets_t findSuitableRDivision(NetworKit::Graph& Graph, int r, int c);

} /* namespace Koala */
