#pragma once

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

using planar_embedding_t = std::unordered_map<NetworKit::node, std::vector<NetworKit::node>>;
using nodeSubsets_t = std::vector<std::vector<NetworKit::node>>;

namespace Koala {

planar_embedding_t findPlanarEmbeding(const NetworKit::Graph& G, bool verbose = false);
NetworKit::Graph makeMaximalPlanar(NetworKit::Graph& G);
NetworKit::Graph convertToMaxDegree3(NetworKit::Graph& G, bool directed = false);
NetworKit::Graph convertDirectedGraphToUndirected(NetworKit::Graph& Graph);
void assert_division(const nodeSubsets_t& division, NetworKit::Graph& Graph);
void printGraph(const NetworKit::Graph& graph);
void print_division(nodeSubsets_t division);

} /* namespace Koala */
