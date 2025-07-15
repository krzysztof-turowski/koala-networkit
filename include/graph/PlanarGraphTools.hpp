#pragma once

#include <networkit/graph/Graph.hpp>

namespace Koala {

namespace PlanarGraphTools {

using planar_embedding_t = std::unordered_map<NetworKit::node, std::vector<NetworKit::node>>;
using node_subsets_t = std::vector<std::vector<NetworKit::node>>;

planar_embedding_t findPlanarEmbeding(const NetworKit::Graph&, bool verbose = false);
NetworKit::Graph makeMaximalPlanar(NetworKit::Graph&);
NetworKit::Graph convertToMaxDegree3(NetworKit::Graph&, bool directed = false);
void assertDivision(const node_subsets_t&, NetworKit::Graph&);

}  // namespace PlanarGraphTools

}  // namespace Koala
