#pragma once

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

using planar_embedding_t = std::unordered_map<NetworKit::node, std::vector<NetworKit::node>>;

namespace Koala {

    planar_embedding_t findPlanarEmbeding(const NetworKit::Graph& G, bool verbose = false);
    NetworKit::Graph makeMaximalPlanar(NetworKit::Graph& G);

} /* namespace Koala */
