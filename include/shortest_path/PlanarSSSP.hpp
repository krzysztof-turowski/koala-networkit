#pragma once

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <set>

namespace Koala {
class PlanarSSSP : public NetworKit::Algorithm {
 public:
    virtual void run() = 0;

    const std::vector<NetworKit::edgeweight> getDistances();
    NetworKit::edgeweight distance(NetworKit::node t);

    explicit PlanarSSSP(
        NetworKit::Graph& graph, NetworKit::node source, NetworKit::node target = NetworKit::none);

 protected:
    NetworKit::node source;
    NetworKit::node target;
    NetworKit::Graph& graph;
    std::vector<NetworKit::edgeweight> distances;
};

} /* namespace Koala */
