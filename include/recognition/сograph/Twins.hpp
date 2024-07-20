#pragma once

#include <vector>
#include <networkit/base/Algorithm.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

class Twins {
 public:
    std::vector<NetworKit::count> used;
    NetworKit::Graph &graph;

    explicit Twins(NetworKit::Graph &Graph);

    bool false_twins(NetworKit::count A, NetworKit::count B, NetworKit::count twins_counter);

    bool true_twins(NetworKit::count A, NetworKit::count B, NetworKit::count twins_counter);

    int twins(NetworKit::count A, NetworKit::count B, NetworKit::count twins_counter);
};

} /* namespace Koala */

