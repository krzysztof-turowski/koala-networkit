#pragma once

#include "PlanarSSSP.hpp"

namespace Koala {

class FredericksonPlanarSSSP : public PlanarSSSP {
 private:
    NetworKit::Graph normal_graph;
    int c = 3;
    int r1 = 100;
    int r2 = 25;

 public:
    FredericksonPlanarSSSP(NetworKit::Graph& graph, NetworKit::node source, NetworKit::node target)
        : PlanarSSSP(graph, source, target) {}

    void run();
};

} /* namespace Koala */
