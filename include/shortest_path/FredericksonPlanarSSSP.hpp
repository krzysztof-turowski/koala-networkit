#pragma once

#include "PlanarSSSP.hpp"

namespace Koala {

class FredericksonPlanarSSSP : public PlanarSSSP {
 private:
    NetworKit::Graph normal_graph;
    int c = 6;
    int r1 = NetworKit::none;
    int r2 = NetworKit::none;

 public:
    FredericksonPlanarSSSP(NetworKit::Graph& graph, NetworKit::node source, NetworKit::node target)
        : PlanarSSSP(graph, source, target) {}

    void setDivisionParameters(int level_1, int level_2);

    void run();
};

} /* namespace Koala */
