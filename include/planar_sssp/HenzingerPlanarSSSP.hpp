#pragma once

#include "PlanarSSSP.hpp"
#include "planar_sssp/HenzingerPriorityQueue.hpp"

using nodeSubsets_t = std::vector<std::vector<NetworKit::node>>;

namespace Koala {

class HenzingerPlanarSSSP : public PlanarSSSP {
   private:
    NetworKit::Graph normal_graph;
    std::vector<int> d;
    HPriorityQueue mainQ;
    std::vector<HPriorityQueue> Q;
    std::vector<std::vector<int>> nodeRegions;
    std::vector<int> isBoundry;
    int numOfRegions;
    int r;

    void initializeQueues(nodeSubsets_t& division);
    void mainThrust();

   public:
    HenzingerPlanarSSSP(NetworKit::Graph& Graph, NetworKit::node source, NetworKit::node target)
        : PlanarSSSP(Graph, source, target) {}

    void run();
};

} /* namespace Koala */
