#pragma once

#include "PlanarSSSP.hpp"

namespace Koala {

    class HenzingerPlanarSSSP : public PlanarSSSP {
    private:
        NetworKit::Graph normal_graph;
        std::vector<int> isBoundry;

    public:
        HenzingerPlanarSSSP(NetworKit::Graph& Graph, NetworKit::node source,
            NetworKit::node target)
            : PlanarSSSP(Graph, source, target) {
        }

        void run();
    };

} /* namespace Koala */
