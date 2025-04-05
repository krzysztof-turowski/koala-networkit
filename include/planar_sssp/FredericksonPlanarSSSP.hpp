#pragma once

#include "PlanarSSSP.hpp"
#include "recognition/CographRecognitionOther.hpp"

namespace Koala
{

    class FredericksonPlanarSSSP : public PlanarSSSP
    {
    private:
        NetworKit::Graph normal_graph;

    public:
        FredericksonPlanarSSSP(NetworKit::Graph &Graph, NetworKit::node source,
                               NetworKit::node target)
            : PlanarSSSP(Graph, source) {}

        void run();
    };

} /* namespace Koala */
