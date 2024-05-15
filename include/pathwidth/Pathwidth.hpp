#pragma once

#include "recognition/CographAlg.hpp"

namespace Koala {
    class Pathwidth {
    public:

        NetworKit::Graph *graph;
        Koala::CographRecognition *recognition;

        Pathwidth(NetworKit::Graph &Graph) {
            graph = &Graph;
        }

        NetworKit::count width = 0;

        void run();

        NetworKit::count pathwidth(NetworKit::count n, NetworKit::count v);

        NetworKit::count SubtreeSize(NetworKit::count n, NetworKit::count v);

    };
}
