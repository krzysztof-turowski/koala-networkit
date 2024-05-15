#pragma once

#include "recognition/CographAlg.hpp"

namespace Koala {
    class MaxClique {
    public:
        Koala::Cotree *cotree;

        MaxClique(Koala::Cotree &CoTree, NetworKit::count N) {
            cotree = &CoTree;
            n = N;
        }

        NetworKit::count size = 0, n = 0;

        void run();

        NetworKit::count recurse_run(NetworKit::count n, NetworKit::count v);

        NetworKit::count BruteForceCliqueSize(NetworKit::Graph &Graph);
    };
}