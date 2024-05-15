#pragma once

#include "recognition/CographAlg.hpp"

namespace Koala {
    class CographIndependentSet {
    public:
        Koala::Cotree *cotree;

        CographIndependentSet(Koala::Cotree &CoTree, NetworKit::count N) {
            cotree = &CoTree;
            n = N;
        }

        std::vector<NetworKit::count> color, number_of_colors;
        NetworKit::count independet_set_size = 0, n;

        void run();

        NetworKit::count recurse_run(NetworKit::count n, NetworKit::count v);

        NetworKit::count BruteForceIndependetSetSize(NetworKit::Graph &Graph);
    };
}
