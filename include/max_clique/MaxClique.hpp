#pragma once

#include "recognition/CographAlg.hpp"

namespace Koala {
    class MaxClique {
    private:
        NetworKit::count recurse_run(NetworKit::count v);
    public:
        Koala::Cotree &cotree;

        MaxClique(Koala::Cotree &CoTree) : cotree(CoTree){

        }

        NetworKit::count size = 0;

        void run();

        NetworKit::count bruteForceCliqueSize(NetworKit::Graph &Graph);
    };
}