#pragma once

#include "recognition/CographAlg.hpp"
#include "IndependentSet.hpp"
namespace Koala {
    class CographIndependentSet {
    private:
        NetworKit::count recurse_run(NetworKit::count v);

        Koala::Cotree &cotree;
    public:
        CographIndependentSet(Koala::Cotree &CoTree) : cotree(CoTree){

        }

        NetworKit::count independet_set_size = 0;

        void run();

        NetworKit::count bruteForceIndependetSetSize(NetworKit::Graph &Graph);
    };
}
