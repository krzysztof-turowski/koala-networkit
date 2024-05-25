#pragma once

#include "recognition/CographAlg.hpp"
#include "IndependentSet.hpp"
namespace Koala {
    class CographIndependentSet:public IndependentSet {
    private:
        std::set<NetworKit::node> recurse_run(NetworKit::count v);

        Koala::Cotree &cotree;

    public:

        CographIndependentSet(NetworKit::Graph &Graph, Koala::Cotree &Cotree) : IndependentSet(Graph), cotree(Cotree)
        {

        };

        void run();

        NetworKit::count bruteForceIndependetSetSize(NetworKit::Graph &Graph);
    };
}
