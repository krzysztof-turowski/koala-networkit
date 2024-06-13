#pragma once

#include "recognition/CographRecognition.hpp"
#include "IndependentSet.hpp"
namespace Koala {
class CographIndependentSet:public IndependentSet {
 private:
    NetworKit::count recurse_run(NetworKit::count v);

    void add_to_set(NetworKit::count v);

    Koala::Cotree &cotree;

    std::vector<NetworKit::count > independent_set_size;

 public:
    CographIndependentSet(NetworKit::Graph &Graph, Koala::Cotree &Cotree) : IndependentSet(Graph), cotree(Cotree) {
    }

    void run();

    NetworKit::count bruteForceIndependetSetSize(NetworKit::Graph &Graph);
};
} /* namespace Koala */

