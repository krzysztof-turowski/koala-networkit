#pragma once

#include "IndependentSet.hpp"

#include "recognition/CographRecognition.hpp"

namespace Koala {

class CographIndependentSet : public IndependentSet {
 private:
    void recurse_run();

    void add_to_set();

    Koala::Cotree &cotree;

    std::vector<NetworKit::count> independent_set_size;

    std::vector<bool> used;

    std::stack<int> st;
 public:
    CographIndependentSet(NetworKit::Graph &Graph, Koala::Cotree &Cotree)
        : IndependentSet(Graph), cotree(Cotree) { }

    void run();

    NetworKit::count bruteForceIndependetSetSize(NetworKit::Graph &Graph);
};

} /* namespace Koala */

