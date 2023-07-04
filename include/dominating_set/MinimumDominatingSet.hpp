#pragma once

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

class MinimumDominatingSet : public NetworKit::Algorithm {
 protected:
    const NetworKit::Graph *G;
    std::vector<bool> dominatingSet;
 public:
    explicit MinimumDominatingSet(const NetworKit::Graph &G);

    const std::vector<bool> &getDominatingSet() {
        assureFinished();
        return dominatingSet;
    }

    bool isDominating(const std::vector<bool> &dominating_set);
    static int dominatingSetSize(const std::vector<bool> &set);

    static int specialCounter1;
    static int specialCounter2;
};

std::vector<bool> &smallerCardinalitySet(std::vector<bool> &lhs, std::vector<bool> &rhs);
