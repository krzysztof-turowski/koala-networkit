#ifndef MINIMUM_DOMINATING_SET_HPP_
#define MINIMUM_DOMINATING_SET_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

class MinimumDominatingSet : public NetworKit::Algorithm {
protected:
    const NetworKit::Graph *G;
    std::vector<bool> dominatingSet;
public:
    MinimumDominatingSet(const NetworKit::Graph &G);

    // void run() override = 0;

    const std::vector<bool> &getDominatingSet() {
        assureFinished();
        return dominatingSet;
    }

    bool isDominating(const std::vector<bool> &dominating_set);
};

std::vector<bool> &smallerCardinalitySet(std::vector<bool> &lhs, std::vector<bool> &rhs);

#endif /* MINIMUM_DOMINATING_SET_HPP_ */
