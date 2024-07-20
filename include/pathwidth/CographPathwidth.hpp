#pragma once

#include "recognition/CographRecognition.hpp"
#include "Pathwidth.hpp"

namespace Koala {

class CographPathwidth : public Pathwidth {
 private:
    void pathwidth();

    void subtree_size();

    std::stack<int> st;

    std::vector<bool> used;

    std::vector<NetworKit::count> path;
 public:
    Koala::Cotree &cotree;

    CographPathwidth(NetworKit::Graph &Graph, Koala::Cotree &CoTree)
        : Pathwidth(Graph), cotree(CoTree) { }

    void run();
};

} /* namespace Koala */
