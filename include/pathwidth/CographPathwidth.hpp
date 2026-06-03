#pragma once

#include <stack>
#include <vector>

#include "pathwidth/Pathwidth.hpp"
#include "structures/Cotree.hpp"

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

    CographPathwidth(NetworKit::Graph &Graph, Koala::Cotree &cotree_ref)
        : Pathwidth(Graph), cotree(cotree_ref) { }

    void run();
};

} /* namespace Koala */
