#pragma once

#include "Clique.hpp"

#include "recognition/CographRecognitionOther.hpp"

namespace Koala {

class CographMaxClique : public MaxClique {
 private:
    void recurse_run();

    void add_to_set();

    Koala::Cotree &cotree;

    std::vector<NetworKit::count> subgraph_clique_size;

    std::vector<bool> used;

    std::stack<int> st;

 public:
    CographMaxClique(NetworKit::Graph &Graph, Koala::Cotree &CoTree)
        : MaxClique(Graph), cotree(CoTree) { }

    void run() override;

    NetworKit::count bruteForceCliqueSize(NetworKit::Graph &Graph);
};

} /* namespace Koala */
