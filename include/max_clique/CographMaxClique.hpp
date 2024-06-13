#pragma once

#include "recognition/CographRecognition.hpp"
#include "MaxClique.hpp"

namespace Koala {
class CographMaxClique : public MaxClique {
 private:
    NetworKit::count  recurse_run(NetworKit::count v);

    void add_to_set(NetworKit::count v);

    Koala::Cotree &cotree;

    std::vector<NetworKit::count > subgraph_clique_size;

 public:

    CographMaxClique(NetworKit::Graph &Graph, Koala::Cotree &CoTree) : MaxClique(Graph), cotree(CoTree) {
    }

    void run() override;

    NetworKit::count bruteForceCliqueSize(NetworKit::Graph &Graph);
};
} /* namespace Koala */
