#pragma once

#include "recognition/CographRecognition.hpp"
#include "MaxClique.hpp"

namespace Koala {
class CographMaxClique : public MaxClique {
 private:
    std::set<NetworKit::node> recurse_run(NetworKit::count v);

 public:
    Koala::Cotree &cotree;

    CographMaxClique(NetworKit::Graph &Graph, Koala::Cotree &CoTree) : MaxClique(Graph), cotree(CoTree) {
    }

    void run() override;

    NetworKit::count bruteForceCliqueSize(NetworKit::Graph &Graph);
};
} /* namespace Koala */
