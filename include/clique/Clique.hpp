#pragma once

#include <optional>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {
class MaxClique : public NetworKit::Algorithm {
 public:
    explicit MaxClique(const NetworKit::Graph &graph);

    std::vector<NetworKit::node> &getMaxClique();

    virtual void run() = 0;

    void check() const;

 protected:
    std::vector<NetworKit::node> max_clique;
    std::optional<NetworKit::Graph> graph;
};

} /* namespace Koala */

