#pragma once

#include <optional>
#include <set>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {
class MaxClique : public NetworKit::Algorithm {
 public:
    explicit MaxClique(const NetworKit::Graph &graph);

    std::set<NetworKit::node> &getMaxCliqueSet();

    virtual void run() = 0;

    void check() const;

 protected:
    std::set<NetworKit::node> max_clique;
    std::optional<NetworKit::Graph> graph;
};

} /* namespace Koala */

