#pragma once
#include "networkit/base/Algorithm.hpp"
#include "networkit/graph/EdgeUtils.hpp"
#include "networkit/graph/Graph.hpp"
#include <vector>

namespace Koala {
class PlanarSeparatorMatching : NetworKit::Algorithm {
  public:
    explicit PlanarSeparatorMatching(NetworKit::Graph &G);
    void run() override;

    std::vector<NetworKit::Edge> matching_set;

  private:
    std::vector<NetworKit::Edge> reduce_procedure(NetworKit::Graph &graph);

    NetworKit::Graph graph;
};
} // namespace Koala
