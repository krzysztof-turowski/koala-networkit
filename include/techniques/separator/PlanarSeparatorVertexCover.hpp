#pragma once
#include "networkit/base/Algorithm.hpp"
#include "networkit/graph/Graph.hpp"

namespace Koala {
struct Bipartite {
    std::vector<int> X;
    std::vector<int> Y;
};

class PlanarSeparatorVertexCover : NetworKit::Algorithm {
  public:
    explicit PlanarSeparatorVertexCover(NetworKit::Graph &G);
    void run() override;
    std::vector<NetworKit::node> getVertexCover() { return vertexCover; };

  private:
    std::vector<NetworKit::node> vertexCover;
    NetworKit::Graph graph;
    NetworKit::Graph residual;
    void prep(std::vector<bool> &U, std::vector<bool> &VC, int n);

    Bipartite bipartite(std::vector<bool> &U, std::vector<int> &deg_residual);
};
} // namespace Koala
