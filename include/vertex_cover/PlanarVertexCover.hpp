#pragma once

#include <vector>
#include "vertex_cover/VertexCover.hpp"
#include <networkit/graph/AdjListGraph.hpp>
namespace Koala {
struct Bipartite {
    std::vector<int> X;
    std::vector<int> Y;
};

class PlanarSeparatorVertexCover : public VertexCover {
public:
    explicit PlanarSeparatorVertexCover(const NetworKit::Graph &G);
    void run() override;

private:
    void prepare(const NetworKit::Graph &G, std::vector<bool> &U, std::vector<bool> &VC, int n);
    Bipartite bipartite(const NetworKit::Graph &G, std::vector<bool> &U,
                        std::vector<int> &deg_residual);

    NetworKit::Graph residual;
};
} // namespace Koala
