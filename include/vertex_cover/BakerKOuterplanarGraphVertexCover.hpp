#pragma once

#include <networkit/graph/AdjListGraph.hpp>

#include "techniques/BakerKOuterplanarGraphScheme.hpp"
#include "vertex_cover/BakerVertexCover.hpp"
#include "vertex_cover/VertexCover.hpp"

namespace Koala {

/**
 * Exact minimum vertex cover algorithm for k-outerplanar graphs.
 */
class BakerKOuterplanarGraphVertexCover final
    : public VertexCover,
      public BakerKOuterplanarGraphScheme<BakerVertexCover> {
 public:
    explicit BakerKOuterplanarGraphVertexCover(
            const NetworKit::Graph &graph)
        : VertexCover(graph),
          BakerKOuterplanarGraphScheme<BakerVertexCover>() { }

    void run() override;
};

}  // namespace Koala
