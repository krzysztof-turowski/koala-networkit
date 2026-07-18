#pragma once

#include <networkit/graph/AdjListGraph.hpp>

#include "techniques/BakerKOuterplanarGraphScheme.hpp"
#include "techniques/BakerPlanarApproximationScheme.hpp"
#include "vertex_cover/BakerVertexCover.hpp"
#include "vertex_cover/VertexCover.hpp"

namespace Koala {

/**
 * Baker PTAS for minimum vertex cover on planar graphs.
 */
class BakerPlanarGraphVertexCover final
    : public VertexCover,
      public BakerPlanarApproximationScheme<
          BakerKOuterplanarGraphScheme, BakerVertexCover> {
 public:
    BakerPlanarGraphVertexCover(
            const NetworKit::Graph &graph, double epsilon)
        : VertexCover(graph),
          BakerPlanarApproximationScheme<
              BakerKOuterplanarGraphScheme, BakerVertexCover>(
                  graph, epsilon) { }

    /**
     * Compute a Baker planar approximation.
     */
    void run() override;
};

}  // namespace Koala
