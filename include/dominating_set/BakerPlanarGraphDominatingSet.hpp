#pragma once

#include <networkit/graph/AdjListGraph.hpp>

#include "dominating_set/BakerDominatingSet.hpp"
#include "dominating_set/DominatingSet.hpp"
#include "techniques/BakerKOuterplanarGraphScheme.hpp"
#include "techniques/BakerPlanarApproximationScheme.hpp"

namespace Koala {

/**
 * Baker PTAS for minimum dominating set on planar graphs.
 */
class BakerPlanarGraphDominatingSet final
    : public DominatingSet,
      public BakerPlanarApproximationScheme<
          BakerKOuterplanarGraphScheme, BakerDominatingSet> {
 public:
    BakerPlanarGraphDominatingSet(
            NetworKit::Graph &graph, double epsilon)
        : DominatingSet(graph),
          BakerPlanarApproximationScheme<
              BakerKOuterplanarGraphScheme, BakerDominatingSet>(
                  graph, epsilon) { }

    /**
     * Compute a Baker planar approximation.
     */
    void run() override;
};

}  // namespace Koala
