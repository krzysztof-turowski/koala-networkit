#pragma once

#include <networkit/graph/AdjListGraph.hpp>

#include "independent_set/BakerIndependentSet.hpp"
#include "independent_set/IndependentSet.hpp"
#include "techniques/BakerKOuterplanarGraphScheme.hpp"
#include "techniques/BakerPlanarApproximationScheme.hpp"

namespace Koala {

/**
 * Baker PTAS for maximum independent set on planar graphs.
 */
class BakerPlanarGraphIndependentSet final
    : public IndependentSet,
      public BakerPlanarApproximationScheme<
          BakerKOuterplanarGraphScheme, BakerIndependentSet> {
 public:
    BakerPlanarGraphIndependentSet(
            const NetworKit::Graph &graph, double epsilon)
        : IndependentSet(graph),
          BakerPlanarApproximationScheme<
              BakerKOuterplanarGraphScheme, BakerIndependentSet>(
                  graph, epsilon) { }

    /**
     * Compute a Baker planar approximation.
     */
    void run() override;
};

}  // namespace Koala
