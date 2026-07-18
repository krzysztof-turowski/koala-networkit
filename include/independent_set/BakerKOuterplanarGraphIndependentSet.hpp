#pragma once

#include <networkit/graph/AdjListGraph.hpp>

#include "independent_set/BakerIndependentSet.hpp"
#include "independent_set/IndependentSet.hpp"
#include "techniques/BakerKOuterplanarGraphScheme.hpp"

namespace Koala {

/**
 * Exact maximum independent set algorithm for k-outerplanar graphs.
 *
 * The planar embedding, levels, fake triangulation edges, Baker forest, and
 * boundary-state dynamic program are constructed when run() is called.
 */
class BakerKOuterplanarGraphIndependentSet final
    : public IndependentSet,
      public BakerKOuterplanarGraphScheme<BakerIndependentSet> {
 public:
    explicit BakerKOuterplanarGraphIndependentSet(const NetworKit::Graph &graph)
        : IndependentSet(graph),
          BakerKOuterplanarGraphScheme<BakerIndependentSet>() { }

    /**
     * Compute a maximum independent set.
     */
    void run() override;
};

}  // namespace Koala
