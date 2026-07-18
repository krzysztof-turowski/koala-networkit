#pragma once

#include <networkit/graph/AdjListGraph.hpp>

#include "dominating_set/BakerDominatingSet.hpp"
#include "dominating_set/DominatingSet.hpp"
#include "techniques/BakerKOuterplanarGraphScheme.hpp"

namespace Koala {

/**
 * Exact minimum dominating set algorithm for k-outerplanar graphs.
 */
class BakerKOuterplanarGraphDominatingSet final
    : public DominatingSet,
      public BakerKOuterplanarGraphScheme<BakerDominatingSet> {
 public:
    explicit BakerKOuterplanarGraphDominatingSet(
            NetworKit::Graph &graph)
        : DominatingSet(graph),
          BakerKOuterplanarGraphScheme<BakerDominatingSet>() { }

    void run() override;
};

}  // namespace Koala
