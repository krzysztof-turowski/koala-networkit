#include "dominating_set/BakerPlanarGraphDominatingSet.hpp"

namespace Koala {

void BakerPlanarGraphDominatingSet::run() {
    dominating_set = BakerPlanarApproximationScheme<
        BakerKOuterplanarGraphScheme, BakerDominatingSet>::solve();
    hasRun = true;
}

}  // namespace Koala
