#include "independent_set/BakerPlanarGraphIndependentSet.hpp"

namespace Koala {

void BakerPlanarGraphIndependentSet::run() {
    independentSet = BakerPlanarApproximationScheme<
        BakerKOuterplanarGraphScheme, BakerIndependentSet>::solve();
    hasRun = true;
}

}  // namespace Koala
