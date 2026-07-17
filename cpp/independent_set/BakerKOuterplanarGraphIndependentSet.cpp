#include "independent_set/BakerKOuterplanarGraphIndependentSet.hpp"

namespace Koala {

void BakerKOuterplanarGraphIndependentSet::run() {
    independentSet =
        BakerKOuterplanarGraphScheme<BakerIndependentSet>::solve(*graph);
    hasRun = true;
}

}  // namespace Koala
