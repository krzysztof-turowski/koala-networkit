#include "dominating_set/BakerKOuterplanarGraphDominatingSet.hpp"

namespace Koala {

void BakerKOuterplanarGraphDominatingSet::run() {
    dominating_set =
        BakerKOuterplanarGraphScheme<BakerDominatingSet>::solve(*graph);
    hasRun = true;
}

}  // namespace Koala
