#include "vertex_cover/BakerPlanarGraphVertexCover.hpp"

namespace Koala {

void BakerPlanarGraphVertexCover::run() {
    vertex_cover = BakerPlanarApproximationScheme<
        BakerKOuterplanarGraphScheme, BakerVertexCover>::solve();
    hasRun = true;
}

}  // namespace Koala
