#include "vertex_cover/BakerKOuterplanarGraphVertexCover.hpp"

namespace Koala {

void BakerKOuterplanarGraphVertexCover::run() {
    vertex_cover =
        BakerKOuterplanarGraphScheme<BakerVertexCover>::solve(*graph);
    hasRun = true;
}

}  // namespace Koala
