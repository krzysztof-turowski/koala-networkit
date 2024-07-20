#include "pathwidth/CographPathwidth.hpp"

namespace Koala {
    Pathwidth::Pathwidth(NetworKit::Graph &Graph)
            : graph(Graph) {}

    NetworKit::count Pathwidth::getPathwidthSize() {
        assureFinished();
        return width;
    }
} /* namespace Koala */
