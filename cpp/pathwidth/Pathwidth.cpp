#include "pathwidth/CographPathwidth.hpp"

namespace Koala {
    Pathwidth::Pathwidth(const NetworKit::Graph &graph)
            : graph(std::make_optional(graph)) { }

    NetworKit::count Pathwidth::getPathwidthSize()
    {
        assureFinished();
        return width;
    }
}