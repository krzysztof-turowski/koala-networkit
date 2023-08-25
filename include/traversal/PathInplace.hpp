/*
 * PathInplace.hpp
 *
 *  Created on: 22.12.2022
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <networkit/graph/Graph.hpp>

namespace Koala {

namespace Traversal {

enum class PathInplaceMode {
    INDUCED_PATH, INDUCED_CYCLE, INDUCED_ODD_HOLE
};

bool NextPathInplace(
        const NetworKit::Graph &graph, NetworKit::count length, std::vector<NetworKit::node> &path,
        PathInplaceMode mode);

}  // namespace Traversal

}  // namespace Koala
