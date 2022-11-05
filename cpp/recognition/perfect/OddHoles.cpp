/*
 * OddHoles.cpp
 *
 *  Created on: 05.11.2022
 *      Author: Adrian Siwiec
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <cassert>

#include <networkit/linkprediction/NeighborhoodUtility.hpp>

#include <recognition/PerfectGraphRecognition.hpp>

#include "../other/oddHoles.h"

namespace Koala {

bool PerfectGraphRecognition::containsOddHole(const NetworKit::Graph &graph) {
    std::vector<NetworKit::node> path;
    return Koala::Traversal::NextPathInplace(
        graph, std::numeric_limits<NetworKit::count>::max(), path, Koala::Traversal::PathInplaceMode::INDUCED_ODD_HOLE);
}

bool PerfectGraphRecognition::containsHole(const NetworKit::Graph &graph, NetworKit::count length) {
    if (length <= 3) {
        return false;
    }
    std::vector<NetworKit::node> path;
    return Koala::Traversal::NextPathInplace(graph, length, path, Koala::Traversal::PathInplaceMode::INDUCED_CYCLE);
}

bool PerfectGraphRecognition::containsT1(const NetworKit::Graph &graph) {
    //TODO(kturowski): temporary check
    Graph G(graph);
    assert(containsHole(graph, 5) == ::containsT1(G));
    if (containsHole(graph, 5)) {
        assert(containsOddHole(graph));
    }
    return containsHole(graph, 5);
}

}
