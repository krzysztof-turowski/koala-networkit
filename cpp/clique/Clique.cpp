#include <graph/GraphTools.hpp>

#include "clique/Clique.hpp"

#include "recognition/CographRecognitionOther.hpp"

namespace Koala {

MaxClique::MaxClique(const NetworKit::Graph &graph)
        : graph(std::make_optional(graph)) {}

std::set<NetworKit::node>& MaxClique::getMaxCliqueSet() {
    assureFinished();
    return max_clique;
}

void MaxClique::check() const {
    assureFinished();
    for (auto x : max_clique) {
        for (auto y : max_clique) {
            assert(x == y || graph->hasEdge(x, y));
        }
    }
}

} /* namespace Koala */
