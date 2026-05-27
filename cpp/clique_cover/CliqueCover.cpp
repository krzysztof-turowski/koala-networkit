#include <graph/GraphTools.hpp>

#include "clique_cover/CliqueCover.hpp"

namespace Koala {

MinCliqueCover::MinCliqueCover(const NetworKit::Graph &graph)
        : graph(std::make_optional(graph)) {}

const std::vector<std::set<NetworKit::node>>& MinCliqueCover::getCliqueCover() const {
    return clique_cover;
}

}