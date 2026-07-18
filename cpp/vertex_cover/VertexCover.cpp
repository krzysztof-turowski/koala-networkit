#include <cassert>

#include "vertex_cover/VertexCover.hpp"

namespace Koala {

VertexCover::VertexCover(const NetworKit::Graph &graph)
    : graph(std::make_optional(graph)) { }

const std::set<NetworKit::node> &VertexCover::getVertexCover() const {
    assureFinished();
    return vertex_cover;
}

void VertexCover::check() const {
    assureFinished();
    graph->forEdges([this](NetworKit::node u, NetworKit::node v) {
        assert(vertex_cover.contains(u) || vertex_cover.contains(v));
    });
}

}  // namespace Koala
