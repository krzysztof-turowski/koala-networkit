#pragma once

#include <optional>
#include <set>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/AdjListGraph.hpp>

namespace Koala {

/**
 * Base class for vertex-cover algorithms.
 */
class VertexCover : public NetworKit::Algorithm {
 public:
    explicit VertexCover(const NetworKit::Graph &graph);

    /**
     * Return the computed vertex cover.
     */
    const std::set<NetworKit::node> &getVertexCover() const;

    /**
     * Execute the vertex-cover algorithm.
     */
    virtual void run() = 0;

    /**
     * Verify that every original edge is covered.
     */
    void check() const;

 protected:
    std::optional<NetworKit::Graph> graph;
    std::set<NetworKit::node> vertex_cover;
};

}  // namespace Koala
