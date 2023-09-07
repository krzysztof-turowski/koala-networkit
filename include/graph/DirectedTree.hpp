/*
 * DirectedTree.hpp
 *
 *  Created on: 23.08.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <optional>

#include <networkit/graph/Graph.hpp>

namespace Koala {

class DirectedTree {
 public:
    explicit DirectedTree(const NetworKit::Graph& G) : tree(G.upperNodeIdBound(), true, true) { }

    NetworKit::Graph& getTree() { return tree; }

    NetworKit::node getRoot() const {
        return tree.upperNodeIdBound() - 1;
    }

    std::optional<NetworKit::node> getParent(NetworKit::node u) const {
        return u != getRoot() ?
            std::make_optional(*(tree.inNeighborRange(u).begin())) : std::nullopt;
    }

    NetworKit::edgeweight edgeWeightToParent(NetworKit::node u) const {
        return u != getRoot() ? tree.weight(*getParent(u), u) : 0;
    }
 private:
    NetworKit::Graph tree;
};

}  // namespace Koala
