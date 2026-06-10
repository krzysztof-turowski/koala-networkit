#pragma once

#include <utility>
#include <vector>

#include <networkit/graph/Graph.hpp>

namespace Koala {
enum class NodeType {
    UNKNOWN,
    LEAF,
    UNION_NODE,
    COMPLEMENT_NODE
};

class Conode {
 public:
    NetworKit::node first_child, next_sibling, previous_sibling, parent;
    NetworKit::count size;
    NodeType type;
    NetworKit::node number;

    Conode(NetworKit::node child, NetworKit::node sibling, NetworKit::node p) noexcept {
        first_child = child;
        next_sibling = sibling;
        previous_sibling = NetworKit::none;
        parent = p;
        type = NodeType::UNKNOWN;
        size = 0;
        number = NetworKit::none;
    }
};

class Cotree {
 private:
    std::vector<Conode> nodes;
    NetworKit::node root;
 public:
    Cotree();

    bool prepared;

    void reserve(NetworKit::count n);

    NetworKit::node add(NodeType type, NetworKit::node number = NetworKit::none);

    void clear();

    void setRoot(NetworKit::node node);

    void addChild(NetworKit::node parent, NetworKit::node child);

    void removeChild(NetworKit::node parent, NetworKit::node child);

    void replaceChild(
            NetworKit::node parent, NetworKit::node old_child, NetworKit::node new_child);

    void moveChildToFront(NetworKit::node parent, NetworKit::node child);

    Conode& getNode(NetworKit::node i) {
        return nodes[i];
    }

    const Conode& getNode(NetworKit::node i) const {
        return nodes[i];
    }

    NetworKit::node getRoot() const;
    NetworKit::node upperNodeIdBound() const;
};
} /* namespace Koala */
