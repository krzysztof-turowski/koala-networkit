#pragma once

#include <graph/GraphTools.hpp>

namespace Koala {
enum class NodeType {
    UNKNOWN,
    LEAF,
    UNION_NODE,
    COMPLEMENT_NODE
};

class Conode {
 public:
    NetworKit::count left_son, right_son, parent, size;
    NodeType type;

    Conode(NetworKit::count l, NetworKit::count r, NetworKit::count p) {
        left_son = l;
        right_son = r;
        parent = p;
        type = NodeType::UNKNOWN;
        size = 0;
    }
};

class Cotree {
 private:
    std::vector<Conode> nodes;
    std::vector<std::pair<std::pair<NetworKit::count, NetworKit::count>, NetworKit::count> > order;
 public:
    NetworKit::Graph *graph;

    explicit Cotree(NetworKit::Graph &Graph);

    bool prepared;

    void buildTree();

    void setOrder(std::vector<std::pair<std::pair<NetworKit::count, NetworKit::count>, NetworKit::count> > a) {
        order = a;
    }

    Conode& getNode(NetworKit::count i) {
        return nodes[i];
    }
};
} /* namespace Koala */
