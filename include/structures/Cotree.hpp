#pragma once

#include <graph/GraphTools.hpp>

namespace Koala {
enum class NodeType {
    UNKNOWN,
    LEAF,
    UNION_NODE,
    COMPLEMENT_NODE
};

class CoNode {
 public:
    NetworKit::count left_son, right_son, parent, size;
    NodeType type;

    CoNode(NetworKit::count l, NetworKit::count r, NetworKit::count p) {
        left_son = l;
        right_son = r;
        parent = p;
        type = NodeType::UNKNOWN;
        size = 0;
    }

};

class Cotree {
 private:
    std::vector<CoNode> nodes;
    std::vector<std::pair<std::pair<NetworKit::count, NetworKit::count>, NetworKit::count> > order;
 public:
    NetworKit::Graph *graph;

    explicit Cotree(NetworKit::Graph &Graph);

    bool prepared;

    void buildTree();

    void setOrder(std::vector<std::pair<std::pair<NetworKit::count, NetworKit::count>, NetworKit::count> > a) {
        order = a;
    }

    CoNode& getNode(NetworKit::count i) {
        return nodes[i];
    }
};
} /* namespace Koala */
