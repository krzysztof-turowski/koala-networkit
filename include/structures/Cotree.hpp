#pragma once

#include <array>
#include <utility>
#include <vector>

#include "graph/GraphTools.hpp"

namespace Koala {
enum class NodeType {
    UNKNOWN,
    LEAF,
    UNION_NODE,
    COMPLEMENT_NODE
};

enum class Marked {
    UNMARKED,
    MARKED,
    MARKED_AND_UNMARKED
};

class Conode {
 public:
    NetworKit::count first_child, next_sibling, previous_sibling, parent, size;
    NodeType type;
    int number;
    Marked marked;
    int md, d;
    bool in_graph;
    std::vector<NetworKit::count> out_edges;
    int number_of_vertices_in_subtree, time_in, time_out;
    std::array<NetworKit::count, 64> get_up;

    Conode(NetworKit::count child, NetworKit::count sibling, NetworKit::count p) {
        first_child = child;
        next_sibling = sibling;
        previous_sibling = NetworKit::none;
        parent = p;
        type = NodeType::UNKNOWN;
        size = 0;
        number = 0;
        marked = Marked::UNMARKED;
        md = 0;
        d = 0;
        in_graph = false;
        number_of_vertices_in_subtree = 0;
        time_in = 0;
        time_out = 0;
        get_up.fill(NetworKit::none);
    }
};

class Cotree {
 private:
    std::vector<Conode> nodes;
    std::vector<std::pair<std::pair<NetworKit::count, NetworKit::count>, NetworKit::count> > order;
    NetworKit::count root;
 public:
    NetworKit::Graph *graph;

    Cotree();

    explicit Cotree(NetworKit::Graph &Graph);

    bool prepared;

    void buildTree();

    void reserve(NetworKit::count n);

    NetworKit::count add(NodeType type, int number);

    void clear();

    void setRoot(NetworKit::count node);

    void addChild(NetworKit::count parent, NetworKit::count child);

    void removeChild(NetworKit::count parent, NetworKit::count child);

    void replaceChild(
            NetworKit::count parent, NetworKit::count old_child, NetworKit::count new_child);

    void moveChildToFront(NetworKit::count parent, NetworKit::count child);

    void unmarkForNewIteration(NetworKit::count node);

    void mark(NetworKit::count node);

    void unmark(NetworKit::count node);

    std::vector<NetworKit::count> removeWereMarked(NetworKit::count node);

    void removeWereNotMarked(NetworKit::count node);

    void setOrder(
            std::vector<std::pair<std::pair<
            NetworKit::count, NetworKit::count>, NetworKit::count> > a) {
        order = a;
    }

    Conode& getNode(NetworKit::count i) {
        return nodes[i];
    }

    NetworKit::count getRoot() const;
    NetworKit::count upperNodeIdBound() const;
};
} /* namespace Koala */
