#pragma once

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

enum class Type {
    ZERO_ONE,
    VERTEX
};

enum class Marked {
    UNMARKED,
    MARKED,
    MARKED_AND_UNMARKED
};

class CoNode {
 public:
    Type type;
    int number;
    Marked marked;
    // d is the current number of children
    // md is the current number of children, which have been both "marked" and "unmarked"
    int md, d;
    bool in_graph;
    CoNode *first_child;
    CoNode *next, *previous;  // in list of children of its parent
    CoNode *parent;
    std::vector<CoNode*> out_edges;  // neighbours of current vertex in G
    int number_of_vertices_in_subtree = 0, time_in = 0, time_out = 0;
    CoNode *get_up[30];

    explicit CoNode(Type type, int number);

    void AddChild(CoNode *x);
    void UnmarkForNewIteration();

    void mark();
    void unmark();

    std::vector<CoNode *> RemoveWereMarked();

    void RemoveWereNotMarked();
};

class CoTree {
 private:
    std::vector<CoNode> save;
 public:
    CoNode *root;

    void ReserveSpace(int n);
    CoNode *Add(Type type, int number);
    void Clear();
};

inline CoNode::CoNode(Type type, int number)
    : type(type), number(number), marked(Marked::UNMARKED), md(0), d(0), in_graph(false),
        first_child(nullptr), next(nullptr), previous(nullptr), parent(nullptr) {}

inline void CoNode::AddChild(CoNode *x) {
    if (first_child == nullptr) {
        first_child = x;
        x->previous = nullptr;
        x->next = nullptr;
    } else {
        first_child->previous = x;
        x->next = first_child;
        x->previous = nullptr;
        first_child = x;
    }
    x->parent = this;
    d++;
}

inline void CoNode::UnmarkForNewIteration() {
    marked = Marked::UNMARKED;
    md = 0;
}

inline void CoNode::mark() {
    marked = Marked::MARKED;
}

inline void CoNode::unmark() {
    marked = Marked::MARKED_AND_UNMARKED;
}

inline std::vector<CoNode*> CoNode::RemoveWereMarked() {
    auto u = first_child;
    std::vector<CoNode*> vec;
    while (u != nullptr) {
        vec.push_back(u);
        d--;
        first_child = u->next;
        if (first_child != nullptr) {
            first_child->previous = nullptr;
        }
        u->previous = nullptr;
        u->next = nullptr;
        u = first_child;
        if (u == nullptr || u->marked != Marked::MARKED_AND_UNMARKED) {
            break;
        }
    }
    return vec;
}

inline void CoNode::RemoveWereNotMarked() {
    auto u = first_child;
    while (u != nullptr && u->marked == Marked::MARKED_AND_UNMARKED) {
        u = u->next;
    }
    while (u != nullptr) {
        d--;
        auto next = u->next;
        auto previous = u->previous;
        if (next != nullptr) {
            next->previous = previous;
        }
        if (previous != nullptr) {
            previous->next = next;
        }
        u->previous = nullptr;
        u->next = nullptr;
        u = next;
    }
}

inline CoNode* CoTree::Add(Type type, int number) {
    return &save.emplace_back(type, number);
}

inline void CoTree::Clear() {
    save.clear();
}

inline void CoTree::ReserveSpace(int n) {
    save.reserve(n);
}
} /* namespace Koala */
