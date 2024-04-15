//
// Created by milana on 15.04.24.
// Copyright 2024 milana
//
#ifndef INCLUDE_RECOGNITION_COTREE_HPP_
#define INCLUDE_RECOGNITION_COTREE_HPP_

#include <vector>

namespace Koala {
    enum class Type {
        ZERO_ONE,
        VERTEX
    };
    enum class Marked {
        UNMARKED,
        MARKED,
        MARKED_AND_UNMARKED
    };
class CoNode{
 public:
        Type type;
        int number;
        Marked marked;
        // d is the current number of children
        // md is the current number of children,
        // which have been both "marked" and "unmarked"
        int md, d;
        bool in_graph;
        CoNode *first_child;
        CoNode *next, *previous;  // in list of children of its parent
        CoNode *parent;
        std::vector<CoNode*> out_edges;  // neighbours of cur vertex in G

        explicit CoNode(Type type, int number);
        void AddChild(CoNode *x);
        void UnmarkForNewIteration();
        void mark();
        void unmark();
        std::vector<CoNode*> RemoveWereMarked();
        void RemoveWereNotMarked();
};

class CoTree {
 private:
        std::vector<CoNode*> save;
 public:
        CoNode *root;

        CoNode* Add(Type type, int number);
        void Clear();
};
}  // namespace Koala

#endif  // INCLUDE_RECOGNITION_COTREE_HPP_
