//
// Created by milana on 30.03.24.
// Copyright 2024 milana
//
#include <list>
#include <graph/GraphTools.hpp>
#include "recognition/CoTree.hpp"

namespace Koala{
    CoNode::CoNode(Type type, int number) : type(type), number(number),
    marked(Marked::UNMARKED), md(0), d(0), in_graph(false),
    first_child(nullptr), next(nullptr), previous(nullptr),
    parent(nullptr) {
    }

    void CoNode::AddChild(CoNode *x) {
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

    void CoNode::UnmarkForNewIteration() {
        marked = Marked::UNMARKED;
        md = 0;
    }

    void CoNode::mark() {
        marked = Marked::MARKED;
    }

    void CoNode::unmark() {
        marked = Marked::MARKED_AND_UNMARKED;
    }

    std::vector<CoNode*> CoNode::RemoveWereMarked() {
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

    void CoNode::RemoveWereNotMarked() {
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

    CoNode* CoTree::Add(Type type, int number) {
        CoNode *x = new CoNode(type, number);
        save.push_back(x);
        return x;
    }

    void CoTree::Clear() {
        for (auto u : save) {
            delete u;
        }
        save.clear();
    }
}  // namespace Koala
