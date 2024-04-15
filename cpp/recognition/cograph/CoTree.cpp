//
// Created by milana on 30.03.24.
//
#include <list>
#include <graph/GraphTools.hpp>

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

    class CoNode {
    public:
        Type type;
        int number;
        Marked marked;
        //d is the current number of children
        //md is the current number of children, which have been both "marked" and "unmarked"
        int md, d;
        bool in_graph;
        CoNode *first_child;
        CoNode *next, *previous;//in list of children of its parent
        CoNode *parent;
        std::vector<CoNode*> out_edges;//neighbours of cur vertex in G

        CoNode(Type type, int number) : type(type), number(number), marked(Marked::UNMARKED), md(0), d(0),
                                        in_graph(false), first_child(nullptr),
                                        next(nullptr), previous(nullptr), parent(nullptr) {

        }

        void AddChild(CoNode *x) {
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

        void UnmarkForNewIteration() {
            marked = Marked::UNMARKED;
            md = 0;
        }

        void mark() {
            marked = Marked::MARKED;
        }

        void unmark() {
            marked = Marked::MARKED_AND_UNMARKED;
        }

        std::vector<CoNode*> RemoveWereMarked() {
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

        void RemoveWereNotMarked() {
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
    };

    class CoTree {
    private:
        std::vector<CoNode*> save;
    public:
        CoNode *root;

        CoNode* Add(Type type, int number) {
            CoNode *x = new CoNode(type, number);
            save.push_back(x);
            return x;
        }

        void Clear() {
            for (auto u: save) {
                delete u;
            }
            save.clear();
        }
    };
}
