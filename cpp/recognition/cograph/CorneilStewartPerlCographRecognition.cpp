/*
 * CorneilStewartPerlCographRecognition.cpp
 *
 *  Created on: 2024
 *      Author: fixikmila
 */
// Copyright 2024 milana

#include <graph/GraphTools.hpp>

#include "recognition/CographRecognition.hpp"
#include "recognition/CoTree.hpp"

namespace Koala {

    void CorneilStewartPerlCographRecognition::
    Unmark() {
        CoNode *u = marked_with_d_equal_to_md.front();
        marked_with_d_equal_to_md.pop();
        u->unmark();
        mark_count--;
        mark_and_unmarked_count++;
        u->md = 0;
        if (u != T.root) {
            auto w = u->parent;
            w->md++;
            if (w->marked == Marked::UNMARKED) {
                mark_count++;
            }
            w->mark();
            if (w->md == w->d) {
                marked_with_d_equal_to_md.push(w);
            }
            auto next = u->next;
            auto previous = u->previous;
            auto head = w->first_child;
            if (previous != nullptr) {
                previous->next = next;
                if (next != nullptr) {
                    next->previous = previous;
                }
                u->previous = nullptr;
                u->next = head;
                head->previous = u;
                w->first_child = u;
            }  // else u is head
        }
    }

    void CorneilStewartPerlCographRecognition::Mark(CoNode *x) {
        mark_count = 0;
        mark_and_unmarked_count = 0;
        mark_ever_count = 0;
        for (auto u : x->out_edges) {
            // !!only neighbours which are already in graph
            if (!(u->in_graph)) {
                continue;
            }
            u->mark();
            mark_ever_count++;
            mark_count++;
            marked_with_d_equal_to_md.push(u);
        }
        while (!marked_with_d_equal_to_md.empty()) {
            Unmark();
        }
        if (mark_count && T.root->d == 1) {
            T.root->mark();
        }
    }

    void ResetAllCoNodes(CoNode *x, int level = 0) {
        x->UnmarkForNewIteration();
        CoNode *y = x->first_child;
        while (y != nullptr) {
            ResetAllCoNodes(y, level + 1);
            y = y->next;
        }
    }

    std::pair<CoNode*, CorneilStewartPerlCographRecognition::State>
    CorneilStewartPerlCographRecognition::FindLowest() {
        auto *y = T.Add(Type::ZERO_ONE, 2);
        if (T.root->marked == Marked::UNMARKED) {
            return {y, CorneilStewartPerlCographRecognition::State::
            GRANDPARENT_IS_NOT_IN_SET};
        }
        if (T.root->md != T.root->d - 1) {
            y = T.root;
        }
        T.root->unmark();
        T.root->md = 0;
        CoNode *w = T.root;
        std::queue<CoNode*> q;
        std::stack<CoNode*> s;
        s.push(T.root);
        while (!s.empty()) {
            auto x = s.top();
            s.pop();
            if (x->marked == Marked::MARKED) {
                q.push(x);
            }
            auto z = x->first_child;
            while (z != nullptr) {
                s.push(z);
                z = z->next;
            }
        }
        while (!q.empty()) {
            CoNode *u = q.front();
            q.pop();
            if (u->marked != Marked::MARKED) {
                continue;
            }
            if (y->number != 2) {  // 1 or 2
                if (y->number == 0) {
                    return {y, CorneilStewartPerlCographRecognition::State::
                    CONTAINS_0_NODE};
                } else {
                    return {y, CorneilStewartPerlCographRecognition::State::
                    EXISTS_1_NODE_NOT_PROPERLY_MARKED};
                }
            }
            CoNode *t;
            if (u->number == 1) {
                if (u->md != u->d - 1) {
                    y = u;
                }
                if (u->parent->marked == Marked::MARKED) {  // 1 or 6
                    if (y->number == 0) {
                        return {y, CorneilStewartPerlCographRecognition::State::
                        CONTAINS_0_NODE};
                    } else {
                        return {y, CorneilStewartPerlCographRecognition::State::
                        WRONG_GRANDPARENT};
                    }
                } else {
                    t = u->parent->parent;
                }
            } else {
                y = u;
                t = u->parent;
            }
            u->unmark();
            u->md = 0;
            while (t != w) {
                if (t == T.root) {  // 4
                    return {y, CorneilStewartPerlCographRecognition::State::
                    NO_ONE_PATH};
                }
                if (t->marked != Marked::MARKED) {  // 3 or 5 or 6
                    if (y->number == 0) {
                        return {y, CorneilStewartPerlCographRecognition::State::
                        WRONG_PARENT};
                    } else {
                        return {y, CorneilStewartPerlCographRecognition::State::
                                WRONG_GRANDPARENT};
                        // if y is alpha, else grandparent not in set
                    }
                }
                if (t->md != t->d - 1) {  // 2
                    return {y, CorneilStewartPerlCographRecognition::State::
                    EXISTS_1_NODE_NOT_PROPERLY_MARKED};
                }
                if (t->parent->marked == Marked::MARKED) {  // 1
                    return {y, CorneilStewartPerlCographRecognition::State::
                    CONTAINS_0_NODE};
                }
                t->unmark();
                t->md = 0;
                t = t->parent->parent;
            }
            w = u;
        }
        return {w, CorneilStewartPerlCographRecognition::State::COGRAPH};
    }

    std::vector<CoNode*> GetWereMarked(CoNode *u) {
        auto x = u->first_child;
        std::vector<CoNode*> a;
        while (x != nullptr && x->marked == Marked::MARKED_AND_UNMARKED) {
            a.push_back(x);
            x = x->next;
        }
        return a;
    }

    CoNode* GetLastFromChildren(CoNode *u) {
        auto x = u->first_child;
        while (x != nullptr && x->marked == Marked::MARKED_AND_UNMARKED) {
            x = x->next;
        }
        return x;
    }

    void CorneilStewartPerlCographRecognition::
    InsertXToCoTree(CoNode *u, CoNode *x) {
        std::vector<CoNode*> a;
        int u_number = u->number;
        a = GetWereMarked(u);
        if ((a.size() == 1 && u_number == 0) ||
        (u->d - a.size() == 1 && u_number == 1)) {
            CoNode *w = a[0];
            if (u_number == 1) {
                w = GetLastFromChildren(u);
            }
            if (w->type == Type::VERTEX) {
                auto *y = T.Add(Type::ZERO_ONE, u_number ^ 1);
                if (u_number == 0) {
                    u->RemoveWereMarked();
                } else {
                    u->RemoveWereNotMarked();
                }
                u->AddChild(y);
                y->AddChild(x);
                y->AddChild(w);
            } else {
                w->AddChild(x);
            }
        } else {
            auto vec = u->RemoveWereMarked();
            auto *y = T.Add(Type::ZERO_ONE, u_number);
            for (auto v : vec) {
                y->AddChild(v);
            }
            if (u_number == 1) {
                auto next = u->next;
                auto previous = u->previous;
                if (previous != nullptr) {
                    previous->next = y;
                }
                if (next != nullptr) {
                    next->previous = y;
                }
                y->previous = previous;
                y->next = next;
                if (previous == nullptr && u->parent != nullptr) {
                    u->parent->first_child = y;
                }
                if (u->parent != nullptr) {
                    y->parent = u->parent;
                } else {
                    T.root = y;
                }
                auto *z = T.Add(Type::ZERO_ONE, 0);
                y->AddChild(z);
                z->AddChild(x);
                z->AddChild(u);
            } else {
                auto *z = T.Add(Type::ZERO_ONE, 1);
                u->AddChild(z);
                z->AddChild(x);
                z->AddChild(y);
            }
        }
    }

    CorneilStewartPerlCographRecognition::State
    CorneilStewartPerlCographRecognition::Recognition() {
        auto *R = T.Add(Type::ZERO_ONE, 1);
        T.root = R;
        G = graph;
        std::vector<NetworKit::node> vertex;
        std::vector<CoNode*> covertex;
        std::map<NetworKit::node, int> pos;
        int count = 0;
        for (auto i : G.nodeRange()) {
            vertex.push_back(i);
            pos[i] = count++;
            auto *C = T.Add(Type::VERTEX, i);
            covertex.push_back(C);
        }
        for (auto i : G.nodeRange()) {
            std::vector<CoNode*> vec;
            for (auto u : G.neighborRange(i)) {
                vec.push_back(covertex[pos[u]]);
            }
            covertex[pos[i]]->out_edges = vec;
        }

        if (count == 0) {
            T.Clear();
            return State::COGRAPH;
        }
        if (count == 1) {
            R->AddChild(covertex[0]);
            T.Clear();
            return State::COGRAPH;
        }
        if (G.hasEdge(vertex[0], vertex[1])) {
            R->AddChild(covertex[0]);
            R->AddChild(covertex[1]);
        } else {
            auto *N = T.Add(Type::ZERO_ONE, 0);
            R->AddChild(N);
            N->AddChild(covertex[0]);
            N->AddChild(covertex[1]);
        }
        covertex[0]->in_graph = true;
        covertex[1]->in_graph = true;

        for (int i = 2; i < count; i++) {
            ResetAllCoNodes(T.root);
            Mark(covertex[i]);
            if (T.root->marked ==
                Marked::MARKED_AND_UNMARKED) {
                // all nodes of T were marked and unmarked <=>
                // R is marked and unmarked
                T.root->AddChild(covertex[i]);
            } else if (mark_ever_count == 0) {
                if (T.root->d == 1) {
                    T.root->first_child->AddChild(covertex[i]);
                } else {
                    auto *R1 = T.Add(Type::ZERO_ONE, 1);
                    auto *R2 = T.Add(Type::ZERO_ONE, 0);
                    R1->AddChild(R2);
                    R2->AddChild(T.root);
                    R2->AddChild(covertex[i]);
                    T.root = R1;
                }
            } else {
                auto [v, state] = FindLowest();
                if (state != State::COGRAPH) {
                    T.Clear();
                    return state;
                }
                InsertXToCoTree(v, covertex[i]);
            }
            covertex[i]->in_graph = true;
        }
        T.Clear();
        return State::COGRAPH;
    }

    bool CorneilStewartPerlCographRecognition::isCograph() const {
        assureFinished();
        return is_cograph == State::COGRAPH;
    }

    CorneilStewartPerlCographRecognition::State
    CorneilStewartPerlCographRecognition::getState() const {
        assureFinished();
        return is_cograph;
    }

    void CorneilStewartPerlCographRecognition::run() {
        hasRun = true;
        is_cograph = Recognition();
        return;
    }

}  // namespace Koala
