#include "recognition/CorneilStewartPerlCographRecognition.hpp"
#include <list>
#include "graph/GraphTools.hpp"
#include "CoTree.cpp"

namespace Koala {
    CoTree *T;
    NetworKit::Graph G;
    int mark_count = 0;
    int mark_and_unmarked_count = 0;
    int mark_ever_count = 0;

    void Unmark(std::queue<CoNode*> &marked_with_d_equal_to_md) {
        CoNode *u = marked_with_d_equal_to_md.front();
        marked_with_d_equal_to_md.pop();
        u->unmark();
        mark_count--;
        mark_and_unmarked_count++;
        u->md = 0;
        if (u != T->root) {
            auto w = u->parent;
            w->md++;
            if (w->marked == Marked::UNMARKED) {
                mark_count++;
            }
            w->mark();
            if (w->md == w->d) {
                marked_with_d_equal_to_md.push(w);
            }
            auto nxt = u->next;
            auto prv = u->prev;
            auto head = w->head_of_list_of_children;
            if (prv != nullptr) {
                prv->next = nxt;
                if (nxt != nullptr) {
                    nxt->prev = prv;
                }
                u->prev = nullptr;
                u->next = head;
                head->prev = u;
                w->head_of_list_of_children = u;
            }//else u is head
        }
    }

    void Mark(CoNode *x) {
        std::queue<CoNode*> marked_with_d_equal_to_md;
        mark_count = 0;
        mark_and_unmarked_count = 0;
        mark_ever_count = 0;
        for (auto u: x->out_edges) {//!!only neighbours which are already in graph
            if (!(u->in_graph)) {
                continue;
            }
            u->mark();
            mark_ever_count++;
            mark_count++;
            marked_with_d_equal_to_md.push(u);
        }
        while (!marked_with_d_equal_to_md.empty()) {
            Unmark(marked_with_d_equal_to_md);
        }
        if (mark_count) {
            if (T->root->d == 1) {
                T->root->mark();
            }
        }

    }

    void ResetAllCoNodes(CoNode *x, int level = 0) {
        x->UnmarkForNewIteration();
        CoNode *y = x->head_of_list_of_children;
        while (y != nullptr) {
            ResetAllCoNodes(y, level + 1);
            y = y->next;
        }
    }

    void MadeQueueOfMarked(CoNode *x, std::queue<CoNode*> &q) {
        if (x->marked == Marked::MARKED) {
            q.push(x);
        }
        auto y = x->head_of_list_of_children;
        while (y != nullptr) {
            MadeQueueOfMarked(y, q);
            y = y->next;
        }
    }

    std::pair<CoNode *, CorneilStewartPerlCographRecognition::State> FindLowest() {
        auto *y = new CoNode(Type::ZERO_ONE, 2);
        T->Add(y);
        CoNode *u, *w, *t;
        if (T->root->marked == Marked::UNMARKED) {
            return {y, CorneilStewartPerlCographRecognition::State::GRANDPARENT_IS_NOT_IN_SET};
        }
        if (T->root->md != T->root->d - 1) {
            y = T->root;
        }
        T->root->unmark();
        T->root->md = 0;
        w = T->root;
        std::queue<CoNode*> q;
        MadeQueueOfMarked(T->root, q);
        while (!q.empty()) {
            u = q.front();
            q.pop();
            if (u->marked != Marked::MARKED) {
                continue;
            }
            if (y->number != 2) {//1 or 2
                if (y->number == 0) {
                    return {y, CorneilStewartPerlCographRecognition::State::CONTAINS_0_NODE};
                } else {
                    return {y, CorneilStewartPerlCographRecognition::State::EXISTS_1_NODE_NOT_PROPERLY_MARKED};
                }
            }
            if (u->number == 1) {
                if (u->md != u->d - 1) {
                    y = u;
                }
                if (u->parent->marked == Marked::MARKED) {//1 or 6
                    if (y->number == 0) {
                        return {y, CorneilStewartPerlCographRecognition::State::CONTAINS_0_NODE};
                    } else {
                        return {y, CorneilStewartPerlCographRecognition::State::WRONG_GRANDPARENT};
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
                if (t == T->root) {//4
                    return {y, CorneilStewartPerlCographRecognition::State::NO_ONE_PATH};
                }
                if (t->marked != Marked::MARKED) {//3 or 5 or 6
                    if (y->number == 0) {
                        return {y, CorneilStewartPerlCographRecognition::State::WRONG_PARENT};
                    } else {
                        return {y,
                                CorneilStewartPerlCographRecognition::State::WRONG_GRANDPARENT};//if y is alpha, else grandparent not in set
                    }
                }
                if (t->md != t->d - 1) {//2
                    return {y, CorneilStewartPerlCographRecognition::State::EXISTS_1_NODE_NOT_PROPERLY_MARKED};
                }
                if (t->parent->marked == Marked::MARKED) {//1
                    return {y, CorneilStewartPerlCographRecognition::State::CONTAINS_0_NODE};
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
        auto x = u->head_of_list_of_children;
        std::vector<CoNode*> a;
        while (x != nullptr && x->marked == Marked::MARKED_AND_UNMARKED) {
            a.push_back(x);
            x = x->next;
        }
        return a;
    }


    CoNode *GetLastFromChildren(CoNode *u) {
        auto x = u->head_of_list_of_children;
        while (x != nullptr && x->marked == Marked::MARKED_AND_UNMARKED) {
            x = x->next;
        }
        return x;
    }

    void InsertXToCoTree(CoNode *u, CoNode *x) {
        std::vector<CoNode*> a;
        int u_number = u->number;
        a = GetWereMarked(u);
        if ((a.size() == 1 && u_number == 0) || (u->d - a.size() == 1 && u_number == 1)) {
            CoNode *w = a[0];
            if (u_number == 1) {
                w = GetLastFromChildren(u);
            }
            if (w->type == Type::VERTEX) {
                auto *y = new CoNode(Type::ZERO_ONE, u_number ^ 1);
                T->Add(y);
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
            auto *y = new CoNode(Type::ZERO_ONE, u_number);
            T->Add(y);
            for (auto v: vec) {
                y->AddChild(v);
            }
            if (u_number == 1) {
                auto nxt = u->next;
                auto prv = u->prev;
                if (prv != nullptr) {
                    prv->next = y;
                }
                if (nxt != nullptr) {
                    nxt->prev = y;
                }
                y->prev = prv;
                y->next = nxt;
                if (prv == nullptr && u->parent != nullptr) {
                    u->parent->head_of_list_of_children = y;
                }
                if (u->parent != nullptr) {
                    y->parent = u->parent;
                } else {
                    T->root = y;
                }
                auto *z = new CoNode(Type::ZERO_ONE, 0);
                T->Add(z);
                y->AddChild(z);
                z->AddChild(x);
                z->AddChild(u);
            } else {
                auto *z = new CoNode(Type::ZERO_ONE, 1);
                T->Add(z);
                u->AddChild(z);
                z->AddChild(x);
                z->AddChild(y);
            }
        }
    }


    CorneilStewartPerlCographRecognition::State CorneilStewartPerlCographRecognition::Cograph_Recognition() {
        auto *R = new CoNode(Type::ZERO_ONE, 1);
        CoTree Tp(R);
        T = &Tp;
        G = graph;
        std::vector<NetworKit::node> vertex;
        std::vector<CoNode*> covertex;
        std::map<NetworKit::node, int> pos;
        int cnt = 0;
        for (auto i: G.nodeRange()) {
            vertex.push_back(i);
            pos[i] = cnt++;
            auto *C = new CoNode(Type::VERTEX, int(i));
            T->Add(C);
            covertex.push_back(C);
        }
        for (auto i: G.nodeRange()) {
            std::vector<CoNode*> vec;
            for (auto u: G.neighborRange(i)) {
                vec.push_back(covertex[pos[u]]);
            }
            covertex[pos[i]]->out_edges = vec;
        }

        if (cnt == 0) {
            T->Clear();
            return State::COGRAPH;
        }
        if (cnt == 1) {
            R->AddChild(covertex[0]);
            T->Clear();
            return State::COGRAPH;
        }
        if (G.hasEdge(vertex[0], vertex[1])) {
            R->AddChild(covertex[0]);
            R->AddChild(covertex[1]);
        } else {
            auto *N = new CoNode(Type::ZERO_ONE, 0);
            T->Add(N);
            R->AddChild(N);
            N->AddChild(covertex[0]);
            N->AddChild(covertex[1]);
        }
        covertex[0]->in_graph = true;
        covertex[1]->in_graph = true;
        CoNode *root;

        for (int i = 2; i < cnt; i++) {
            root = T->root;
            ResetAllCoNodes(root);
            Mark(covertex[i]);
            if (root->marked ==
                Marked::MARKED_AND_UNMARKED) {//all nodes of T were marked and unmarked <=> R is marked and unmarked
                root->AddChild(covertex[i]);
            } else if (mark_ever_count == 0) {
                if (root->d == 1) {
                    (root->head_of_list_of_children)->AddChild(covertex[i]);
                } else {
                    auto *R1 = new CoNode(Type::ZERO_ONE, 1);
                    auto *R2 = new CoNode(Type::ZERO_ONE, 0);
                    T->Add(R1);
                    T->Add(R2);
                    R1->AddChild(R2);
                    R2->AddChild(root);
                    R2->AddChild(covertex[i]);
                    T->root = R1;
                }
            } else {
                auto [v, state] = FindLowest();
                if (state != State::COGRAPH) {
                    T->Clear();
                    return state;
                }
                InsertXToCoTree(v, covertex[i]);
            }
            covertex[i]->in_graph = true;
        }
        T->Clear();
        return State::COGRAPH;
    }

    CorneilStewartPerlCographRecognition::CorneilStewartPerlCographRecognition(
        NetworKit::Graph &graph) : graph(graph), is_cograph(State::UNKNOWN) {

    }

    bool CorneilStewartPerlCographRecognition::isCograph() const {
        assureFinished();
        return is_cograph == State::COGRAPH;
    }

    CorneilStewartPerlCographRecognition::State CorneilStewartPerlCographRecognition::getState() const {
        assureFinished();
        return is_cograph;
    }



    void CorneilStewartPerlCographRecognition::run() {
        hasRun = true;
        if (is_cograph != State::UNKNOWN) {
            return;
        }
        is_cograph = Cograph_Recognition();
        if (is_cograph != State::UNKNOWN) {
            return;
        }
        is_cograph = State::COGRAPH;
    }

}
