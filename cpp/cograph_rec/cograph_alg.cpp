#include <list>
#include <set>

#include <unordered_map>
#include <map>

#include <graph/GraphTools.hpp>
#include <cograph_rec/cograph_alg.hpp>

#include "cograph/lists.cpp"
#include "cograph/Twins.cpp"
#include "cograph/MaxClique.cpp"
#include "cograph/pathwidth.cpp"

namespace Koala {

    CographRecognition::CographRecognition(NetworKit::Graph &Graph) {
        original_graph = &Graph;
        graph = new NetworKit::Graph(Graph);
        status = CographRecognition::State::UNKNOWN;
        prepared = 0;
    }

    bool CographRecognition::IsCograph() const {
        return (status == State::COGRAPH);
    }

    CographRecognition::State CographRecognition::GetState() const {
        return status;
    }

    void CographRecognition::Clear() {
        prepared = 0;
        num_of_parts = 0;
        l = nullptr;
        r = nullptr;
        first_part = nullptr;
        last_part = nullptr;
        origin = nullptr;
        order.clear();
        left_son.clear();
        right_son.clear();
        parent.clear();
        type.clear();
        size.clear();
        for (auto v: E) {
            delete v;
        }
        for (auto v: unused_parts) {
            delete v;
        }
        for (auto v: garbage) {
            delete v;
        }
        E.clear();
        unused_parts.clear();
        garbage.clear();
        status = CographRecognition::State::UNKNOWN;
    }


    void CographRecognition::run() {

        Clear();
        num_of_nodes = graph->numberOfNodes();
        long long i, twins_counter = 0;

        part *P = new part();

        first_part = P;
        last_part = P;

        for (i = 0; i < num_of_nodes; i++) {
            element *E1 = new element();
            E.push_back(E1);
            E[i]->num = i;
            AddElementToPart(P, E[i]);
            used.push_back(0);
        }

        origin = E[0];
        num_of_parts = 1;

        while (num_of_parts < num_of_nodes) {
            part *H = origin->my_part;
            if (H->first != H->last) {
                r = origin->my_part;
                l = origin->my_part;

                if (H->pivot == nullptr) {
                    H->pivot = H->first;
                }
                part *n_part = new part();


                for (auto v: graph->neighborRange(origin->num)) {
                    if (E[v]->my_part == H) {
                        EraseElement(E[v]);
                        AddElementToPart(n_part, E[v]);
                    }
                }

                if (n_part->first != nullptr) {
                    AddPart(H, n_part);
                    num_of_parts++;
                    unused_parts.push_back(n_part);
                }

                if (H->first != H->last) {
                    part *o_part = new part();
                    EraseElement(origin);
                    AddElementToPart(o_part, origin);
                    o_part->pivot = origin;
                    AddPart(H, o_part);
                    num_of_parts++;
                    unused_parts.push_back(H);

                    r = origin->my_part;
                    l = origin->my_part;
                }

                Lcheck();
                Rcheck();


            }

            while (unused_parts.size() > 0) {

                H = unused_parts.front();

                unused_parts.pop_front();
                if (H != nullptr && H->first != nullptr) {

                    if (H->pivot == nullptr) {
                        H->pivot = H->first;
                    }

                    for (auto v: graph->neighborRange(H->pivot->num)) {
                        part *v_part = E[v]->my_part;
                        if (v_part != H) {
                            v_part->amount++;
                        }
                    }

                    for (auto v: graph->neighborRange(H->pivot->num)) {
                        part *v_part = E[v]->my_part;
                        if (v_part != H && v_part->size > v_part->amount) {
                            if (v_part->division == nullptr) {
                                if (v_part->first == v_part->last) {
                                    continue;
                                }
                                part *n_part = new part();
                                v_part->division = n_part;
                                AddPart(v_part, n_part);
                                num_of_parts++;
                            }
                            if (v_part->pivot == E[v]) {
                                unused_parts.push_back(v_part);
                                v_part->next->pivot = E[v];
                            }
                            EraseElement(E[v]);
                            v_part->amount--;
                            AddElementToPart(v_part->next, E[v]);
                        }
                    }

                    for (auto v: graph->neighborRange(H->pivot->num)) {
                        part *v_part = E[v]->my_part;
                        if (v_part->amount > 0) {
                            v_part->amount = 0;
                        }
                        if (v_part != H && v_part != first_part && v_part->previous->division != nullptr) {
                            part *previous_v_part = v_part->previous;
                            previous_v_part->division = nullptr;
                            previous_v_part->amount = 0;

                            if (previous_v_part->first == nullptr) {
                                ErasePart(previous_v_part);
                                num_of_parts--;
                            }

                            if (v_part->pivot == nullptr) {
                                unused_parts.push_back(E[v]->my_part);
                            }

                        }
                    }
                }

                Lcheck();
                Rcheck();

            }


            element *lpivot = l->pivot, *rpivot = r->pivot;

            if (r == origin->my_part) {
                origin = lpivot;
            } else {
                if (l == origin->my_part) {
                    origin = rpivot;
                } else {
                    if (graph->hasEdge(lpivot->num, rpivot->num)) {
                        origin = lpivot;
                    } else {
                        origin = rpivot;
                    }
                }
            }

        }


        part *X_0 = new part(), *X_N = new part(), *Z;
        element *Y_0 = new element(), *Y_N = new element();
        Y_0->num = -1;
        Y_N->num = -1;
        AddElementToPart(X_0, Y_0);
        AddElementToPart(X_N, Y_N);


        first_part->previous = X_0;
        last_part->next = X_N;
        X_0->next = first_part;
        X_N->previous = last_part;
        first_part = X_0;
        last_part = X_N;


        Z = first_part->next;
        long long flag = 0;
        StartTest(num_of_nodes);
        while (Z != last_part) {
            int twin = Twins(Z->first->num, Z->previous->first->num, graph, twins_counter);
            if (twin < 2) {
                order.push_back({{Z->previous->first->num, Z->first->num}, twin});
                graph->removeNode(Z->previous->first->num);
                ErasePart(Z->previous);
                flag++;
            } else {
                twin = Twins(Z->first->num, Z->previous->first->num, graph, twins_counter);
                if (twin < 2) {
                    Z = Z->next;
                    order.push_back({{Z->previous->first->num, Z->first->num}, twin});
                    graph->removeNode(Z->previous->first->num);
                    ErasePart(Z->previous);
                    flag++;
                } else {
                    Z = Z->next;
                }
            }
        }

        delete graph;
        prepared = 1;
        if (flag == num_of_nodes - 1) {
            status = CographRecognition::State::COGRAPH;
            order.push_back({{first_part->next->pivot->num, -1}, 3});
        } else {

            status = CographRecognition::State::IS_NOT_COGRAPH;
        }
    }

    void CographRecognition::BuildTree() {
        reverse(order.begin(), order.end());
        long long n = num_of_nodes;
        long long i;

        for (i = 0; i < 2 * n; i++) {
            left_son.push_back(0);
            right_son.push_back(0);
            parent.push_back(0);
            type.push_back(0);
            size.push_back(0);
        }

        left_son[n] = order[0].first.first;
        right_son[n] = -1;
        parent[n] = -1;
        type[n] = 0;

        left_son[order[0].first.first] = -1;
        right_son[order[0].first.first] = -1;
        type[order[0].first.first] = 2;
        parent[order[0].first.first] = n;

        for (i = 1; i < n; i++) {
            left_son[order[i].first.first] = -1;
            right_son[order[i].first.first] = -1;
            type[order[i].first.first] = 2;
            parent[order[i].first.first] = n + i; //&

            long long p = parent[order[i].first.second];
            if (left_son[p] == order[i].first.second) {
                left_son[p] = n + i;
            } else {
                right_son[p] = n + i;
            }

            parent[order[i].first.second] = n + i;

            left_son[n + i] = order[i].first.first;
            right_son[n + i] = order[i].first.first;
            parent[n + i] = p;
            type[n + i] = order[i].second;
        }
        prepared = 2;
    }

}  // namespace Koala
