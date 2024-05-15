#include "recognition/CographAlg.hpp"
#include "cograph/Twins.cpp"

namespace Koala {
    CographRecognition::CographRecognition(NetworKit::Graph &Graph) {
        original_graph = &Graph;
        graph = new NetworKit::Graph(Graph);
        status = CographRecognition::State::UNKNOWN;
        cotree = new Cotree(Graph);
        permutation = new FactorizingPermutation();
    }

    CographRecognition::~CographRecognition() {
        for (auto v: unused_parts) {
            delete v;
        }
        delete cotree;
        delete graph;
    }

    bool CographRecognition::IsCograph() const {
        return (status == State::COGRAPH);
    }

    CographRecognition::State CographRecognition::GetState() const {
        return status;
    }

    void CographRecognition::Clear() {
        num_of_parts = 0;
        order.clear();
        for (auto v: unused_parts) {
            delete v;
        }
        unused_parts.clear();
    }

    void CographRecognition::run() {
        Clear();
        num_of_nodes = graph->numberOfNodes();
        long long i, twins_counter = 0;
        part *P = new part();

        permutation->first_part = P;
        permutation->last_part = P;

        for (i = 0; i < num_of_nodes; i++) {
            element *E1 = new element();
            permutation->E.push_back(E1);
            permutation->E[i]->num = i;
            permutation->AddElementToPart(P, permutation->E[i]);

        }
        permutation->origin = permutation->E[0];
        num_of_parts = 1;
        while (num_of_parts < num_of_nodes) {
            part *H = permutation->origin->my_part;
            if (H->first != H->last) {
                permutation->r = permutation->origin->my_part;
                permutation->l = permutation->origin->my_part;

                if (H->pivot == nullptr) {
                    H->pivot = H->first;
                }
                part *n_part = new part();

                for (auto v: graph->neighborRange(permutation->origin->num)) {
                    if (permutation->E[v]->my_part == H) {
                        permutation->EraseElement(permutation->E[v]);
                        permutation->AddElementToPart(n_part, permutation->E[v]);
                    }
                }

                if (n_part->first != nullptr) {
                    permutation->AddPart(H, n_part);
                    num_of_parts++;
                    unused_parts.push_back(n_part);
                }

                if (H->first != H->last) {
                    part *o_part = new part();
                    permutation->EraseElement(permutation->origin);
                    permutation->AddElementToPart(o_part, permutation->origin);
                    o_part->pivot = permutation->origin;
                    permutation->AddPart(H, o_part);
                    num_of_parts++;
                    unused_parts.push_back(H);
                    permutation->r = permutation->origin->my_part;
                    permutation->l = permutation->origin->my_part;
                }

                permutation->Lcheck();
                permutation->Rcheck();
            }

            while (unused_parts.size() > 0) {
                H = unused_parts.front();
                unused_parts.pop_front();
                if (H != nullptr && H->first != nullptr) {
                    if (H->pivot == nullptr) {
                        H->pivot = H->first;
                    }

                    for (auto v: graph->neighborRange(H->pivot->num)) {
                        part *v_part = permutation->E[v]->my_part;
                        if (v_part != H) {
                            v_part->amount++;
                        }
                    }

                    for (auto v: graph->neighborRange(H->pivot->num)) {
                        part *v_part = permutation->E[v]->my_part;
                        if (v_part != H && v_part->size > v_part->amount) {
                            if (v_part->division == nullptr) {
                                if (v_part->first == v_part->last) {
                                    continue;
                                }
                                part *n_part = new part();
                                v_part->division = n_part;
                                permutation->AddPart(v_part, n_part);
                                num_of_parts++;
                            }
                            if (v_part->pivot == permutation->E[v]) {
                                unused_parts.push_back(v_part);
                                v_part->next->pivot = permutation->E[v];
                            }
                            permutation->EraseElement(permutation->E[v]);
                            v_part->amount--;
                            permutation->AddElementToPart(v_part->next, permutation->E[v]);
                        }
                    }

                    for (auto v: graph->neighborRange(H->pivot->num)) {
                        part *v_part = permutation->E[v]->my_part;
                        if (v_part->amount > 0) {
                            v_part->amount = 0;
                        }
                        if (v_part != H && v_part != permutation->first_part && v_part->previous->division != nullptr) {
                            part *previous_v_part = v_part->previous;
                            previous_v_part->division = nullptr;
                            previous_v_part->amount = 0;

                            if (previous_v_part->first == nullptr) {
                                permutation->ErasePart(previous_v_part);
                                num_of_parts--;
                            }

                            if (v_part->pivot == nullptr) {
                                unused_parts.push_back(permutation->E[v]->my_part);
                            }
                        }
                    }
                }
                permutation->Lcheck();
                permutation->Rcheck();
            }

            element *lpivot = permutation->l->pivot, *rpivot = permutation->r->pivot;

            if (permutation->r == permutation->origin->my_part) {
                permutation->origin = lpivot;
            } else {
                if (permutation->l == permutation->origin->my_part) {
                    permutation->origin = rpivot;
                } else {
                    if (graph->hasEdge(lpivot->num, rpivot->num)) {
                        permutation->origin = lpivot;
                    } else {
                        permutation->origin = rpivot;
                    }
                }
            }
        }
        part *X_0 = new part(), *X_N = new part(), *Z;
        element *Y_0 = new element(), *Y_N = new element();
        Y_0->num = -1;
        Y_N->num = -1;
        permutation->AddElementToPart(X_0, Y_0);
        permutation->AddElementToPart(X_N, Y_N);

        permutation->first_part->previous = X_0;
        permutation->last_part->next = X_N;
        X_0->next = permutation->first_part;
        X_N->previous = permutation->last_part;
        permutation->first_part = X_0;
        permutation->last_part = X_N;

        Z = permutation->first_part->next;
        long long flag = 0;
        StartTest(num_of_nodes);
        while (Z != permutation->last_part) {
            int twin = Twins(Z->first->num, Z->previous->first->num, graph, twins_counter);
            twins_counter += 4;
            if (twin < 2) {
                order.push_back({{Z->previous->first->num, Z->first->num}, twin});
                graph->removeNode(Z->previous->first->num);
                permutation->ErasePart(Z->previous);
                flag++;
            } else {
                twin = Twins(Z->first->num, Z->previous->first->num, graph, twins_counter);
                twins_counter += 4;
                if (twin < 2) {
                    Z = Z->next;
                    order.push_back({{Z->previous->first->num, Z->first->num}, twin});
                    graph->removeNode(Z->previous->first->num);
                    permutation->ErasePart(Z->previous);
                    flag++;
                } else {
                    Z = Z->next;
                }
            }
        }
        if (flag == num_of_nodes - 1) {
            status = CographRecognition::State::COGRAPH;
            order.push_back({{permutation->first_part->next->pivot->num, -1}, 3});
        } else {

            status = CographRecognition::State::NOT_COGRAPH;
        }
        cotree->SetOrder(order);
    }


}
