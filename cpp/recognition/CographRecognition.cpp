#include "recognition/CographRecognition.hpp"

namespace Koala {
    CographRecognition::CographRecognition(NetworKit::Graph &Graph) :
            cotree(Cotree(Graph)), graph(Graph), T(Twins(graph)) {
        original_graph = Graph;
        status = CographRecognition::State::UNKNOWN;
    }

    bool CographRecognition::isCograph() const {
        return (status == State::COGRAPH);
    }

    CographRecognition::State CographRecognition::getState() const {
        return status;
    }

    void CographRecognition::Clear() {
        num_of_parts = 0;
        order.clear();
        unused_parts.clear();
    }

    void CographRecognition::run() {
        Clear();
        hasRun = true;
        num_of_nodes = graph.numberOfNodes();

        NetworKit::count twins_counter = 0;


        permutation.E.resize(graph.numberOfNodes(), element());

        permutation.P.resize(graph.numberOfNodes() + 2, part());

        permutation.first_part = &(permutation.P[0]);
        permutation.last_part = &permutation.P[0];

        for (const auto &u : graph.nodeRange()) {
            permutation.E[u].num = u;
            permutation.addElementToPart(permutation.P[0], permutation.E[u]);
        }
        permutation.origin = &permutation.E[0];
        num_of_parts = 1;

        while (num_of_parts < num_of_nodes) {
            H = permutation.origin->my_part;
            if (H->first != H->last) {
                permutation.r = permutation.origin->my_part;
                permutation.l = permutation.origin->my_part;

                if (H->pivot == nullptr) {
                    H->pivot = H->first;
                }

                NetworKit::count n_part = permutation.newPart();


                for (auto v : graph.neighborRange(permutation.origin->num)) {
                    if (permutation.E[v].my_part == H) {
                        permutation.eraseElement(permutation.E[v]);
                        permutation.addElementToPart(permutation.P[n_part], permutation.E[v]);
                    }
                }

                if (permutation.P[n_part].first != nullptr) {
                    permutation.AddPart(*H, permutation.P[n_part]);
                    num_of_parts++;
                    unused_parts.push_back(&permutation.P[n_part]);
                }

                if (H->first != H->last) {
                    NetworKit::count o_part = permutation.newPart();

                    permutation.eraseElement(*permutation.origin);
                    permutation.addElementToPart(permutation.P[o_part], *permutation.origin);
                    permutation.P[o_part].pivot = permutation.origin;

                    permutation.AddPart(*H, permutation.P[o_part]);

                    num_of_parts++;
                    unused_parts.push_back(H);
                    permutation.r = permutation.origin->my_part;
                    permutation.l = permutation.origin->my_part;
                }


                permutation.lCheck();

                permutation.rCheck();
            }


            while (unused_parts.size() > 0) {
                if (unused_parts.front() != nullptr && unused_parts.front()->first != nullptr) {
                    H = unused_parts.front();
                    unused_parts.pop_front();
                    if (H->pivot == nullptr) {
                        H->pivot = H->first;
                    }

                    for (auto v : graph.neighborRange(H->pivot->num)) {
                        part *v_part = permutation.E[v].my_part;
                        if (v_part != H) {
                            v_part->amount++;
                        }
                    }

                    for (auto v : graph.neighborRange(H->pivot->num)) {
                        part *v_part = permutation.E[v].my_part;
                        if (v_part != H && v_part->size > v_part->amount) {
                            if (v_part->division == nullptr) {
                                if (v_part->first == v_part->last) {
                                    continue;
                                }
                                NetworKit::count n_part = permutation.newPart();
                                v_part->division = &permutation.P[n_part];
                                permutation.AddPart(*v_part, permutation.P[n_part]);
                                num_of_parts++;
                            }
                            if (v_part->pivot == &permutation.E[v]) {
                                unused_parts.push_back(v_part);
                                v_part->next->pivot = &permutation.E[v];
                            }
                            permutation.eraseElement(permutation.E[v]);
                            v_part->amount--;
                            permutation.addElementToPart(*v_part->next, permutation.E[v]);
                        }
                    }

                    for (auto v : graph.neighborRange(H->pivot->num)) {
                        part *v_part = permutation.E[v].my_part;
                        if (v_part->amount > 0) {
                            v_part->amount = 0;
                        }
                        if (v_part != H && v_part != permutation.first_part && v_part->previous->division != nullptr) {
                            part *previous_v_part = v_part->previous;
                            previous_v_part->division = nullptr;
                            previous_v_part->amount = 0;

                            if (previous_v_part->first == nullptr) {
                                permutation.erasePart(*previous_v_part);
                                num_of_parts--;
                            }

                            if (v_part->pivot == nullptr) {
                                unused_parts.push_back(permutation.E[v].my_part);
                            }
                        }
                    }
                } else {
                    unused_parts.pop_front();
                }
                permutation.lCheck();
                permutation.rCheck();
            }

            element *lpivot = permutation.l->pivot, *rpivot = permutation.r->pivot;

            if (permutation.r == permutation.origin->my_part) {
                permutation.origin = lpivot;
            } else {
                if (permutation.l == permutation.origin->my_part) {
                    permutation.origin = rpivot;
                } else {
                    if (graph.hasEdge(lpivot->num, rpivot->num)) {
                        permutation.origin = lpivot;
                    } else {
                        permutation.origin = rpivot;
                    }
                }
            }
        }


        part X_0 = part(), X_N = part(), *Z;
        element Y_0 = element(), Y_N = element();
        Y_0.num = NetworKit::none;
        Y_N.num = NetworKit::none;
        permutation.addElementToPart(X_0, Y_0);
        permutation.addElementToPart(X_N, Y_N);

        permutation.first_part->previous = &X_0;
        permutation.last_part->next = &X_N;
        X_0.next = permutation.first_part;
        X_N.previous = permutation.last_part;
        permutation.first_part = &X_0;
        permutation.last_part = &X_N;

        Z = permutation.first_part->next;
        NetworKit::count flag = 0;


        while (Z != permutation.last_part) {
            int twin = T.twins(Z->first->num, Z->previous->first->num, twins_counter);
            twins_counter += 4;
            if (twin < 2) {
                order.push_back({{Z->previous->first->num, Z->first->num}, twin});
                graph.removeNode(Z->previous->first->num);
                permutation.erasePart(*Z->previous);
                flag++;
            } else {
                twin = T.twins(Z->first->num, Z->next->first->num, twins_counter);
                twins_counter += 4;
                if (twin < 2) {
                    Z = Z->next;

                    order.push_back({{Z->previous->first->num, Z->first->num}, twin});


                    graph.removeNode(Z->previous->first->num);

                    permutation.erasePart(*(Z->previous));

                    flag++;
                } else {
                    Z = Z->next;
                }
            }
        }
        if (flag == num_of_nodes - 1) {
            status = CographRecognition::State::COGRAPH;
            order.push_back({{permutation.first_part->next->pivot->num, -1}, 3});
            cotree.setOrder(order);
            cotree.buildTree();
        } else {
            status = CographRecognition::State::NOT_COGRAPH;
        }
    }


} /* namespace Koala */

