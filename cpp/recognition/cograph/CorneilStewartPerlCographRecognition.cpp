/*
 * CorneilStewartPerlCographRecognition.cpp
 *
 *  Created on: 06.03.2024
 *      Author: fixikmila
 */

#include <queue>
#include <utility>
#include <vector>

#include "graph/GraphTools.hpp"
#include "recognition/CographRecognition.hpp"
#include "structures/Cotree.hpp"

namespace Koala {

void CorneilStewartPerlCographRecognition::run() {
    hasRun = true;

    T.clear();
    if (graph.numberOfNodes() == 0) {
        is_cograph = State::COGRAPH;
        return;
    }

    const auto nodes = 3 * graph.upperNodeIdBound() + 1;
    status.assign(nodes, Marked::UNMARKED);
    md.assign(nodes, 0);
    T.reserve(nodes);
    auto R = T.add(NodeType::COMPLEMENT_NODE);
    T.setRoot(R);
    std::vector<NetworKit::node> vertices(graph.nodeRange().begin(), graph.nodeRange().end());
    std::vector<NetworKit::node> covertex(graph.upperNodeIdBound(), NetworKit::none);

    if (vertices.size() == 1) {
        T.addChild(R, T.add(NodeType::LEAF, vertices[0]));
        is_cograph = State::COGRAPH;
        return;
    }
    auto cv0 = covertex[vertices[0]] = T.add(NodeType::LEAF, vertices[0]);
    auto cv1 = covertex[vertices[1]] = T.add(NodeType::LEAF, vertices[1]);
    if (graph.hasEdge(vertices[0], vertices[1])) {
        T.addChild(R, cv0);
        T.addChild(R, cv1);
    } else {
        auto N = T.add(NodeType::UNION_NODE);
        T.addChild(R, N);
        T.addChild(N, cv0);
        T.addChild(N, cv1);
    }

    for (NetworKit::index i = 2; i < vertices.size(); ++i) {
        auto u = vertices[i], cu = covertex[u] = T.add(NodeType::LEAF, u);
        std::vector<NetworKit::node> processed;
        processed.reserve(graph.degree(u));
        for (auto v : graph.neighborRange(u)) {
            if (v < u) {  // already inserted vertices
                processed.push_back(covertex[v]);
            }
        }

        auto marked = mark(processed);
        if (status[T.getRoot()] == Marked::MARKED_AND_UNMARKED) {
            // all nodes of T were marked and unmarked <=> R is marked and unmarked
            T.addChild(T.getRoot(), cu);
        } else if (processed.empty()) {
            if (T.getNode(T.getRoot()).size == 1) {
                T.addChild(T.getNode(T.getRoot()).first_child, cu);
            } else {
                auto R1 = T.add(NodeType::COMPLEMENT_NODE);
                auto R2 = T.add(NodeType::UNION_NODE);
                T.addChild(R1, R2);
                T.addChild(R2, T.getRoot());
                T.addChild(R2, cu);
                T.setRoot(R1);
            }
        } else {
            auto low = find_lowest(marked);
            if (low == NetworKit::none) {
                is_cograph = State::NOT_COGRAPH;
                touched.clear();
                return;
            }
            attach_to_cotree(low, cu);
        }
        for (auto v : touched) {
            status[v] = Marked::UNMARKED, md[v] = 0;
        }
        touched.clear();
    }
    is_cograph = State::COGRAPH;
}

std::vector<NetworKit::node> CorneilStewartPerlCographRecognition::mark(
        const std::vector<NetworKit::node> &processed) {
    std::vector<NetworKit::node> marked;
    NetworKit::count mark_count = 0;
    std::queue<NetworKit::node> ready;
    auto mark_node = [&](NetworKit::node u) {
        if (status[u] != Marked::MARKED) {
            marked.push_back(u), status[u] = Marked::MARKED, mark_count++;
        }
        touched.push_back(u);
    };
    for (auto u : processed) {
        mark_node(u);
        if (md[u] == T.getNode(u).size) {
            ready.push(u);
        }
    }
    while (!ready.empty()) {
        NetworKit::node u = ready.front();
        ready.pop();
        if (status[u] != Marked::MARKED || md[u] != T.getNode(u).size) {
            continue;
        }
        status[u] = Marked::MARKED_AND_UNMARKED, md[u] = 0, mark_count--;
        touched.push_back(u);
        if (u != T.getRoot()) {
            auto w = T.getNode(u).parent;
            md[w]++;
            mark_node(w);
            if (md[w] == T.getNode(w).size) {
                ready.push(w);
            }
            T.moveChildToFront(w, u);
        }
    }
    if (mark_count > 0 && T.getNode(T.getRoot()).size == 1
            && status[T.getRoot()] != Marked::MARKED) {
        mark_node(T.getRoot());
    }
    return marked;
}

NetworKit::node CorneilStewartPerlCographRecognition::find_lowest(
        const std::vector<NetworKit::node> &marked) {
    NetworKit::node y = NetworKit::none;
    if (status[T.getRoot()] == Marked::UNMARKED) {
        return NetworKit::none;  // GRANDPARENT_IS_NOT_IN_SET
    }
    auto finish = [&](NetworKit::node u) {
        status[u] = Marked::MARKED_AND_UNMARKED, md[u] = 0;
        if (u != T.getRoot()) {
            T.moveChildToFront(T.getNode(u).parent, u);
        }
    };
    if (md[T.getRoot()] + 1 != T.getNode(T.getRoot()).size) {
        y = T.getRoot();
    }
    finish(T.getRoot());
    NetworKit::node w = T.getRoot();
    for (auto u : marked) {
        if (status[u] != Marked::MARKED) {
            continue;
        }
        if (y != NetworKit::none) {
            if (T.getNode(y).type == NodeType::UNION_NODE) {
                return NetworKit::none;  // case 1: CONTAINS_0_NODE
            } else {
                return NetworKit::none;  // case 2: EXISTS_1_NODE_NOT_PROPERLY_MARKED
            }
        }
        NetworKit::node t;
        if (T.getNode(u).type == NodeType::COMPLEMENT_NODE) {
            if (md[u] + 1 != T.getNode(u).size) {
                y = u;
            }
            if (status[T.getNode(u).parent] == Marked::MARKED) {
                if (y == NetworKit::none || T.getNode(y).type == NodeType::UNION_NODE) {
                    return NetworKit::none;  // case 1: CONTAINS_0_NODE
                } else {
                    return NetworKit::none;  // case 6: WRONG_GRANDPARENT
                }
            }
            t = T.getNode(T.getNode(u).parent).parent;
        } else {
            y = u;
            t = T.getNode(u).parent;
        }
        finish(u);
        while (t != w) {
            if (t == T.getRoot()) {
                return NetworKit::none;  // case 4: NO_ONE_PATH
            }
            if (status[t] != Marked::MARKED) {
                if (y == NetworKit::none || T.getNode(y).type == NodeType::UNION_NODE) {
                    return NetworKit::none;  // cases 3 and 5: WRONG_PARENT
                } else {
                    return NetworKit::none;  // case 6: WRONG_GRANDPARENT
                }
            }
            if (md[t] != T.getNode(t).size - 1) {
                return NetworKit::none;  // case 2: EXISTS_1_NODE_NOT_PROPERLY_MARKED
            }
            if (status[T.getNode(t).parent] == Marked::MARKED) {
                return NetworKit::none;  // case 1: CONTAINS_0_NODE
            }
            finish(t);
            t = T.getNode(T.getNode(t).parent).parent;
        }
        w = u;
    }
    return w;
}

void CorneilStewartPerlCographRecognition::attach_to_cotree(NetworKit::node u, NetworKit::node x) {
    std::vector<NetworKit::node> A;
    auto child = T.getNode(u).first_child;
    while (child != NetworKit::none && status[child] == Marked::MARKED_AND_UNMARKED) {
        A.push_back(child);
        child = T.getNode(child).next_sibling;
    }

    const bool u_is_zero = T.getNode(u).type == NodeType::UNION_NODE;
    if ((A.size() == 1 && u_is_zero) || (T.getNode(u).size == A.size() + 1 && !u_is_zero)) {
        NetworKit::node w = u_is_zero ? A[0] : child;
        if (T.getNode(w).type == NodeType::LEAF) {
            auto y = T.add(u_is_zero ? NodeType::COMPLEMENT_NODE : NodeType::UNION_NODE);
            T.removeChild(u, w);
            T.addChild(u, y);
            T.addChild(y, x);
            T.addChild(y, w);
        } else {
            T.addChild(w, x);
        }
    } else {
        auto y = T.add(T.getNode(u).type);
        for (auto v : A) {
            T.removeChild(u, v);
            T.addChild(y, v);
        }
        if (!u_is_zero) {
            auto parent = T.getNode(u).parent;
            if (parent != NetworKit::none) {
                T.replaceChild(parent, u, y);
            } else {
                T.setRoot(y);
            }
            auto z = T.add(NodeType::UNION_NODE);
            T.addChild(y, z);
            T.addChild(z, x);
            T.addChild(z, u);
        } else {
            auto z = T.add(NodeType::COMPLEMENT_NODE);
            T.addChild(u, z);
            T.addChild(z, x);
            T.addChild(z, y);
        }
    }
}

}  // namespace Koala
