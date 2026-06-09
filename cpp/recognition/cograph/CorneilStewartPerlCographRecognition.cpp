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
    is_cograph = recognition();
}

void CorneilStewartPerlCographRecognition::unmark(
        NetworKit::node u, std::queue<NetworKit::node> &ready,
        std::vector<NetworKit::node> &marked_nodes, std::vector<NetworKit::node> &touched) {
    marked[u] = Marked::MARKED_AND_UNMARKED;
    touched.push_back(u);
    mark_count--;
    md[u] = 0;
    if (u != T.getRoot()) {
        auto w = T.getNode(u).parent;
        md[w]++;
        touched.push_back(w);
        if (marked[w] != Marked::MARKED) {
            mark_count++;
            marked_nodes.push_back(w);
        }
        marked[w] = Marked::MARKED;
        if (md[w] == T.getNode(w).d) {
            ready.push(w);
        }
        T.moveChildToFront(w, u);
    }
}

void CorneilStewartPerlCographRecognition::mark(
        NetworKit::node x, const std::vector<NetworKit::node> &covertex,
        const std::vector<bool> &inserted, std::vector<NetworKit::node> &marked_nodes,
        std::vector<NetworKit::node> &touched) {
    mark_count = 0;
    mark_ever_count = 0;
    std::queue<NetworKit::node> ready;
    for (auto v : graph.neighborRange(T.getNode(x).number)) {
        if (!inserted[v]) {
            continue;
        }
        auto u = covertex[v];
        if (marked[u] != Marked::MARKED) {
            mark_count++;
            marked_nodes.push_back(u);
        }
        marked[u] = Marked::MARKED;
        touched.push_back(u);
        mark_ever_count++;
        if (md[u] == T.getNode(u).d) {
            ready.push(u);
        }
    }
    while (!ready.empty()) {
        NetworKit::node u = ready.front();
        ready.pop();
        if (marked[u] == Marked::MARKED && md[u] == T.getNode(u).d) {
            unmark(u, ready, marked_nodes, touched);
        }
    }
    if (mark_count && T.getNode(T.getRoot()).d == 1) {
        if (marked[T.getRoot()] != Marked::MARKED) {
            mark_count++;
            marked_nodes.push_back(T.getRoot());
            touched.push_back(T.getRoot());
        }
        marked[T.getRoot()] = Marked::MARKED;
    }
}

std::pair<NetworKit::node, CorneilStewartPerlCographRecognition::State>
CorneilStewartPerlCographRecognition::find_lowest(
        const std::vector<NetworKit::node> &marked_nodes) {
    NetworKit::node y = NetworKit::none;
    if (marked[T.getRoot()] == Marked::UNMARKED) {
        return {y, CorneilStewartPerlCographRecognition::State::GRANDPARENT_IS_NOT_IN_SET};
    }
    if (md[T.getRoot()] != T.getNode(T.getRoot()).d - 1) {
        y = T.getRoot();
    }
    marked[T.getRoot()] = Marked::MARKED_AND_UNMARKED;
    md[T.getRoot()] = 0;
    NetworKit::node w = T.getRoot();
    for (auto u : marked_nodes) {
        if (marked[u] != Marked::MARKED) {
            continue;
        }
        if (y != NetworKit::none) {  // 1 or 2
            if (T.getNode(y).number == 0) {
                return {y, CorneilStewartPerlCographRecognition::State::CONTAINS_0_NODE};
            } else {
                return {y, CorneilStewartPerlCographRecognition::State::
                EXISTS_1_NODE_NOT_PROPERLY_MARKED};
            }
        }
        NetworKit::node t;
        if (T.getNode(u).number == 1) {
            if (md[u] != T.getNode(u).d - 1) {
                y = u;
            }
            if (marked[T.getNode(u).parent] == Marked::MARKED) {  // 1 or 6
                if (y == NetworKit::none || T.getNode(y).number == 0) {
                    return {y, CorneilStewartPerlCographRecognition::State::CONTAINS_0_NODE};
                } else {
                    return {y, CorneilStewartPerlCographRecognition::State::WRONG_GRANDPARENT};
                }
            } else {
                t = T.getNode(T.getNode(u).parent).parent;
            }
        } else {
            y = u;
            t = T.getNode(u).parent;
        }
        marked[u] = Marked::MARKED_AND_UNMARKED;
        md[u] = 0;
        while (t != w) {
            if (t == T.getRoot()) {  // 4
                return {y, CorneilStewartPerlCographRecognition::State::NO_ONE_PATH};
            }
            if (marked[t] != Marked::MARKED) {  // 3 or 5 or 6
                if (y == NetworKit::none || T.getNode(y).number == 0) {
                    return {y, CorneilStewartPerlCographRecognition::State::WRONG_PARENT};
                } else {
                    return {y, CorneilStewartPerlCographRecognition::State::WRONG_GRANDPARENT};
                    // if y is alpha, else grandparent not in set
                }
            }
            if (md[t] != T.getNode(t).d - 1) {  // 2
                return {y, CorneilStewartPerlCographRecognition::State::
                EXISTS_1_NODE_NOT_PROPERLY_MARKED};
            }
            if (marked[T.getNode(t).parent] == Marked::MARKED) {  // 1
                return {y, CorneilStewartPerlCographRecognition::State::CONTAINS_0_NODE};
            }
            marked[t] = Marked::MARKED_AND_UNMARKED;
            md[t] = 0;
            t = T.getNode(T.getNode(t).parent).parent;
        }
        w = u;
    }
    return {w, CorneilStewartPerlCographRecognition::State::COGRAPH};
}

std::vector<NetworKit::node> get_marked(
        Cotree &T, NetworKit::node u, std::vector<Marked> &marked) {
    auto x = T.getNode(u).first_child;
    std::vector<NetworKit::node> a;
    while (x != NetworKit::none && marked[x] == Marked::MARKED_AND_UNMARKED) {
        a.push_back(x);
        x = T.getNode(x).next_sibling;
    }
    return a;
}

std::vector<NetworKit::node> CorneilStewartPerlCographRecognition::remove_marked(
        Cotree &T, NetworKit::node node, std::vector<Marked> &marked) {
    auto child = T.getNode(node).first_child;
    std::vector<NetworKit::node> removed;
    while (child != NetworKit::none && marked[child] == Marked::MARKED_AND_UNMARKED) {
        removed.push_back(child);
        auto next = T.getNode(child).next_sibling;
        T.removeChild(node, child);
        child = next;
    }
    return removed;
}

void remove_not_marked(Cotree &T, NetworKit::node node, std::vector<Marked> &marked) {
    auto child = T.getNode(node).first_child;
    while (child != NetworKit::none && marked[child] == Marked::MARKED_AND_UNMARKED) {
        child = T.getNode(child).next_sibling;
    }
    while (child != NetworKit::none) {
        auto next = T.getNode(child).next_sibling;
        T.removeChild(node, child);
        child = next;
    }
}

void CorneilStewartPerlCographRecognition::insert_to_cotree(
        NetworKit::node u, NetworKit::node x) {
    NetworKit::node u_number = T.getNode(u).number;
    std::vector<NetworKit::node> a = get_marked(T, u, marked);
    if ((a.size() == 1 && u_number == 0)
            || (T.getNode(u).d - static_cast<int>(a.size()) == 1 && u_number == 1)) {
        NetworKit::node w = a[0];
        if (u_number == 1) {
            w = get_marked(T, u, marked).back();
        }
        if (T.getNode(w).type == NodeType::LEAF) {
            auto y = T.add(
                u_number == 0 ? NodeType::COMPLEMENT_NODE : NodeType::UNION_NODE, u_number ^ 1);
            if (u_number == 0) {
                remove_marked(T, u, marked);
            } else {
                remove_not_marked(T, u, marked);
            }
            T.addChild(u, y);
            T.addChild(y, x);
            T.addChild(y, w);
        } else {
            T.addChild(w, x);
        }
    } else {
        auto y = T.add(u_number == 0 ? NodeType::UNION_NODE : NodeType::COMPLEMENT_NODE, u_number);
        for (auto v : remove_marked(T, u, marked)) {
            T.addChild(y, v);
        }
        if (u_number == 1) {
            auto parent = T.getNode(u).parent;
            if (parent != NetworKit::none) {
                T.replaceChild(parent, u, y);
            } else {
                T.setRoot(y);
            }
            auto z = T.add(NodeType::UNION_NODE, 0);
            T.addChild(y, z);
            T.addChild(z, x);
            T.addChild(z, u);
        } else {
            auto z = T.add(NodeType::COMPLEMENT_NODE, 1);
            T.addChild(u, z);
            T.addChild(z, x);
            T.addChild(z, y);
        }
    }
}

CorneilStewartPerlCographRecognition::State CorneilStewartPerlCographRecognition::recognition() {
    T.clear();
    const auto max_cotree_node_bound = 3 * graph.upperNodeIdBound();
    marked.assign(max_cotree_node_bound, Marked::UNMARKED);
    md.assign(max_cotree_node_bound, 0);
    T.reserve(max_cotree_node_bound);
    auto R = T.add(NodeType::COMPLEMENT_NODE, 1);
    T.setRoot(R);
    std::vector<NetworKit::node> vertex;
    std::vector<NetworKit::node> covertex(graph.upperNodeIdBound(), NetworKit::none);
    std::vector<bool> inserted(graph.upperNodeIdBound(), false);
    int count = 0;
    for (auto i : graph.nodeRange()) {
        vertex.push_back(i);
        covertex[i] = T.add(NodeType::LEAF, i);
        count++;
    }

    if (count == 0) {
        return State::COGRAPH;
    }
    if (count == 1) {
        T.addChild(R, covertex[vertex[0]]);
        return State::COGRAPH;
    }
    if (graph.hasEdge(vertex[0], vertex[1])) {
        T.addChild(R, covertex[vertex[0]]);
        T.addChild(R, covertex[vertex[1]]);
    } else {
        auto N = T.add(NodeType::UNION_NODE, 0);
        T.addChild(R, N);
        T.addChild(N, covertex[vertex[0]]);
        T.addChild(N, covertex[vertex[1]]);
    }
    inserted[vertex[0]] = true;
    inserted[vertex[1]] = true;

    for (int i = 2; i < count; i++) {
        std::vector<NetworKit::node> marked_nodes;
        std::vector<NetworKit::node> touched;
        mark(covertex[vertex[i]], covertex, inserted, marked_nodes, touched);
        if (marked[T.getRoot()] == Marked::MARKED_AND_UNMARKED) {
            // all nodes of T were marked and unmarked <=>
            // R is marked and unmarked
            T.addChild(T.getRoot(), covertex[vertex[i]]);
        } else if (mark_ever_count == 0) {
            if (T.getNode(T.getRoot()).d == 1) {
                T.addChild(T.getNode(T.getRoot()).first_child, covertex[vertex[i]]);
            } else {
                auto R1 = T.add(NodeType::COMPLEMENT_NODE, 1);
                auto R2 = T.add(NodeType::UNION_NODE, 0);
                T.addChild(R1, R2);
                T.addChild(R2, T.getRoot());
                T.addChild(R2, covertex[vertex[i]]);
                T.setRoot(R1);
            }
        } else {
            auto [v, state] = find_lowest(marked_nodes);
            if (state != State::COGRAPH) {
                return state;
            }
            insert_to_cotree(v, covertex[vertex[i]]);
        }
        for (auto u : touched) {
            marked[u] = Marked::UNMARKED;
            md[u] = 0;
        }
        inserted[vertex[i]] = true;
    }
    return State::COGRAPH;
}

}  // namespace Koala
