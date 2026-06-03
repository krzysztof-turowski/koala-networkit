/*
 * CorneilStewartPerlCographRecognition.cpp
 *
 *  Created on: 06.03.2024
 *      Author: fixikmila
 */

#include <map>
#include <queue>
#include <stack>
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

void CorneilStewartPerlCographRecognition::unmark() {
    NetworKit::count u = marked_with_d_equal_to_md.front();
    marked_with_d_equal_to_md.pop();
    T.unmark(u);
    mark_count--;
    mark_and_unmarked_count++;
    T.getNode(u).md = 0;
    if (u != T.getRoot()) {
        auto w = T.getNode(u).parent;
        T.getNode(w).md++;
        if (T.getNode(w).marked == Marked::UNMARKED) {
            mark_count++;
        }
        T.mark(w);
        if (T.getNode(w).md == T.getNode(w).d) {
            marked_with_d_equal_to_md.push(w);
        }
        T.moveChildToFront(w, u);
    }
}

void CorneilStewartPerlCographRecognition::mark(NetworKit::count x) {
    mark_count = 0;
    mark_and_unmarked_count = 0;
    mark_ever_count = 0;
    for (auto u : T.getNode(x).out_edges) {
        // !!only neighbours which are already in graph
        if (!(T.getNode(u).in_graph)) {
            continue;
        }
        T.mark(u);
        mark_ever_count++;
        mark_count++;
        marked_with_d_equal_to_md.push(u);
    }
    while (!marked_with_d_equal_to_md.empty()) {
        unmark();
    }
    if (mark_count && T.getNode(T.getRoot()).d == 1) {
        T.mark(T.getRoot());
    }
}

void reset_all_conodes(Cotree &T, NetworKit::count x) {
    T.unmarkForNewIteration(x);
    auto y = T.getNode(x).first_child;
    while (y != NetworKit::none) {
        reset_all_conodes(T, y);
        y = T.getNode(y).next_sibling;
    }
}

std::pair<NetworKit::count, CorneilStewartPerlCographRecognition::State>
CorneilStewartPerlCographRecognition::find_lowest() {
    NetworKit::count y = NetworKit::none;
    if (T.getNode(T.getRoot()).marked == Marked::UNMARKED) {
        return {y, CorneilStewartPerlCographRecognition::State::GRANDPARENT_IS_NOT_IN_SET};
    }
    if (T.getNode(T.getRoot()).md != T.getNode(T.getRoot()).d - 1) {
        y = T.getRoot();
    }
    T.unmark(T.getRoot());
    T.getNode(T.getRoot()).md = 0;
    NetworKit::count w = T.getRoot();
    std::queue<NetworKit::count> q;
    std::stack<NetworKit::count> s;
    s.push(T.getRoot());
    while (!s.empty()) {
        auto x = s.top();
        s.pop();
        if (T.getNode(x).marked == Marked::MARKED) {
            q.push(x);
        }
        auto z = T.getNode(x).first_child;
        while (z != NetworKit::none) {
            s.push(z);
            z = T.getNode(z).next_sibling;
        }
    }
    while (!q.empty()) {
        NetworKit::count u = q.front();
        q.pop();
        if (T.getNode(u).marked != Marked::MARKED) {
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
        NetworKit::count t;
        if (T.getNode(u).number == 1) {
            if (T.getNode(u).md != T.getNode(u).d - 1) {
                y = u;
            }
            if (T.getNode(T.getNode(u).parent).marked == Marked::MARKED) {  // 1 or 6
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
        T.unmark(u);
        T.getNode(u).md = 0;
        while (t != w) {
            if (t == T.getRoot()) {  // 4
                return {y, CorneilStewartPerlCographRecognition::State::NO_ONE_PATH};
            }
            if (T.getNode(t).marked != Marked::MARKED) {  // 3 or 5 or 6
                if (y == NetworKit::none || T.getNode(y).number == 0) {
                    return {y, CorneilStewartPerlCographRecognition::State::WRONG_PARENT};
                } else {
                    return {y, CorneilStewartPerlCographRecognition::State::WRONG_GRANDPARENT};
                    // if y is alpha, else grandparent not in set
                }
            }
            if (T.getNode(t).md != T.getNode(t).d - 1) {  // 2
                return {y, CorneilStewartPerlCographRecognition::State::
                EXISTS_1_NODE_NOT_PROPERLY_MARKED};
            }
            if (T.getNode(T.getNode(t).parent).marked == Marked::MARKED) {  // 1
                return {y, CorneilStewartPerlCographRecognition::State::CONTAINS_0_NODE};
            }
            T.unmark(t);
            T.getNode(t).md = 0;
            t = T.getNode(T.getNode(t).parent).parent;
        }
        w = u;
    }
    return {w, CorneilStewartPerlCographRecognition::State::COGRAPH};
}

std::vector<NetworKit::count> get_were_marked(Cotree &T, NetworKit::count u) {
    auto x = T.getNode(u).first_child;
    std::vector<NetworKit::count> a;
    while (x != NetworKit::none && T.getNode(x).marked == Marked::MARKED_AND_UNMARKED) {
        a.push_back(x);
        x = T.getNode(x).next_sibling;
    }
    return a;
}

NetworKit::count get_last_from_children(Cotree &T, NetworKit::count u) {
    auto x = T.getNode(u).first_child;
    while (x != NetworKit::none && T.getNode(x).marked == Marked::MARKED_AND_UNMARKED) {
        x = T.getNode(x).next_sibling;
    }
    return x;
}

void CorneilStewartPerlCographRecognition::insert_x_to_cotree(
        NetworKit::count u, NetworKit::count x) {
    std::vector<NetworKit::count> a;
    int u_number = T.getNode(u).number;
    a = get_were_marked(T, u);
    if ((a.size() == 1 && u_number == 0)
            || (T.getNode(u).d - static_cast<int>(a.size()) == 1 && u_number == 1)) {
        NetworKit::count w = a[0];
        if (u_number == 1) {
            w = get_last_from_children(T, u);
        }
        if (T.getNode(w).type == NodeType::LEAF) {
            auto y = T.add(
                u_number == 0 ? NodeType::COMPLEMENT_NODE : NodeType::UNION_NODE, u_number ^ 1);
            if (u_number == 0) {
                T.removeWereMarked(u);
            } else {
                T.removeWereNotMarked(u);
            }
            T.addChild(u, y);
            T.addChild(y, x);
            T.addChild(y, w);
        } else {
            T.addChild(w, x);
        }
    } else {
        auto vec = T.removeWereMarked(u);
        auto y = T.add(u_number == 0 ? NodeType::UNION_NODE : NodeType::COMPLEMENT_NODE, u_number);
        for (auto v : vec) {
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
    T.reserve(3 * graph.numberOfNodes());
    auto R = T.add(NodeType::COMPLEMENT_NODE, 1);
    T.setRoot(R);
    std::vector<NetworKit::node> vertex;
    std::vector<NetworKit::count> covertex;
    std::map<NetworKit::node, int> pos;
    int count = 0;
    for (auto i : graph.nodeRange()) {
        vertex.push_back(i);
        pos[i] = count++;
        auto C = T.add(NodeType::LEAF, static_cast<int>(i));
        covertex.push_back(C);
    }
    for (auto i : graph.nodeRange()) {
        std::vector<NetworKit::count> vec;
        for (auto u : graph.neighborRange(i)) {
            vec.push_back(covertex[pos[u]]);
        }
        T.getNode(covertex[pos[i]]).out_edges = vec;
    }

    if (count == 0) {
        T.clear();
        return State::COGRAPH;
    }
    if (count == 1) {
        T.addChild(R, covertex[0]);
        T.clear();
        return State::COGRAPH;
    }
    if (graph.hasEdge(vertex[0], vertex[1])) {
        T.addChild(R, covertex[0]);
        T.addChild(R, covertex[1]);
    } else {
        auto N = T.add(NodeType::UNION_NODE, 0);
        T.addChild(R, N);
        T.addChild(N, covertex[0]);
        T.addChild(N, covertex[1]);
    }
    T.getNode(covertex[0]).in_graph = true;
    T.getNode(covertex[1]).in_graph = true;

    for (int i = 2; i < count; i++) {
        reset_all_conodes(T, T.getRoot());
        mark(covertex[i]);
        if (T.getNode(T.getRoot()).marked == Marked::MARKED_AND_UNMARKED) {
            // all nodes of T were marked and unmarked <=>
            // R is marked and unmarked
            T.addChild(T.getRoot(), covertex[i]);
        } else if (mark_ever_count == 0) {
            if (T.getNode(T.getRoot()).d == 1) {
                T.addChild(T.getNode(T.getRoot()).first_child, covertex[i]);
            } else {
                auto R1 = T.add(NodeType::COMPLEMENT_NODE, 1);
                auto R2 = T.add(NodeType::UNION_NODE, 0);
                T.addChild(R1, R2);
                T.addChild(R2, T.getRoot());
                T.addChild(R2, covertex[i]);
                T.setRoot(R1);
            }
        } else {
            auto [v, state] = find_lowest();
            if (state != State::COGRAPH) {
                T.clear();
                return state;
            }
            insert_x_to_cotree(v, covertex[i]);
        }
        T.getNode(covertex[i]).in_graph = true;
    }
    T.clear();
    return State::COGRAPH;
}

}  // namespace Koala
