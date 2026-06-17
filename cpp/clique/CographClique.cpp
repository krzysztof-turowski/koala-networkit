#include <algorithm>
#include <vector>

#include "graph/GraphTools.hpp"

#include "clique/CographClique.hpp"
#include "recognition/CographRecognition.hpp"

namespace Koala {

void CographMaxClique::recurse_run() {
    while (!st.empty()) {
        NetworKit::node v = st.top();
        Conode &V = cotree.getNode(v);
        if (used[v] == false) {
            used[v] = true;
            for (auto child = V.first_child; child != NetworKit::none;
                    child = cotree.getNode(child).next_sibling) {
                st.push(child);
            }
        } else {
            st.pop();
            if (V.type == NodeType::LEAF) {
                subgraph_clique_size[v] = 1;
            } else if (V.type == NodeType::UNION_NODE) {
                for (auto child = V.first_child; child != NetworKit::none;
                        child = cotree.getNode(child).next_sibling) {
                    subgraph_clique_size[v] = std::max(
                        subgraph_clique_size[v], subgraph_clique_size[child]);
                }
            } else {
                for (auto child = V.first_child; child != NetworKit::none;
                        child = cotree.getNode(child).next_sibling) {
                    subgraph_clique_size[v] += subgraph_clique_size[child];
                }
            }
        }
    }
}

void CographMaxClique::add_to_set() {
    while (!st.empty()) {
        NetworKit::node v = st.top();
        st.pop();
        Conode &V = cotree.getNode(v);
        if (V.type == NodeType::LEAF) {
            max_clique.insert(v);
        } else if (V.type == NodeType::UNION_NODE) {
            NetworKit::node best = NetworKit::none;
            for (auto child = V.first_child; child != NetworKit::none;
                    child = cotree.getNode(child).next_sibling) {
                if (best == NetworKit::none
                        || subgraph_clique_size[child] > subgraph_clique_size[best]) {
                    best = child;
                }
            }
            if (best != NetworKit::none) {
                st.push(best);
            }
        } else {
            for (auto child = V.first_child; child != NetworKit::none;
                    child = cotree.getNode(child).next_sibling) {
                st.push(child);
            }
        }
    }
}

void CographMaxClique::run() {
    hasRun = true;
    subgraph_clique_size.resize(cotree.upperNodeIdBound());
    used.resize(cotree.upperNodeIdBound());
    st.push(cotree.getRoot());
    recurse_run();
    used.resize(cotree.upperNodeIdBound());
    st.push(cotree.getRoot());
    add_to_set();
}

NetworKit::count CographMaxClique::bruteForceCliqueSize(NetworKit::Graph &Graph) {
    NetworKit::count n = Graph.numberOfNodes(), ans = 0, flag;
    std::vector<NetworKit::count> st(31), clique_nodes;
    st[0] = 1;
    for (int i = 1; i <= 30; i++) {
        st[i] = st[i - 1] * 2;
    }
    for (NetworKit::count mask = 1; mask < st[n]; mask++) {
        clique_nodes.clear();
        for (int i = 0; i <= 30; i++) {
            if ((st[i] & mask) > 0) {
                clique_nodes.push_back(i);
            }
        }
        flag = 0;

        for (std::size_t i = 0; i < clique_nodes.size(); i++) {
            for (std::size_t j = i + 1; j < clique_nodes.size(); j++) {
                if (Graph.hasEdge(clique_nodes[i], clique_nodes[j]) == false) {
                    flag = 1;
                    break;
                }
            }
            if (flag == 1) {
                break;
            }
        }
        if (flag == 0 && clique_nodes.size() > ans) {
            ans = clique_nodes.size();
        }
    }
    return ans;
}
} /* namespace Koala */
