#include <graph/GraphTools.hpp>

#include "clique/CographClique.hpp"
#include "recognition/CographRecognitionOther.hpp"

namespace Koala {

void CographMaxClique::recurse_run() {
    while (!st.empty()) {
        int v = st.top();
        CoNode &V = cotree.getNode(v);
        if (used[v] == false) {
            used[v] = true;
            if (V.left_son != NetworKit::none) {
                st.push(V.left_son);
            }

            if (V.right_son != NetworKit::none) {
                st.push(V.right_son);
            }
        } else {
            st.pop();
            if (V.type == NodeType::LEAF) {
                subgraph_clique_size[v] = 1;
            } else {
                NetworKit::count l = 0, r = 0;
                if (V.left_son != NetworKit::none) {
                    l = subgraph_clique_size[V.left_son];
                }

                if (V.right_son != NetworKit::none) {
                    r = subgraph_clique_size[V.right_son];
                }

                if (V.type == NodeType::UNION_NODE) {
                    subgraph_clique_size[v] = std::max(l, r);
                } else {
                    subgraph_clique_size[v] = l + r;
                }
            }
        }
    }
}

void CographMaxClique::add_to_set() {
    while (!st.empty()) {
        int v = st.top();
        st.pop();
        CoNode &V = cotree.getNode(v);
        if (V.type == NodeType::LEAF) {
            max_clique.insert(v);
        } else {
            NetworKit::count l = 0, r = 0;
            if (V.left_son != NetworKit::none) {
                l = subgraph_clique_size[V.left_son];
            }
            if (V.right_son != NetworKit::none) {
                r = subgraph_clique_size[V.right_son];
            }
            if (V.type == NodeType::UNION_NODE) {
                if (l >= r) {
                    st.push(V.left_son);
                } else {
                    st.push(V.right_son);
                }
            } else {
                if (V.left_son != NetworKit::none) {
                    st.push(V.left_son);
                }
                if (V.right_son != NetworKit::none) {
                    st.push(V.right_son);
                }
            }
        }
    }
}

void CographMaxClique::run() {
    hasRun = true;
    NetworKit::count n = cotree.graph->numberOfNodes();
    subgraph_clique_size.resize(2 * n + 1);
    used.resize(2 * n + 1);
    st.push(n);
    recurse_run();
    used.resize(2 * n + 1);
    st.push(n);
    add_to_set();
}

NetworKit::count CographMaxClique::bruteForceCliqueSize(NetworKit::Graph &Graph) {
    NetworKit::count n = Graph.numberOfNodes(), ans = 0, flag;
    std::vector<NetworKit::count> st(31), clique_nodes;
    st[0] = 1;
    for (int i = 1; i <= 30; i++) {
        st[i] = st[i - 1] * 2;
    }
    for (int mask = 1; mask < st[n]; mask++) {
        clique_nodes.clear();
        for (int i = 0; i <= 30; i++) {
            if ((st[i] & mask) > 0) {
                clique_nodes.push_back(i);
            }
        }
        flag = 0;

        for (int i = 0; i < clique_nodes.size(); i++) {
            for (int j = i + 1; j < clique_nodes.size(); j++) {
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
