#include "pathwidth/CographPathwidth.hpp"

namespace Koala {

void CographPathwidth::subtree_size() {
    while (!st.empty()) {
        int v = st.top();
        Conode &V = cotree.getNode(v);
        if (used[v] == false) {
            used[v] = true;
            if (V.left_son != NetworKit::none) {
                st.push(V.left_son);
            }
            if (V.right_son != NetworKit::none) {
                st.push(V.right_son);
            }
        } else {
            V.size = 1;
            st.pop();
            if (V.left_son != NetworKit::none) {
                V.size += cotree.getNode(V.left_son).size;
            }

            if (V.right_son != NetworKit::none) {
                V.size += cotree.getNode(V.right_son).size;
            }
        }
    }
}

void CographPathwidth::pathwidth() {
    while (!st.empty()) {
        int v = st.top();
        Conode &V = cotree.getNode(v);
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
                path[v] = 1;
            } else {
                NetworKit::count l = 0, r = 0;
                if (V.left_son != NetworKit::none) {
                    l = path[V.left_son];
                }
                if (V.right_son != NetworKit::none) {
                    r = path[V.right_son];
                }
                if (V.type == NodeType::UNION_NODE) {
                    path[v] = std::max(l, r);
                } else {
                    NetworKit::count pathwidth = 0;
                    if (V.left_son != NetworKit::none) {
                        pathwidth = std::max(pathwidth, r + cotree.getNode(V.left_son).size);
                    }
                    if (V.right_son != NetworKit::none) {
                        pathwidth = std::max(pathwidth, l + cotree.getNode(V.right_son).size);
                    }
                    path[v] = pathwidth;
                }
            }
        }
    }
}

void CographPathwidth::run() {
    hasRun = true;
    int n = cotree.graph->numberOfNodes();
    used.resize(2 * n + 1, false);
    path.resize(2 * n + 1, 0);
    st.push(n);
    subtree_size();
    for (int i = 0; i < 2 * n + 1; i++) {
        used[i] = false;
    }
    st.push(n);
    pathwidth();
    width = path[n];
}

} /* namespace Koala */
