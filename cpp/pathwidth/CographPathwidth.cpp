#include <algorithm>

#include <pathwidth/CographPathwidth.hpp>

namespace Koala {

void CographPathwidth::subtree_size() {
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
            subtree_sizes[v] = 1;
            for (auto child = V.first_child; child != NetworKit::none;
                    child = cotree.getNode(child).next_sibling) {
                subtree_sizes[v] += subtree_sizes[child];
            }
        }
    }
}

void CographPathwidth::pathwidth() {
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
                path[v] = 1;
            } else if (V.type == NodeType::UNION_NODE) {
                for (auto child = V.first_child; child != NetworKit::none;
                        child = cotree.getNode(child).next_sibling) {
                    path[v] = std::max(path[v], path[child]);
                }
            } else {
                for (auto child = V.first_child; child != NetworKit::none;
                        child = cotree.getNode(child).next_sibling) {
                    path[v] = std::max(
                        path[v], path[child] + subtree_sizes[v] - subtree_sizes[child] - 1);
                }
            }
        }
    }
}

void CographPathwidth::run() {
    hasRun = true;
    used.resize(cotree.upperNodeIdBound(), false);
    path.resize(cotree.upperNodeIdBound(), 0);
    subtree_sizes.resize(cotree.upperNodeIdBound(), 0);
    st.push(cotree.getRoot());
    subtree_size();
    for (NetworKit::node i = 0; i < cotree.upperNodeIdBound(); i++) {
        used[i] = false;
    }
    st.push(cotree.getRoot());
    pathwidth();
    width = path[cotree.getRoot()];
}

} /* namespace Koala */
