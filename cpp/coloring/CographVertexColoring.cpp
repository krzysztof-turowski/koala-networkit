#include <algorithm>

#include <coloring/CographVertexColoring.hpp>

namespace Koala {

void CographVertexColoring::subtree_colors() {
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
                color[v] = 0;
                number_of_colors[v] = 1;
            } else if (V.type == NodeType::COMPLEMENT_NODE) {
                NetworKit::count color_offset = 0;
                for (auto child = V.first_child; child != NetworKit::none;
                        child = cotree.getNode(child).next_sibling) {
                    color[child] += color_offset;
                    color_offset += number_of_colors[child];
                    number_of_colors[v] += number_of_colors[child];
                }
            } else {
                for (auto child = V.first_child; child != NetworKit::none;
                        child = cotree.getNode(child).next_sibling) {
                    number_of_colors[v] = std::max(number_of_colors[v], number_of_colors[child]);
                }
            }
        }
    }
}

void CographVertexColoring::end_of_coloring() {
    while (!st.empty()) {
        NetworKit::node v = st.top();
        Conode &V = cotree.getNode(v);
        if (used[v] == false) {
            used[v] = true;
            for (auto child = V.first_child; child != NetworKit::none;
                    child = cotree.getNode(child).next_sibling) {
                color[child] += color[v];
                st.push(child);
            }
        } else {
            st.pop();
        }
    }
}

void CographVertexColoring::run() {
    hasRun = true;
    color.assign(cotree.upperNodeIdBound(), 0);
    number_of_colors.assign(cotree.upperNodeIdBound(), 0);
    used.assign(cotree.upperNodeIdBound(), false);
    st.push(cotree.getRoot());
    subtree_colors();
    used.assign(cotree.upperNodeIdBound(), false);
    st.push(cotree.getRoot());
    end_of_coloring();
    for (const auto &u : graph->nodeRange()) {
        colors[u] = color[u];
    }
}

bool CographVertexColoring::checkColoring() {
    auto &colors = getColoring();
    for (const auto &u : graph->nodeRange()) {
        for (const auto &v : graph->neighborRange(u)) {
            if (colors.at(u) == colors.at(v)) {
                return false;
            }
        }
    }
    return true;
}

} /* namespace Koala */
