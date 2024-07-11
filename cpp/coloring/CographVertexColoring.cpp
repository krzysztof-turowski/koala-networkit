#include "coloring/CographVertexColoring.hpp"

namespace Koala {
    void CographVertexColoring::subtree_colors() {
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
                    color[v] = 0;
                    number_of_colors[v] = 1;
                } else {
                    if (V.type == NodeType::COMPLEMENT_NODE) {
                        if (V.left_son != NetworKit::none) {
                            number_of_colors[v] += number_of_colors[V.left_son];
                        }

                        if (V.right_son != NetworKit::none) {
                            color[V.right_son] += number_of_colors[V.left_son];
                            number_of_colors[v] += number_of_colors[V.right_son];
                        }
                    } else {
                        if (V.left_son != NetworKit::none) {
                            number_of_colors[v] = number_of_colors[V.left_son];
                        }

                        if (V.right_son != NetworKit::none &&
                            number_of_colors[V.right_son] >
                            number_of_colors[V.left_son]) {
                            number_of_colors[v] = number_of_colors[V.right_son];
                        }
                    }
                }
            }
        }
    }

    void CographVertexColoring::end_of_coloring() {
        while (!st.empty()) {
            int v = st.top();
            CoNode &V = cotree.getNode(v);
            if (used[v] == false) {
                used[v] = true;
                if (V.left_son != NetworKit::none) {
                    color[V.left_son] += color[v];
                    st.push(V.left_son);
                }

                if (V.right_son != NetworKit::none) {
                    color[V.right_son] += color[v];
                    st.push(V.right_son);
                }
            } else {
                st.pop();
            }
        }
    }

    void CographVertexColoring::run() {
        hasRun = true;
        NetworKit::count n = graph->numberOfNodes();
        color.resize(2 * n + 1, 0);
        number_of_colors.resize(2 * n + 1, 0);
        used.resize(2 * n + 1, false);
        st.push(n);
        subtree_colors();
        used.resize(2 * n + 1, false);
        st.push(n);
        end_of_coloring();
        std::cout << number_of_colors[n] << std::endl;
        for (const auto &u: graph->nodeRange()) {
            colors[u] = color[u];
        }
    }

    bool CographVertexColoring::checkColoring() {
        auto &coloring = getColoring();
        for (const auto &u: graph->nodeRange()) {
            for (const auto &v: graph->neighborRange(u)) {
                if (colors[u] == colors[v]) {
                    return false;
                }
            }
        }
        return true;
    }
} /* namespace Koala */
