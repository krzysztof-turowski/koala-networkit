#include "coloring/CographVertexColoring.hpp"

namespace Koala {
    void CographVertexColoring::subtree_colors(NetworKit::count v) {
        CoNode V(0,0,0);
        V=cotree.getNode(v);
        if (V.left_son != NetworKit::none) {
            subtree_colors(V.left_son);
        }

        if (V.right_son != NetworKit::none) {
            subtree_colors(V.right_son);
        }

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

    void CographVertexColoring::end_of_coloring(NetworKit::count v) {
        CoNode V(0,0,0);
        V=cotree.getNode(v);
        if (V.left_son != NetworKit::none) {
            color[V.left_son] += color[v];
            end_of_coloring(V.left_son);
        }

        if (V.right_son != NetworKit::none) {
            color[V.right_son] += color[v];
            end_of_coloring(V.right_son);
        }
    }

    void CographVertexColoring::run() {
        hasRun = true;
        NetworKit::count n = graph->numberOfNodes();
        if (!cotree.prepared) {
            cotree.buildTree();
        }
        color.resize(2 * n, 0);
        number_of_colors.resize(2 * n, 0);
        subtree_colors(n);
        end_of_coloring(n);
        for (const auto &u : graph->nodeRange())
        {
            colors[u] = color[u];
        }
    }

    bool CographVertexColoring::checkColoring()
    {
        auto &coloring = getColoring();
        for (const auto &u : graph->nodeRange())
        {
            for (const auto &v : graph->neighborRange(u))
            {
                if (colors[u] == colors[v])
                {
                    return false;
                }
            }
        }
        return true;
    }
}