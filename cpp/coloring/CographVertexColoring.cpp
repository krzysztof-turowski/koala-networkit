#include "coloring/CographVertexColoring.hpp"

namespace Koala {
    void CographVertexColoring::SubtreeColors(NetworKit::count v) {
        if (recognition->cotree->nodes[v].left_son != -1) {
            SubtreeColors(recognition->cotree->nodes[v].left_son);
        }

        if (recognition->cotree->nodes[v].right_son != -1) {
            SubtreeColors(recognition->cotree->nodes[v].right_son);
        }

        if (recognition->cotree->nodes[v].type == 2) {
            color[v] = 0;
            number_of_colors[v] = 1;
        } else {
            if (recognition->cotree->nodes[v].type == 1) {
                if (recognition->cotree->nodes[v].left_son != -1) {
                    number_of_colors[v] += number_of_colors[recognition->cotree->nodes[v].left_son];
                }

                if (recognition->cotree->nodes[v].right_son != -1) {
                    color[recognition->cotree->nodes[v].right_son] += number_of_colors[recognition->cotree->nodes[v].left_son];
                    number_of_colors[v] += number_of_colors[recognition->cotree->nodes[v].right_son];
                }
            } else {
                if (recognition->cotree->nodes[v].left_son != -1) {
                    number_of_colors[v] = number_of_colors[recognition->cotree->nodes[v].left_son];
                }

                if (recognition->cotree->nodes[v].right_son != -1 &&
                    number_of_colors[recognition->cotree->nodes[v].right_son] >
                    number_of_colors[recognition->cotree->nodes[v].left_son]) {
                    number_of_colors[v] = number_of_colors[recognition->cotree->nodes[v].right_son];
                }
            }

        }
    }

    void CographVertexColoring::EndOfColoring(NetworKit::count v) {
        if (recognition->cotree->nodes[v].left_son != -1) {
            color[recognition->cotree->nodes[v].left_son] += color[v];
            EndOfColoring(recognition->cotree->nodes[v].left_son);
        }

        if (recognition->cotree->nodes[v].right_son != -1) {
            color[recognition->cotree->nodes[v].right_son] += color[v];
            EndOfColoring(recognition->cotree->nodes[v].right_son);
        }
    }

    void CographVertexColoring::run() {
        NetworKit::count n = graph->numberOfNodes();
        recognition->run();
        if (recognition->cotree->prepared == false) {
            recognition->cotree->BuildTree();
        }
        color.resize(2 * n, 0);
        number_of_colors.resize(2 * n, 0);
        SubtreeColors(n);
        EndOfColoring(n);
        for (int i = 0; i < n; i++) {
            colors[i] = color[i];
        }
    }

    NetworKit::count CographVertexColoring::GetColor(NetworKit::count i) {
        if (recognition->cotree->prepared == false) {
            recognition->cotree->BuildTree();
        }
        return color[i];
    }

    bool CographVertexColoring::CheckColoring() {
        NetworKit::count n = graph->numberOfNodes(), i, j;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                if (graph->hasEdge(i, j) && GetColor(i) == GetColor(j)) {
                    return false;
                }
            }
        }
        return true;
    }
}