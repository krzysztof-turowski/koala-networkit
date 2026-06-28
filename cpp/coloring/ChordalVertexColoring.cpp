/*
 * ChordalVertexColoring.cpp
 *
 *  Created on: 2026-06-26
 *      Author: Mateusz Przebieracz
 */

#include "coloring/ChordalVertexColoring.hpp"

#include <vector>

namespace Koala {

ChordalVertexColoring::ChordalVertexColoring(
    const NetworKit::Graph &graph, const PerfectEliminationOrdering &peo)
    : VertexColoring(graph), peo(peo) {
}

void ChordalVertexColoring::run() {
    hasRun = true;
    NetworKit::count n = graph->numberOfNodes();

    if (n == 0) return;

    if (!peo.has_value()) {
        MaximumCardinalitySearchChordalGraphRecognition recognizer(*graph);
        recognizer.run();
        peo = recognizer.getPEO();
    }

    std::vector<uint8_t> used_colors(n + 1, false);
    std::vector<NetworKit::count> color_arr(graph->upperNodeIdBound(), 0);
    std::vector<NetworKit::node> n_plus;

    for (NetworKit::index i = n; i >= 1; i--) {
        NetworKit::node v = peo->alpha_inv[i];
        n_plus.clear();

        for (auto w : graph->neighborRange(v)) {
            if (peo->alpha[w] > i) {
                n_plus.push_back(w);
                used_colors[color_arr[w]] = true;
            }
        }

        int c = 1;
        while (used_colors[c]) {
            c++;
        }
        color_arr[v] = c;

        for (auto w : n_plus) {
            used_colors[color_arr[w]] = false;
        }
    }

    colors.clear();
    for (auto v : graph->nodeRange()) {
        colors.emplace_hint(colors.end(), v, color_arr[v]);
    }
}

} /* namespace Koala */
