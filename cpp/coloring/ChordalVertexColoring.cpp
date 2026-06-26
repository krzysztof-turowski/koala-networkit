#include "coloring/ChordalVertexColoring.hpp"

#include <vector>

namespace Koala {

ChordalVertexColoring::ChordalVertexColoring(NetworKit::Graph &graph, const PerfectEliminationOrdering &peo)
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
    colors.clear();

    for (NetworKit::index i = n; i >= 1; i--) {
        NetworKit::node v = peo->alpha_inv[i];
        std::vector<NetworKit::node> n_plus;

        for (auto w : graph->neighborRange(v)) {
            if (peo->alpha[w] > i) {
                n_plus.push_back(w);
                used_colors[colors[w]] = true;
            }
        }

        int c = 1;
        while (used_colors[c]) {
            c++;
        }
        colors[v] = c;

        for (auto w : n_plus) {
            used_colors[colors[w]] = false;
        }
    }
}

} /* namespace Koala */