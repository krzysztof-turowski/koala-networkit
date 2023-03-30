/*
 * PerfectGraphColoring.cpp
 *
 *  Created on: 30.03.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <coloring/PerfectGraphColoring.hpp>

#include "perfect/color.h"

namespace Koala {

void PerfectGraphColoring::run() {
    vec<vec<int>> M;
    M.resize(graph->numberOfNodes());
    graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
        M[u].push_back(v), M[v].push_back(u);
    });
    Graph G(M);
    auto C = color(G);
    for (int i = 0; i < C.size(); i++) {
        std::cout << "VERTEX " << i << " COLOR " << C[i] << std::endl;
    }
    hasRun = true;
}

} /* namespace Koala */
