/*
 * DimacsGraphWriter.cpp
 *
 *  Created on: 27.10.2021
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <fstream>

#include <networkit/auxiliary/Enforce.hpp>

#include <io/DimacsGraphWriter.hpp>

namespace Koala {

void DimacsGraphWriter::write(const NetworKit::Graph &G, const std::string &path) {
    std::ofstream graphFile(path);
    Aux::enforceOpened(graphFile);

    graphFile << "p edge " << G.numberOfNodes() << ' ' << G.numberOfEdges() << std::endl;
    std::string edge_type = G.isDirected() ? "a" : "e";
    G.forEdges([&](NetworKit::node u, NetworKit::node v) {
        graphFile << edge_type << ' ' << u + 1 << ' ' << v + 1 << std::endl;
    });
}

} /* namespace Koala */
