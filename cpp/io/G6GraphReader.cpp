/*
 * G6GraphReader.cpp
 *
 *  Created on: 06.10.2021
 *      Author: Krzysztof Turowski
 */

#include <fstream>

#include <networkit/auxiliary/Enforce.hpp>

#include <io/G6GraphReader.hpp>

namespace Koala {

NetworKit::Graph G6GraphReader::read(const std::string &path) {
    std::ifstream graphFile(path);
    Aux::enforceOpened(graphFile);
    std::string line;
    std::getline(graphFile, line);
    return readline(line);
}

NetworKit::Graph G6GraphReader::readline(const std::string &line) {
    auto it = line.cbegin();

    const char LOW = 0x3f, HIGH = 0x7e, MASK = 0x20;
    const int SHORT_N = 1, MEDIUM_N = 2, LONG_N = 6, LENGTH = 6;
    int nodes_length = SHORT_N;
    if (*it >= HIGH) {
        nodes_length = MEDIUM_N, ++it;
        if (*it >= HIGH) {
          nodes_length = LONG_N, ++it;
        }
    }

    NetworKit::index nodes = 0;
    for (int i = 0; i < nodes_length; i++, ++it) {
      nodes = (nodes << LENGTH) | (*it - LOW);
    }

    NetworKit::Graph graph(nodes, false, false);
    char mask = 0, bits = 0;
    for (NetworKit::index v = 1; v < nodes; v++) {
        for (NetworKit::index u = 0; u < v; u++, mask >>= 1) {
            if (!mask) {
                assert(*it >= LOW && *it <= HIGH);
                bits = *it - LOW, mask = MASK, ++it;
            }
            if (bits & mask) {
                graph.addEdge(u, v);
            }
        }
    }
    graph.shrinkToFit();
    return graph;
}

} /* namespace Koala */
