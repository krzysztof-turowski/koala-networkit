/*
 * D6GraphReader.cpp
 *
 *  Created on: 27.10.2021
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <fstream>

#include <networkit/auxiliary/Enforce.hpp>

#include <io/D6GraphReader.hpp>

namespace Koala {

NetworKit::Graph D6GraphReader::read(std::string_view path) {
    std::ifstream graphFile(std::string{path});
    Aux::enforceOpened(graphFile);
    std::string line;
    std::getline(graphFile, line);
    return readline(line);
}

NetworKit::Graph D6GraphReader::readline(const std::string &line) {
    auto it = line.cbegin();
    assert(*it == '&');
    ++it;

    const char LOW = 0x3f, HIGH = 0x7e, MASK = 0x20;
    const int SHORT_N = 1, MEDIUM_N = 3, LONG_N = 6, LENGTH = 6;
    NetworKit::count nodes_length = SHORT_N;
    if (*it >= HIGH) {
        nodes_length = MEDIUM_N, ++it;
        if (*it >= HIGH) {
          nodes_length = LONG_N, ++it;
        }
    }

    NetworKit::count nodes = 0;
    for (NetworKit::count i = 0; i < nodes_length; i++, ++it) {
      nodes = (nodes << LENGTH) | (*it - LOW);
    }

    NetworKit::Graph graph(nodes, false, true);
    char mask = 0, bits = 0;
    for (NetworKit::node u = 0; u < nodes; u++) {
        for (NetworKit::node v = 0; v < nodes; v++, mask >>= 1) {
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
