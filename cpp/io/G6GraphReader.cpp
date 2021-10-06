/*
 * G6GraphReader.cpp
 *
 *  Created on: 06.10.2021
 *      Author: Krzysztof Turowski
 */

#include <fstream>

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/io/G6GraphReader.hpp>

namespace Koala {

Graph G6GraphReader::read(const std::string &path) {
    std::ifstream graphFile(path);
    Aux::enforceOpened(graphFile);
    std::string line;
    std::getline(graphFile, line);
    return readline(line);
}

Graph G6GraphReader::readline(const std::string &line) {
    auto it = line.cbegin();

    auto scanGraphOrder = [&it]() -> index {
        const char LOW = 0x3f, HIGH = 0x7e;
        const int SHORT_N = 2, LONG_N = 6, LENGTH = 6;
        index nodes = 0;
        if (*it < HIGH) {
            nodes = *it - LOW;
            return nodes;
        }
        ++it;
        if (*it < HIGH) {
            for (int i = 0; i < SHORT_N; i++, ++it) {
              nodes = (nodes << LENGTH) | (*it - LOW);
            }
            return nodes;
        }
        ++it;
        for (int i = 0; i < LONG_N; i++, ++it) {
          nodes = (nodes << LENGTH) | (*it - LOW);
        }
        return nodes;
    };

    auto scanNode = [&](int &v, char &mask) -> char {
        const char LOW = 0x3f, HIGH = 0x7e, MASK = 0x20;
        for (int u = 0; u < v; u++, mask >>= 1) {
            if (!mask) {
                ++it;
                assert(*it >= LOW && *it <= HIGH);
                *it -= LOW;
                mask = MASK;
            }
            if (*it & mask) {
                graph.addEdge(u, v);
            }
        }
        return mask;
    };

    index nodes = scanGraphOrder();
    Graph graph(nodes, false, false);
    char mask = 0;
    for (int v = 1; v < nodes; v++) {
        mask = scanNode(v, mask);
    }
    graph.shrinkToFit();
    return graph;
}

} /* namespace Koala */
