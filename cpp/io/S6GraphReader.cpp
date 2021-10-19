/*
 * S6GraphReader.cpp
 *
 *  Created on: 06.10.2021
 *      Author: Krzysztof Turowski
 */

#include <fstream>

#include <networkit/auxiliary/Enforce.hpp>

#include <io/S6GraphReader.hpp>

namespace Koala {

constexpr index log2(index n) {
  return n > 1 ? 1 + log2(n >> 1) : 0;
}

NetworKit::Graph S6GraphReader::read(const std::string &path) {
    std::ifstream graphFile(path);
    Aux::enforceOpened(graphFile);
    std::string line;
    std::getline(graphFile, line);
    return readline(line);
}

NetworKit::Graph S6GraphReader::readline(const std::string &line) {
    auto it = line.cbegin();
    assert(*it == ':');
    ++it;

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

    auto getBitsGroup = [&]() {
        const char LOW = 0x3f, HIGH = 0x7e;
        const int LENGTH = 6;
        ++it;
        assert(*it >= LOW && *it <= HIGH);
        bits = (bits << LENGTH) + (*it - LOW);
        length += LENGTH;
    };

    index nodes = scanGraphOrder();
    Graph graph(nodes, false, false);

    char MASKS[] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40 };
    index u = 0, v = 0, lognodes = log2(nodes), length = 0;
    while (true) {
        if (length == 0) {
            getBitsGroup();
        }
        --length, v += (bits & MASKS[length]), bits ^= MASKS[length];
        while (length < lognodes) {
            getBitsGroup();
        }
        length -= lognodes, u = bits >> length, bits ^= (u << length);
        if (u >= nodes) {
            break;
        }
        if (u > v) {
            v = u;
        } else {
            graph.addEdge(u, v);
        }
    }
    graph.shrinkToFit();
    return graph;
}

} /* namespace Koala */
