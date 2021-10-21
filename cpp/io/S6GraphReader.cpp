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

constexpr NetworKit::count log2(NetworKit::count n) {
    return n > 1 ? 1 + log2(n >> 1) : n;
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

    const char LOW = 0x3f, HIGH = 0x7e;
    const int SHORT_N = 1, MEDIUM_N = 2, LONG_N = 6, LENGTH = 6;
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

    NetworKit::Graph graph(nodes, false, false);
    const char MASKS[] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40 };
    NetworKit::count lognodes = log2(nodes - 1), length = 0, bits = 0;
    NetworKit::node u = 0, v = 0;
    while (true) {
        if (length == 0) {
            if (it == line.cend()) {
                break;
            }
            assert(*it >= LOW && *it <= HIGH);
            bits = (bits << LENGTH) | (*it - LOW);
            length += LENGTH, ++it;
        }
        if (bits & MASKS[--length]) {
            v++, bits &= ~MASKS[length];
        }
        while (length < lognodes) {
            if (it == line.cend()) {
                break;
            }
            assert(*it >= LOW && *it <= HIGH);
            bits = (bits << LENGTH) | (*it - LOW);
            length += LENGTH, ++it;
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
