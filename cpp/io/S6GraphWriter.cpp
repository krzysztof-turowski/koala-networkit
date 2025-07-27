/*
 * S6GraphWriter.cpp
 *
 *  Created on: 20.10.2021
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <fstream>

#include <networkit/auxiliary/Enforce.hpp>

#include <io/S6GraphWriter.hpp>

namespace Koala {

constexpr NetworKit::count log2(NetworKit::count n) {
    return n > 1 ? 1 + log2(n >> 1) : n;
}

void S6GraphWriter::write(const NetworKit::Graph &G, std::string_view path) {
    std::ofstream graphFile(std::string{path});
    Aux::enforceOpened(graphFile);
    std::string s6String = writeline(G);
    graphFile << s6String << std::endl;
}

std::string S6GraphWriter::writeline(const NetworKit::Graph &G) {
    std::string output = ":";

    const char LOW = 0x3f, HIGH = 0x7e;
    const NetworKit::count SHORT_N = 1, MEDIUM_N = 2, LONG_N = 6, LENGTH = 6;
    const NetworKit::count SHORT_NODES = 63, LONG_NODES = 258048;
    NetworKit::count nodes_length = SHORT_N;
    if (G.numberOfNodes() >= SHORT_NODES) {
        output.push_back(HIGH), nodes_length = MEDIUM_N;
        if (G.numberOfNodes() >= LONG_NODES) {
            output.push_back(HIGH), nodes_length = LONG_N;
        }
    }

    NetworKit::count nodes = G.numberOfNodes();
    std::string nodes_string(nodes_length, LOW);
    for (int i = nodes_length - 1; i >= 0; i--) {
        nodes_string[i] += (nodes & LOW), nodes >>= LENGTH;
    }
    output += nodes_string;

    NetworKit::count lognodes = log2(G.numberOfNodes() - 1) + 1, flag = 1 << (lognodes - 1);
    NetworKit::count bits = 0, length = 0;
    NetworKit::index v_previous = 0;
    const int MASKS[] = {
        0x0, 0x1, 0x3, 0x7, 0xf, 0x1f, 0x3f, 0x7f, 0xff, 0x1ff, 0x3ff, 0x7ff, 0xfff };
    for (const NetworKit::node v : G.nodeRange()) {
        for (const NetworKit::node u : G.neighborRange(v)) {
            if (u > v) {
                continue;
            }
            if (v == v_previous) {
                bits = (bits << lognodes) | u, length += lognodes;
            } else if (v == v_previous + 1) {
                bits = (bits << lognodes) | flag | u, length += lognodes;
            } else {
                bits = (bits << lognodes) | flag | v, length += lognodes;
                while (length >= LENGTH) {
                    length -= LENGTH, output.push_back((bits >> length) + LOW);
                    bits &= MASKS[length];
                }
                bits = (bits << lognodes) | u, length += lognodes;
            }
            v_previous = v;
            while (length >= LENGTH) {
                length -= LENGTH, output.push_back((bits >> length) + LOW);
                bits &= MASKS[length];
            }
        }
    }
    if (length > 0) {
        int special_case = static_cast<int>(
            G.upperNodeIdBound() == flag && v_previous == G.upperNodeIdBound() - 2);
        char padding = (bits << (LENGTH - length)) + MASKS[LENGTH - length - special_case];
        output.push_back(padding + LOW);
    }
    return output;
}

} /* namespace Koala */
