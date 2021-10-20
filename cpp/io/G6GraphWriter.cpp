/*
 * G6GraphWriter.cpp
 *
 *  Created on: 20.10.2021
 *      Author: Krzysztof Turowski
 */

#include <fstream>

#include <networkit/auxiliary/Enforce.hpp>

#include <io/G6GraphWriter.hpp>

namespace Koala {

void G6GraphWriter::write(const NetworKit::Graph &G, const std::string &path) {
    std::ofstream graphFile(path);
    Aux::enforceOpened(graphFile);
    std::string g6String = writeline(G);
    graphFile << g6String << std::endl;
}

std::string G6GraphWriter::writeline(const NetworKit::Graph &G) {
    std::string output;

    const char LOW = 0x3f, HIGH = 0x7e;
    const NetworKit::count SHORT_N = 1, MEDIUM_N = 2, LONG_N = 6, LENGTH = 6;
    const NetworKit::count SHORT_NODES = 63, LONG_NODES = 258048;
    NetworKit::count nodes_length = SHORT_N;
    if (nodes_length >= SHORT_NODES) {
        output.push_back(HIGH), nodes_length = MEDIUM_N;
        if (nodes_length >= LONG_NODES) {
            output.push_back(HIGH), nodes_length = LONG_N;
        }
    }

    NetworKit::count nodes = G.numberOfNodes();
    std::string nodes_string(nodes_length, LOW);
    for (int i = nodes_length - 1; i >= 0; i--) {
        nodes_string[i] += (nodes & LOW), nodes >>= LENGTH;
    }
    output += nodes_string;

    NetworKit::count start = output.size(), index = 0, shift = 0;
    int length = (G.numberOfNodes() * (G.numberOfNodes() - 1) / 2) / LENGTH + 1;
    output.append(length, 0x0);
    const char MASKS[] = { 0x20, 0x10, 0x08, 0x04, 0x02, 0x01 };
    for (const NetworKit::index v : G.nodeRange()) {
        for (const NetworKit::index u : G.neighborRange(v)) {
            if (u < v) {
                int position = shift + u;
                output[start + position / LENGTH] |= MASKS[position % LENGTH];
            }
        }
        shift += index, index++;
    }
    for (int i = start; i < output.size(); i++) {
        output[i] += LOW;
    }
    return output;
}

} /* namespace Koala */