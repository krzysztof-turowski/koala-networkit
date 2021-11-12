/*
 * D6GraphWriter.cpp
 *
 *  Created on: 27.10.2021
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <fstream>

#include <networkit/auxiliary/Enforce.hpp>

#include <io/D6GraphWriter.hpp>

namespace Koala {

void D6GraphWriter::write(const NetworKit::Graph &G, const std::string &path) {
    std::ofstream graphFile(path);
    Aux::enforceOpened(graphFile);
    std::string d6String = writeline(G);
    graphFile << d6String << std::endl;
}

std::string D6GraphWriter::writeline(const NetworKit::Graph &G) {
    std::string output = "&";

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

    NetworKit::count start = output.size(), shift = 0;
    int length = (G.numberOfNodes() * G.numberOfNodes()) / LENGTH + 1;
    output.append(length, 0x0);
    const char MASKS[] = { 0x20, 0x10, 0x08, 0x04, 0x02, 0x01 };
    for (const NetworKit::node u : G.nodeRange()) {
        for (const NetworKit::node v : G.neighborRange(u)) {
            int position = shift + v;
            output[start + position / LENGTH] |= MASKS[position % LENGTH];
        }
        shift += G.numberOfNodes();
    }
    for (unsigned i = start; i < output.size(); i++) {
        output[i] += LOW;
    }
    return output;
}

} /* namespace Koala */
