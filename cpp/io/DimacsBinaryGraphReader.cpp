/*
 * DimacsBinaryGraphReader.cpp
 *
 *  Created on: 27.10.2021
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <fstream>
#include <sstream>

#include <networkit/auxiliary/Enforce.hpp>

#include <io/DimacsBinaryGraphReader.hpp>

namespace Koala {

NetworKit::Graph DimacsBinaryGraphReader::read(const std::string &path) {
    std::ifstream graphFile(path, std::ifstream::binary);
    Aux::enforceOpened(graphFile);

    int preamble_size = 0;
    graphFile >> preamble_size;
    const int MAX = 2048;
    graphFile.ignore(MAX, '\n');
    std::string preamble(preamble_size, 0);
    graphFile.read(preamble.data(), preamble_size + 1);
    std::stringstream preambleFile(preamble);

    NetworKit::Graph graph(0, false, false);
    char command = 0;
    std::string format;
    NetworKit::count nodes = 0, edges = 0;
    while (true) {
        preambleFile >> command;
        if (preambleFile.eof()) {
            break;
        }
        switch (command) {
            case 'p':
                preambleFile >> format >> nodes >> edges;
                graph.addNodes(nodes);
                break;
            case 'c':
                break;
            default:
                throw std::runtime_error("Unknown line type");
        }
        preambleFile.ignore(MAX, '\n');
    }

    const unsigned char MASKS[] = { 0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01 };
    std::string data((nodes >> 3) + 1, 0);
    for (NetworKit::node u = 0; u < nodes; u++) {
        graphFile.read(data.data(), (u >> 3) + 1);
        for (NetworKit::node v = 0; v <= u; v++) {
            unsigned char mask = MASKS[v & 7];
            if ((data[v >> 3] & mask) == mask) {
                graph.addEdge(u, v);
            }
        }
    }
    graph.shrinkToFit();
    return graph;
}

} /* namespace Koala */
