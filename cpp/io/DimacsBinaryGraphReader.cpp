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

NetworKit::Graph DimacsBinaryGraphReader::read(std::string_view path) {
    std::ifstream graph_file(std::string{path}, std::ios::binary);
    Aux::enforceOpened(graph_file);

    int preamble_size = 0;
    graph_file >> preamble_size;
    const int MAX = 2048;
    graph_file.ignore(MAX, '\n');
    std::string preamble(preamble_size, 0);
    graph_file.read(preamble.data(), preamble_size);
    std::stringstream preamble_file(preamble);

    NetworKit::Graph graph(0, false, false);
    char command = 0;
    std::string format;
    NetworKit::count nodes = 0, edges = 0;
    while (true) {
        preamble_file >> command;
        if (preamble_file.eof()) {
            break;
        }
        switch (command) {
            case 'p':
                preamble_file >> format >> nodes >> edges;
                graph.addNodes(nodes);
                break;
            case 'c':
                break;
            default:
                throw std::runtime_error("Unknown line type");
        }
        preamble_file.ignore(MAX, '\n');
    }

    const unsigned char MASKS[] = { 0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01 };
    std::string data((nodes >> 3) + 1, 0);
    for (NetworKit::node u = 0; u < nodes; u++) {
        graph_file.read(data.data(), (u >> 3) + 1);
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
