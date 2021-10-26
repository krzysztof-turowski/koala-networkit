/*
 * DimacsGraphReader.cpp
 *
 *  Created on: 21.10.2021
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <fstream>
#include <iostream>

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/graph/GraphTools.hpp>

#include <io/DimacsGraphReader.hpp>

namespace Koala {

NetworKit::Graph DimacsGraphReader::read(const std::string &path) {
    std::ifstream graphFile(path);
    Aux::enforceOpened(graphFile);

    NetworKit::Graph graph(0, false, true);
    const int MAX = 2048;
    char command = 0;
    std::string format;
    NetworKit::count nodes = 0, edges = 0;
    NetworKit::node u = 0, v = 0;
    while (true) {
        graphFile >> command;
        if (graphFile.eof()) {
          break;
        }
        switch (command) {
            case 'p':
                graphFile >> format >> nodes >> edges;
                graph.addNodes(nodes);
                break;
            case 'c':
                break;
            case 'a':
                graphFile >> u >> v;
                graph.addEdge(u - 1, v - 1);
                break;
            case 'e':
                if (graph.isDirected()) {
                  graph = NetworKit::GraphTools::toUndirected(graph);
                }
                graphFile >> u >> v;
                graph.addEdge(u - 1, v - 1);
                break;
            case 'n':
                throw std::runtime_error("Labels not implemented");
            default:
                throw std::runtime_error("Unknown line type");
        }
        graphFile.ignore(MAX, '\n');
    }
    graph.shrinkToFit();
    return graph;
}

} /* namespace Koala */
