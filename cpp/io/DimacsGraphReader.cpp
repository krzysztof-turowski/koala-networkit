/*
 * DimacsGraphReader.cpp
 *
 *  Created on: 21.10.2021
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <fstream>
#include <map>
#include <tuple>

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/graph/GraphTools.hpp>

#include <io/DimacsGraphReader.hpp>

namespace Koala {

enum class Format {
    edge,
    min,
    max,
    sp,
    mat
};

std::map<std::string, Format> convert {
    {"edge", Format::edge}, {"min", Format::min}, {"max", Format::max}, {"sp", Format::sp},
    {"mat", Format::mat}
};

NetworKit::Graph create_graph(const std::string &format) {
    NetworKit::Graph graph;
    switch (convert[format]) {
        case Format::edge:
            graph = NetworKit::Graph(0, false, false);
            break;
        case Format::min:
            graph = NetworKit::Graph(0, true, true, true);
            break;
        case Format::max:
        case Format::sp:
            graph = NetworKit::Graph(0, true, true);
            break;
        default:
            throw std::runtime_error("Format not supported");
    }
    return graph;
}


void read_edge(std::ifstream &graphFile, NetworKit::Graph &graph,
    std::map<NetworKit::edgeid, long long> &costs, const std::string &format) {
    NetworKit::node u = 0, v = 0;
    NetworKit::edgeweight w = 0;
    NetworKit::edgeid id;
    int cost, _;

    switch (convert[format]) {
        case Format::edge:
            graphFile >> u >> v;
            graph.addEdge(u - 1, v - 1);
            break;
        case Format::min: 
        graphFile >> u >> v >> _ >> w >> cost;
            id = graph.upperEdgeIdBound();
            graph.addEdge(u - 1, v - 1, w);
            costs[id] = cost;
            break;
        case Format::max:
        case Format::sp:
            graphFile >> u >> v >> w;
            graph.increaseWeight(u - 1, v - 1, w);
            graph.increaseWeight(v - 1, u - 1, 0);
            break;
        default:
            throw std::runtime_error("Format not supported");
    }
}

void read_node(std::ifstream &graphFile, NetworKit::Graph &graph,
        std::map<NetworKit::node, long long> &b,
        NetworKit::node &s, 
        NetworKit::node &t, 
        const std::string &format) {
    NetworKit::node v = NetworKit::none;
    std::string label;
    int supply;

    switch (convert[format]) {
        case Format::min: 
            graphFile >> v >> supply;
            b[v] = supply;
            break;
        default:
            graphFile >> v >> label;
            if (label == "s") {
                s = v - 1;
                break;
            }
            if (label == "t") {
                t = v - 1;
                break;
            }
            throw std::runtime_error("Unknown label");
            
    }
}


NetworKit::Graph DimacsGraphReader::read(std::string_view path) {
    return std::get<0>(read_all(std::string{path}));
}

std::tuple<NetworKit::Graph, 
        std::map<NetworKit::edgeid, long long>,
        std::map<NetworKit::node, long long>, 
        NetworKit::node, NetworKit::node> DimacsGraphReader::read_all_mcf(
        const std::string &path) {
    std::ifstream graphFile(path);
    Aux::enforceOpened(graphFile);

    NetworKit::Graph graph;
    const int MAX = 2048;
    char command = 0;
    std::string format;
    NetworKit::count nodes = 0, edges = 0;
    NetworKit::node s, t;
    s = t = NetworKit::none;
    std::map<NetworKit::node, long long> b;
    std::map<NetworKit::edgeid, long long> costs;

    while (true) {
        graphFile >> command;
        if (graphFile.eof()) {
            break;
        }
        switch (command) {
            case 'p':
                graphFile >> format >> nodes >> edges;
                graph = create_graph(format);
                graph.addNodes(nodes);
                break;
            case 'c':
                break;
            case 'a':
                read_edge(graphFile, graph, costs, format);
                break;
            case 'e':
                if (graph.isDirected()) {
                    graph = NetworKit::GraphTools::toUndirected(graph);
                }
                read_edge(graphFile, graph, costs, format);
                break;
            case 'n':
                read_node(graphFile, graph, b, s, t, format);
                break;
            default:
                throw std::runtime_error("Unknown line type");
        }
        graphFile.ignore(MAX, '\n');
    }
    graph.shrinkToFit();
    return std::make_tuple(graph, costs, b, s, t);
}

std::tuple<NetworKit::Graph, NetworKit::node, NetworKit::node> DimacsGraphReader::read_all(
        const std::string &path) {
    auto [graph, costs, b, s, t] = read_all_mcf(path);
    return { graph, s, t }; 
}

std::tuple<NetworKit::Graph,
        std::map<NetworKit::edgeid, long long>,
        std::map<NetworKit::node, long long>> 
        DimacsGraphReader::read_minimum_cost_flow(const std::string &path) {
    auto [graph, costs, b, s, t] = read_all_mcf(path);
    return { graph, costs, b }; 
}

} /* namespace Koala */
