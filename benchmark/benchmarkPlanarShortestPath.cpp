#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/distance/Dijkstra.hpp>

#include <benchmark/utils.hpp>
#include <graph/GraphTools.hpp>
#include <shortest_path/PlanarSSSP.hpp>

template <typename T>
std::vector<NetworKit::edgeweight> run_algorithm(NetworKit::Graph& G) {
    auto algorithm = T(G, 0, 0);
    algorithm.run();
    return algorithm.getDistances();
}

enum class Algorithm : uint32_t { ALL, FREDERICKSON, HENZINGER };

std::map<std::string, Algorithm> ALGORITHM = {
    {"all", Algorithm::ALL},
    {"Frederickson", Algorithm::FREDERICKSON},
    {"Henzinger", Algorithm::HENZINGER}
};

void choose_algorithm(NetworKit::Graph &G, NetworKit::Graph &G_directed, Algorithm algorithm) {
    std::vector<NetworKit::edgeweight> F, H, D;
    switch (algorithm) {
    case Algorithm::ALL:
        D = run_algorithm<NetworKit::Dijkstra>(G);
        F = run_algorithm<Koala::FredericksonPlanarSSSP>(G);
        H = run_algorithm<Koala::HenzingerPlanarSSSP>(G_directed);
        for (auto node : G.nodeRange()) {
            assert(std::set({D[node], F[node], H[node]}).size() == 1);
        }
        break;
    case Algorithm::FREDERICKSON:
        run_algorithm<Koala::FredericksonPlanarSSSP>(G);
        break;
    case Algorithm::HENZINGER:
        run_algorithm<Koala::HenzingerPlanarSSSP>(G_directed);
        break;
    default:
        throw std::logic_error("Unhandled algorithm");
    }
}

void run_test(NetworKit::Graph G_directed, Algorithm algorithm) {
    auto G = Koala::GraphTools::convertDirectedGraphToUndirected(G_directed, true);
    G_directed = Koala::GraphTools::convertUndirectedGraphToDirected(G, true);
    NetworKit::ConnectedComponents connected_components(G);
    connected_components.run();
    if (connected_components.getComponents().size() > 1) {
        throw std::runtime_error("Graph is not connected");
    }
    choose_algorithm(G, G_directed, algorithm);
}

int main(int argc, const char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <algorithm> <file>" << std::endl;
        return 1;
    }
    std::string algorithm_name(argv[1]);
    auto algorithm = ALGORITHM.find(algorithm_name);
    if (algorithm == ALGORITHM.end()) {
        throw std::invalid_argument("Unknown algorithm: " + algorithm_name);
    }
    Koala::Benchmark::executeForEachGraph(
        argv[2], [&](const std::string &label, NetworKit::Graph G) {
        std::cout << label << " " << std::flush;
        run_test(std::move(G), algorithm->second);
        std::cout << std::endl;
    });
    return 0;
}
