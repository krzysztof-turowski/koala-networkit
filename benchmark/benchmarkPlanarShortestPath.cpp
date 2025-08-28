#include <cassert>
#include <iostream>
#include <map>
#include <set>

#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/distance/Dijkstra.hpp>

#include <graph/GraphTools.hpp>
#include <io/DimacsGraphReader.hpp>
#include <io/G6GraphReader.hpp>
#include <shortest_path/FredericksonPlanarSSSP.hpp>
#include <shortest_path/HenzingerPlanarSSSP.hpp>
#include <shortest_path/PlanarSSSP.hpp>

template <typename T>
std::vector<NetworKit::edgeweight> run_algorithm(NetworKit::Graph& G) {
    auto algorithm = T(G, 0, 0);
    algorithm.run();
    return algorithm.getDistances();
}

std::map<std::string, int> ALGORITHM = {{"all", 0}, {"Frederickson", 1}, {"Henzinger", 2}};

void run_g6_tests(const std::string& path, const std::string& algorithm) {
    std::fstream file(path, std::fstream::in);
    std::map<int, int> classification;
    while (true) {
        std::string line;
        file >> line;
        if (!file.good()) {
            break;
        }
        auto G_directed = Koala::G6GraphReader().readline(line);
        auto G = Koala::GraphTools::convertDirectedGraphToUndirected(G_directed, true);
        G_directed = Koala::GraphTools::convertUndirectedGraphToDirected(G, true);

        NetworKit::ConnectedComponents connected_components(G);
        connected_components.run();
        if (connected_components.getComponents().size() > 1) {
            std::cout << "graph is not connected" << std::endl;
            throw std::runtime_error("graph is not connected");
        }

        std::vector<NetworKit::edgeweight> F, H, D;
        std::cout << line << " " << std::flush;
        switch (ALGORITHM[algorithm]) {
            case 0:
                D = run_algorithm<NetworKit::Dijkstra>(G);
                F = run_algorithm<Koala::FredericksonPlanarSSSP>(G);
                H = run_algorithm<Koala::HenzingerPlanarSSSP>(G_directed);
                for (auto node : G.nodeRange()) {
                    assert(std::set({D[node], F[node], H[node]}).size() == 1);
                }
                F.clear();
                H.clear();
                D.clear();
                break;
            case 1:
                run_algorithm<Koala::FredericksonPlanarSSSP>(G);
                break;
            case 2:
                run_algorithm<Koala::HenzingerPlanarSSSP>(G_directed);
                break;
        }
        std::cout << std::endl;
    }
}

void run_dimacs_tests(const std::string& path, const std::string& algorithm) {
    auto G_directed = Koala::DimacsGraphReader().read(path);
    auto G = Koala::GraphTools::convertDirectedGraphToUndirected(G_directed, true);
    G_directed = Koala::GraphTools::convertUndirectedGraphToDirected(G, true);

    std::cout << path << " " << std::flush;
    auto F = run_algorithm<Koala::FredericksonPlanarSSSP>(G);
    auto H = run_algorithm<Koala::HenzingerPlanarSSSP>(G_directed);
    auto D = run_algorithm<NetworKit::Dijkstra>(G);
    for (auto node : G.nodeRange()) {
        assert(std::set({D[node], F[node], H[node]}).size() == 1);
    }
    std::cout << std::endl;
    return;
}

int main(int argc, const char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <algorithm> <file>" << std::endl;
        return 1;
    }
    std::string path(argv[2]);
    auto position = path.find_last_of(".");
    if (path.substr(position + 1) == "g6") {
        run_g6_tests(path, std::string(argv[1]));
    } else if (path.substr(position + 1) == "gr") {
        run_dimacs_tests(path, std::string(argv[1]));
    } else {
        std::cerr << "File type not supported: " << path << std::endl;
    }
    return 0;
}
