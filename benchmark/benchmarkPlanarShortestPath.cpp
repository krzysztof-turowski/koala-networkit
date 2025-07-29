#include <cassert>
#include <io/DimacsGraphReader.hpp>
#include <io/G6GraphReader.hpp>
#include <iostream>
#include <map>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <set>
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
        auto G = NetworKit::Graph(G_directed.numberOfNodes(), true, false);
        G_directed.forEdges([&](NetworKit::node u, NetworKit::node v) {
            if (!G.hasEdge(u, v) && !G.hasEdge(v, u)) {
                G.addEdge(u, v);
            }
        });
        auto G_directed_both_ways = NetworKit::Graph(G.numberOfNodes(), true, true);
        G.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
            G_directed_both_ways.addEdge(u, v, w);
            G_directed_both_ways.addEdge(v, u, w);
        });
        NetworKit::ConnectedComponents cc(G);
        cc.run();
        if (cc.getComponents().size() > 1) {
            std::cout << "graph is not connected" << std::endl;
            throw std::runtime_error("graph not connected");
            continue;
        }
        std::cout << std::endl;
        std::vector<NetworKit::edgeweight> F, H, D;
        std::cout << line << " " << std::flush;
        switch (ALGORITHM[algorithm]) {
            case 0:
                F = run_algorithm<Koala::FredericksonPlanarSSSP>(G);
                H = run_algorithm<Koala::HenzingerPlanarSSSP>(G_directed_both_ways);
                D = run_algorithm<NetworKit::Dijkstra>(G);
                for (auto node : G.nodeRange()) {
                    assert(H[node] == D[node] && F[node] == D[node]);
                }
                F.clear();
                H.clear();
                D.clear();
                break;
            case 1:
                run_algorithm<Koala::FredericksonPlanarSSSP>(G);
                break;
            case 2:
                run_algorithm<Koala::HenzingerPlanarSSSP>(G_directed_both_ways);
                break;
        }
        std::cout << std::endl;
    }
}

void run_dimacs_tests(const std::string& path, const std::string& algorithm) {
    auto G_directed = Koala::DimacsGraphReader().read(path);
    auto G = NetworKit::Graph(G_directed.numberOfNodes(), true, false);
    G_directed.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        if (!G.hasEdge(u, v) && !G.hasEdge(v, u) && w > 0) {
            G.addEdge(u, v, w);
        }
    });
    auto G_directed_both_ways = NetworKit::Graph(G.numberOfNodes(), true, true);
    G.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        G_directed_both_ways.addEdge(u, v, w);
        G_directed_both_ways.addEdge(v, u, w);
    });
    std::cout << path << " " << std::flush;
    auto F = run_algorithm<Koala::FredericksonPlanarSSSP>(G);
    auto H = run_algorithm<Koala::HenzingerPlanarSSSP>(G_directed_both_ways);
    auto D = run_algorithm<NetworKit::Dijkstra>(G);
    for (auto node : G.nodeRange()) {
        assert(F[node] == H[node] && H[node] == D[node]);
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
