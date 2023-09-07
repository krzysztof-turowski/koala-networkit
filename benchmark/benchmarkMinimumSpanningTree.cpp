#include <cassert>
#include <iostream>
#include <map>
#include <set>

#include <networkit/graph/GraphTools.hpp>

#include <io/G6GraphReader.hpp>
#include <io/DimacsGraphReader.hpp>
#include <mst/MinimumSpanningTree.hpp>

template <typename T>
NetworKit::edgeweight run_algorithm(NetworKit::Graph &G) {
    auto algorithm = T(G);
    algorithm.run();
    auto &spanning_tree = algorithm.getForest();
    std::cout << spanning_tree.totalEdgeWeight() << " " << std::flush;
    algorithm.check();
    return spanning_tree.totalEdgeWeight();
}

std::map<std::string, int> ALGORITHM = {
    { "exact", 0 },
    { "Kruskal", 1 }, { "Prim", 2 }, { "Boruvka", 3 }, { "KKT", 4 }
};

void run_g6_tests(const std::string &path, const std::string &algorithm) {
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
        std::set<NetworKit::edgeweight> T;
        std::cout << line << " " << std::flush;
        switch (ALGORITHM[algorithm]) {
        case 0:
            T.insert(run_algorithm<Koala::KruskalMinimumSpanningTree>(G));
            T.insert(run_algorithm<Koala::PrimMinimumSpanningTree>(G));
            T.insert(run_algorithm<Koala::BoruvkaMinimumSpanningTree>(G));
            for (int i = 0; i < 5; i++) {
                T.insert(run_algorithm<Koala::KargerKleinTarjanMinimumSpanningTree>(G));
            }
            assert(T.size() == 1);
            break;
        case 1:
            run_algorithm<Koala::KruskalMinimumSpanningTree>(G);
            break;
        case 2:
            run_algorithm<Koala::PrimMinimumSpanningTree>(G);
            break;
        case 3:
            run_algorithm<Koala::BoruvkaMinimumSpanningTree>(G);
            break;
        }
        std::cout << std::endl;
    }
}

void run_dimacs_tests(const std::string &path, const std::string &algorithm) {
    auto G_directed = Koala::DimacsGraphReader().read(path);
    auto G = NetworKit::Graph(G_directed.numberOfNodes(), true, false);
    G_directed.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        if (!G.hasEdge(u, v) && !G.hasEdge(v, u) && w > 0) {
            G.addEdge(u, v, w);
        }
    });
    std::cout << path << " " << std::flush;
    std::set<NetworKit::edgeweight> T;
    T.insert(run_algorithm<Koala::KruskalMinimumSpanningTree>(G));
    T.insert(run_algorithm<Koala::PrimMinimumSpanningTree>(G));
    T.insert(run_algorithm<Koala::BoruvkaMinimumSpanningTree>(G));
    for (int i = 0; i < 5; i++) {
        T.insert(run_algorithm<Koala::KargerKleinTarjanMinimumSpanningTree>(G));
    }
    assert(T.size() == 1);
    std::cout << std::endl;
    return;
}

int main(int argc, const char *argv[]) {
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
