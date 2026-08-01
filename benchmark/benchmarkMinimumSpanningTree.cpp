#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>

#include <networkit/components/ConnectedComponents.hpp>

#include <graph/GraphTools.hpp>
#include <io/DimacsGraphReader.hpp>
#include <io/G6GraphReader.hpp>
#include <mst/MinimumSpanningTree.hpp>

template <typename T>
NetworKit::edgeweight run_algorithm(NetworKit::Graph &G) {
    auto algorithm = T(G);
    algorithm.run();
    auto &spanning_tree = algorithm.getForest();
    std::cout << spanning_tree.totalEdgeWeight() << " " << std::flush;
    return spanning_tree.totalEdgeWeight();
}

template <typename T>
NetworKit::edgeweight run_algorithm(NetworKit::Graph &G, float eps) {
    auto algorithm = T(G);
    int max_weight = 0xFF;
    algorithm.run(max_weight, eps);
    std::cout << algorithm.getTreeWeight() << " " << std::flush;
    return algorithm.getTreeWeight();
}

enum class Algorithm : uint32_t {
    EXACT = 1,
    KRUSKAL,
    PRIM,
    BORUVKA,
    KKT,
    CRT,
    CHAZELLE = 100,
};

std::map<std::string, Algorithm> ALGORITHM = {
    { "exact", Algorithm::EXACT},
    { "Kruskal", Algorithm::KRUSKAL },
    { "Prim", Algorithm::PRIM },
    { "Boruvka", Algorithm::BORUVKA },
    { "KKT", Algorithm::KKT },
    { "Chazelle", Algorithm::CHAZELLE },
    { "CRT", Algorithm::CRT },
};

double EPS = 0.1;

NetworKit::edgeweight run_exact(NetworKit::Graph &G) {
    std::set<NetworKit::edgeweight> T;
    T.insert(run_algorithm<Koala::KruskalMinimumSpanningTree>(G));
    T.insert(run_algorithm<Koala::PrimMinimumSpanningTree>(G));
    T.insert(run_algorithm<Koala::BoruvkaMinimumSpanningTree>(G));
    for (int i = 0; i < 5; i++) {
        T.insert(run_algorithm<Koala::KargerKleinTarjanMinimumSpanningTree>(G));
    }
    T.insert(run_algorithm<Koala::Chazelle2000MinimumSpanningTree>(G));
    assert(T.size() == 1);

    if (Koala::GraphTools::hasDistinctIntegerWeights(G)) {
        for (int i = 0; i < 5; i++) {
            auto result = run_algorithm<Koala::ChazelleRubinfeldTrevisanMinimumSpanningTree>(
                G, EPS);
            assert((1 - EPS) * (*T.begin()) <= result);
            assert((1 + EPS) * (*T.begin()) >= result);
        }
    }
    return *T.begin();
}

NetworKit::edgeweight choose_algorithm(NetworKit::Graph &G, const std::string &algorithm_name) {
    switch (ALGORITHM[algorithm_name]) {
    case Algorithm::EXACT:
        return run_exact(G);
    case Algorithm::KRUSKAL:
        return run_algorithm<Koala::KruskalMinimumSpanningTree>(G);
    case Algorithm::PRIM:
        return run_algorithm<Koala::PrimMinimumSpanningTree>(G);
    case Algorithm::BORUVKA:
        return run_algorithm<Koala::BoruvkaMinimumSpanningTree>(G);
    case Algorithm::KKT:
        return run_algorithm<Koala::KargerKleinTarjanMinimumSpanningTree>(G);
    case Algorithm::CHAZELLE:
        return run_algorithm<Koala::Chazelle2000MinimumSpanningTree>(G);
    case Algorithm::CRT: {
        if (Koala::GraphTools::hasDistinctIntegerWeights(G)) {
            return run_algorithm<Koala::ChazelleRubinfeldTrevisanMinimumSpanningTree>(G, EPS);
        } else {
            throw std::invalid_argument(
                "Algorithm " + algorithm_name + " requires distinct edge weights");
        }
    }
    default:
        throw std::logic_error("Unhandled algorithm");
    }
}

void run_g6_tests(const std::string &path, const std::string &algorithm_name) {
    std::fstream file(path, std::fstream::in);
    while (true) {
        std::string line;
        file >> line;
        if (!file.good()) {
            break;
        }
        auto G = Koala::G6GraphReader().readline(line);
        G = Koala::GraphTools::convertDirectedGraphToUndirected(G, true);
        std::cout << line << " " << std::flush;
        choose_algorithm(G, algorithm_name);
        std::cout << std::endl;
    }
}

void run_dimacs_tests(const std::string &path, const std::string &algorithm_name) {
    auto G = Koala::DimacsGraphReader().read(path);
    G = Koala::GraphTools::convertDirectedGraphToUndirected(G, true);

    bool is_connected = Koala::GraphTools::isConnected(G);
    if (!is_connected) {
        throw std::invalid_argument("Graph is not connected");
    }

    std::cout << path << " " << std::flush;
    choose_algorithm(G, algorithm_name);
    std::cout << std::endl;
}

int main(int argc, const char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <algorithm> <file>" << std::endl;
        return 1;
    }
    std::string algorithm_name(argv[1]);
    std::string path(argv[2]);

    auto algorithm = ALGORITHM.find(algorithm_name);
    if (algorithm == ALGORITHM.end()) {
        throw std::invalid_argument("Unknown algorithm: " + algorithm_name);
    }

    auto position = path.find_last_of(".");
    if (path.substr(position + 1) == "g6") {
        run_g6_tests(path, algorithm_name);
    } else if (path.substr(position + 1) == "gr") {
        run_dimacs_tests(path, algorithm_name);
    } else {
        throw std::invalid_argument("File type not supported: " + path);
    }
    return 0;
}
