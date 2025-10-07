#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <iomanip>
#include <algorithm>

#include <networkit/components/ConnectedComponents.hpp>

#include <graph/GraphTools.hpp>
#include <io/G6GraphReader.hpp>
#include <io/DimacsGraphReader.hpp>
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
    int max_w = 0xFF;
    algorithm.run(max_w, eps);
    std::cout << algorithm.getTreeWeight() << " " << std::flush;
    return algorithm.getTreeWeight();
}

enum class Algorithm : uint32_t{
    EXACT = 0,
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
    { "CRT", Algorithm::CRT },
    { "Chazelle", Algorithm::CHAZELLE },
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
        auto G = Koala::GraphTools::convertDirectedGraphToUndirected(G_directed, true);
        std::set<NetworKit::edgeweight> T;
        std::cout << line << " " << std::flush;
        switch (ALGORITHM[algorithm]) {
        case Algorithm::EXACT:
            T.insert(run_algorithm<Koala::KruskalMinimumSpanningTree>(G));
            T.insert(run_algorithm<Koala::PrimMinimumSpanningTree>(G));
            T.insert(run_algorithm<Koala::BoruvkaMinimumSpanningTree>(G));
            for (int i = 0; i < 5; i++) {
                T.insert(run_algorithm<Koala::KargerKleinTarjanMinimumSpanningTree>(G));
            }
            assert(T.size() == 1);
            break;
        case Algorithm::KRUSKAL:
            run_algorithm<Koala::KruskalMinimumSpanningTree>(G);
            break;
        case Algorithm::PRIM:
            run_algorithm<Koala::PrimMinimumSpanningTree>(G);
            break;
        case Algorithm::BORUVKA:
            run_algorithm<Koala::BoruvkaMinimumSpanningTree>(G);
            break;
        case Algorithm::CHAZELLE:
            run_algorithm<Koala::Chazelle2000MinimumSpanningTree>(G);
            break;
        case Algorithm::CRT:
            run_algorithm<Koala::ChazelleRubinfeldTrevisanMinimumSpanningTree>(G, 0.1);
        }
        std::cout << std::endl;
    }
}

void run_dimacs_tests(const std::string &path, const std::string &algorithm) {
    auto G_directed = Koala::DimacsGraphReader().read(path);
    auto G_distinct = NetworKit::Graph(G_directed.numberOfNodes(), true, false);
    double max_ew = 0;
    G_directed.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight ew){
        max_ew = std::max(max_ew, ew);
    });
    std::set<double> ews;
    G_directed.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        if (ews.contains(w)) {
            max_ew += 1;
            w = max_ew;
        }
        ews.insert(w);
        if (!G_distinct.hasEdge(u, v) && !G_distinct.hasEdge(v, u) && w > 0) {
            G_distinct.addEdge(u, v, w);
        }
    });

    // Ensure that the graph is connected
    auto connected_components = NetworKit::ConnectedComponents(G_distinct);
    connected_components.run();
    auto components = connected_components.getComponents();
    for (auto i = 1; i < connected_components.numberOfComponents(); i++) {
        G_distinct.addEdge(components[0][0], components[i][0], ++max_ew);
    }

    std::cout << path << " " << std::flush;
    std::set<NetworKit::edgeweight> T;
    T.insert(run_algorithm<Koala::KruskalMinimumSpanningTree>(G_distinct));
    T.insert(run_algorithm<Koala::PrimMinimumSpanningTree>(G_distinct));
    T.insert(run_algorithm<Koala::BoruvkaMinimumSpanningTree>(G_distinct));
    T.insert(run_algorithm<Koala::Chazelle2000MinimumSpanningTree>(G_distinct));
    for (int i = 0; i < 5; i++) {
        T.insert(run_algorithm<Koala::KargerKleinTarjanMinimumSpanningTree>(G_distinct));
    }
    assert(T.size() == 1);

    float eps = 0.1;
    for (int i = 0; i < 5; i++) {
        auto w = run_algorithm<Koala::ChazelleRubinfeldTrevisanMinimumSpanningTree>(
            G_distinct, eps);
        assert((1 - eps) * (*T.begin()) <= w);
        assert((1 + eps) * (*T.begin()) >= w);
    }
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
