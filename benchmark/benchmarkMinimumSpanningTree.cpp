#include <cassert>
#include <iostream>
#include <map>
#include <set>

#include <networkit/graph/GraphTools.hpp>

#include <io/G6GraphReader.hpp>
#include <mst/MinimumSpanningTree.hpp>

#include <chrono>
#include <iomanip>

template <typename T>
int run_algorithm(NetworKit::Graph &G, bool print = true) {
    auto start = std::chrono::high_resolution_clock::now();
    auto algorithm = T(G);
    algorithm.run();
    auto stop = std::chrono::high_resolution_clock::now();
    auto &spanning_tree = algorithm.getForest();
    if (print) {
      std::cout << std::endl;
      std::cout << std::fixed << std::setw(10) << std::setprecision(3) << spanning_tree.numberOfEdges() << " " << spanning_tree.totalEdgeWeight() << " " << std::flush;
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
      std::cout << std::fixed << std::setw(10) << std::setprecision(3) << duration / 1000.0 << " ";
    }
    return spanning_tree.totalEdgeWeight();
}

std::map<std::string, int> ALGORITHM = {
    { "all", 0 },
    { "Kruskal", 1 }, { "Prim", 2 }, { "Boruvka", 3 }, { "KKT", 4 }
};

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <algorithm>" << std::endl;
        return 1;
    }
    while (true) {
        std::string line;
        std::cin >> line;
        if (!std::cin.good()) {
            break;
        }
        NetworKit::Graph G = Koala::G6GraphReader().readline(line);
        NetworKit::GraphTools::randomizeWeights(G);
        std::set<int> T;
        std::cout << line << " " << std::flush;
        switch (ALGORITHM[std::string(argv[1])]) {
        case 0:
            T.insert(run_algorithm<Koala::KruskalMinimumSpanningTree>(G));
            T.insert(run_algorithm<Koala::PrimMinimumSpanningTree>(G));
            T.insert(run_algorithm<Koala::BoruvkaMinimumSpanningTree>(G));
            T.insert(run_algorithm<Koala::KargerKleinTarjanMinimumSpanningTree>(G));
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
        case 4:
            run_algorithm<Koala::KargerKleinTarjanMinimumSpanningTree>(G);
            break;
        }
        std::cout << std::endl;
    }
    return 0;
}
