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

int x = 0xFF;
int D = 200;
float eps = 0.1;

template <typename T>
NetworKit::edgeweight run_algorithm(NetworKit::Graph &G) {
    auto algorithm = T(G);
    algorithm.run();
    auto &spanning_tree = algorithm.getForest();
    std::cout << spanning_tree.totalEdgeWeight() << " " << std::flush;
    return spanning_tree.totalEdgeWeight();
}

template<>
NetworKit::edgeweight
    run_algorithm<Koala::ChazelleRubinfeldTrevisanMinimumSpanningTree>(NetworKit::Graph &G) {
    auto algorithm = Koala::ChazelleRubinfeldTrevisanMinimumSpanningTree(G);
    int max_w = x;
    algorithm.run(max_w, eps);
    std::cout << algorithm.getTreeWeight() << " " << std::flush;
    return algorithm.getTreeWeight();
}

namespace {
enum class Algorithm : uint32_t{
    EXACT = 0,
    KRUSKAL,
    PRIM,
    BORUVKA,
    KKT,
    CRT,
    CHAZELLE = 100,
};
}  // namespace

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
        case Algorithm::CRT:
            run_algorithm<Koala::ChazelleRubinfeldTrevisanMinimumSpanningTree>(G);
            break;
        case Algorithm::CHAZELLE:
            run_algorithm<Koala::Chazelle2000MinimumSpanningTree>(G);
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

    // Ensure that the graph is connected...
    auto connected_components = NetworKit::ConnectedComponents(G_distinct);
    connected_components.run();
    auto components = connected_components.getComponents();
    for (auto i = 1; i < connected_components.numberOfComponents(); i++) {
        G_distinct.addEdge(
            components[0][0], components[i][0], ++max_ew);
    }

    std::cout << path << " " << std::flush;
    std::set<NetworKit::edgeweight> T;
    T.insert(run_algorithm<Koala::KruskalMinimumSpanningTree>(G_distinct));
    T.insert(run_algorithm<Koala::PrimMinimumSpanningTree>(G_distinct));
    T.insert(run_algorithm<Koala::BoruvkaMinimumSpanningTree>(G_distinct));
    T.insert(run_algorithm<Koala::Chazelle2000MinimumSpanningTree>(G_distinct));

    for (int i = 0; i < 5; i++) {
        run_algorithm<Koala::ChazelleRubinfeldTrevisanMinimumSpanningTree>(G_distinct);
    }

    assert(T.size() == 1);

    for (int xx : std::vector{0xFF, 0x7F, 0x3F, 0x1F,  0xF}) {
        x = xx;
        auto GNow = NetworKit::Graph(G_directed.numberOfNodes(), true, false);
        int count1 = 0;
        int countMax = 0;
        G_directed.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
            if (!GNow.hasEdge(u, v) && !GNow.hasEdge(v, u) && w > 0) {
                int new_w = std::min(static_cast<int>(w) / D + 1, x);
                if (new_w == 1) count1 += 1;
                if (new_w == x) countMax += 1;
                GNow.addEdge(u, v, new_w);
            }
        });

        std::cout << "#x:" << x << std::endl;
        for (int i = 0; i < 5; ++i) {
            std::cout << "#Boruvka" << std::endl;
            run_algorithm<Koala::BoruvkaMinimumSpanningTree>(GNow);
        }
        for (float ee : std::vector{.01, 0.02, 0.03, 0.05, 0.1}) {
            eps = ee;
            std :: cout << std::setprecision(3) << "#eps:" << eps << '\n';
            for (int i = 0; i < 5; ++i) {
                std::cout << "#CRT" << std::endl;
                run_algorithm<Koala::ChazelleRubinfeldTrevisanMinimumSpanningTree>(GNow);
            }
        }
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
