#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <iomanip>

#include <networkit/graph/GraphTools.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <io/G6GraphReader.hpp>
#include <io/DimacsGraphReader.hpp>
#include <mst/MinimumSpanningTree.hpp>

int x = 0xFF; 
int D = 200; 
float eps = 0.1; 

template <typename T>
std::pair<NetworKit::edgeweight, NetworKit::Graph> run_algorithm(NetworKit::Graph &G) {
    auto algorithm = T(G);

    auto start = std::chrono::high_resolution_clock::now();
    algorithm.run();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    auto &spanning_tree = algorithm.getForest();
    std::cout << std::fixed << std::setprecision(0) << duration.count() << ' ' << spanning_tree.totalEdgeWeight() << '\n';
    // std::cout << "BRUUUH 0" << std::endl;
    // algorithm.check();
    // std::cout << "BRUUUH 1" << std::endl;
    return {spanning_tree.totalEdgeWeight(), spanning_tree};
}

template<>
std::pair<NetworKit::edgeweight, NetworKit::Graph> run_algorithm<Koala::ChazelleRubinfeldTrevisanMinimumSpanningTree>(NetworKit::Graph &G) {
    auto algorithm = Koala::ChazelleRubinfeldTrevisanMinimumSpanningTree(G);

    int max_w = x;
    // G.forEdges([&max_w](NetworKit::node, NetworKit::node, NetworKit::edgeweight ew, NetworKit::edgeid) {
    //     max_w = std::max(max_w, static_cast<int>(ew));
    // });

    auto start = std::chrono::high_resolution_clock::now();
    algorithm.run(max_w, eps);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << duration.count() << ' ' << algorithm.getTreeWeight() << '\n';
    return {algorithm.getTreeWeight(), NetworKit::Graph()};
}
namespace {
    enum class Algorithm {
        EXACT = 0,
        KRUSKAL,
        PRIM,
        BORUVKA,
        KKT,
        CRT,
        CHAZELLE,
    };
}

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
        auto G = NetworKit::Graph(G_directed.numberOfNodes(), true, false);
        G_directed.forEdges([&](NetworKit::node u, NetworKit::node v) {
            if (!G.hasEdge(u, v) && !G.hasEdge(v, u)) {
                G.addEdge(u, v);
            }
        });
        std::set<NetworKit::edgeweight> T;
        // std::cout << line << " " << std::flush;
        switch (ALGORITHM[algorithm]) {
        case Algorithm::EXACT:
            T.insert(run_algorithm<Koala::KruskalMinimumSpanningTree>(G).first);
            T.insert(run_algorithm<Koala::PrimMinimumSpanningTree>(G).first);
            T.insert(run_algorithm<Koala::BoruvkaMinimumSpanningTree>(G).first);
            for (int i = 0; i < 5; i++) {
                T.insert(run_algorithm<Koala::KargerKleinTarjanMinimumSpanningTree>(G).first);
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
    auto G = NetworKit::Graph(G_directed.numberOfNodes(), true, false);
    // int count1 = 0;
    // int countMax = 0;
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
        if (!G.hasEdge(u, v) && !G.hasEdge(v, u) && w > 0) {
            // int new_w = std::min((int)w / D + 1, x);
            // if (new_w == 1) count1 += 1;
            // if (new_w == x) countMax += 1;
            G.addEdge(u, v, w);
        }
    });

    // Ensure that the graph is connected...
    auto connected_components = NetworKit::ConnectedComponents(G);
    connected_components.run();
    auto components = connected_components.getComponents();
    for (auto i = 1; i < connected_components.numberOfComponents(); i++) {
        G.addEdge(
            components[0][0], components[i][0], D);
    }

    // std::cout << "Num 1 " << count1 << '\n'; 
    // std::cout << "Num X " << countMax << '\n'; 
    // std::cout << path << " " << std::flush;
    std::set<NetworKit::edgeweight> T;
    // std::cout << "BEFORE KRUSKAL" << std::endl;
    // T.insert(run_algorithm<Koala::KruskalMinimumSpanningTree>(G));
    // std::cout << "Run Kruskal" << std::endl;
    // std::cout << "Run Kruskal" << std::endl;
    // T.insert(run_algorithm<Koala::PrimMinimumSpanningTree>(G));
    // std::cout << "RUN BORUVKA" << std::endl;
    
    
    // std::cout << path << " " << std::flush;
    // std::cout << "time weight\n";
    // for (int xx : std::vector{0xFF, 0x7F, 0xF}) {
    //     x = xx;
    //     auto GNow = NetworKit::Graph(G_directed.numberOfNodes(), true, false);
    //     int count1 = 0;
    //     int countMax = 0;
    //     G_directed.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
    //         if (!GNow.hasEdge(u, v) && !GNow.hasEdge(v, u) && w > 0) {
    //             int new_w = std::min((int)w / D + 1, x);
    //             if (new_w == 1) count1 += 1;
    //             if (new_w == x) countMax += 1;
    //             GNow.addEdge(u, v, new_w);
    //         }
    //     });
    //     std::cout << "Num 1 " << count1 << '\n'; 
    //     std::cout << "Num X " << countMax << '\n'; 
        
    //     for (float ee : std::vector{.01, .05, 0.1}) {
    //         eps = ee;
    //         std :: cout << '\n';
    //         std :: cout << "x: " << x << ", eps: " << eps << '\n';
    //         for (int i = 0; i < 5; ++i) {
    //             T.insert(run_algorithm<Koala::BoruvkaMinimumSpanningTree>(GNow));
    //             // T.insert(run_algorithm<Koala::ChazelleRubinfeldTrevisanMinimumSpanningTree>(GNow));
    //         }
    //     }
    // }
    // std::cout << "Original graph" << '\n';
    // for (int i = 0; i < 5; i++) {
    //     T.insert(run_algorithm<Koala::BoruvkaMinimumSpanningTree>(G));
    //     // T.insert(run_algorithm<Koala::ChazelleRubinfeldTrevisanMinimumSpanningTree>(G));
    // }

    
    std::cout<< "Running chazelle now!" << std::endl;
    auto [chazelle, chazelleMST] = run_algorithm<Koala::Chazelle2000MinimumSpanningTree>(G);
    std::cout << "Finished Chazelle" << std::endl;
    auto [boruvka, boruvkaMST] = run_algorithm<Koala::BoruvkaMinimumSpanningTree>(G);
    std::cout << "Finished Boruvka" << std::endl;
    T.insert(chazelle);
    T.insert(boruvka);

    using NetworKit::node;
    using NetworKit::edgeweight;
    std::set<std::tuple<node, node, edgeweight>> boruvkaEdges;
    boruvkaMST.forEdges([&](node u, node v, edgeweight ew){
        auto [uu, vv] = std::minmax(u, v);
        boruvkaEdges.insert({uu, vv, ew});
    });


    NetworKit::Graph g(3, true);
    g.addEdge(1, 2, 2);
    NetworKit::Graph g2 = NetworKit::GraphTools::subgraphFromNodes(g, {1, 2});
    g2.forNodes([](node u){std::cout << u << ' ';}); std::cout <<std::endl;
    g2.forEdges([](node u, node v){std::cout << u << '-' << v << ' ';}); std::cout <<std::endl;

    // boruvkaEdges.forEdges([&](node u, node v, edgeweight ew){
    //     auto [uu, vv] = std::minmax(u, v);
    //     assert(chazelleMST.contains({uu, vv, ew}));
    // });
    // chazelleMST.forEdges([&](node u, node v, edgeweight ew){
    //     auto [uu, vv] = std::minmax(u, v);
    //     assert(boruvkaEdges.contains({uu, vv, ew}));
    // });

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
