#include <cassert>
#include <iostream>
#include <map>

#include <io/G6GraphReader.hpp>
#include <io/DimacsGraphReader.hpp>
#include <independent_set/ExactIndependentSet.hpp>
#include <independent_set/PlanarGraphIndependentSet.hpp>

template <typename T>
int run_algorithm(NetworKit::Graph &G) {
    auto algorithm = T(G);
    algorithm.run();
    auto &independent_set = algorithm.getIndependentSet();
    std::cout << independent_set.size() << " " << std::flush;
    return independent_set.size();
}

template <typename T>
int run_algorithm(NetworKit::Graph &G, double epsilon) {
    auto algorithm = T(G, epsilon);
    algorithm.run();
    auto &independent_set = algorithm.getIndependentSet();
    std::cout << independent_set.size() << " " << std::flush;
    return independent_set.size();
}

const std::map<std::string, int> ALGORITHM = {
    { "exact", 0 },
    { "MIS1", 1 }, { "MIS2", 2 }, { "MIS3", 3 }, { "MIS4", 4 }, { "MIS5", 5 },
    { "MeasureAndConquer", 6 },
    { "baker", 10 }, { "bodlaender", 11 }
};

const double epsilon = 0.5;

void run_g6_tests(const std::string &path, const std::string &algorithm) {
    std::fstream file(path, std::fstream::in);
    std::map<int, int> classification;
    while (true) {
        std::string line;
        file >> line;
        if (!file.good()) {
            break;
        }
        NetworKit::Graph G = Koala::G6GraphReader().readline(line);
        std::set<int> I;
        std::cout << line << " " << std::flush;
        switch (ALGORITHM[algorithm]) {
        case 0:
            I.insert(run_algorithm<Koala::Mis1IndependentSet>(G));
            I.insert(run_algorithm<Koala::Mis2IndependentSet>(G));
            I.insert(run_algorithm<Koala::Mis3IndependentSet>(G));
            I.insert(run_algorithm<Koala::Mis4IndependentSet>(G));
            I.insert(run_algorithm<Koala::Mis5IndependentSet>(G));
            I.insert(run_algorithm<Koala::MeasureAndConquerIndependentSet>(G));
            I.insert(run_algorithm<Koala::BakerPlanarGraphIndependentSet>(G, epsilon));
            I.insert(run_algorithm<Koala::BodlaenderPlanarGraphIndependentSet>(G, epsilon));
            assert(I.size() == 1);
            classification[*I.begin()]++;
            break;
        case 1:
            run_algorithm<Koala::Mis1IndependentSet>(G);
            break;
        case 2:
            run_algorithm<Koala::Mis2IndependentSet>(G);
            break;
        case 3:
            run_algorithm<Koala::Mis3IndependentSet>(G);
            break;
        case 4:
            run_algorithm<Koala::Mis4IndependentSet>(G);
            break;
        case 5:
            run_algorithm<Koala::Mis5IndependentSet>(G);
            break;
        case 6:
            run_algorithm<Koala::MeasureAndConquerIndependentSet>(G);
            break;
        case 10:
            run_algorithm<Koala::BakerPlanarGraphIndependentSet>(G, epsilon);
            break;
        case 11:
            run_algorithm<Koala::BodlaenderPlanarGraphIndependentSet>(G, epsilon);
            break;
        }
        std::cout << std::endl;
    }
    std::cout << "List of graphs counted by solution size:" << std::endl;
    for (const auto &[k, v] : classification) {
        std::cout << k << ": " << v << std::endl;
    }
}

void run_dimacs_tests(const std::string &path, const std::string &algorithm) {
    auto G_directed = Koala::DimacsGraphReader().read(path);
    NetworKit::Graph G = NetworKit::Graph(G_directed.numberOfNodes(), true, false);
    G_directed.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        if (!G.hasEdge(u, v) && !G.hasEdge(v, u) && w > 0) {
            G.addEdge(u, v);
        }
    });
    std::cout << path << " " << std::flush;
    std::set<int> I;
    // I.insert(run_algorithm<Koala::Mis1IndependentSet>(G));
    // I.insert(run_algorithm<Koala::Mis2IndependentSet>(G));
    // I.insert(run_algorithm<Koala::Mis3IndependentSet>(G));
    // I.insert(run_algorithm<Koala::Mis4IndependentSet>(G));
    // I.insert(run_algorithm<Koala::Mis5IndependentSet>(G));
    // I.insert(run_algorithm<Koala::MeasureAndConquerIndependentSet>(G));
    I.insert(run_algorithm<Koala::BakerPlanarGraphIndependentSet>(G, epsilon));
    I.insert(run_algorithm<Koala::BodlaenderPlanarGraphIndependentSet>(G, epsilon));
    // assert(I.size() == 1);
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
