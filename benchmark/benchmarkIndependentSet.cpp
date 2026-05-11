#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>

#include <io/DimacsGraphReader.hpp>
#include <io/G6GraphReader.hpp>

#include <independent_set/IndependentSet.hpp>
#include <independent_set/ExactIndependentSet.hpp>
#include <independent_set/PlanarGraphIndependentSet.hpp>
#include <independent_set/CographIndependentSet.hpp>
#include <recognition/CographRecognition.hpp>

template<typename T, typename... Args>
std::size_t run_algorithm(NetworKit::Graph &G, Args&&... args) {
    auto algorithm = T(G, std::forward<Args>(args)...);
    algorithm.run();
    auto independent_set = algorithm.getIndependentSet();
    std::cout << independent_set.size() << " " << std::flush;
    return independent_set.size();
}

template<typename T>
requires std::is_same_v<T, Koala::CographIndependentSet>
std::size_t run_algorithm(NetworKit::Graph &G) {
    auto recognition = Koala::HabibPaulCographRecognition(G);
    recognition.run();
    if (!recognition.isCograph()) {
        throw std::logic_error("Graph is not a cograph");
    }
    auto algorithm = T(G, recognition.cotree);
    algorithm.run();
    auto independent_set = algorithm.getIndependentSet();
    std::cout << independent_set.size() << " " << std::flush;
    return independent_set.size();
}

enum class Algorithm {
    Exact,
    MIS1, MIS2, MIS3, MIS4, MIS5, MeasureAndConquer,
    Baker, Bodlaender,
    Cograph
};

const std::map<std::string, Algorithm> ALGORITHM = {
    { "exact", Algorithm::Exact },
    { "MIS1", Algorithm::MIS1 }, { "MIS2", Algorithm::MIS2 }, { "MIS3", Algorithm::MIS3 },
    { "MIS4", Algorithm::MIS4 }, { "MIS5", Algorithm::MIS5 },
    { "M&C", Algorithm::MeasureAndConquer },
    { "Baker", Algorithm::Baker }, { "Bodlaender", Algorithm::Bodlaender },
    { "cograph", Algorithm::Cograph}
};

std::size_t run_exact_algorithms(NetworKit::Graph &G) {
    std::set<std::size_t> I;
    I.insert(run_algorithm<Koala::Mis1IndependentSet>(G));
    I.insert(run_algorithm<Koala::Mis2IndependentSet>(G));
    I.insert(run_algorithm<Koala::Mis3IndependentSet>(G));
    I.insert(run_algorithm<Koala::Mis4IndependentSet>(G));
    I.insert(run_algorithm<Koala::Mis5IndependentSet>(G));
    I.insert(run_algorithm<Koala::MeasureAndConquerIndependentSet>(G));
    assert(I.size() == 1);
    return *I.begin();
}

std::size_t run(NetworKit::Graph &G, Algorithm algorithm) {
    switch (algorithm) {
    case Algorithm::Exact:
        return run_exact_algorithms(G);
    case Algorithm::MIS1:
        return run_algorithm<Koala::Mis1IndependentSet>(G);
    case Algorithm::MIS2:
        return run_algorithm<Koala::Mis2IndependentSet>(G);
    case Algorithm::MIS3:
        return run_algorithm<Koala::Mis3IndependentSet>(G);
    case Algorithm::MIS4:
        return run_algorithm<Koala::Mis4IndependentSet>(G);
    case Algorithm::MIS5:
        return run_algorithm<Koala::Mis5IndependentSet>(G);
    case Algorithm::MeasureAndConquer:
        return run_algorithm<Koala::MeasureAndConquerIndependentSet>(G);
    case Algorithm::Baker: {
        double epsilon = 1.0 / (G.numberOfNodes() + 1);
        return run_algorithm<Koala::BakerPlanarGraphIndependentSet>(G, epsilon);
    }
    case Algorithm::Bodlaender: {
        double epsilon = 1.0 / (G.numberOfNodes() + 1);
        return run_algorithm<Koala::BodlaenderPlanarGraphIndependentSet>(G, epsilon);
    }
    case Algorithm::Cograph:
        return run_algorithm<Koala::CographIndependentSet>(G);
    default:
        throw std::logic_error("Algorithm not implemented");
    }
}

void run_g6_tests(const std::string &path, Algorithm algorithm) {
    std::fstream file(path, std::fstream::in);
    std::map<std::size_t, std::size_t> classification;
    while (true) {
        std::string line;
        file >> line;
        if (!file.good()) {
            break;
        }
        NetworKit::Graph G = Koala::G6GraphReader().readline(line);
        std::cout << line << " " << std::flush;
        ++classification[run(G, algorithm)];
        std::cout << std::endl;
    }
    std::cout << "List of graphs counted by solution size:" << std::endl;
    for (const auto &[k, v] : classification) {
        std::cout << k << ": " << v << std::endl;
    }
}

void run_dimacs_tests(const std::string &path, Algorithm algorithm) {
    auto G_directed = Koala::DimacsGraphReader().read(path);
    NetworKit::Graph G = NetworKit::Graph(G_directed.numberOfNodes(), true, false);
    G_directed.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        if (!G.hasEdge(u, v) && !G.hasEdge(v, u) && w > 0) {
            G.addEdge(u, v);
        }
    });
    std::cout << path << " " << std::flush;
    run(G, algorithm);
    std::cout << std::endl;
    return;
}

int main(int argc, const char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <algorithm> <file>" << std::endl;
        return 1;
    }
    std::string algorithm(argv[1]), path(argv[2]);
    if (!ALGORITHM.count(algorithm)) {
      std::cerr << "Algorithm not supported: " << path << std::endl;
    } else if (path.ends_with(".g6")) {
        run_g6_tests(path, ALGORITHM[algorithm]);
    } else if (path.ends_with(".gr")) {
        run_dimacs_tests(path, ALGORITHM[algorithm]);
    } else {
        std::cerr << "File type not supported: " << path << std::endl;
    }
    return 0;
}
