#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <string>

#include <io/DimacsGraphReader.hpp>
#include <io/G6GraphReader.hpp>

#include <independent_set/IndependentSet.hpp>
#include <independent_set/CographIndependentSet.hpp>
#include <recognition/CographRecognition.hpp>

template<typename T, typename... Args>
int run_algorithm(NetworKit::Graph &G, Args&&... args) {
    auto algorithm = T(G, std::forward<Args>(args)...);
    algorithm.run();
    auto independent_set = algorithm.getIndependentSet();
    std::cout << independent_set.size() << " " << std::flush;
    return independent_set.size();
}

template<typename T>
requires std::is_same_v<T, Koala::CographIndependentSet>
int run_algorithm(NetworKit::Graph &G) {
    auto recognition = Koala::HabibPaulCographRecognition(G);
    recognition.run();
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
        std::cout << line << " " << std::flush;
        switch (ALGORITHM[algorithm]) {
        case Algorithm::Exact: {
            std::set<int> I;
            double epsilon = 1 / (G.numberOfNodes() + 1);
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
        }
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
        case 20:
            run_algorithm<Koala::CographIndependentSet>(G);
            break;
        default:
            throw std::logic_error("Unrecognized algorithm");
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
    switch (ALGORITHM[algorithm]) {
    case 0: {
        I.insert(run_algorithm<Koala::Mis1IndependentSet>(G));
        I.insert(run_algorithm<Koala::Mis2IndependentSet>(G));
        I.insert(run_algorithm<Koala::Mis3IndependentSet>(G));
        I.insert(run_algorithm<Koala::Mis4IndependentSet>(G));
        I.insert(run_algorithm<Koala::Mis5IndependentSet>(G));
        I.insert(run_algorithm<Koala::MeasureAndConquerIndependentSet>(G));
        I.insert(run_algorithm<Koala::BakerPlanarGraphIndependentSet>(G, epsilon));
        I.insert(run_algorithm<Koala::BodlaenderPlanarGraphIndependentSet>(G, epsilon));
        assert(I.size() == 1);
        break;
    }
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
        run_algorithm<Koala::BakerPlanarGraphIndependentSet>(G);
        break;
    case 11:
        run_algorithm<Koala::BodlaenderPlanarGraphIndependentSet>(G);
        break;
    case 20:
        run_algorithm<Koala::CographIndependentSet>(G);
        break;
    default:
        throw std::logic_error("Unrecognized algorithm");
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
    if (path.ends_with(".g6")) {
        run_g6_tests(path, std::string(argv[1]));
    } else if (path.ends_with(".gr")) {
        run_dimacs_tests(path, std::string(argv[1]));
    } else {
        std::cerr << "File type not supported: " << path << std::endl;
    }
    return 0;
}
