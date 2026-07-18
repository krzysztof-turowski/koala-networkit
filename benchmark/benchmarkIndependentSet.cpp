#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <benchmark/utils.hpp>
#include <independent_set/ChordalIndependentSet.hpp>
#include <independent_set/IndependentSet.hpp>
#include <independent_set/BakerKOuterplanarGraphIndependentSet.hpp>
#include <independent_set/BakerPlanarGraphIndependentSet.hpp>
#include <independent_set/CographIndependentSet.hpp>
#include <recognition/ChordalGraphRecognition.hpp>
#include <recognition/CographRecognition.hpp>

template<typename T>
int run_algorithm(NetworKit::Graph &G) {
    std::vector<NetworKit::node> independent_set;
    if constexpr (std::is_same_v<T, Koala::CographIndependentSet>) {
        auto recognition = Koala::HabibPaulCographRecognition(G);
        recognition.run();
        auto algorithm = T(G, recognition.cotree);
        algorithm.run();
        independent_set = algorithm.getIndependentSet();
    } else if constexpr (std::is_same_v<T, Koala::ChordalIndependentSet>) {
        auto recognition = Koala::MaximumCardinalitySearchChordalGraphRecognition(G);
        recognition.run();
        auto algorithm = T(G, recognition.getPEO());
        algorithm.run();
        independent_set = algorithm.getIndependentSet();
    } else {
        auto algorithm = T(G);
        algorithm.run();
        independent_set = algorithm.getIndependentSet();
    }

    std::cout << independent_set.size() << " " << std::flush;
    return independent_set.size();
}

template <typename T>
int run_approximation_algorithm(NetworKit::Graph &graph, double epsilon) {
    auto algorithm = T(graph, epsilon);
    algorithm.run();
    algorithm.check();
    const auto &independent_set = algorithm.getIndependentSet();
    std::cout << independent_set.size() << " " << std::flush;
    return independent_set.size();
}

enum class Algorithm : uint32_t {
    EXACT, BRUTE_FORCE, MIS1, MIS2, MIS3, MIS4, MIS5,
    MEASURE_AND_CONQUER, COGRAPH, KOUTERPLANAR, PLANAR, CHORDAL
};

std::map<std::string, Algorithm> ALGORITHM = {
    { "exact", Algorithm::EXACT },
    { "bruteforce", Algorithm::BRUTE_FORCE },
    { "MIS1", Algorithm::MIS1 }, { "MIS2", Algorithm::MIS2 }, { "MIS3", Algorithm::MIS3 },
    { "MIS4", Algorithm::MIS4 }, { "MIS5", Algorithm::MIS5 },
    { "MeasureAndConquer", Algorithm::MEASURE_AND_CONQUER },
    { "cograph", Algorithm::COGRAPH },
    { "k-outerplanar", Algorithm::KOUTERPLANAR },
    { "planar", Algorithm::PLANAR },
    { "Baker", Algorithm::PLANAR },
    { "chordal", Algorithm::CHORDAL }
};

NetworKit::Graph normalize_graph(NetworKit::Graph &G_directed) {
    NetworKit::Graph G = NetworKit::Graph(G_directed.numberOfNodes(), true, false);
    G_directed.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        if (!G.hasEdge(u, v) && !G.hasEdge(v, u) && w > 0) {
            G.addEdge(u, v);
        }
    });
    return G;
}

void choose_algorithm(
        NetworKit::Graph &G, Algorithm algorithm, double epsilon) {
    std::set<int> independent_sets;
    switch (algorithm) {
    case Algorithm::EXACT:
        independent_sets.insert(run_algorithm<Koala::BruteForceIndependentSet>(G));
        independent_sets.insert(run_algorithm<Koala::Mis1IndependentSet>(G));
        independent_sets.insert(run_algorithm<Koala::Mis2IndependentSet>(G));
        independent_sets.insert(run_algorithm<Koala::Mis3IndependentSet>(G));
        independent_sets.insert(run_algorithm<Koala::Mis4IndependentSet>(G));
        independent_sets.insert(run_algorithm<Koala::Mis5IndependentSet>(G));
        independent_sets.insert(run_algorithm<Koala::MeasureAndConquerIndependentSet>(G));
        assert(independent_sets.size() == 1);
        break;
    case Algorithm::BRUTE_FORCE:
        run_algorithm<Koala::BruteForceIndependentSet>(G);
        break;
    case Algorithm::MIS1:
        run_algorithm<Koala::Mis1IndependentSet>(G);
        break;
    case Algorithm::MIS2:
        run_algorithm<Koala::Mis2IndependentSet>(G);
        break;
    case Algorithm::MIS3:
        run_algorithm<Koala::Mis3IndependentSet>(G);
        break;
    case Algorithm::MIS4:
        run_algorithm<Koala::Mis4IndependentSet>(G);
        break;
    case Algorithm::MIS5:
        run_algorithm<Koala::Mis5IndependentSet>(G);
        break;
    case Algorithm::MEASURE_AND_CONQUER:
        run_algorithm<Koala::MeasureAndConquerIndependentSet>(G);
        break;
    case Algorithm::COGRAPH:
        run_algorithm<Koala::CographIndependentSet>(G);
        break;
    case Algorithm::KOUTERPLANAR:
        run_algorithm<Koala::BakerKOuterplanarGraphIndependentSet>(G);
        break;
    case Algorithm::PLANAR:
        run_approximation_algorithm<
            Koala::BakerPlanarGraphIndependentSet>(G, epsilon);
        break;
    case Algorithm::CHORDAL:
        run_algorithm<Koala::ChordalIndependentSet>(G);
        break;
    default:
        throw std::invalid_argument("Unknown algorithm");
    }
}

int main(int argc, const char *argv[]) {
    if (argc < 3 || argc > 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <algorithm> <file> [epsilon]" << std::endl;
        return 1;
    }
    std::string algorithm_name(argv[1]);
    auto algorithm = ALGORITHM.find(algorithm_name);
    if (algorithm == ALGORITHM.end()) {
        throw std::invalid_argument("Unknown algorithm: " + algorithm_name);
    }
    const bool is_approximation = algorithm->second == Algorithm::PLANAR;
    if (argc != (is_approximation ? 4 : 3)) {
        throw std::invalid_argument(
            "The planar Baker algorithm requires an epsilon argument");
    }
    const double epsilon = is_approximation ? std::stod(argv[3]) : 0.0;
    Koala::Benchmark::executeForEachGraph(
        argv[2], [&](const std::string &label, NetworKit::Graph input_graph) {
        auto G = normalize_graph(input_graph);
        std::cout << label << " " << std::flush;
        choose_algorithm(G, algorithm->second, epsilon);
        std::cout << std::endl;
    });
    return 0;
}
