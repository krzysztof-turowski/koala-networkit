#include <cassert>
#include <exception>
#include <iostream>
#include <map>
#include <set>
#include <string>

#include <benchmark/utils.hpp>
#include <matching/MaximumMatching.hpp>
#include <matching/gaussian_matching/BipartiteGaussianMatching.hpp>
#include <matching/gaussian_matching/GeneralGaussianMatching.hpp>
#include <matching/gaussian_matching/NaiveGaussianMatching.hpp>

enum class Algorithm : uint32_t {
    ALL, MICALI_VAZIRANI, NAIVE_GAUSSIAN, GAUSSIAN, EDMONDS, GABOW
};

std::map<std::string, Algorithm> ALGORITHM = {
    { "all", Algorithm::ALL },
    { "MicaliVazirani", Algorithm::MICALI_VAZIRANI },
    { "NaiveGaussian", Algorithm::NAIVE_GAUSSIAN },
    { "Gaussian", Algorithm::GAUSSIAN },
    { "Edmonds", Algorithm::EDMONDS },
    { "Gabow", Algorithm::GABOW },
};

template <typename Algorithm>
NetworKit::count run_algorithm(NetworKit::Graph &G) {
    auto algorithm = Algorithm(G);
    algorithm.run();
    if constexpr (std::same_as<Algorithm, Koala::MicaliVaziraniMatching>) {
        auto matching = algorithm.getMatching();
        NetworKit::count size = std::count_if(
            matching.begin(), matching.end(),
            [](auto pair) { return pair.second != NetworKit::none; });
        std::cout << size / 2 << " " << std::flush;
        return size / 2;
    }
    if constexpr (std::same_as<Algorithm, Koala::GeneralGaussianMatching>) {
        std::cout << algorithm.getMatching().size() << " " << std::flush;
        return algorithm.getMatching().size();
    }
    std::cout << algorithm.getMatching().size() / 2 << " " << std::flush;
    return algorithm.getMatching().size() / 2;
}

void run_test(NetworKit::Graph &G, Algorithm algorithm) {
    std::set<NetworKit::count> T;
    switch (algorithm) {
    case Algorithm::ALL:
        T.insert(run_algorithm<Koala::MicaliVaziraniMatching>(G));
        T.insert(run_algorithm<Koala::GeneralGaussianMatching>(G));
        T.insert(run_algorithm<Koala::NaiveGaussianMatching>(G));
        T.insert(run_algorithm<Koala::EdmondsMaximumMatching>(G));
        T.insert(run_algorithm<Koala::GabowMaximumMatching>(G));
        assert(T.size() == 1);
        break;
    case Algorithm::MICALI_VAZIRANI:
        T.insert(run_algorithm<Koala::MicaliVaziraniMatching>(G));
        break;
    case Algorithm::NAIVE_GAUSSIAN:
        T.insert(run_algorithm<Koala::NaiveGaussianMatching>(G));
        break;
    case Algorithm::GAUSSIAN:
        T.insert(run_algorithm<Koala::GeneralGaussianMatching>(G));
        break;
    case Algorithm::EDMONDS:
        T.insert(run_algorithm<Koala::EdmondsMaximumMatching>(G));
        break;
    case Algorithm::GABOW:
        T.insert(run_algorithm<Koala::GabowMaximumMatching>(G));
        break;
    }
    std::cout << std::endl;
}

int main(int argc, const char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <algorithm> <file>" << std::endl;
        return 1;
    }

    std::string algorithm_name(argv[1]);
    auto algorithm = ALGORITHM.find(algorithm_name);
    if (algorithm == ALGORITHM.end()) {
        throw std::invalid_argument("Unknown algorithm: " + algorithm_name);
    }
    Koala::Benchmark::executeForEachGraph(
        argv[2], [&](const std::string &label, NetworKit::Graph G) {
        G.indexEdges(true);
        std::cout << label << " " << std::flush;
        run_test(G, algorithm->second);
    });
    return 0;
}
