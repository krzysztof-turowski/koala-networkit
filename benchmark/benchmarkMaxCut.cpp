#include <cassert>
#include <chrono>
#include <iostream>
#include <map>

#include <benchmark/utils.hpp>

#include <max_cut/BranchAndBoundMaxCut.hpp>
#include <max_cut/GoemansWilliamsonMaxCut.hpp>
#include <max_cut/NaiveMaxCut.hpp>
#include <max_cut/RankTwoRelaxationMaxCut.hpp>

template <typename T>
std::pair<double, double> run_algorithm(NetworKit::Graph &G) {
    auto start = std::chrono::high_resolution_clock::now();

    auto algorithm = T(G);
    algorithm.run();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    double cutValue = algorithm.getMaxCutValue();
    return {duration.count(), cutValue};
}

enum class Algorithm : uint32_t {
    NAIVE, BRANCH_AND_BOUND, RANK_TWO_RELAXATION, GOEMANS_WILLIAMSON
};

std::map<std::string, Algorithm> ALGORITHM = {
    { "naive", Algorithm::NAIVE },
    { "branchAndBound", Algorithm::BRANCH_AND_BOUND },
    { "rankTwoRelaxation", Algorithm::RANK_TWO_RELAXATION },
    { "goemansWilliamson", Algorithm::GOEMANS_WILLIAMSON }
};

void process_graph(NetworKit::Graph &G, const std::string &algorithm_name, Algorithm algorithm) {
    std::pair<double, double> result;

    switch (algorithm) {
    case Algorithm::NAIVE:
        result = run_algorithm<Koala::NaiveMaxCut>(G);
        break;
    case Algorithm::BRANCH_AND_BOUND:
        result = run_algorithm<Koala::BranchAndBoundMaxCut>(G);
        break;
    case Algorithm::RANK_TWO_RELAXATION:
        result = run_algorithm<Koala::RankTwoRelaxationMaxCut>(G);
        break;
    case Algorithm::GOEMANS_WILLIAMSON:
        result = run_algorithm<Koala::GoemansWilliamsonMaxCut>(G);
        break;
    }

    double duration = result.first;
    double cutValue = result.second;

    std::cout << "Algorithm: " << algorithm_name << " "
              << "Time: " << duration << "s "
              << "Cut_Value: " << cutValue << std::endl;
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
        std::cout << label << " " << std::flush;
        process_graph(G, algorithm_name, algorithm->second);
    });
    return 0;
}
