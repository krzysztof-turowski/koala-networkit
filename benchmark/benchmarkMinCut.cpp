#include <cassert>
#include <chrono>
#include <iostream>
#include <map>

#include <benchmark/utils.hpp>
#include <flow/MaximumFlow.hpp>
#include <flow/KingRaoTarjanMaximumFlow.hpp>
#include <min_cut/HaoOrlinMinCut.hpp>
#include <min_cut/KargerMinCut.hpp>
#include <min_cut/KargerSteinMinCut.hpp>
#include <min_cut/StoerWagnerMinCut.hpp>

template <typename T>
std::pair<double, double> run_algorithm(NetworKit::Graph &G) {
    auto start = std::chrono::high_resolution_clock::now();

    auto algorithm = T(G);
    algorithm.run();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    double cutValue = algorithm.getMinCutValue();
    return {duration.count(), cutValue};
}

enum class Algorithm : uint32_t { HAO_ORLIN, KARGER, KARGER_STEIN, STOER_WAGNER };

std::map<std::string, Algorithm> ALGORITHM = {
    { "haoOrlin", Algorithm::HAO_ORLIN },
    { "karger", Algorithm::KARGER },
    { "kargerStein", Algorithm::KARGER_STEIN },
    { "stoerWagner", Algorithm::STOER_WAGNER }
};

void process_graph(NetworKit::Graph &G, const std::string &algorithm_name, Algorithm algorithm) {
    std::pair<double, double> result;

    switch (algorithm) {
    case Algorithm::HAO_ORLIN:
        result = run_algorithm<Koala::HaoOrlinMinCut<Koala::KingRaoTarjanMaximumFlow>>(G);
        break;
    case Algorithm::KARGER:
        result = run_algorithm<Koala::KargerMinCut>(G);
        break;
    case Algorithm::KARGER_STEIN:
        result = run_algorithm<Koala::KargerSteinMinCut>(G);
        break;
    case Algorithm::STOER_WAGNER:
        result = run_algorithm<Koala::StoerWagnerMinCut>(G);
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
