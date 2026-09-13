#include <cassert>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>

#include <networkit/distance/APSP.hpp>

#include <benchmark/utils.hpp>
#include <shortest_path/Aingworth.hpp>
#include <shortest_path/Seidel.hpp>
#include <shortest_path/UnweightedShoshanZwick.hpp>
#include <shortest_path/Spira.hpp>

template <typename T>
std::pair<double, double> run_algorithm(NetworKit::Graph &G) {
    auto start = std::chrono::high_resolution_clock::now();

    auto algorithm = T(G);
    algorithm.run();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    return {duration.count(), static_cast<double>(algorithm.getDiameter())};
}

enum class Algorithm : uint32_t {
    ALL, SEIDEL, SHOSHAN_ZWICK, AINGWORTH, SPIRA
};

std::map<std::string, Algorithm> ALGORITHM = {
    { "all", Algorithm::ALL },
    { "Seidel", Algorithm::SEIDEL },
    { "UnweightedShoshanZwick", Algorithm::SHOSHAN_ZWICK },
    { "Aingworth", Algorithm::AINGWORTH },
    { "Spira", Algorithm::SPIRA }
};

void check_against_reference(NetworKit::Graph &G) {
    NetworKit::APSP reference(G);
    reference.run();
    Koala::SeidelAPSP seidel(G);
    seidel.run();
    Koala::UnweightedShoshanZwickAPSP shoshan_zwick(G);
    shoshan_zwick.run();
    Koala::AingworthAPSP aingworth(G);
    aingworth.run();
    for (auto u : G.nodeRange()) {
        for (auto v : G.nodeRange()) {
            int exact = static_cast<int>(reference.getDistance(u, v));
            assert(seidel.getDistance(u, v) == exact);
            assert(shoshan_zwick.getDistance(u, v) == exact);
            assert(aingworth.getDistance(u, v) >= exact
                && aingworth.getDistance(u, v) <= exact + 2);
        }
    }
}

void process_graph(NetworKit::Graph &G, const std::string &algorithm_name, Algorithm algorithm) {
    if (algorithm == Algorithm::ALL) {
        check_against_reference(G);
        return;
    }

    std::pair<double, double> result;
    switch (algorithm) {
    case Algorithm::SEIDEL:
        result = run_algorithm<Koala::SeidelAPSP>(G);
        break;
    case Algorithm::SHOSHAN_ZWICK:
        result = run_algorithm<Koala::UnweightedShoshanZwickAPSP>(G);
        break;
    case Algorithm::AINGWORTH:
        result = run_algorithm<Koala::AingworthAPSP>(G);
        break;
    case Algorithm::SPIRA:
        result = run_algorithm<Koala::SpiraAPSP>(G);
        break;
    default:
        throw std::logic_error("Unhandled algorithm");
    }

    std::cout << "Algorithm: " << algorithm_name << " "
              << "Time: " << result.first << "s "
              << "Diameter: " << result.second;
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
        std::cout << std::endl;
    });
    return 0;
}
