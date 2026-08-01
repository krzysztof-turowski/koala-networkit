#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>

#include <benchmark/utils.hpp>
#include "recognition/ChordalGraphRecognition.hpp"

template <typename T>
int run_algorithm(NetworKit::Graph &G, bool verbose = false) {
    auto algorithm = T(G);
    algorithm.run();
    if (verbose) {
        std::cout << static_cast<int>(algorithm.getState()) << " " << std::flush;
    }
    return static_cast<int>(algorithm.getState());
}

enum class Algorithm : uint32_t { ALL, MCS, LEXBFS };

std::map<std::string, Algorithm> ALGORITHM = {
    { "all", Algorithm::ALL },
    { "MCS", Algorithm::MCS },
    { "LexBFS", Algorithm::LEXBFS }
};

int choose_algorithm(NetworKit::Graph &G, Algorithm algorithm) {
    switch (algorithm) {
    case Algorithm::ALL: {
        std::set<int> states;
        states.insert(run_algorithm<Koala::MaximumCardinalitySearchChordalGraphRecognition>(
            G, true));
        states.insert(run_algorithm<Koala::LexBFSChordalGraphRecognition>(G, true));
        assert(states.size() == 1);
        return *states.begin();
    }
    case Algorithm::MCS:
        return run_algorithm<Koala::MaximumCardinalitySearchChordalGraphRecognition>(G);
    case Algorithm::LEXBFS:
        return run_algorithm<Koala::LexBFSChordalGraphRecognition>(G);
    default:
        throw std::logic_error("Unhandled algorithm");
    }
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
    std::string types[] = { "UNKNOWN", "CHORDAL", "NOT_CHORDAL" };
    std::map<int, int> classification;
    Koala::Benchmark::executeForEachGraph(
        argv[2], [&](const std::string &label, NetworKit::Graph G) {
        std::cout << label << " " << std::flush;
        classification[choose_algorithm(G, algorithm->second)]++;
        std::cout << std::endl;
    });
    for (const auto &[state, count] : classification) {
        std::cout << types[state] << ": " << count << std::endl;
    }
    return 0;
}
