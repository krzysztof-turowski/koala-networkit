#include <iostream>
#include <map>
#include <set>
#include <string>

#include <benchmark/utils.hpp>
#include "recognition/CographRecognition.hpp"

std::string types[] = {
    "UNKNOWN",
    "COGRAPH",
    "NOT_COGRAPH",
    "CONTAINS_0_NODE",
    "EXISTS_1_NODE_NOT_PROPERLY_MARKED",
    "GRANDPARENT_IS_NOT_IN_SET",
    "NO_ONE_PATH",
    "WRONG_PARENT",
    "WRONG_GRANDPARENT"
};

template <typename T>
int run_algorithm(NetworKit::Graph &G, bool verbose = false) {
    auto algorithm = T(G);
    algorithm.run();
    if (verbose) {
        std::cout << static_cast<int>(algorithm.getState()) << " " << std::flush;
    }
    if (algorithm.isCograph()) {
        algorithm.check();
    }
    return static_cast<int>(algorithm.getState());
}

enum class Algorithm : uint32_t { ALL, BCHP, CSP, DAHLHAUS, HABIB_PAUL };

std::map<std::string, Algorithm> ALGORITHM = {
    { "all", Algorithm::ALL },
    { "BCHP", Algorithm::BCHP },
    { "CSP", Algorithm::CSP },
    { "Dahlhaus", Algorithm::DAHLHAUS },
    { "HabibPaul", Algorithm::HABIB_PAUL }
};

int choose_algorithm(NetworKit::Graph &G, Algorithm algorithm) {
    switch (algorithm) {
    case Algorithm::ALL: {
        std::set<int> states;
        states.insert(run_algorithm<Koala::BretscherCorneilHabibPaulCographRecognition>(G, true));
        states.insert(
            run_algorithm<Koala::CorneilStewartPerlCographRecognition>(G, true) != 1 ? 2 : 1);
        states.insert(run_algorithm<Koala::DahlhausCographRecognition>(G, true));
        states.insert(run_algorithm<Koala::HabibPaulCographRecognition>(G, true));
        assert(states.size() == 1);
        return *states.begin();
    }
    case Algorithm::BCHP:
        return run_algorithm<Koala::BretscherCorneilHabibPaulCographRecognition>(G);
    case Algorithm::CSP:
        return run_algorithm<Koala::CorneilStewartPerlCographRecognition>(G);
    case Algorithm::DAHLHAUS:
        return run_algorithm<Koala::DahlhausCographRecognition>(G);
    case Algorithm::HABIB_PAUL:
        return run_algorithm<Koala::HabibPaulCographRecognition>(G);
    }
    throw std::logic_error("Unhandled algorithm");
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
