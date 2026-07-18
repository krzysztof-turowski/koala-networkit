#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>

#include <benchmark/utils.hpp>
#include <dominating_set/BakerKOuterplanarGraphDominatingSet.hpp>
#include <dominating_set/BakerPlanarGraphDominatingSet.hpp>
#include <dominating_set/ExactDominatingSet.hpp>
#include <set_cover/BranchAndReduceSetCover.hpp>

template <typename T>
int run_algorithm(NetworKit::Graph &G) {
    auto algorithm = T(G);
    algorithm.run();
    auto &dominating_set = algorithm.getDominatingSet();
    std::cout << dominating_set.size() << " " << std::flush;
    algorithm.check();
    return dominating_set.size();
}

template <typename T>
int run_approximation_algorithm(
        NetworKit::Graph &graph, double epsilon) {
    auto algorithm = T(graph, epsilon);
    algorithm.run();
    const auto &dominating_set = algorithm.getDominatingSet();
    std::cout << dominating_set.size() << " " << std::flush;
    algorithm.check();
    return dominating_set.size();
}

enum class Algorithm : uint32_t {
    EXACT, FKW, SCHIERMEYER, GRANDONI, FGK, ROOIJ, KOUTERPLANAR, PLANAR
};

std::map<std::string, Algorithm> ALGORITHM = {
    { "exact", Algorithm::EXACT },
    { "FKW", Algorithm::FKW },
    { "Schiermeyer", Algorithm::SCHIERMEYER },
    { "Grandoni", Algorithm::GRANDONI },
    { "FGK", Algorithm::FGK },
    { "Rooij", Algorithm::ROOIJ },
    { "k-outerplanar", Algorithm::KOUTERPLANAR },
    { "planar", Algorithm::PLANAR },
    { "Baker", Algorithm::PLANAR }
};

void choose_algorithm(
        NetworKit::Graph &G, Algorithm algorithm, double epsilon) {
    std::set<int> dominating_sets;
    switch (algorithm) {
    case Algorithm::EXACT:
        dominating_sets.insert(run_algorithm<Koala::FominKratschWoegingerDominatingSet>(G));
        dominating_sets.insert(run_algorithm<Koala::SchiermeyerDominatingSet>(G));
        dominating_sets.insert(
            run_algorithm<Koala::BranchAndReduceDominatingSet<Koala::GrandoniSetCover>>(G));
        dominating_sets.insert(run_algorithm<
            Koala::BranchAndReduceDominatingSet<Koala::FominGrandoniKratschSetCover>>(G));
        dominating_sets.insert(run_algorithm<
            Koala::BranchAndReduceDominatingSet<Koala::RooijBodlaenderSetCover>>(G));
        assert(dominating_sets.size() == 1);
        break;
    case Algorithm::FKW:
        run_algorithm<Koala::FominKratschWoegingerDominatingSet>(G);
        break;
    case Algorithm::SCHIERMEYER:
        run_algorithm<Koala::SchiermeyerDominatingSet>(G);
        break;
    case Algorithm::GRANDONI:
        run_algorithm<Koala::BranchAndReduceDominatingSet<Koala::GrandoniSetCover>>(G);
        break;
    case Algorithm::FGK:
        run_algorithm<Koala::BranchAndReduceDominatingSet<Koala::FominGrandoniKratschSetCover>>(G);
        break;
    case Algorithm::ROOIJ:
        run_algorithm<Koala::BranchAndReduceDominatingSet<Koala::RooijBodlaenderSetCover>>(G);
        break;
    case Algorithm::KOUTERPLANAR:
        run_algorithm<Koala::BakerKOuterplanarGraphDominatingSet>(G);
        break;
    case Algorithm::PLANAR:
        run_approximation_algorithm<
            Koala::BakerPlanarGraphDominatingSet>(G, epsilon);
        break;
    default:
        throw std::invalid_argument("Unknown algorithm");
    }
}

int main(int argc, char **argv) {
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
        argv[2], [&](const std::string &label, NetworKit::Graph G) {
        std::cout << label << " " << std::flush;
        choose_algorithm(G, algorithm->second, epsilon);
        std::cout << std::endl;
    });
    return 0;
}
