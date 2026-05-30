#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <string>

#include <benchmark/utils.hpp>
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

enum class Algorithm : uint32_t { EXACT, FKW, SCHIERMEYER, GRANDONI, FGK, ROOIJ };

std::map<std::string, Algorithm> ALGORITHM = {
    { "exact", Algorithm::EXACT },
    { "FKW", Algorithm::FKW },
    { "Schiermeyer", Algorithm::SCHIERMEYER },
    { "Grandoni", Algorithm::GRANDONI },
    { "FGK", Algorithm::FGK },
    { "Rooij", Algorithm::ROOIJ }
};

void choose_algorithm(NetworKit::Graph &G, Algorithm algorithm) {
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
    }
}

int main(int argc, char **argv) {
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
        choose_algorithm(G, algorithm->second);
        std::cout << std::endl;
    });
    return 0;
}
