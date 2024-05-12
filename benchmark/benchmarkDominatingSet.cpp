#include <cassert>
#include <iostream>
#include <map>

#include <dominating_set/ExactDominatingSet.hpp>
#include <io/G6GraphReader.hpp>
#include <set_cover/BranchAndReduceSetCover.hpp>

template <typename T>
int run_algorithm(NetworKit::Graph &G, bool print = false) {
    auto algorithm = T(G);
    algorithm.run();
    auto &dominating_set = algorithm.getDominatingSet();
    std::cout << dominating_set.size() << " " << std::flush;
    algorithm.check();
    return dominating_set.size();
}

std::map<std::string, int> ALGORITHM = {
    { "exact", 0 },
    { "FKW", 1 }, { "Schiermeyer", 2 }, { "Grandoni", 3 }, { "FGK", 4 }, { "Rooij", 5 }
};

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <algorithm>" << std::endl;
        return 1;
    }
    while (true) {
        std::string line;
        std::cin >> line;
        if (!std::cin.good()) {
            break;
        }
        NetworKit::Graph G = Koala::G6GraphReader().readline(line);
        std::set<int> D;
        std::cout << line << " " << std::flush;
        switch (ALGORITHM[std::string(argv[1])]) {
        case 0:
            D.insert(run_algorithm<Koala::FominKratschWoegingerDominatingSet>(G));
            D.insert(run_algorithm<Koala::SchiermeyerDominatingSet>(G));
            D.insert(run_algorithm<
                Koala::BranchAndReduceDominatingSet<Koala::GrandoniSetCover>>(G));
            D.insert(run_algorithm<
                Koala::BranchAndReduceDominatingSet<Koala::FominGrandoniKratschSetCover>>(G));
            D.insert(run_algorithm<
                Koala::BranchAndReduceDominatingSet<Koala::RooijBodlaenderSetCover>>(G));
            assert(D.size() == 1);
            break;
        case 1:
            run_algorithm<Koala::FominKratschWoegingerDominatingSet>(G);
            break;
        case 2:
            run_algorithm<Koala::SchiermeyerDominatingSet>(G);
            break;
        case 3:
            run_algorithm<Koala::BranchAndReduceDominatingSet<Koala::GrandoniSetCover>>(G);
            break;
        case 4:
            run_algorithm<
                Koala::BranchAndReduceDominatingSet<Koala::FominGrandoniKratschSetCover>>(G);
            break;
        case 5:
            run_algorithm<Koala::BranchAndReduceDominatingSet<Koala::RooijBodlaenderSetCover>>(G);
            break;
        }
        std::cout << std::endl;
    }
    return 0;
}
