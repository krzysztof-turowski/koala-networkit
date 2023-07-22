#include <cassert>
#include <iostream>
#include <map>

#include <dominating_set/ExactDominatingSet.hpp>
#include <io/G6GraphReader.hpp>
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
        std::cout << line << " " << G.numberOfEdges() << " ";

        std::set<int> D;
        switch (std::stoi(argv[1])) {
        case 0:
            D.insert(run_algorithm<Koala::FominKratschWoegingerDominatingSet>(G));
            D.insert(run_algorithm<Koala::SchiermeyerDominatingSet>(G));
            D.insert(run_algorithm<
                Koala::BranchAndReduceDominatingSet<Koala::RooijBodlaenderSetCover>>(G));
            D.insert(run_algorithm<
                Koala::BranchAndReduceDominatingSet<Koala::FominGrandoniKratschSetCover>>(G));
            D.insert(run_algorithm<
                Koala::BranchAndReduceDominatingSet<Koala::GrandoniSetCover>>(G));
            break;
        case 1:
            D.insert(run_algorithm<Koala::FominKratschWoegingerDominatingSet>(G));
            break;
        case 2:
            D.insert(run_algorithm<Koala::SchiermeyerDominatingSet>(G));
            break;
        case 3:
            D.insert(run_algorithm<
                Koala::BranchAndReduceDominatingSet<Koala::RooijBodlaenderSetCover>>(G));
            break;
        case 4:
            D.insert(run_algorithm<
                Koala::BranchAndReduceDominatingSet<Koala::FominGrandoniKratschSetCover>>(G));
            break;
        case 5:
            D.insert(run_algorithm<
                Koala::BranchAndReduceDominatingSet<Koala::GrandoniSetCover>>(G));
            break;
        }
        std::cout << std::endl;
        assert(D.size() == 1);
    }
    return 0;
}
