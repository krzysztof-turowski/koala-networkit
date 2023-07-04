#include <cassert>
#include <iostream>
#include <map>

#include <dominating_set/BranchAndReduceMDS.hpp>
#include <dominating_set/SchiermeyerMDS.hpp>
#include <dominating_set/FominKratschWoeginger.hpp>
#include <dominating_set/ExhaustiveMDS.hpp>
#include <io/G6GraphReader.hpp>
#include <set_cover/BranchAndReduceMSC.hpp>

template <typename T>
int run_algorithm(NetworKit::Graph &G) {
    auto algorithm = T(G);
    algorithm.run();
    auto dominating_set = algorithm.getDominatingSet();
    assert(algorithm.isDominating(dominating_set));
    int size = std::count(dominating_set.begin(), dominating_set.end(), true);
    std::cout << size << " ";
    return size;
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
            D.insert(run_algorithm<FominKratschWoegingerMDS>(G));
            D.insert(run_algorithm<SchiermeyerMDS>(G));
            D.insert(run_algorithm<BranchAndReduceMDS<RooijBodlaenderMSC>>(G));
            D.insert(run_algorithm<BranchAndReduceMDS<FominGrandoniKratschMSC>>(G));
            D.insert(run_algorithm<BranchAndReduceMDS<GrandoniMSC>>(G));
            break;
        case 1:
            D.insert(run_algorithm<FominKratschWoegingerMDS>(G));
            break;
        case 2:
            D.insert(run_algorithm<SchiermeyerMDS>(G));
            break;
        case 3:
            D.insert(run_algorithm<BranchAndReduceMDS<RooijBodlaenderMSC>>(G));
            break;
        case 4:
            D.insert(run_algorithm<BranchAndReduceMDS<FominGrandoniKratschMSC>>(G));
            break;
        case 5:
            D.insert(run_algorithm<BranchAndReduceMDS<GrandoniMSC>>(G));
            break;
        }
        assert(D.size() == 1);

        std::cout << std::endl;
    }
    return 0;
}
