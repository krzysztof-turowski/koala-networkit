#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <fstream>

#include <networkit/graph/GraphTools.hpp>
#include <io/G6GraphReader.hpp>
#include <io/DimacsGraphReader.hpp>

#include <max_cut/NaiveMaxCut.hpp>
#include <max_cut/BranchAndBoundMaxCut.hpp>
#include <max_cut/RankTwoRelaxationMaxCut.hpp>
#include <max_cut/GoemansWilliamsonMaxCut.hpp>

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

std::map<std::string, int> ALGORITHM = {
    { "naive", 0 },
    { "branchAndBound", 1 },
    { "rankTwoRelaxation", 2 },
    { "goemansWilliamson", 3 }
};

void process_graph(NetworKit::Graph &G, const std::string &algorithm) {
    std::pair<double, double> result;

    auto it = ALGORITHM.find(algorithm);
    if (it == ALGORITHM.end()) {
        std::cerr << "Unknown algorithm: " << algorithm << std::endl;
        return;
    }

    switch (it->second) {
    case 0:
        result = run_algorithm<Koala::NaiveMaxCut>(G);
        break;
    case 1:
        result = run_algorithm<Koala::BranchAndBoundMaxCut>(G);
        break;
    case 2:
        result = run_algorithm<Koala::RankTwoRelaxationMaxCut>(G);
        break;
    case 3:
        result = run_algorithm<Koala::GoemansWilliamsonMaxCut>(G);
        break;
    default:
        std::cerr << "Unknown algorithm: " << algorithm << std::endl;
        return;
    }

    double duration = result.first;
    double cutValue = result.second;

    std::cout << "Algorithm: " << algorithm << " "
              << "Time: " << duration << "s "
              << "Cut_Value: " << cutValue << std::endl;
}

void run_g6_tests(const std::string &algorithm) {
    while (true) {
        std::string line;
        std::cin >> line;
        if (!std::cin.good()) {
            break;
        }
        NetworKit::Graph G = Koala::G6GraphReader().read(line);
        std::cout << line << " " << std::flush;

        process_graph(G, algorithm);
    }
}

void run_dimacs_tests(const std::string &path, const std::string &algorithm) {
    auto G = Koala::DimacsGraphReader().read(path);

    std::cout << path << " " << std::flush;

    process_graph(G, algorithm);
}

int main(int argc, const char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <algorithm> <file>" << std::endl;
        return 1;
    }

    std::string algorithm(argv[1]);
    std::string path(argv[2]);
    auto position = path.find_last_of(".");
    std::string ext = path.substr(position + 1);

    if (ext == "g6") {
        run_g6_tests(algorithm);
    } else if (ext == "gr") {
        run_dimacs_tests(path, algorithm);
    } else {
        std::cerr << "File type not supported: " << path << std::endl;
        return 1;
    }

    return 0;
}
