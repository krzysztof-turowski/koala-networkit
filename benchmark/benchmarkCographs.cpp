//
// Created by milana on 19.03.24.
//
#include <io/G6GraphReader.hpp>

#include <iostream>
#include <map>
#include <string>

#include "recognition/CographRecognition.hpp"

std::string types[] = {
        "UNKNOWN",
        "COGRAPH",
        "CONTAINS_0_NODE",
        "EXISTS_1_NODE_NOT_PROPERLY_MARKED",
        "GRANDPARENT_IS_NOT_IN_SET",
        "NO_ONE_PATH",
        "WRONG_PARENT",
        "WRONG_GRANDPARENT"
};

std::string default_types[] = {
        "NOT_COGRAPH",
        "COGRAPH"
};

std::map<std::string, int> classification;

void run_algorithm(NetworKit::Graph &G, int number) {
    if (number == 0) {
        auto algorithm = Koala::CorneilStewartPerlCographRecognition(G);
        algorithm.run();
        classification[types[static_cast<int>(algorithm.getState())]]++;
    } else if (number == 1) {
        auto algorithm = Koala::BretscherCorneilHabibPaulCographRecognition(G);
        algorithm.run();
        classification[default_types[algorithm.isCograph()]]++;
    }
}

std::map<std::string, int> ALGORITHM = {
        { "decomposition", 0 },
        { "lexbfs", 1 }
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
        run_algorithm(G, std::stoi(argv[1]));
    }
    for (const auto &[k, v] : classification) {
        std::cout << k << ": " << v << std::endl;
    }
    return 0;
}
