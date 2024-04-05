//
// Created by milana on 19.03.24.
//
#include <iostream>
#include <map>
#include "io/G6GraphReader.hpp"
#include "recognition/CorneilStewartPerlCographRecognition.h"

int main() {
    std::map<Koala::CorneilStewartPerlCographRecognition::State, int> classification;
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
    while (true) {
        std::string line;
        std::cin >> line;
        if (!std::cin.good()) {
            break;
        }
        NetworKit::Graph G = Koala::G6GraphReader().readline(line);
        auto recognize = Koala::CorneilStewartPerlCographRecognition(G);
        recognize.run();
        classification[recognize.getState()]++;
    }
    for (const auto &[k, v]: classification) {
        std::cout << types[static_cast<int>(k)] << ": " << v << std::endl;
    }
    return 0;
}
