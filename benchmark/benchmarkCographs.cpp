//
// Created by milana on 19.03.24.
//
#include <cassert>
#include <iostream>
#include <map>

#include <io/G6GraphReader.hpp>
#include <cograph_recognition/CographRecognition.hpp>

int main() {
    std::map<Koala::CographRecognition::State, int> classification;
    std::string types[] = {
            "UNKNOWN",
            "COMPLEMENT_REDUCIBLE",
            "COND_1",
            "COND_2",
            "COND_3",
            "COND_4",
            "COND_5",
            "COND_6"
    };

    while (true) {
        std::string line;
        std::cin >> line;
        if (!std::cin.good()) {
            break;
        }
        NetworKit::Graph G = Koala::G6GraphReader().readline(line);
        auto recognize = Koala::CographRecognition(G);
        recognize.run();
        classification[recognize.getState()]++;
        recognize.check();
    }
    for (const auto &[k, v] : classification) {
        std::cout << types[static_cast<int>(k)] << ": " << v << std::endl;
    }
    return 0;
}
