#include <iostream>
#include <map>

#include <io/G6GraphReader.hpp>

#include "recognition/ChordalGraphRecognition.hpp"

int main() {
    std::map<Koala::ChordalGraphRecognition::State, int> classification;
    std::string types[] = {
        "UNKNOWN",
        "CHORDAL",
        "NOT_CHORDAL"
    };

    while (true) {
        std::string line;
        std::cin >> line;
        if (!std::cin.good()) {
            break;
        }

        NetworKit::Graph G = Koala::G6GraphReader().readline(line);
        auto recognize = Koala::LexBFSChordalGraphRecognition(G);
        recognize.run();
        classification[recognize.getState()]++;
    }

    for (const auto &[k, v] : classification) {
        std::cout << types[static_cast<int>(k)] << ": " << v << std::endl;
    }
    return 0;
}
