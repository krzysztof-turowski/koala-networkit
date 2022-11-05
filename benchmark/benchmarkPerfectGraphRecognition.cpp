#include <cassert>
#include <iostream>
#include <map>

#include <io/G6GraphReader.hpp>
#include <recognition/PerfectGraphRecognition.hpp>

int main() {
    std::map<Koala::PerfectGraphRecognition::State, int> classification;
    std::string types[] = {
        "UNKNOWN",
        "PERFECT",
        "HAS_JEWEL",
        "HAS_PYRAMID",
        "HAS_T1",
        "HAS_T2",
        "HAS_T3",
        "HAS_NEAR_CLEANER_ODD_HOLE"
    };

    while (true) {
        std::string line;
        std::cin >> line;
        if (!std::cin.good()) {
            break;
        }
        NetworKit::Graph G = Koala::G6GraphReader().readline(line);
        auto recognize = Koala::PerfectGraphRecognition(G);
        recognize.run();
        classification[recognize.getState()]++;
        switch (recognize.getState()) {
            case Koala::PerfectGraphRecognition::State::HAS_T1:
              std::cout << "T1: " << line << std::endl;
              break;
            case Koala::PerfectGraphRecognition::State::HAS_T2:
              std::cout << "T2: " << line << std::endl;
              break;
            case Koala::PerfectGraphRecognition::State::HAS_T3:
              std::cout << "T3: " << line << std::endl;
              break;
        }
        recognize.check();
    }
    for (const auto &[k, v] : classification) {
        std::cout << types[static_cast<int>(k)] << ": " << v << std::endl;
    }
    return 0;
}
