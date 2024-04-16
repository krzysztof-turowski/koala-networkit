#include <cassert>
#include <iostream>
#include <map>

#include <coloring/PerfectGraphVertexColoring.hpp>
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
        std::ifstream fin("../../input/");
        std::string line;
        std::cin >> line;
        if (!std::cin.good()) {
            break;
        }
        NetworKit::Graph G = Koala::G6GraphReader().readline(line);
        auto recognize = Koala::PerfectGraphRecognition(G);
        recognize.run();
        classification[recognize.getState()]++;
        recognize.check();
        if (recognize.getState() == Koala::PerfectGraphRecognition::State::PERFECT) {
            auto color = Koala::PerfectGraphVertexColoring(G);
            color.run();
            color.check();
        }
    }
    for (const auto &[k, v] : classification) {
        std::cout << types[static_cast<int>(k)] << ": " << v << std::endl;
    }
    return 0;
}
