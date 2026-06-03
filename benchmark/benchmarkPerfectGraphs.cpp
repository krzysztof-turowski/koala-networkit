#include <cassert>
#include <iostream>
#include <map>
#include <string>

#include <benchmark/utils.hpp>
#include <coloring/PerfectGraphVertexColoring.hpp>
#include <recognition/PerfectGraphRecognition.hpp>

int main(int argc, const char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <file>" << std::endl;
        return 1;
    }
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

    Koala::Benchmark::executeForEachGraph(argv[1], [&](const std::string &, NetworKit::Graph G) {
        auto recognize = Koala::PerfectGraphRecognition(G);
        recognize.run();
        classification[recognize.getState()]++;
        recognize.check();
        if (recognize.getState() == Koala::PerfectGraphRecognition::State::PERFECT) {
            auto color = Koala::PerfectGraphVertexColoring(G);
            color.run();
            color.check();
        }
    });
    for (const auto &[k, v] : classification) {
        std::cout << types[static_cast<int>(k)] << ": " << v << std::endl;
    }
    return 0;
}
