#include <cassert>
#include <iostream>
#include <map>

#include <independent_set/PlanarGraphIndependentSet.hpp>
#include <io/G6GraphReader.hpp>
#include <recognition/PlanarGraphRecognition.hpp>

int main() {
    std::map<bool, int> classification;
    std::string types[] = {
        "NON_PLANAR",
        "PLANAR"
    };

    while (true) {
        std::string line;
        std::cin >> line;
        if (!std::cin.good()) {
            break;
        }
        NetworKit::Graph G = Koala::G6GraphReader().readline(line);
        auto recognize = Koala::PlanarGraphRecognition(G);
        recognize.run();
        classification[recognize.isPlanar()]++;
        if (recognize.isPlanar()) {
            auto independent_set = Koala::BakerPlanarGraphIndependentSet(G, 1.0);
            independent_set.run();
        }
    }
    for (const auto &[k, v] : classification) {
        std::cout << types[k] << ": " << v << std::endl;
    }
    return 0;
}
