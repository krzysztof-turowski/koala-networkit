#include <sys/stat.h>
#include <sys/types.h>

#include <cassert>
#include <fstream>
#include <io/G6GraphReader.hpp>
#include <iostream>
#include <map>
#include <recognition/planar/PlanarGraphRecognition.hpp>

int main() {
    std::map<int, int> classification;
    std::string types[] = {"NOT_PLANAR", "PLANAR"};
    while (true) {
        std::string line;
        std::cin >> line;
        if (!std::cin.good()) {
            break;
        }

        NetworKit::Graph G = Koala::G6GraphReader().readline(line);

       //auto recognize = Koala::HopcroftTarjan(G, false);
       // recognize.run();
        auto recognize = Koala::BoyerMyrvold(G, false);
        recognize.run();
        classification[recognize.isPlanar() == Koala::PlanarGraphRecognition::State::PLANAR ? 1 : 0]++;
    }
    for (const auto &[k, v] : classification) {
        std::cout << types[static_cast<int>(k)] << ": " << v << std::endl;
    }
    return 0;
}
