#include <iostream>
#include <map>
#include <cassert>

#include <io/G6GraphReader.hpp>
#include "recognition/CographAlg.hpp"
#include "coloring/CographVertexColoring.hpp"
#include "independent_set/CographIndependentSet.hpp"
#include "max_clique/MaxClique.hpp"

int main() {
    std::map<Koala::CographRecognition::State, int> classification;
    std::string types[] = {
            "UNKNOWN",
            "COGRAPH",
            "NOT_COGRAPH"
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

        if (recognize.GetState() == Koala::CographRecognition::State::COGRAPH) {
            auto IndependentSet = Koala::CographIndependentSet(*(recognize.cotree), G.numberOfNodes());
            IndependentSet.run();

            assert(("Wrong independent set", IndependentSet.independet_set_size ==
                                             IndependentSet.BruteForceIndependetSetSize(G)));

            auto Clique = Koala::MaxClique(*(recognize.cotree), G.numberOfNodes());
            Clique.run();

            assert(("Wrong max clique size", Clique.size == Clique.BruteForceCliqueSize(G)));

            auto Coloring = Koala::CographVertexColoring(G);
            Coloring.run();

            assert(("Wrong coloring", Coloring.CheckColoring()));
        }
        classification[recognize.GetState()]++;
    }

    for (const auto &[k, v]: classification) {
        std::cout << types[static_cast<int>(k)] << ": " << v << std::endl;
    }
    return 0;
}
