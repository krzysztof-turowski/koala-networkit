#include <iostream>
#include <map>
#include <cassert>

#include <io/G6GraphReader.hpp>
#include "recognition/CographAlg.hpp"
#include "coloring/CographVertexColoring.hpp"
#include "independent_set/CographIndependentSet.hpp"
#include "max_clique/CographMaxClique.hpp"

int main() {
    std::map<Koala::CographRecognition::State, int> classification;
    std::string types[] = {
            "UNKNOWN",
            "COGRAPH",
            "NOT_COGRAPH"
    };
    std::ifstream fin("../../input/cographConnected9.g6");
    //std::ifstream fin("../../input/graph9c.g6");
    int kol=0;
    while (true) {
        std::string line;
        fin >> line;
        kol++;
        if (!fin.good()) {
           break;
        }

        NetworKit::Graph G = Koala::G6GraphReader().readline(line);

        auto recognize = Koala::CographRecognition(G);

        recognize.run();

        if (recognize.getState() == Koala::CographRecognition::State::COGRAPH) {
            Koala::CographIndependentSet IndependentSet =Koala::CographIndependentSet(G,recognize.cotree);
            IndependentSet.run();
            IndependentSet.check();
            assert(("Wrong independent set", IndependentSet.getIndependentSet().size() ==IndependentSet.bruteForceIndependetSetSize(G)));

            auto Clique = Koala::CographMaxClique(G,recognize.cotree);

            Clique.run();
            Clique.check();
            assert(("Wrong max clique size", Clique.getMaxCliqueSet().size()==Clique.bruteForceCliqueSize(G)));

            auto Coloring = Koala::CographVertexColoring(G,recognize.cotree);
            Coloring.run();

            assert(("Wrong coloring", Coloring.checkColoring()));
        }
        classification[recognize.getState()]++;
    }

        for (const auto &[k, v]: classification) {
            std::cout << types[static_cast<int>(k)] << ": " << v << std::endl;
        }
        return 0;
}
