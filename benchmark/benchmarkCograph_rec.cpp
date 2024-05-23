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
    std::ifstream fin("../../input/cographConnected10.g6");
    //std::ifstream fin("../../input/graph9c.g6");
    int kol=0;
    while (true) {
        std::string line;
        fin >> line;
        kol++;
        if (!fin.good()) {
           break;
        }
        //line="Ds_";

        NetworKit::Graph G = Koala::G6GraphReader().readline(line);

        auto recognize = Koala::CographRecognition(G);

        //std::cout<<"line="<<line<<std::endl;

        recognize.run();

        if (recognize.getState() == Koala::CographRecognition::State::COGRAPH) {
            auto IndependentSet = Koala::CographIndependentSet((recognize.cotree));
            IndependentSet.run();

            assert(("Wrong independent set", IndependentSet.independet_set_size ==
                                             IndependentSet.bruteForceIndependetSetSize(G)));

            auto Clique = Koala::MaxClique((recognize.cotree));
            Clique.run();

            assert(("Wrong max clique size", Clique.size == Clique.bruteForceCliqueSize(G)));

            auto Coloring = Koala::CographVertexColoring(G);
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
