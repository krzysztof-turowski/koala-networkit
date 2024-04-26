#include <cassert>
#include <iostream>
#include <map>
#include "vector"

#include <io/G6GraphReader.hpp>
#include "cograph_rec/CographAlg.hpp"

int main() {
    std::map<Koala::CographRecognition::State, int> classification;
    std::string types[] = {
        "UNKNOWN",
        "COGRAPH",
        "IS_NOT_COGRAPH"
    };
    std::ifstream fin("../../input/cographConnected11.g6");
    //std::ifstream fin("../../input/graph9c.g6");

    //cographConnected
    //grapc
    int kol=0;

    while (true)
    {
        std::string line;
        fin >> line;

        if (!fin.good())
        {
            break;
        }
        NetworKit::Graph G = Koala::G6GraphReader().readline(line);

        auto recognize = Koala::CographRecognition(G);
        kol++;
        recognize.run();


        //std::cout<<recognize.cotree->GetMaxIndependetSetSize()<<std::endl;
        if(recognize.GetState()==Koala::CographRecognition::State::COGRAPH && recognize.cotree->GetMaxIndependetSetSize()!=recognize.cotree->BruteForceIndependetSetSize())
        {
            std::cout<<"error1"<<std::endl;
        }
        if(recognize.GetState()==Koala::CographRecognition::State::COGRAPH && recognize.cotree->GetMaxCliqueSize()!=recognize.cotree->BruteForceCliqueSize())
        {
            std::cout<<"error2"<<std::endl;
        }

        classification[recognize.GetState()]++;

    }

    for (const auto &[k, v] : classification) {
        std::cout << types[static_cast<int>(k)] << ": " << v << std::endl;
    }
    return 0;
}
