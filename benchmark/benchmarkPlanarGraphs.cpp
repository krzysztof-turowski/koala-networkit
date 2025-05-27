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
    int c = 0;
    while (true) {
        std::string line;
        std::cin >> line;
        if (!std::cin.good()) {
            break;
        }

        NetworKit::Graph G1 = Koala::G6GraphReader().readline(line);
        NetworKit::Graph G2 = Koala::G6GraphReader().readline(line);

       /*
        NetworKit::Graph G(6, false);
        G.addEdge(0,3); G.addEdge(0,5); G.addEdge(0,4);
        G.addEdge(1,3); G.addEdge(1,5); G.addEdge(1,4);
        G.addEdge(2,3); G.addEdge(2,5); G.addEdge(2,4);
        */
        //Koala::PlanarityTester tester(G);
        //std::vector<std::vector<NetworKit::node>> embedding;

       

        Koala::LeftRightPlanarity tester;
        bool result = tester.isPlanar(G1);
        //std::cout << result;
        //if (result)
           // std::cout << "Graph is PLANAR.\n";
        //else
           // std::cout << "Graph is NON-PLANAR.\n";
    
     //   if recognize.isPlanar() == Koala::PlanarGraphRecognition::State::PLANAR && 
        auto recognize = Koala::HopcroftTarjan(G2, false);
        recognize.run();
        //auto recognize = Koala::BoyerMyrvold(G, false);
       // recognize.run();
         //classification[recognize.isPlanar() == Koala::PlanarGraphRecognition::State::PLANAR ? 1 : 0]++;
          classification[result]++;
    
        if( (recognize.isPlanar() == Koala::PlanarGraphRecognition::State::PLANAR) != result) {
           // std::cout << result << " " << (recognize.isPlanar() == Koala::PlanarGraphRecognition::State::PLANAR) << "\n\n\n";
           if(recognize.isPlanar() == Koala::PlanarGraphRecognition::State::PLANAR) {
            std::cout << "SHOULD BE PLANAR, BUT ALGO SAID NOT PLANAR\n";
           } else {
            std::cout << "SHOULD BE NOT PLANAR, BUT ALGO SAID PLANAR\n";
           }
           for (auto i : G1.nodeRange()){
                for (auto u : G1.neighborRange(i)) {
                    if (i < u)
                    std::cout << i << " " << u << '\n';
                }
            }
            std::cout << "\n\n";
           // break;
        }
        // break;
       
    }
    for (const auto &[k, v] : classification) {
        std::cout << types[static_cast<int>(k)] << ": " << v << std::endl;
    }
    return 0;
}
