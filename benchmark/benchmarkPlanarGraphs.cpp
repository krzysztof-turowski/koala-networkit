#include <iostream>
#include <map>
#include <io/G6GraphReader.hpp>
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

        Koala::LeftRightPlanarity tester;
        bool result = tester.isPlanar(G1);

        auto recognize = Koala::HopcroftTarjan(G2, false);
        recognize.run();

        classification[result]++;
    
        if( (recognize.isPlanar() == Koala::PlanarGraphRecognition::State::PLANAR) != result) {
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
        }
       
    }
    for (const auto &[k, v] : classification) {
        std::cout << types[static_cast<int>(k)] << ": " << v << std::endl;
    }
    return 0;
}
