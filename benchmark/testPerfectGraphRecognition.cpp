#include <cassert>
#include <iostream>

#include <io/G6GraphReader.hpp>
#include <recognition/PerfectGraphRecognition.hpp>

int main() {
    int total = 0, perfect = 0;
    while (true) {
        std::string line;
        std::cin >> line;
        if (!std::cin.good()) {
            break;
        }
        NetworKit::Graph G = Koala::G6GraphReader().readline(line);
        auto recognize = Koala::PerfectGraphRecognition(G);
        recognize.run();
        if (recognize.isPerfect()) {
            perfect++;
        }
        total++;
    }
    std::cout << perfect << "/" << total << " graphs in this suite are perfect" << std::endl;
    return 0;
}
