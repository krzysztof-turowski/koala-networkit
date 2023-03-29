#include <cassert>
#include <iostream>
#include <map>

#include <io/G6GraphReader.hpp>
#include <independent_set/SimpleIndependentSet.hpp>

int main() {
    std::map<int, int> classification;

    while (true) {
        std::string line;
        std::cin >> line;
        if (!std::cin.good()) {
            break;
        }
        NetworKit::Graph G = Koala::G6GraphReader().readline(line);
        auto recognize = Koala::Mis1IndependentSet(G);
        recognize.run();
        int count = 0;
        for (const auto &[_, belongs] : recognize.getIndependentSet()) {
          count += belongs;
        }
        classification[count]++;
    }
    for (const auto &[k, v] : classification) {
        std::cout << "SIZE " << k << ": " << v << std::endl;
    }
    return 0;
}
