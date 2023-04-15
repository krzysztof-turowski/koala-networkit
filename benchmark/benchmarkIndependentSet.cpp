#include <cassert>
#include <iostream>
#include <map>

#include <io/G6GraphReader.hpp>
#include <independent_set/SimpleIndependentSet.hpp>

template <typename T>
void benchmark() {
    std::map<int, int> classification;

    while (true) {
        std::string line;
        std::cin >> line;
        if (!std::cin.good()) {
            break;
        }
        NetworKit::Graph G = Koala::G6GraphReader().readline(line);
        auto recognize = T(G);
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
}

int main(int argc, const char *argv[]) {
    switch (std::stoi(argv[1])) {
        case 0:
            benchmark<Koala::BruteForceIndependentSet>();
            break;
        case 1:
            benchmark<Koala::Mis1IndependentSet>();
            break;
        case 2:
            benchmark<Koala::Mis2IndependentSet>();
            break;
        case 3:
            benchmark<Koala::Mis3IndependentSet>();
            break;
        case 4:
            benchmark<Koala::Mis4IndependentSet>();
            break;
        case 5:
            benchmark<Koala::Mis5IndependentSet>();
            break;
        case 6:
            benchmark<Koala::MeasureAndConquerIndependentSet>();
            break;
    }
    return 0;
}
