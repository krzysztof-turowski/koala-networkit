#include <cassert>
#include <iostream>
#include <map>

#include <coloring/EnumerationVertexColoring.hpp>
#include <io/G6GraphReader.hpp>

template <typename T> void benchmark() {
    std::map<int, int> chromatic_number;

    int i = 0;

    while (true) {
        std::string line;
        std::cin >> line;
        if (!std::cin.good()) {
            break;
        }
        NetworKit::Graph G = Koala::G6GraphReader().readline(line);
        auto algorithm = T(G);
        algorithm.run();
        auto colors = algorithm.getColoring();
        int max_color = 0;
        G.forEdges([&](NetworKit::node u, NetworKit::node v) { assert(colors[u] != colors[v]); });
        for (const auto& [v, c] : colors) {
            max_color = std::max(max_color, c);
        }
        chromatic_number[i++] = max_color;
    }
    for (const auto& [k, v] : chromatic_number) {
        std::cout << "GRAPH " << k << ": " << v << std::endl;
    }
}

int main(int argc, const char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <algorithm>" << std::endl;
        return 1;
    }
    switch (std::stoi(argv[1])) {
    case 0: benchmark<Koala::BrownsOrdinaryEnumerationVertexColoring>(); break;
    case 1: benchmark<Koala::ChristofidesEnumerationVertexColoring>(); break;
    case 2: benchmark<Koala::BrelazEnumerationVertexColoring>(); break;
    case 3: benchmark<Koala::KormanEnumerationVertexColoring>(); break;
    }
    return 0;
}
