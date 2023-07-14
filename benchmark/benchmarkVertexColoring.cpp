#include <cassert>
#include <iostream>
#include <map>

#include <coloring/ExactVertexColoring.hpp>
#include <io/G6GraphReader.hpp>

template <typename T>
int run_algorithm(NetworKit::Graph &G) {
    auto algorithm = T(G);
    algorithm.run();
    auto colors = algorithm.getColoring();
    G.forEdges([&](NetworKit::node u, NetworKit::node v) { assert(colors[u] != colors[v]); });
    int max_color = 0;
    for (const auto& [v, c] : colors) {
        max_color = std::max(max_color, c);
    }
    std::cout << max_color << " ";
    return max_color;
}

int main() {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <algorithm>" << std::endl;
        return 1;
    }
    while (true) {
        std::string line;
        std::cin >> line;
        if (!std::cin.good()) {
            break;
        }
        NetworKit::Graph G = Koala::G6GraphReader().readline(line);
        std::cout << line << " " << G.numberOfEdges() << " ";

        std::set<int> C;
        switch (std::stoi(argv[1])) {
        case 0:
            C.insert(run_algorithm<Koala::BrownsOrdinaryEnumerationVertexColoring>(G));
            C.insert(run_algorithm<Koala::ChristofidesEnumerationVertexColoring>(G));
            C.insert(run_algorithm<Koala::BrelazEnumerationVertexColoring>(G));
            C.insert(run_algorithm<Koala::KormanEnumerationVertexColoring>(G));
            break;
        case 1:
            C.insert(run_algorithm<Koala::BrownsOrdinaryEnumerationVertexColoring>(G));
            break;
        case 2:
            C.insert(run_algorithm<Koala::ChristofidesEnumerationVertexColoring>(G));
            break;
        case 3:
            C.insert(run_algorithm<Koala::BrelazEnumerationVertexColoring>(G));
            break;
        case 4:
            C.insert(run_algorithm<Koala::KormanEnumerationVertexColoring>(G));
            break;
        }
        assert(C.size() == 1);

        std::cout << std::endl;
    }
    return 0;
}
