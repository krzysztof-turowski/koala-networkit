#include <cassert>
#include <iostream>
#include <map>

#include <coloring/EnumerationVertexColoring.hpp>
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
    while (true) {
        std::string line;
        std::cin >> line;
        if (!std::cin.good()) {
            break;
        }
        NetworKit::Graph G = Koala::G6GraphReader().readline(line);
        std::cout << line << " " << G.numberOfEdges() << " ";

        std::set<int> C;
        C.insert(run_algorithm<Koala::BrownsOrdinaryEnumerationVertexColoring>(G));
        C.insert(run_algorithm<Koala::ChristofidesEnumerationVertexColoring>(G));
        C.insert(run_algorithm<Koala::BrelazEnumerationVertexColoring>(G));
        C.insert(run_algorithm<Koala::KormanEnumerationVertexColoring>(G));
        assert(C.size() == 1);

        std::cout << std::endl;
    }
    return 0;
}
