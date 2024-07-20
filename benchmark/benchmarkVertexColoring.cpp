#include <cassert>
#include <iostream>
#include <map>

#include <coloring/ExactVertexColoring.hpp>
#include <coloring/GreedyVertexColoring.hpp>
#include <coloring/PerfectGraphVertexColoring.hpp>
#include <coloring/CographVertexColoring.hpp>
#include <io/G6GraphReader.hpp>

#include <recognition/CographRecognitionOther.hpp>

template<typename T>
int run_algorithm(NetworKit::Graph &G) {
    std::map<NetworKit::node, int> colors;
    if constexpr (std::is_same_v<T, Koala::CographVertexColoring>) {
        auto recognition = Koala::HabibPaulCographRecognition(G);
        recognition.run();
        auto algorithm = T(G, recognition.cotree);
        algorithm.run();
        colors = algorithm.getColoring();
    } else {
        auto algorithm = T(G);
        algorithm.run();
        colors = algorithm.getColoring();
    }

    G.forEdges([&](NetworKit::node u, NetworKit::node v) { assert(colors[u] != colors[v]); });
    int max_color = 0;
    for (const auto &[v, c] : colors) {
        max_color = std::max(max_color, c);
    }
    std::cout << max_color << " " << std::flush;
    return max_color;
}

std::map<std::string, int> ALGORITHM = {
    { "exact", 0 },
    { "RS", 1 }, { "LF", 2 }, { "SL", 3 }, { "SLF", 4 }, { "GIS", 5 },
    { "Brown", 10 }, { "Christofides", 11 }, { "Brelaz", 12 }, { "Korman", 13 },
    { "perfect", 20 }, { "cograph", 30}
};

int main(int argc, char **argv) {
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
        std::set<int> C;
        std::cout << line << " " << std::flush;
        switch (ALGORITHM[std::string(argv[1])]) {
            case 0:
                C.insert(run_algorithm<Koala::BrownEnumerationVertexColoring>(G));
                C.insert(run_algorithm<Koala::ChristofidesEnumerationVertexColoring>(G));
                C.insert(run_algorithm<Koala::BrelazEnumerationVertexColoring>(G));
                C.insert(run_algorithm<Koala::KormanEnumerationVertexColoring>(G));
                assert(C.size() == 1);
                break;
            case 1:
                run_algorithm<Koala::RandomSequentialVertexColoring>(G);
                break;
            case 2:
                run_algorithm<Koala::LargestFirstVertexColoring>(G);
                break;
            case 3:
                run_algorithm<Koala::SmallestLastVertexColoring>(G);
                break;
            case 4:
                run_algorithm<Koala::SaturatedLargestFirstVertexColoring>(G);
                break;
            case 5:
                run_algorithm<Koala::GreedyIndependentSetVertexColoring>(G);
                break;
            case 10:
                run_algorithm<Koala::BrownEnumerationVertexColoring>(G);
                break;
            case 11:
                run_algorithm<Koala::ChristofidesEnumerationVertexColoring>(G);
                break;
            case 12:
                run_algorithm<Koala::BrelazEnumerationVertexColoring>(G);
                break;
            case 13:
                run_algorithm<Koala::KormanEnumerationVertexColoring>(G);
                break;
            case 20:
                run_algorithm<Koala::PerfectGraphVertexColoring>(G);
                break;
            case 30:
                run_algorithm<Koala::CographVertexColoring>(G);
        }
        std::cout << std::endl;
    }
    return 0;
}
