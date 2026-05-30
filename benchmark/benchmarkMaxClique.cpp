#include <iostream>
#include <set>
#include <string>

#include <benchmark/utils.hpp>
#include <clique/Clique.hpp>
#include <clique/CographClique.hpp>
#include <recognition/CographRecognition.hpp>

template<typename T>
int run_algorithm(NetworKit::Graph &G) {
    std::set<NetworKit::node> max_clique;
    if constexpr (std::is_same_v<T, Koala::CographMaxClique>) {
        auto recognition = Koala::HabibPaulCographRecognition(G);
        recognition.run();
        Koala::CographMaxClique algorithm = T(G, recognition.cotree);
        algorithm.run();
        max_clique = algorithm.getMaxCliqueSet();
    } else {
        auto algorithm = T(G);
        algorithm.run();
        max_clique = algorithm.getMaxCliqueSet();
    }

    std::cout << max_clique.size() << " " << std::flush;
    return max_clique.size();
}

int main(int argc, const char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <file>" << std::endl;
        return 1;
    }
    Koala::Benchmark::executeForEachGraph(
        argv[1], [](const std::string &label, NetworKit::Graph G) {
        std::cout << label << " " << std::flush;
        run_algorithm<Koala::CographMaxClique>(G);
        std::cout << std::endl;
    });
    return 0;
}
