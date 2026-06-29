#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <benchmark/utils.hpp>
#include <clique/ChordalClique.hpp>
#include <clique/CographClique.hpp>
#include <recognition/ChordalGraphRecognition.hpp>
#include <recognition/CographRecognition.hpp>

template<typename T>
int run_algorithm(NetworKit::Graph &G) {
    std::vector<NetworKit::node> max_clique;
    if constexpr (std::is_same_v<T, Koala::CographMaxClique>) {
        auto recognition = Koala::HabibPaulCographRecognition(G);
        recognition.run();
        Koala::CographMaxClique algorithm = T(G, recognition.cotree);
        algorithm.run();
        max_clique = algorithm.getMaxClique();
    } else if constexpr (std::is_same_v<T, Koala::ChordalMaxClique>) {
        auto recognition = Koala::MaximumCardinalitySearchChordalGraphRecognition(G);
        recognition.run();
        if (!recognition.isChordal()) {
            throw std::invalid_argument("Graph is not chordal");
        }
        auto algorithm = T(G, recognition.getPEO());
        algorithm.run();
        max_clique = algorithm.getMaxClique();
    } else {
        auto algorithm = T(G);
        algorithm.run();
        max_clique = algorithm.getMaxClique();
    }

    std::cout << max_clique.size() << " " << std::flush;
    return max_clique.size();
}

enum class Algorithm : uint32_t { COGRAPH, CHORDAL };

std::map<std::string, Algorithm> ALGORITHM = {
    { "cograph", Algorithm::COGRAPH },
    { "chordal", Algorithm::CHORDAL }
};

void choose_algorithm(NetworKit::Graph &G, Algorithm algorithm) {
    switch (algorithm) {
    case Algorithm::COGRAPH:
        run_algorithm<Koala::CographMaxClique>(G);
        break;
    case Algorithm::CHORDAL:
        run_algorithm<Koala::ChordalMaxClique>(G);
        break;
    }
}

int main(int argc, const char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <algorithm> <file>" << std::endl;
        return 1;
    }
    std::string algorithm_name(argv[1]);
    auto algorithm = ALGORITHM.find(algorithm_name);
    if (algorithm == ALGORITHM.end()) {
        throw std::invalid_argument("Unknown algorithm: " + algorithm_name);
    }
    Koala::Benchmark::executeForEachGraph(
        argv[2], [&](const std::string &label, NetworKit::Graph G) {
        std::cout << label << " " << std::flush;
        choose_algorithm(G, algorithm->second);
        std::cout << std::endl;
    });
    return 0;
}
