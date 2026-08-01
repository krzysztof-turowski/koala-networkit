#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <benchmark/utils.hpp>
#include <clique_cover/ChordalCliqueCover.hpp>
#include <recognition/ChordalGraphRecognition.hpp>

template<typename T>
int run_algorithm(NetworKit::Graph &G) {
    std::vector<std::vector<NetworKit::node>> cover;
    if constexpr (std::is_same_v<T, Koala::ChordalMinCliqueCover>) {
        auto recognition = Koala::MaximumCardinalitySearchChordalGraphRecognition(G);
        recognition.run();
        if (!recognition.isChordal()) {
            throw std::invalid_argument("Graph is not chordal");
        }
        auto algorithm = T(G, recognition.getPEO());
        algorithm.run();
        cover = algorithm.getCliqueCover();
    }

    std::cout << cover.size() << " " << std::flush;
    return cover.size();
}

enum class Algorithm : uint32_t { CHORDAL };

std::map<std::string, Algorithm> ALGORITHM = {
    { "chordal", Algorithm::CHORDAL }
};

void choose_algorithm(NetworKit::Graph &G, Algorithm algorithm) {
    switch (algorithm) {
    case Algorithm::CHORDAL:
        run_algorithm<Koala::ChordalMinCliqueCover>(G);
        break;
    default:
        throw std::logic_error("Unhandled algorithm");
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
