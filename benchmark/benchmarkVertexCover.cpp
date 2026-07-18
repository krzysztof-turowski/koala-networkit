#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

#include "benchmark/utils.hpp"
#include "vertex_cover/BakerKOuterplanarGraphVertexCover.hpp"

template <typename T>
int run_algorithm(NetworKit::Graph &graph) {
    auto algorithm = T(graph);
    algorithm.run();
    const auto &vertex_cover = algorithm.getVertexCover();
    std::cout << vertex_cover.size() << " " << std::flush;
    algorithm.check();
    return vertex_cover.size();
}

enum class Algorithm : std::uint32_t {
    KOUTERPLANAR
};

std::map<std::string, Algorithm> ALGORITHM = {
    {"k-outerplanar", Algorithm::KOUTERPLANAR}
};

void choose_algorithm(NetworKit::Graph &graph, Algorithm algorithm) {
    switch (algorithm) {
        case Algorithm::KOUTERPLANAR:
            run_algorithm<Koala::BakerKOuterplanarGraphVertexCover>(graph);
            break;
        default:
            throw std::invalid_argument("Unknown algorithm");
    }
}

int main(int argc, const char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <algorithm> <file>"
                  << std::endl;
        return 1;
    }
    const std::string algorithm_name(argv[1]);
    const auto algorithm = ALGORITHM.find(algorithm_name);
    if (algorithm == ALGORITHM.end()) {
        throw std::invalid_argument(
            "Unknown algorithm: " + algorithm_name);
    }
    Koala::Benchmark::executeForEachGraph(
        argv[2], [&](const std::string &label, NetworKit::Graph graph) {
            std::cout << label << " " << std::flush;
            choose_algorithm(graph, algorithm->second);
            std::cout << std::endl;
        });
    return 0;
}
