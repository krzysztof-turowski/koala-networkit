#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

#include "benchmark/utils.hpp"
#include "vertex_cover/BakerKOuterplanarGraphVertexCover.hpp"
#include "vertex_cover/BakerPlanarGraphVertexCover.hpp"

template <typename T>
int run_algorithm(NetworKit::Graph &graph) {
    auto algorithm = T(graph);
    algorithm.run();
    const auto &vertex_cover = algorithm.getVertexCover();
    std::cout << vertex_cover.size() << " " << std::flush;
    algorithm.check();
    return vertex_cover.size();
}

template <typename T>
int run_approximation_algorithm(
        NetworKit::Graph &graph, double epsilon) {
    auto algorithm = T(graph, epsilon);
    algorithm.run();
    const auto &vertex_cover = algorithm.getVertexCover();
    std::cout << vertex_cover.size() << " " << std::flush;
    algorithm.check();
    return vertex_cover.size();
}

enum class Algorithm : std::uint32_t {
    KOUTERPLANAR,
    PLANAR
};

std::map<std::string, Algorithm> ALGORITHM = {
    {"k-outerplanar", Algorithm::KOUTERPLANAR},
    {"planar", Algorithm::PLANAR},
    {"Baker", Algorithm::PLANAR}
};

void choose_algorithm(
        NetworKit::Graph &graph, Algorithm algorithm, double epsilon) {
    switch (algorithm) {
        case Algorithm::KOUTERPLANAR:
            run_algorithm<Koala::BakerKOuterplanarGraphVertexCover>(graph);
            break;
        case Algorithm::PLANAR:
            run_approximation_algorithm<
                Koala::BakerPlanarGraphVertexCover>(graph, epsilon);
            break;
        default:
            throw std::invalid_argument("Unknown algorithm");
    }
}

int main(int argc, const char *argv[]) {
    if (argc < 3 || argc > 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <algorithm> <file> [epsilon]"
                  << std::endl;
        return 1;
    }
    const std::string algorithm_name(argv[1]);
    const auto algorithm = ALGORITHM.find(algorithm_name);
    if (algorithm == ALGORITHM.end()) {
        throw std::invalid_argument(
            "Unknown algorithm: " + algorithm_name);
    }
    const bool is_approximation = algorithm->second == Algorithm::PLANAR;
    if (argc != (is_approximation ? 4 : 3)) {
        throw std::invalid_argument(
            "The planar Baker algorithm requires an epsilon argument");
    }
    const double epsilon = is_approximation ? std::stod(argv[3]) : 0.0;
    Koala::Benchmark::executeForEachGraph(
        argv[2], [&](const std::string &label, NetworKit::Graph graph) {
            std::cout << label << " " << std::flush;
            choose_algorithm(graph, algorithm->second, epsilon);
            std::cout << std::endl;
        });
    return 0;
}
