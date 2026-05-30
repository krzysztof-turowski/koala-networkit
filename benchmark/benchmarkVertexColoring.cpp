#include <concepts>
#include <cassert>

#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <string>

#include <benchmark/utils.hpp>
#include <coloring/CographVertexColoring.hpp>
#include <coloring/ExactVertexColoring.hpp>
#include <coloring/GreedyVertexColoring.hpp>
#include <coloring/PerfectGraphVertexColoring.hpp>
#include <recognition/CographRecognition.hpp>
#include <recognition/PerfectGraphRecognition.hpp>

template<typename T>
NetworKit::count run_algorithm(NetworKit::Graph &G) {
    auto algorithm = T(G);
    algorithm.run();
    auto colors = algorithm.getColoring();

    algorithm.check();
    auto max_color = algorithm.getMaximumColor();
    std::cout << max_color << " " << std::flush;
    return max_color;
}

template<typename T>
requires std::derived_from<T, Koala::CographVertexColoring>
NetworKit::count run_algorithm(NetworKit::Graph &G) {
    auto recognition = Koala::HabibPaulCographRecognition(G);
    recognition.run();

    if (!recognition.isCograph()) {
        throw std::invalid_argument("Graph is not a cograph");
    }
    auto algorithm = T(G, recognition.cotree);
    algorithm.run();
    auto colors = algorithm.getColoring();

    algorithm.check();
    auto max_color = algorithm.getMaximumColor();
    std::cout << max_color << " " << std::flush;
    return max_color;
}

template<typename T>
requires std::derived_from<T, Koala::PerfectGraphVertexColoring>
NetworKit::count run_algorithm(NetworKit::Graph &G) {
    auto recognition = Koala::PerfectGraphRecognition(G);
    recognition.run();

    if (recognition.getState() != Koala::PerfectGraphRecognition::State::PERFECT) {
        throw std::invalid_argument("Graph is not perfect");
    }
    auto algorithm = T(G);
    algorithm.run();
    auto colors = algorithm.getColoring();

    algorithm.check();
    auto max_color = algorithm.getMaximumColor();
    std::cout << max_color << " " << std::flush;
    return max_color;
}

enum class Algorithm : uint32_t {
    EXACT = 1,
    RS, LF, SL, SLF, GIS,
    BROWN = 10, CHRISTOFIDES, BRELAZ, KORMAN,
    PERFECT = 100, COGRAPH
};

std::map<std::string, Algorithm> ALGORITHM = {
    { "exact", Algorithm::EXACT },
    { "RS", Algorithm::RS },
    { "LF", Algorithm::LF },
    { "SL", Algorithm::SL },
    { "SLF", Algorithm::SLF },
    { "GIS", Algorithm::GIS },
    { "Brown", Algorithm::BROWN },
    { "Christofides", Algorithm::CHRISTOFIDES },
    { "Brelaz", Algorithm::BRELAZ },
    { "Korman", Algorithm::KORMAN },
    { "perfect", Algorithm::PERFECT },
    { "cograph", Algorithm::COGRAPH }
};

void choose_algorithm(NetworKit::Graph &G, Algorithm algorithm) {
    std::set<int> colors;
    switch (algorithm) {
    case Algorithm::EXACT:
        colors.insert(run_algorithm<Koala::BrownEnumerationVertexColoring>(G));
        colors.insert(run_algorithm<Koala::ChristofidesEnumerationVertexColoring>(G));
        colors.insert(run_algorithm<Koala::BrelazEnumerationVertexColoring>(G));
        colors.insert(run_algorithm<Koala::KormanEnumerationVertexColoring>(G));
        assert(colors.size() == 1);
        break;
    case Algorithm::RS:
        run_algorithm<Koala::RandomSequentialVertexColoring>(G);
        break;
    case Algorithm::LF:
        run_algorithm<Koala::LargestFirstVertexColoring>(G);
        break;
    case Algorithm::SL:
        run_algorithm<Koala::SmallestLastVertexColoring>(G);
        break;
    case Algorithm::SLF:
        run_algorithm<Koala::SaturatedLargestFirstVertexColoring>(G);
        break;
    case Algorithm::GIS:
        run_algorithm<Koala::GreedyIndependentSetVertexColoring>(G);
        break;
    case Algorithm::BROWN:
        run_algorithm<Koala::BrownEnumerationVertexColoring>(G);
        break;
    case Algorithm::CHRISTOFIDES:
        run_algorithm<Koala::ChristofidesEnumerationVertexColoring>(G);
        break;
    case Algorithm::BRELAZ:
        run_algorithm<Koala::BrelazEnumerationVertexColoring>(G);
        break;
    case Algorithm::KORMAN:
        run_algorithm<Koala::KormanEnumerationVertexColoring>(G);
        break;
    case Algorithm::PERFECT:
        run_algorithm<Koala::PerfectGraphVertexColoring>(G);
        break;
    case Algorithm::COGRAPH:
        run_algorithm<Koala::CographVertexColoring>(G);
        break;
    }
}

int main(int argc, char **argv) {
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
