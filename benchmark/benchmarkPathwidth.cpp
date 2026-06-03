#include <iostream>
#include <string>

#include <benchmark/utils.hpp>
#include "pathwidth/CographPathwidth.hpp"
#include "recognition/CographRecognition.hpp"

template<typename T>
int run_algorithm(NetworKit::Graph &G) {
    NetworKit::count pathwidth_size;
    if constexpr (std::is_same_v<T, Koala::CographPathwidth>) {
        auto recognition = Koala::HabibPaulCographRecognition(G);
        recognition.run();
        Koala::CographPathwidth algorithm = T(G, recognition.cotree);
        algorithm.run();
        pathwidth_size = algorithm.getPathwidthSize();
    } else {
        auto algorithm = T(G);
        algorithm.run();
        pathwidth_size = algorithm.getPathwidthSize();
    }

    std::cout << pathwidth_size << " " << std::flush;
    return pathwidth_size;
}

int main(int argc, const char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <file>" << std::endl;
        return 1;
    }
    Koala::Benchmark::executeForEachGraph(
        argv[1], [](const std::string &label, NetworKit::Graph G) {
        std::cout << label << " " << std::flush;
        run_algorithm<Koala::CographPathwidth>(G);
        std::cout << std::endl;
    });
    return 0;
}
