#include <cassert>
#include <filesystem>
#include <iostream>

#include <matching/MaximumMatching.hpp>
#include <io/DimacsGraphReader.hpp>

int main(int argc, char **argv) {   
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " algorithm filename" << std::endl;
        return 0;
    }
    std::string algorithm(argv[1]);
    std::string path(argv[2]);
    if (!std::filesystem::exists(path)) {
        std::cout << "File " << path << " does not exist" << std::endl;
        return 0;
    }
    auto G = Koala::DimacsGraphReader().read(path);
    G.indexEdges(true);
    Koala::MaximumMatching * maximum_matching;
    if (algorithm == "edmonds") {
        maximum_matching = new Koala::EdmondsMaximumMatching(G);
    } else if (algorithm == "gabow") {
        maximum_matching = new Koala::GabowMaximumMatching(G);
    } else if (algorithm == "micali") {
        maximum_matching = new Koala::MicaliGabowMaximumMatching(G);
    } else {
        std::cout << "Unknown algorithm: " << algorithm << std::endl;
    }
    maximum_matching->run();
    auto matching = maximum_matching->getMatching();
    NetworKit::edgeweight weight = 0.0;
    for (auto [u, v] : matching) {
        std::cout << u << " " << v << std::endl;
        if (v != NetworKit::none)
            weight += G.weight(u, v);
    }
    std::cout << weight / 2 << std::endl;

    delete maximum_matching;
    return 0;
}
