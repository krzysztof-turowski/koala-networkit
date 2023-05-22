#include <cassert>
#include <filesystem>
#include <iostream>

#include <matching/MaximumMatching.hpp>
#include <io/DimacsGraphReader.hpp>

int main(int argc, char **argv) {   
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " filename" << std::endl;
        return 0;
    }
    std::string path(argv[1]);
    if (!std::filesystem::exists(path)) {
        std::cout << "File " << path << " does not exist" << std::endl;
        return 0;
    }
    auto G = Koala::DimacsGraphReader().read(path);
    G.indexEdges(true);
    auto maximum_matching = Koala::GabowMaximumMatching(G);
    maximum_matching.run();
    auto matching = maximum_matching.getMatching();
    NetworKit::edgeweight weight = 0.0;
    G.forNodes([&G, &matching, &weight] (NetworKit::node v) {
        if (matching[v] != NetworKit::none && v < matching[v]) {
            weight += G.weight(v, matching[v]);
        }
        std::cout << v << " " << matching[v] << std::endl;
    });
    std::cout << weight << std::endl;
    return 0;
}
