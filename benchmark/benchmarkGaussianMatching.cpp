#include <cassert>
#include <filesystem>
#include <iostream>

#include <matching/GaussianMatching.hpp>
#include <io/DimacsGraphReader.hpp>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " filename" << std::endl;
        return 1;
    }
    std::string path(argv[1]);
    if (!std::filesystem::exists(path)) {
        std::cout << "File " << path << " does not exist" << std::endl;
        return 1;
    }

    auto G = Koala::DimacsGraphReader().read(path);
    G.indexEdges(true);
    Koala::GaussianMatching Gmatching(G, true);
    Gmatching.run();

    auto matching = Gmatching.getMatching();
    int size = matching.size();
    if (size >= G.numberOfNodes()) {
        std::cout << "Perfect match (" << size << ")" << std::endl;
    } else {
        std::cout << "Maximum match (" << size << ")" << std::endl;
    }
    
    for (auto [u, v] : matching) {
        std::cout << u << " " << v << std::endl;
    }

    return 0;
}