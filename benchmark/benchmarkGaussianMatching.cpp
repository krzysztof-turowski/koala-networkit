#include <cassert>
#include <filesystem>
#include <iostream>
#include <ctime>
#include <random>

#include <matching/gaussian_matching/GeneralGaussianMatching.hpp>
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
    int n = G.numberOfNodes();

    srand( atoi(argv[2]) );
    Koala::GeneralGaussianMatching gen(G);
    gen.run();

    auto M = gen.getMatching();
    if (M.size() != n / 2) {
        std::cout << "Matching found size is incorrect: " << M.size() << "(instead of " << n/2 << ")\n";
        return 1;
    }

    std::vector<int> counts(n, 0);
    for (auto [u, v] : M) {
        counts[u]++;
        counts[v]++;

        if (!G.hasEdge(u, v)) {
            std::cout << "Match (" << u << ", " << v << ") is not a correct edge\n";
            return 2;
        }
    }

    for (int v = 0; v < n; ++v) {
        if (counts[v] != 1) {
            std::cout << "Vertex " << v << "was matched " << counts[v] << "times\n";
            return 3;
        }
    }

    return 0;
}