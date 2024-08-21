#include <cassert>
#include <filesystem>
#include <iostream>
#include <chrono>

#include <matching/MaximumMatching.hpp>
#include <io/DimacsGraphReader.hpp>

int calc_matching_weight(
        const NetworKit::Graph& G, const std::map<NetworKit::node, NetworKit::node> matching) {
    int weight = 0;
    for (auto [u, v] : matching) {
        if (v != NetworKit::none)
            weight += static_cast<int>(G.weight(u, v));
    }
    return weight / 2;
}

template<typename Algorithm>
void test_algorithm(NetworKit::Graph& G, bool perfect) {
    auto start = std::chrono::high_resolution_clock::now();
    Algorithm algo(G, perfect);
    algo.run();
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed_seconds = end - start;

    auto matching = algo.getMatching();
    auto weight = calc_matching_weight(G, matching);
    
    std::cout << "Time: " << elapsed_seconds.count() << "s\n";
    std::cout << weight << std::endl;
}

bool check_matching_perfect(const std::map<NetworKit::node, NetworKit::node>& matching) {
    for (auto [v, u] : matching)
        if (u == NetworKit::none)
            return false;
    return true;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " algorithm filename [--pefrect]" << std::endl;
        return 1;
    }
    std::string algorithm(argv[1]);
    std::string path(argv[2]);
    if (!std::filesystem::exists(path)) {
        std::cout << "File " << path << " does not exist" << std::endl;
        return 1;
    }
    bool perfect = false;
    if (argc > 3 && std::string(argv[3]) == "--perfect") {
        perfect = true;
    } 
    auto G = Koala::DimacsGraphReader().read(path);
    G.indexEdges(true);

    if (perfect) {
        Koala::MicaliVaziraniMatching mcm(G);
        mcm.run();
        if (!check_matching_perfect(mcm.getMatching())) {
            std::cout << "The provided graph does not contain a perfect matching" << std::endl;
            return 1;
        }
    }

    if (algorithm == "edmonds") {
        test_algorithm<Koala::EdmondsMaximumMatching>(G, perfect);
    } else if (algorithm == "gabow") {
        test_algorithm<Koala::GabowMaximumMatching>(G, perfect);
    } else if (algorithm == "micali") {
        test_algorithm<Koala::GalilMicaliGabowMaximumMatching>(G, perfect);
    } else if (algorithm == "scaling") {
        test_algorithm<Koala::GabowScalingMatching>(G, perfect);
    } else {
        std::cout << "Unknown algorithm: " << algorithm << std::endl;
        return 1;
    }

    return 0;
}
