#include <cassert>
#include <exception>
#include <filesystem>
#include <iostream>
#include <map>

#include <io/DimacsGraphReader.hpp>
#include <matching/MaximumMatching.hpp>
#include <matching/gaussian_matching/BipartiteGaussianMatching.hpp>
#include <matching/gaussian_matching/GeneralGaussianMatching.hpp>
#include <matching/gaussian_matching/NaiveGaussianMatching.hpp>

std::map<std::string, int> ALGORITHM = {
    { "all", 0 },
    { "MicaliVazirani", 1 },
    { "NaiveGaussian", 2 },
    { "Gaussian", 3 },
    { "Edmonds", 4 },
    { "Gabow", 5 },
};

template <typename Algorithm>
NetworKit::count run_algorithm(NetworKit::Graph &G) {
    auto algorithm = Algorithm(G);
    algorithm.run();
    if constexpr (std::same_as<Algorithm, Koala::MicaliVaziraniMatching>) {
        auto matching = algorithm.getMatching();
        NetworKit::count size = std::count_if(
            matching.begin(), matching.end(),
            [](auto pair) { return pair.second != NetworKit::none; });
        std::cout << size / 2 << " " << std::flush;
        return size / 2;
    }
    if constexpr (std::same_as<Algorithm, Koala::GeneralGaussianMatching>) {
        std::cout << algorithm.getMatching().size() << " " << std::flush;
        return algorithm.getMatching().size();
    }
    std::cout << algorithm.getMatching().size() / 2 << " " << std::flush;
    return algorithm.getMatching().size() / 2;
}

void run_test(NetworKit::Graph &G, const std::string &algorithm) {
    std::set<NetworKit::count> T;
    switch (ALGORITHM[algorithm]) {
    case 0:
        T.insert(run_algorithm<Koala::MicaliVaziraniMatching>(G));
        T.insert(run_algorithm<Koala::GeneralGaussianMatching>(G));
        T.insert(run_algorithm<Koala::NaiveGaussianMatching>(G));
        T.insert(run_algorithm<Koala::EdmondsMaximumMatching>(G));
        T.insert(run_algorithm<Koala::GabowMaximumMatching>(G));
        assert(T.size() == 1);
        break;
    case 1:
        T.insert(run_algorithm<Koala::MicaliVaziraniMatching>(G));
        break;
    case 2:
        T.insert(run_algorithm<Koala::GeneralGaussianMatching>(G));
        break;
    case 3:
        T.insert(run_algorithm<Koala::NaiveGaussianMatching>(G));
        break;
    case 4:
        T.insert(run_algorithm<Koala::EdmondsMaximumMatching>(G));
        break;
    case 5:
        T.insert(run_algorithm<Koala::GabowMaximumMatching>(G));
        break;
    default:
        std::cout << "Unknown algorithm: " << algorithm << std::endl;
        throw std::exception();
    }
    std::cout << std::endl;
}

void run_dimacs_tests(const std::string &path, const std::string &algorithm) {
    auto G = Koala::DimacsGraphReader().read(path);
    G.indexEdges(true);
    std::cout << path << " " << std::flush;
    run_test(G, algorithm);
}

int main(int argc, const char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <algorithm> <file>" << std::endl;
        return 1;
    }

    std::string algorithm(argv[1]), path(argv[2]);
    if (!std::filesystem::exists(path)) {
        std::cerr << "File " << path << " does not exist" << std::endl;
        return 1;
    }
    if (std::filesystem::is_directory(path)) {
        std::cerr << path << " is a directory" << std::endl;
        return 1;
    }

    if (path.substr(path.find_last_of(".") + 1) == "max") {
        run_dimacs_tests(path, algorithm);
    } else {
        std::cerr << "File type not supported: " << path << std::endl;
    }
    return 0;
}
