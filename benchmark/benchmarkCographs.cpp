#include <io/G6GraphReader.hpp>

#include <iostream>
#include <map>
#include <string>

#include "recognition/CographRecognition.hpp"

std::string types[] = {
    "UNKNOWN",
    "COGRAPH",
    "NOT_COGRAPH",
    "CONTAINS_0_NODE",
    "EXISTS_1_NODE_NOT_PROPERLY_MARKED",
    "GRANDPARENT_IS_NOT_IN_SET",
    "NO_ONE_PATH",
    "WRONG_PARENT",
    "WRONG_GRANDPARENT"
};

template <typename T>
int run_algorithm(NetworKit::Graph &G, bool verbose = false) {
    auto algorithm = T(G);
    algorithm.run();
    if (verbose) {
        std::cout << static_cast<int>(algorithm.getState()) << " " << std::flush;
    }
    if (algorithm.isCograph()) {
        algorithm.check();
    }
    return static_cast<int>(algorithm.getState());
}

std::map<std::string, int> ALGORITHM = {
    { "all", 0 }, { "BCHP", 1 }, { "CSP", 2 }, { "Dahlhaus", 3 }
};

void run_g6_tests(const std::string &path, const std::string &algorithm) {
    std::fstream file(path, std::fstream::in);
    std::map<int, int> classification;
    while (true) {
        std::string line;
        file >> line;
        if (!file.good()) {
            break;
        }
        auto G = Koala::G6GraphReader().readline(line);
        std::set<int> T;
        switch (ALGORITHM[algorithm]) {
        case 0:
            std::cout << line << " " << std::flush;
            T.insert(run_algorithm<Koala::BretscherCorneilHabibPaulCographRecognition>(G, true));
            T.insert(
                run_algorithm<Koala::CorneilStewartPerlCographRecognition>(G, true) != 1 ? 2 : 1);
            T.insert(run_algorithm<Koala::DahlhausCographRecognition>(G, true));
            std::cout << std::endl;
            assert(T.size() == 1);
            break;
        case 1:
            classification[run_algorithm<Koala::BretscherCorneilHabibPaulCographRecognition>(G)]++;
            break;
        case 2:
            classification[run_algorithm<Koala::CorneilStewartPerlCographRecognition>(G)]++;
            break;
        case 3:
            classification[run_algorithm<Koala::DahlhausCographRecognition>(G)]++;
            break;
        }
    }
    if (ALGORITHM[algorithm]) {
        for (const auto &[k, v] : classification) {
            std::cout << types[k] << ": " << v << std::endl;
        }
    }
}

int main(int argc, const char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <algorithm> <file>" << std::endl;
        return 1;
    }
    std::string path(argv[2]);
    auto extension = path.substr(path.find_last_of(".") + 1);
    if (extension == "g6") {
        run_g6_tests(path, std::string(argv[1]));
    } else {
        std::cerr << "File type not supported: " << path << std::endl;
    }
    return 0;
}
