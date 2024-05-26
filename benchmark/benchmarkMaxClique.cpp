#include <iostream>

#include <io/G6GraphReader.hpp>
#include "max_clique/MaxClique.hpp"
#include "max_clique/CographMaxClique.hpp"

template<typename T>
int run_algorithm(NetworKit::Graph &G) {
    std::set<NetworKit::node> max_clique;
    if constexpr (std::is_same_v<T, Koala::CographMaxClique>) {
        auto recognition = Koala::CographRecognition(G);
        recognition.run();
        Koala::CographMaxClique algorithm = T(G, recognition.cotree);
        algorithm.run();
        max_clique = algorithm.getMaxCliqueSet();
    } else {
        auto algorithm = T(G);
        algorithm.run();
        max_clique = algorithm.getMaxCliqueSet();
    }

    std::cout << max_clique.size() << " " << std::flush;
    return max_clique.size();
}

std::map<std::string, int> ALGORITHM = {
        {"cograph", 0}
};


int main(int argc, const char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <algorithm> <file>" << std::endl;
        return 1;
    }
    std::string path(argv[2]);
    auto position = path.find_last_of(".");
    if (path.substr(position + 1) == "g6") {
        std::fstream file(path, std::fstream::in);
        std::map<int, int> classification;
        while (true) {
            std::string line;
            file >> line;
            if (!file.good()) {
                break;
            }
            NetworKit::Graph G = Koala::G6GraphReader().readline(line);
            std::cout << line << " " << std::flush;
            switch (ALGORITHM[std::string(argv[1])]) {
                case 0:
                    run_algorithm<Koala::CographMaxClique>(G);
                    break;
            }
            std::cout << std::endl;
        }
    } else {
        std::cerr << "File type not supported: " << path << std::endl;
    }
    return 0;
}
