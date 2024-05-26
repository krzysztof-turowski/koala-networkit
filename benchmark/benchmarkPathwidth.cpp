#include <iostream>

#include <io/G6GraphReader.hpp>

#include "pathwidth/CographPathwidth.hpp"

template<typename T>
int run_algorithm(NetworKit::Graph &G) {

    NetworKit::count pathwidth_size;
    if constexpr (std::is_same_v<T, Koala::CographPathwidth>) {
        auto recognition = Koala::CographRecognition(G);
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
                    run_algorithm<Koala::CographPathwidth>(G);
                    break;
            }
            std::cout << std::endl;
        }
    } else {
        std::cerr << "File type not supported: " << path << std::endl;
    }
    return 0;
}
