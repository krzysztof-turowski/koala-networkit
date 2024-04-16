#include <cassert>
#include <iostream>
#include <map>
#include <fstream>
#include <sys/stat.h>
#include <io/G6GraphReader.hpp>
#include <recognition/planar/HopcroftTarjan.hpp>
#include <sys/types.h>

int main() {
   // mkdir("folder", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    freopen("./input/planar_conn.9.g6", "r", stdin);
    std::map<int, int> classification;
    std::string types[] = {
        "NOT_PLANAR",
        "PLANAR"
    };
    int i = 0;
    while (true) {
        i++;
        std::cout << i << '\n';
        std::string line;
        std::cin >> line;
        if (!std::cin.good()) {
            break;
        }

        NetworKit::Graph G = Koala::G6GraphReader().readline(line);
       
        auto recognize = Koala::HopcroftTarjan(G, false);
        recognize.run();
        classification[recognize.isPlanar() ? 1 : 0]++;
    }
    for (const auto &[k, v] : classification) {
        std::cout << types[static_cast<int>(k)] << ": " << v << std::endl;
    }
    return 0;
}
