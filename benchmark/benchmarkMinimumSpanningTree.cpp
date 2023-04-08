#include <cassert>
#include <iostream>

#include <networkit/graph/KruskalMSF.hpp>

#include <io/G6GraphReader.hpp>

int main() {
    while (true) {
        std::string line;
        std::cin >> line;
        if (!std::cin.good()) {
            break;
        }
        NetworKit::Graph G = Koala::G6GraphReader().readline(line);
    }
    return 0;
}
