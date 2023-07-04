#include <cassert>
#include <filesystem>
#include <iostream>
#include <map>

#include <flow/MaximumFlow.hpp>
#include <io/DimacsGraphReader.hpp>

int main(int argc, char **argv) {
    if (argc == 1) {
        std::cout << "Filename is empty" << std::endl;
        return 0;
    }
    std::string path(argv[1]);
    if (!std::filesystem::exists(path)) {
        std::cout << "File " << path << " does not exist" << std::endl;
        return 0;
    }
    auto [G, s, t] = Koala::DimacsGraphReader().read_all(path);
    auto maximum_flow = Koala::KingRaoTarjanMaximumFlow(G, s, t);
    maximum_flow.run();
    std::cout << path << " has maximum flow of size " << maximum_flow.getFlowSize() << std::endl;
    return 0;
}
