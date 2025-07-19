#include <cassert>
#include <filesystem>
#include <iostream>
#include <map>

#include <electric_flow/ElectricFlow.hpp>
#include <io/DimacsGraphReader.hpp>

int main(int argc, char **argv) {
    if (argc <= 1) {
        std::cout << "Filename is empty" << std::endl;
        return 1;
    }

    std::string path(argv[1]);
    if (!std::filesystem::exists(path)) {
        std::cout << "File " << path << " does not exist" << std::endl;
        return 1;
    }

    auto [G, s, t] = Koala::DimacsGraphReader().read_all(path);
    auto maximum_flow = Koala::ElectricFlow(G, s, t);
    maximum_flow.run();
    std::cout << path << " has maximum flow of size " << maximum_flow.getMaxFlow() << std::endl;
    return 0;
}
