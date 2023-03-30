#include <cassert>
#include <iostream>
#include <map>

#include <io/DimacsGraphReader.hpp>
#include <flow/MaximumFlow.hpp>

int main(int argc, char **argv) {
    std::string path(argv[1]);
    const auto &[G, s, t] = Koala::DimacsGraphReader().read_all(path);
    auto maximum_flow = Koala::KingRaoTarjanMaximumFlow(G, s, t);
    maximum_flow.run();
    std::cout << path << " has maximum flow of size " << maximum_flow.getFlowSize() << std::endl;
    return 0;
}
