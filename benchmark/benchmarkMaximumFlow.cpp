#include <cassert>
#include <exception>
#include <filesystem>
#include <iostream>
#include <map>

#include <flow/MaximumFlow.hpp>
#include <flow/PushRelabel.hpp>
#include <flow/BoykovKolmogorovFlow.hpp>
#include <flow/KingRaoTarjanMaximumFlow.hpp>
#include <flow/MalhotraKumarMaheshwariFlow.hpp>
#include <flow/electrical_flow/ElectricalFlow.hpp>
#include <graph/GraphTools.hpp>
#include <io/DimacsGraphReader.hpp>

std::map<std::string, int> ALGORITHM = {
    { "all", 0 },
    { "PushRelabel", 1 },
    { "BK", 2 },
    { "MKM", 3 },
    { "KRT", 4 },
    { "ElectricalFlow", 5 }
};

template <typename Algorithm>
NetworKit::edgeweight run_algorithm(NetworKit::Graph &G, NetworKit::node s, NetworKit::node t) {
    auto flow = Algorithm(G, s, t);
    flow.run();
    std::cout << flow.getFlowSize() << " " << std::flush;
    return flow.getFlowSize();
}

void run_test(
        NetworKit::Graph &G, NetworKit::node s, NetworKit::node t, const std::string &algorithm) {
    std::set<NetworKit::edgeweight> T;
    switch (ALGORITHM[algorithm]) {
    case 0:
        T.insert(run_algorithm<Koala::PushRelabel>(G, s, t));
        T.insert(run_algorithm<Koala::BoykovKolmogorovFlow>(G, s, t));
        T.insert(run_algorithm<Koala::MalhotraKumarMaheshwariFlow>(G, s, t));
        T.insert(run_algorithm<Koala::KingRaoTarjanMaximumFlow>(G, s, t));
        T.insert(run_algorithm<Koala::ElectricalFlow>(G, s, t));
        assert(T.size() == 1);
        break;
    case 1:
        T.insert(run_algorithm<Koala::PushRelabel>(G, s, t));
        break;
    case 2:
        T.insert(run_algorithm<Koala::BoykovKolmogorovFlow>(G, s, t));
        break;
    case 3:
        T.insert(run_algorithm<Koala::MalhotraKumarMaheshwariFlow>(G, s, t));
        break;
    case 4:
        T.insert(run_algorithm<Koala::KingRaoTarjanMaximumFlow>(G, s, t));
        break;
    case 5:
        T.insert(run_algorithm<Koala::ElectricalFlow>(G, s, t));
        break;
    default:
        std::cout << "Unknown algorithm: " << algorithm << std::endl;
        throw std::exception();
    }
    std::cout << std::endl;
}

void run_dimacs_tests(const std::string &path, const std::string &algorithm) {
    auto [G, s, t] = Koala::DimacsGraphReader().read_all(path);
    G = Koala::GraphTools::convertDirectedGraphToUndirected(G, true);
    G = Koala::GraphTools::convertUndirectedGraphToDirected(G, true);
    std::cout << path << " " << std::flush;
    run_test(G, s, t, algorithm);
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
