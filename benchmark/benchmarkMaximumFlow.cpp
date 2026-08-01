#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>

#include <benchmark/utils.hpp>
#include <flow/BoykovKolmogorovMaximumFlow.hpp>
#include <flow/ElectricalFlow.hpp>
#include <flow/KingRaoTarjanMaximumFlow.hpp>
#include <flow/MalhotraKumarMaheshwariFlow.hpp>
#include <flow/MaximumFlow.hpp>
#include <flow/GoldbergTarjanPushRelabelMaximumFlow.hpp>
#include <graph/GraphTools.hpp>
#include <io/DimacsGraphReader.hpp>

enum class Algorithm : uint32_t { ALL, PUSH_RELABEL, BK, MKM, KRT, ELECTRICAL_FLOW };

std::map<std::string, Algorithm> ALGORITHM = {
    { "all", Algorithm::ALL },
    { "GoldbergTarjan", Algorithm::PUSH_RELABEL },
    { "BK", Algorithm::BK },
    { "MKM", Algorithm::MKM },
    { "KRT", Algorithm::KRT },
    { "ElectricalFlow", Algorithm::ELECTRICAL_FLOW }
};

template <typename Algorithm>
NetworKit::edgeweight run_algorithm(NetworKit::Graph &G, NetworKit::node s, NetworKit::node t) {
    auto flow = Algorithm(G, s, t);
    flow.run();
    std::cout << flow.getFlowSize() << " " << std::flush;
    return flow.getFlowSize();
}

void run_test(
        NetworKit::Graph &G, NetworKit::node s, NetworKit::node t, Algorithm algorithm) {
    std::set<NetworKit::edgeweight> T;
    switch (algorithm) {
    case Algorithm::ALL:
        T.insert(run_algorithm<Koala::GoldbergTarjanPushRelabelMaximumFlow>(G, s, t));
        T.insert(run_algorithm<Koala::BoykovKolmogorovMaximumFlow>(G, s, t));
        T.insert(run_algorithm<Koala::MalhotraKumarMaheshwariFlow>(G, s, t));
        T.insert(run_algorithm<Koala::KingRaoTarjanMaximumFlow>(G, s, t));
        T.insert(run_algorithm<Koala::ElectricalFlow>(G, s, t));
        assert(T.size() == 1);
        break;
    case Algorithm::PUSH_RELABEL:
        T.insert(run_algorithm<Koala::GoldbergTarjanPushRelabelMaximumFlow>(G, s, t));
        break;
    case Algorithm::BK:
        T.insert(run_algorithm<Koala::BoykovKolmogorovMaximumFlow>(G, s, t));
        break;
    case Algorithm::MKM:
        T.insert(run_algorithm<Koala::MalhotraKumarMaheshwariFlow>(G, s, t));
        break;
    case Algorithm::KRT:
        T.insert(run_algorithm<Koala::KingRaoTarjanMaximumFlow>(G, s, t));
        break;
    case Algorithm::ELECTRICAL_FLOW:
        T.insert(run_algorithm<Koala::ElectricalFlow>(G, s, t));
        break;
    default:
        throw std::logic_error("Unhandled algorithm");
    }
    std::cout << std::endl;
}

void run_g6_tests(const std::string &path, Algorithm algorithm) {
    std::fstream file(path, std::fstream::in);
    std::string line;
    while (file >> line) {
        auto G = Koala::G6GraphReader().readline(line);
        G = Koala::GraphTools::convertUndirectedGraphToDirected(G, true);
        std::cout << line << " " << std::flush;
        run_test(G, 0, G.upperNodeIdBound() - 1, algorithm);
    }
}

void run_dimacs_tests(const std::string &path, Algorithm algorithm) {
    auto [G, s, t] = Koala::DimacsGraphReader().read_all(path);
    std::cout << path << " " << std::flush;
    run_test(G, s, t, algorithm);
}

int main(int argc, const char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <algorithm> <file>" << std::endl;
        return 1;
    }

    std::string algorithm_name(argv[1]), path(argv[2]);
    auto algorithm = ALGORITHM.find(algorithm_name);
    if (algorithm == ALGORITHM.end()) {
        throw std::invalid_argument("Unknown algorithm: " + algorithm_name);
    }

    if (Koala::Benchmark::hasExtension(path, "g6")) {
        run_g6_tests(path, algorithm->second);
    } else if (Koala::Benchmark::hasExtension(path, "gr")
            || Koala::Benchmark::hasExtension(path, "max")) {
        run_dimacs_tests(path, algorithm->second);
    } else {
        throw std::invalid_argument("File type not supported: " + path);
    }
    return 0;
}
