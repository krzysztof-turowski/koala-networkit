#include <cassert>
#include <filesystem>
#include <iostream>
#include <map>
#include <chrono>

#include <electric_flow/ElectricFlow.hpp>
#include <io/DimacsGraphReader.hpp>
#include <networkit/flow/EdmondsKarp.hpp>

using namespace std;
using namespace chrono;

struct FlowBenchmark {
    FlowBenchmark(string name, int maxFlow, steady_clock::time_point begin, steady_clock::time_point end) :
        name(name),
        maxFlow(maxFlow),
        duration(duration_cast<microseconds>(end - begin).count()) {}

    string name;
    int maxFlow;
    uint64_t duration;

    friend std::ostream& operator<<( std::ostream& os, const FlowBenchmark& fb ) {
        os << "Benchmark for " << fb.name << "algorithm:\n";
        os << "\tMaximum flow:\t" << fb.maxFlow << '\n';
        os << "\tTime:\t" << fb.duration << "μs\n";
        return os;
    }
};

FlowBenchmark runEdmondsKarp(const NetworKit::Graph& G, int s, int t) {
    auto begin = steady_clock::now();
    NetworKit::EdmondsKarp algo(G, s, t);
    algo.run();
    auto maxFlow = algo.getMaxFlow();
    auto end = steady_clock::now();
    return { "EdmondsKarp", maxFlow, begin, end };
}

FlowBenchmark runElectricFlow(const NetworKit::Graph& G, int s, int t) {
    auto begin = steady_clock::now();
    Koala::ElectricFlow algo(G, s, t);
    algo.run();
    auto maxFlow = algo.getMaxFlow();
    auto end = steady_clock::now();
    return { "ElectricFlow", maxFlow, begin, end };
}

int main(int argc, char** argv) {
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
    G.indexEdges();
    cout << runElectricFlow(G, s, t) << '\n';
    cout << runEdmondsKarp(G, s, t) << '\n';
    return 0;
}
