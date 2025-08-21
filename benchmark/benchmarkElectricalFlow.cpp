#include <cassert>
#include <filesystem>
#include <iostream>
#include <map>
#include <chrono>
#include <functional>
#include <sstream>

#include <io/DimacsGraphReader.hpp>

#include <networkit/flow/EdmondsKarp.hpp>

#include <flow/electrical_flow/ElectricalFlow.hpp>
#include <flow/BoykovKolmogorovFlow.hpp>
#include <flow/KingRaoTarjanMaximumFlow.hpp>
#include <flow/PushRelabel.hpp>

using namespace std;
using namespace chrono;
using namespace NetworKit;

struct FlowBenchmark {
    FlowBenchmark(string name, int maxFlow, steady_clock::time_point begin, steady_clock::time_point end) :
        name(name),
        maxFlow(maxFlow),
        duration(duration_cast<microseconds>(end - begin).count()) {}

    string name;
    int maxFlow;
    uint64_t duration;

    friend std::ostream& operator<<(std::ostream& os, const FlowBenchmark& fb) {
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

FlowBenchmark runElectricalFlow(const NetworKit::Graph& G, int s, int t) {
    auto begin = steady_clock::now();
    Koala::ElectricalFlow algo(G, s, t);
    algo.run();
    auto maxFlow = algo.getMaxFlow();
    auto end = steady_clock::now();
    return { "ElectricalFlow", maxFlow, begin, end };
}

FlowBenchmark runFractionalElectricalFlow(const NetworKit::Graph& G, int s, int t) {
    auto begin = steady_clock::now();
    Koala::ElectricalFlow algo(G, s, t, false);
    algo.run();
    auto maxFlow = algo.getMaxFlow();
    auto end = steady_clock::now();
    return { "ElectricalFlow w/o rounding", maxFlow, begin, end };
}

FlowBenchmark runBoykovKolmogorovFlow(NetworKit::Graph& G, int s, int t) {
    auto begin = steady_clock::now();
    Koala::BoykovKolmogorovFlow algo(G, s, t);
    algo.run();
    auto maxFlow = algo.getFlowSize();
    auto end = steady_clock::now();
    return { "BoykovKolmogorovFlow", maxFlow, begin, end };
}

FlowBenchmark runKingRaoTarjanMaximumFlow(NetworKit::Graph& G, int s, int t) {
    auto begin = steady_clock::now();
    Koala::KingRaoTarjanMaximumFlow algo(G, s, t);
    algo.run();
    auto maxFlow = algo.getFlowSize();
    auto end = steady_clock::now();
    return { "KingRaoTarjanFlow", maxFlow, begin, end };
}

FlowBenchmark runPushRelabel(NetworKit::Graph& G, int s, int t) {
    auto begin = steady_clock::now();
    Koala::PushRelabel algo(G, s, t);
    algo.run();
    auto maxFlow = algo.getFlowSize();
    auto end = steady_clock::now();
    return { "PushRelabel", maxFlow, begin, end };
}

tuple<Graph, node, node> readGraph(const string& filepath) {
    if (!filesystem::exists(filepath)) {
        std::stringstream ss;
        ss << "File " << filepath << " does not exist";
        throw ss.str();
    }
    auto [G, s, t] = Koala::DimacsGraphReader().read_all(filepath);
    G.indexEdges();
    return { G, s, t };
}

void benchmarkAverage(function<FlowBenchmark(Graph&, node s, node t)> runAlgorithm, vector<vector<string>> testCases) {
    vector<uint64_t> benchmark;
    string name;
    for (int i = 0; i < testCases.size(); ++i) {
        auto testCase = testCases[i];
        uint64_t totalTime = 0;
        cerr << "Running tc " << i << "...\n";
        for (auto testPath : testCase) {
            auto [G, s, t] = readGraph(testPath);
            G.indexEdges();
            cerr << "\tRunning " << testPath << "...";
            auto result = runAlgorithm(G, s, t);
            cerr << "\tdone\n";
            totalTime += result.duration;
            name = result.name;
        }
        cerr << "done\n";
        benchmark.push_back(totalTime / testCase.size());
    }

    cout << name << '\n';
    for (auto time : benchmark) {
        cout << '\t' << time << "μs\n";
    }
}

int main() {
    tuple<int, int, int> sizes[] = { { 100, 2000, 5 }, { 500, 10000, 2 } };
    int K = 10;

    vector<vector<string>> testCases;
    for (auto [N, M, U] : sizes) {
        vector<string> testCase;
        for (int i = 0; i < K; ++i) {
            std::stringstream ss;
            ss << "./input/test_" << N << "_" << M << "_" << U << "_" << i << ".flow";
            testCase.push_back(ss.str());
        }
        testCases.push_back(testCase);
    }

    benchmarkAverage(runElectricalFlow, testCases);
    benchmarkAverage(runFractionalElectricalFlow, testCases);
    benchmarkAverage(runKingRaoTarjanMaximumFlow, testCases);
    benchmarkAverage(runBoykovKolmogorovFlow, testCases);
    benchmarkAverage(runPushRelabel, testCases);
    benchmarkAverage(runEdmondsKarp, testCases);

    return 0;
}
