// #include <iostream>
#include <dominatingset/RooijBodlaenderMDS.hpp>
#include <dominatingset/BranchAndReduceMDS.hpp>
#include <dominatingset/BranchAndReduceMSC.hpp>
#include <dominatingset/SchiermeyerMDS.hpp>
#include <dominatingset/FominKratschWoeginger.hpp>
#include <dominatingset/ExhaustiveMDS.hpp>
#include <chrono>
#include <io/G6GraphReader.hpp>
#include <fstream>
#include<algorithm>
#include<iterator>


int dsSize(const std::vector<bool> &set) {
    return std::count(set.begin(), set.end(), true);
}

// void printGraph(const NetworKit::Graph &graph) {
//     std::cout << "graph\n";
//     graph.forNodes([&graph](NetworKit::node u) {
//         std::cout << u << " has neighbors [";
//         graph.forEdgesOf(u, [u](NetworKit::node v) {
//             std::cout << v << " ";
//         });
//         std::cout << "]\n";
//     });
// }

template <typename A>
void test_algo_g6(const std::string inputFileName) {
    std::ifstream inputStream(inputFileName);
    assert(inputStream.good());    
    std::string g6Graph;
    int processed = 0;
    while (std::getline(inputStream, g6Graph)) {
        Koala::G6GraphReader reader;
        NetworKit::Graph graph = reader.readline(g6Graph);
        A algo(graph);
        ExhaustiveMDS exhaustive(graph);

        algo.run();
        exhaustive.run();

        assert(algo.isDominating(algo.getDominatingSet()));
        assert(exhaustive.isDominating(exhaustive.getDominatingSet()));
        assert(dsSize(algo.getDominatingSet()) == dsSize(exhaustive.getDominatingSet()));

        processed++;
        if ((processed % 100'000) == 0) {
            std::cout << processed << " graphs processed\n";
        }
    }
}

template <typename A, typename B>
void testPowerOfTwoAdjacency(bool exactlySame) {
    int dim = 6;
    int nNodes = 1 << dim;
    NetworKit::Graph graph(nNodes);
    for (int i = 0; i < dim; i++) {
        int diff = 1 << i;
        for (int v = 0; v < nNodes - diff; v++) {
            graph.addEdge(v, v + diff);
        }
    }

    A algorithmA(graph);
    B algorithmB(graph);

    {
        auto start = std::chrono::high_resolution_clock::now();
        algorithmA.run();
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << duration << "\n";
    }
    
    {
        auto start = std::chrono::high_resolution_clock::now();
        algorithmB.run();
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << duration << "\n";
    }
    assert(algorithmA.isDominating(algorithmA.getDominatingSet()));
    assert(algorithmB.isDominating(algorithmB.getDominatingSet()));
    assert(dsSize(algorithmA.getDominatingSet()) == dsSize(algorithmB.getDominatingSet()));
    if (exactlySame) {
        for (int i = 0; i < graph.numberOfNodes(); i++) {
            assert(algorithmA.getDominatingSet().at(i) == algorithmB.getDominatingSet().at(i));
        }
    }
}

void testReduceCalls() {
    int dim = 6;
    int nNodes = 1 << dim;
    NetworKit::Graph graph(nNodes);
    for (int i = 0; i < dim; i++) {
        int diff = 1 << i;
        for (int v = 0; v < nNodes - diff; v++) {
            graph.addEdge(v, v + diff);
        }
    }

    RooijBodlaenderMDS algorithmA(graph);
    BranchAndReduceMDS<RooijBodlaenderMSC> algorithmB(graph);

    {
        auto start = std::chrono::high_resolution_clock::now();
        algorithmA.run();
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << duration << "\n";
    }
    
    {
        auto start = std::chrono::high_resolution_clock::now();
        algorithmB.run();
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << duration << "\n";
    }
    assert(algorithmA.isDominating(algorithmA.getDominatingSet()));
    assert(algorithmB.isDominating(algorithmB.getDominatingSet()));
    assert(dsSize(algorithmA.getDominatingSet()) == dsSize(algorithmB.getDominatingSet()));

    for (int i = 0; i < graph.numberOfNodes(); i++) {
        assert(algorithmA.getDominatingSet().at(i) == algorithmB.getDominatingSet().at(i));
    }
    std::cout << MinimumDominatingSet::specialCounter1 << " vs " << MinimumDominatingSet::specialCounter2 << " calls\n";
    assert(MinimumDominatingSet::specialCounter1 == MinimumDominatingSet::specialCounter2);
}

// https://users.cecs.anu.edu.au/~bdm/data/graphs.html
int main() {
    auto fileNames = {"input/graph2.g6", "input/graph3.g6", "input/graph4.g6", "input/graph5.g6", "input/graph6.g6", "input/graph7.g6", "input/graph8.g6", "input/graph9.g6", "input/graph10.g6"};
    for (auto& name : fileNames) {
        std::cout << "run " << name << "\n";
        // test_algo_g6<FominKratschWoegingerMDS>(name);
        // test_algo_g6<SchiermeyerMDS>(name);
        test_algo_g6<FominKratschWoegingerMDS>(name);
    }

    // testPowerOfTwoAdjacency<RooijBodlaenderMDS, BranchAndReduceMDS<RooijBodlaenderMSC>>(true);
    // testReduceCalls();
    // testPowerOfTwoAdjacency<BranchAndReduceMDS<FominGrandoniKratschMSC>, BranchAndReduceMDS<RooijBodlaenderMSC>>(false);
    // testPowerOfTwoAdjacency<BranchAndReduceMDS<GrandoniMSC>, BranchAndReduceMDS<RooijBodlaenderMSC>>(false);
}


// void testSpecific() {
//     int nNodes = 10;
//     NetworKit::Graph graph(nNodes);
//     graph.addEdge(5, 0);
//     graph.addEdge(6, 1);
//     graph.addEdge(6, 2);
//     graph.addEdge(7, 1);
//     graph.addEdge(7, 3);
//     graph.addEdge(8, 0);
//     graph.addEdge(8, 2);
//     graph.addEdge(8, 4);
//     graph.addEdge(9, 3);
//     graph.addEdge(9, 4);
//     graph.addEdge(9, 5);
//     graph.addEdge(9, 8);
//     RooijBodlaenderMDS algorithm(graph);
//     // ExhaustiveMDS algorithm(graph);
//     auto start = std::chrono::high_resolution_clock::now();
//     algorithm.run();
//     auto end = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
//     std::cout << duration << "\n";
//     std::vector<bool> dominating_set = algorithm.getDominatingSet();
//     int ds_size = 0;
//     for (bool is_in_dominating_set : dominating_set) {
//         if (is_in_dominating_set) {
//             ds_size++;
//         }
//         std::cout << is_in_dominating_set;
//     } std::cout << " ";
//     std::cout << "with dominating set size = " << ds_size << "\n";
//     std::cout << (algorithm.isDominating(dominating_set) ? "correct\n" : "incorrect\n");
// }
