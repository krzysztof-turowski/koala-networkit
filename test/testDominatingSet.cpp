// #include <iostream>
#include <dominatingset/RooijBodlaenderMDS.hpp>
#include <dominatingset/ExhaustiveMDS.hpp>
// #include <chrono>
#include <io/G6GraphReader.hpp>
#include <fstream>

int dsSize(const std::vector<bool> &set) {
    return std::count(set.begin(), set.end(), true);
}

void test_graph_g6(const std::string inputFileName) {
    std::ifstream inputStream(inputFileName);
    assert(inputStream.good());
    std::string g6Graph;
    while (std::getline(inputStream, g6Graph)) {
        Koala::G6GraphReader reader;
        NetworKit::Graph graph = reader.readline(g6Graph);
        RooijBodlaenderMDS rooij(graph);
        ExhaustiveMDS exhaustive(graph);

        rooij.run();
        exhaustive.run();

        assert(rooij.isDominating(rooij.getDominatingSet()));
        assert(rooij.isDominating(rooij.getDominatingSet()));
        assert(dsSize(rooij.getDominatingSet()) == dsSize(exhaustive.getDominatingSet()));
    }
}

// https://users.cecs.anu.edu.au/~bdm/data/graphs.html
int main() {
    test_graph_g6("input/graph2.g6");
    test_graph_g6("input/graph3.g6");
    test_graph_g6("input/graph4.g6");
    test_graph_g6("input/graph5.g6");
    test_graph_g6("input/graph6.g6");
    test_graph_g6("input/graph7.g6");
    test_graph_g6("input/graph8.g6");
    test_graph_g6("input/graph9.g6");
    test_graph_g6("input/graph10.g6");
}

// void testPowerOfTwoAdjacency() {
//     int dim = 6;
//     int nNodes = 1 << dim;
//     NetworKit::Graph graph(nNodes);
//     for (int i = 0; i < dim; i++) {
//         int diff = 1 << i;
//         for (int v = 0; v < nNodes - diff; v++) {
//             graph.addEdge(v, v + diff);
//         }
//     }

//     RooijBodlaenderMDS algorithm(graph);
//     // ExhaustiveDominatingSet algorithm(graph);
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
